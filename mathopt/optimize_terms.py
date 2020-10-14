# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2020 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from __future__ import division

from sympy import *
from math import isclose
from sympy.core import Add, Mul, Number
from mathopt.addition_chain import minimum_addition_chain_multi_heuristic

__all__ = ['replace_inv', 'replace_power_sqrts', 'horner_expr',
           'optimize_expression_for_var', 'optimize_expression',
           'recursive_find_power', 'make_pow_sym', 'replace_intpowers',
           'replace_fracpowers',
           'integer_chain_symbolic_path', 'simplify_powers_as_fractions',
           'singleton_variables_inline']

def replace_inv(expr, var, assignments=None, expressions=None):
    '''Accepts and expression, and replaces a specified variable and replaces
    it by its inverse where ever its inverse is used.
    
    Cases where replacement happens:
        
    >>> x, y = symbols('x, y')
    >>> replace_inv(x + 1/x, x)[0]
    x + x_inv
    
    >>> replace_inv(sin(x) + sin(1)/(5*x**1*y), x)[0]
    x_inv*sin(1)/(5*y) + sin(x)
    
    >>> tau, delta = symbols('tau, delta')
    >>> expr = 0.000233594806142*delta**11*tau**3.25*exp(-delta**2)
    >>> replace_inv(expr, tau)[0]
    0.000233594806142*delta**11*tau**3.25*exp(-delta**2)

    Case where replacement wasn't happening because of a bracket:
        
    >>> tau, delta = symbols('tau, delta')
    >>> expr = 0.16*delta*tau**(3/5)*exp(-delta) +0.23*delta/tau**(67/100) - 0.008*delta**4/tau**(4/5)
    >>> expr = simplify_powers_as_fractions(expr, tau)
    >>> expr = simplify_powers_as_fractions(expr, delta)
    >>> replace_inv(expr, tau)[0]
    -0.008*delta**4*tau_inv**(4/5) + 0.16*delta*tau**(3/5)*exp(-delta) + 0.23*delta*tau_inv**(67/100)
        
    Cases where replacement does not happen
    
    >>> replace_inv(sin(x) + 1/sin(x), x)[0]
    sin(x) + 1/sin(x)
    
    >>> tau, tau_inv, delta = symbols('tau, tau_inv, delta')
    >>> expr = tau_inv*(tau_inv*(tau_inv*(4.20549538e-5 - 1.8171582e-7*tau_inv) + 0.000158860716) + 2.490888032)
    >>> replace_inv(expr, delta)[0]
    tau_inv*(tau_inv*(tau_inv*(4.20549538e-5 - 1.8171582e-7*tau_inv) + 0.000158860716) + 2.490888032)
    '''
    if assignments is None:
        assignments = []
    if expressions is None:
        expressions = []
    new = 0
    find_pow_inv = str(var)+'**-'
    find_inv = str(var)
    var_inv = symbols(var.name + '_inv') # Make it even if we don't need it

    def change_term(arg):
        numer, denom = fraction(arg)
        str_denom = str(denom)
        if find_pow_inv in str_denom or find_inv in str_denom: #and not '(' in str_denom:
            coeff, power = denom.as_coeff_exponent(var)
            arg = arg.replace(1/var**power, var_inv**power)
        return arg

    if isinstance(expr, Add):
        for arg in expr.args:
            new += change_term(arg)
    elif isinstance(expr, Mul):
        new = 1
        for arg in expr.args:
            new *= replace_inv(arg, var)[0]
    elif isinstance(expr, (Number, Pow, Function, Symbol)) or 1:
        new = expr
    else:
        new = 0
    
    if var_inv in new.free_symbols:
        assignments.append(var_inv)
        expressions.append(1.0/var)
    return new, assignments, expressions


def replace_power_sqrts(expr, var):
    '''
    >>> x, y = symbols('x, y')
    >>> replace_power_sqrts(x**Rational(3,2)*y, x)
    (x*xrt2*y, [xrt2], [sqrt(x)])
    >>> replace_power_sqrts(x**Rational(7,2)*y, x)
    (x**3*xrt2*y, [xrt2], [sqrt(x)])
    >>> replace_power_sqrts(x**35.5*y, x)
    (x**35*xrt2*y, [xrt2], [sqrt(x)])
    >>> replace_power_sqrts(x**-35.5*y, x)
    (xrt2*y/x**36, [xrt2], [sqrt(x)])
    >>> replace_power_sqrts(x**-.5*y, x)
    (xrt2inv*y, [xrt2inv], [1/sqrt(x)])
    >>> replace_power_sqrts(x**35.25*y, x)
    (x**35*xrt4*y, [xrt2, xrt4], [sqrt(x), sqrt(xrt2)])
    >>> replace_power_sqrts(x**.25*y, x)
    (xrt4*y, [xrt2, xrt4], [sqrt(x), sqrt(xrt2)])
    >>> replace_power_sqrts(x**.75*y, x)
    (xrt2*xrt4*y, [xrt2, xrt4], [sqrt(x), sqrt(xrt2)])
    
    >>> replace_power_sqrts(x**-.25*y, x)
    (xrt4inv*y, [xrt2, xrt4, xrt4inv], [sqrt(x), sqrt(xrt2), 1/xrt4inv])
    >>> replace_power_sqrts(x**-.75*y, x)
    (xrt34inv*y, [xrt2, xrt4, xrt34inv], [sqrt(x), sqrt(xrt2), 1/(xrt2*xrt4)])
    
    >>> replace_power_sqrts(x**1.5*y+ x**1.5, x)
    (x*xrt2*y + x*xrt2, [xrt2], [sqrt(x)])

    
    Case where replacement was not happening because of depth
    >>> delta, tau, tau_inv = symbols('delta, tau, tau_inv')
    >>> expr = delta*(0.1*delta**10*tau**(5/4)*exp(-delta**2) - 0.03*delta**5*tau_inv**(3/4)*exp(-delta))
    >>> replace_power_sqrts(expr, tau)
    (delta*(0.1*delta**10*tau*taurt4*exp(-delta**2) - 0.03*delta**5*tau_inv**0.75*exp(-delta)), [taurt2, taurt4], [sqrt(tau), sqrt(taurt2)])
    '''
    assignments = []
    expressions = []
    new = 0
    def newvar(var, root, suffix=''):
        name = var.name + 'rt' + str(root) + suffix
        sym = symbols(name)
        return sym
        
    def change_term(arg, assignments, expressions):
        factor, power = arg.as_coeff_exponent(var)
        if isinstance(power, Number):
            
            pow_float_rem = float(power %1)
            is05 = isclose(pow_float_rem, 0.5, rel_tol=1e-12)
            is025 = (power == -0.25 or isclose(pow_float_rem, 0.25, rel_tol=1e-12)) and not power == -0.75
            is075 = power == -0.75 or isclose(pow_float_rem, 0.75, rel_tol=1e-12)
            if is05:
                if power == -0.5:
                    new_power = 0
                else:
                    new_power = int(power - .5)
            elif is025:
                if power == -0.25:
                    new_power = 0
                else:
                    new_power = int(power - .25)
            elif is075:
                if power == -0.75:
                    new_power = 0
                else:
                    new_power = int(power - .75)
                
            if is05 or is025 or is075:
                if power == -0.5:
                    # Removing power completely
                    rtvar = newvar(var, 2, 'inv')
                    rtexpr = 1/sqrt(var)
                else:
                    rtvar = newvar(var, 2)
                    rtexpr = sqrt(var)
                if rtvar not in assignments:
                    assignments.append(rtvar)
                    expressions.append(rtexpr)

            if is025 or is075:
                rtexpr = sqrt(rtvar)
                rtvar = newvar(var, 4)
                if rtvar not in assignments:
                    assignments.append(rtvar)
                    expressions.append(rtexpr)
            if is025:
                if power == -0.25:
                    rtvar = newvar(var, 4, 'inv')
                    rtexpr = 1/rtvar
                    if rtvar not in assignments:
                        assignments.append(rtvar)
                        expressions.append(rtexpr)
                else:
                    pass
            elif is075:
                if power == -0.75:
                    rtexpr = 1/(rtvar*assignments[-2])
                    rtvar = newvar(var, 34, 'inv')
                    if rtvar not in assignments:
                        assignments.append(rtvar)
                        expressions.append(rtexpr)
                else:
                    rtvar = rtvar*assignments[-2]
                
            if is05 or is025 or is075:
                arg = factor*rtvar*var**(new_power)
        return arg
    
    if isinstance(expr, Add):
        for arg in expr.args:
            new += change_term(arg, assignments, expressions)
    elif isinstance(expr, Pow):
        new = change_term(expr, assignments, expressions)
    elif isinstance(expr, Mul):
        #new = change_term(expr)
        new = 1
        for arg in expr.args:
            to_mul, temp_assign, temp_expr = replace_power_sqrts(arg, var)
            new *= to_mul
            assignments += temp_assign
            expressions += temp_expr
    elif isinstance(expr, (Number, Pow, Function, Symbol)) or 1:
        return expr, [], []
    return new, assignments, expressions

def recursive_find_power(expr, var, powers=None, selector=lambda x: True):
    '''Recursively find all powers of `var` in `expr`. Optionally, a selection
    criteria such as only finding integers an be applied.

    Does not return 0 or 1 obviously.
    
    >>> x, y = symbols('x, y')
    >>> test = x**3*log(x**2)*sin(x**20)*y**5*exp(log(sin(x**8))) + y**3*x**15
    >>> list(sorted(list(recursive_find_power(test, x))))
    [2, 3, 8, 15, 20]
    >>> list(sorted(list(recursive_find_power(test, x, selector=lambda x: x > 3))))
    [8, 15, 20]
    >>> list(sorted(list(recursive_find_power(test,y))))
    [3, 5]
    >>> test = x**3.1*log(x**2.2)*sin(x**20.5)*y**5*exp(log(sin(x**8))) + y**3*x**15
    >>> list(sorted(list(recursive_find_power(test, x))))
    [2.20000000000000, 3.10000000000000, 8, 15, 20.5000000000000]
    >>> list(sorted(list(recursive_find_power(test, x, selector=lambda x: int(x) == x))))
    [8, 15]
    '''
    if powers is None:
        powers = set([])
    for arg in expr.args:
        coeff, exponent = arg.as_coeff_exponent(var)
        if isinstance(exponent, Number) and exponent != 0 and exponent != 1 and selector(exponent):
            powers.add(exponent)
        else:
            recursive_find_power(arg, var, powers, selector)
    return powers

def simplify_powers_as_fractions(expr, var, max_denom=1000):
    '''Takes an expression and replaces
    
    >>> x, y = symbols('x, y')
    >>> expr = x**.15 + x**.2 + x**.33 + x**.35 + x**.8 + x**1.01 + x**1.6
    >>> simplify_powers_as_fractions(expr, x)
    x**(101/100) + x**(33/100) + x**(7/20) + x**(3/20) + x**(8/5) + x**(4/5) + x**(1/5)
    >>> expr = x**.15
    >>> simplify_powers_as_fractions(expr, x)
    x**(3/20)
    
    >>> simplify_powers_as_fractions(x**2.15*sin(x**3.22)*y+y*x**20, x)
    x**(43/20)*y*sin(x**(161/50)) + x**20*y
    '''
    def change_term(arg):
        coeff, exponent = arg.as_coeff_exponent(var)
        if isinstance(exponent, Number) and exponent != 0 and exponent != 1 :
            exponent_simplified = nsimplify(exponent)
            if exponent_simplified.denominator() <= max_denom:
                return arg.replace(var**exponent, var**exponent_simplified)
            return arg
        else:
            if isinstance(arg, Mul):
                base = 1
                for a in arg.args:
                    base *= simplify_powers_as_fractions(a, var, max_denom)
                return base
            return arg


    if isinstance(expr, Add):
        base = 0
        for i in range(len(expr.args)):
            base += change_term(expr.args[i])
        return base
    elif isinstance(expr, Mul) or isinstance(expr, Pow):
        return change_term(expr)
    elif isinstance(expr, Function):
        return type(expr)(*(simplify_powers_as_fractions(v, var, max_denom) for v in expr.args))
    else:
        return expr

#def convert_numbers_to_floats(expr):
#    '''
#    
#    >>> x, y = symbols('x, y')
#    >>> expr = Rational("1.5")
#    >>> convert_numbers_to_floats(expr)
#    1.5
#    >>> expr = sin(Rational("1.5")*x)
#    '''
#    def change_term(arg):
#        if isinstance(arg, Number):
#            return float(arg)
#        else:
#            return arg
#
#
#    if isinstance(expr, Add):
#        base = 0
#        for i in range(len(expr.args)):
#            base += change_term(expr.args[i])
#        return base
#    elif isinstance(expr, Mul) or isinstance(expr, Pow):
#        return change_term(expr)
#    elif isinstance(expr, Function):
#        return type(expr)(*(convert_numbers_to_floats(v) for v in expr.args))
#    elif isinstance(expr, Number):
#        return float(expr)
#    else:
#        return expr

def horner_expr(expr, var):
    '''Basic wrapper around sympy's horner which does not raise an exception if
    there is nothing to do.
    
    >>> x = symbols('x')
    >>> horner_expr(x**3 + x**2 + x**1 + x, x)
    x*(x*(x + 1) + 2)
    
    Case where horner's method does not work:
    
    >>> horner_expr(x**3 + x**2 + x**1 + x + 1/x, x)
    x**3 + x**2 + 2*x + 1/x
    '''
    try:
        expr = horner(expr, var)
    except Exception as e:
        pass
    return expr

def make_pow_sym(var, power, suffix=''):
    '''Create a new symbol for a specified symbol.
    
    >>> x = symbols('x')
    >>> make_pow_sym(x, 100)
    x100
    >>> make_pow_sym(x, 1)
    x
    '''
    if power == 1:
        name = var.name + suffix
    else:
        name = var.name + str(power) + suffix
    sym = symbols(name)
    return sym

def integer_chain_symbolic_path(chain, var, suffix=''):
    '''Returns a tuple of assignments, expressions which can be joined together
    to calculate all of the necessary powers for an operation.
    
    Although this function returns UnevaluatedExprs, they can be simplified
    removed with the simplify() function. 
    
    >>> x = symbols('x')
    >>> chain = [[2], [2, 3], [2, 3, 5, 10, 13], [2, 3, 5, 10, 20]]
    >>> integer_chain_symbolic_path(chain, x)
    ([x2, x3, x5, x10, x13, x20], [x*x, x*x2, x2*x3, x5*x5, x3*x10, x10*x10])
    '''
    assignments = []
    expressions = []
    for l in chain:
        for i, v in enumerate(l):
            to_add_asign = make_pow_sym(var, v, suffix)
            if i == 0:
                assert v == 2
                to_add_expr = UnevaluatedExpr(make_pow_sym(var, 1, suffix))*make_pow_sym(var, 1, suffix)
            else:
                prev = l[i-1]
                delta = v-l[i-1]
                to_add_expr = UnevaluatedExpr(make_pow_sym(var, prev, suffix))*make_pow_sym(var, delta, suffix)
            
            if to_add_asign not in assignments:
                assignments.append(to_add_asign)
                expressions.append(to_add_expr)
    return assignments, expressions
    
def replace_intpowers(expr, var):
    '''
    >>> x, y = symbols('x, y')
    >>> test = x**2*sin(x**3)*y+y*x**20
    >>> replace_intpowers(test, x)[0]
    x2*y*sin(x3) + x20*y
    >>> replace_intpowers(test, x)[1]
    [x2, x3, x5, x10, x20]
    >>> replace_intpowers(test, x)[2]    
    [x*x, x*x2, x2*x3, x5*x5, x10*x10]
    >>> replace_intpowers(test, y)[0]
    x**20*y + x**2*y*sin(x**3)
    '''
    powers = list(sorted(list(recursive_find_power(expr, var, selector=lambda x: int(x) == x))))
    powers_int = [int(i) for i in powers]
    chain_length, chain = minimum_addition_chain_multi_heuristic(powers_int, small_chain_length=0)
    assignments, expressions = integer_chain_symbolic_path(chain, var)
    replacement_vars = [make_pow_sym(var, p) for p in powers]
    for power, replacement in zip(powers[::-1], replacement_vars[::-1]):
        # iterate from highest to low
        expr = expr.replace(var**power, replacement)
    return expr, assignments, expressions

def replace_fracpowers(expr, var):
    '''
    >>> x, y = symbols('x, y')
    >>> test = x**2.15*sin(x**3.22)*y+y*x**20
    >>> test = simplify_powers_as_fractions(test, x)
    >>> replace_fracpowers(test, x)
    (x**20*y + x215_100*y*sin(x322_100), [x_100, x2_100, x4_100, x5_100, x10_100, x20_100, x40_100, x45_100, x85_100, x170_100, x215_100, x80_100, x160_100, x320_100, x322_100], [x**0.01, x_100*x_100, x2_100*x2_100, x_100*x4_100, x5_100*x5_100, x10_100*x10_100, x20_100*x20_100, x5_100*x40_100, x40_100*x45_100, x85_100*x85_100, x45_100*x170_100, x40_100*x40_100, x80_100*x80_100, x160_100*x160_100, x2_100*x320_100])

    >>> tau, delta = symbols('tau, delta')
    >>> test = - 0.042053322884200002*delta**4*tau**0.200000000000000011 + 0.0349008431981999989*delta**4*tau**0.349999999999999978 
    >>> test = simplify_powers_as_fractions(test, tau)
    >>> test = simplify_powers_as_fractions(test, delta)
    >>> replace_fracpowers(test, tau)
    (-0.0420533228842*delta**4*tau4_20 + 0.0349008431982*delta**4*tau7_20, [tau_20, tau2_20, tau4_20, tau6_20, tau7_20], [tau**0.05, tau_20*tau_20, tau2_20*tau2_20, tau2_20*tau4_20, tau_20*tau6_20])

    Test case of one power
    
    >>> tau, delta, tau_inv, delta_inv = symbols('tau, delta, tau_inv, delta_inv')
    >>> expr = 0.16*delta*tau**(3/5)*exp(-delta) +0.23*delta/tau**(67/100) - 0.008*delta**4/tau**(4/5)
    >>> expr = simplify_powers_as_fractions(expr, tau)
    >>> expr = simplify_powers_as_fractions(expr, delta)
    >>> expr = replace_inv(expr, tau)[0]
    >>> replace_fracpowers(expr, tau)
    (-0.008*delta**4*tau_inv**(4/5) + 0.16*delta*tau**(3/5)*exp(-delta) + 0.23*delta*tau_inv**(67/100), [], [])
    >>> replace_fracpowers(expr, tau_inv)
    (-0.008*delta**4*tau_inv80_100 + 0.16*delta*tau**(3/5)*exp(-delta) + 0.23*delta*tau_inv67_100, [tau_inv_100, tau_inv2_100, tau_inv4_100, tau_inv8_100, tau_inv16_100, tau_inv32_100, tau_inv64_100, tau_inv66_100, tau_inv67_100, tau_inv80_100], [tau_inv**0.01, tau_inv_100*tau_inv_100, tau_inv2_100*tau_inv2_100, tau_inv4_100*tau_inv4_100, tau_inv8_100*tau_inv8_100, tau_inv16_100*tau_inv16_100, tau_inv32_100*tau_inv32_100, tau_inv2_100*tau_inv64_100, tau_inv_100*tau_inv66_100, tau_inv16_100*tau_inv64_100])

    >>> expr = - 1.6*delta*tau**(1/100) - 0.5*delta*tau**(13/5)*exp(-delta**2) - 0.16*delta*tau**(3/5)*exp(-delta) 
    >>> expr = simplify_powers_as_fractions(expr, tau)
    >>> replace_fracpowers(expr, tau)
    (-0.5*delta*tau260_100*exp(-delta**2) - 0.16*delta*tau60_100*exp(-delta) - 1.6*delta*tau_100, [tau_100, tau2_100, tau4_100, tau8_100, tau16_100, tau20_100, tau40_100, tau60_100, tau32_100, tau64_100, tau128_100, tau256_100, tau260_100], [tau**0.01, tau_100*tau_100, tau2_100*tau2_100, tau4_100*tau4_100, tau8_100*tau8_100, tau4_100*tau16_100, tau20_100*tau20_100, tau20_100*tau40_100, tau16_100*tau16_100, tau32_100*tau32_100, tau64_100*tau64_100, tau128_100*tau128_100, tau4_100*tau256_100])

    '''
    fractional_powers = recursive_find_power(expr, var, selector=lambda x: int(x) != x and abs(x)%.25 != 0)
    if not fractional_powers or len(fractional_powers) == 1:
        return expr, [], []
    fractional_powers = list(sorted(list(fractional_powers)))
    base_power = gcd(fractional_powers)
    powers_int = [int(i/base_power) for i in fractional_powers]
    powers_int = [i for i in powers_int if i != 1] # Remove the base_power if it appears as it is handled separately
    prefix = '_' + str(base_power.numerator()) 
    suffix = '_' + str(base_power.denominator())
    var_suffix = symbols(var.name + suffix)
    
    chain_length, chain = minimum_addition_chain_multi_heuristic(powers_int, small_chain_length=0)
    assignments, expressions = integer_chain_symbolic_path(chain, var, suffix)
    replacement_vars = [make_pow_sym(var, p, suffix) for p in powers_int]
    
#    subs = {var**power: replacement for power, replacement in zip(fractional_powers[::-1], replacement_vars[::-1])}
#    subs = {var**power: replacement for power, replacement in zip(fractional_powers, replacement_vars)}
#    expr = expr.subs(subs, simultaneous=True)
    for power, replacement in zip(fractional_powers[::-1], replacement_vars[::-1]):
#        iterate from highest to low
        expr = expr.replace(var**power, replacement)

    # Handle the case the base power is in there already
    expr = expr.replace(var**base_power, var_suffix)
        
    assignments.insert(0, var_suffix)
    expressions.insert(0, var**float(base_power))
    return expr, assignments, expressions

    
def singleton_variables_inline(assignments, expressions, expr):
    '''Replaces variables which are used only once by putting them right in
    the final expression, so they are never stored.
    
    >>> delta2, delta4, delta8, delta10, taurt2, taurt4, delta, tau = symbols('delta2, delta4, delta8, delta10, taurt2, taurt4, delta, tau')
    >>> assignments = [delta2, delta4, delta8, delta10, taurt2, taurt4]
    >>> expressions = [delta*delta, delta2*delta2, delta4*delta4, delta2*delta8, sqrt(tau), sqrt(taurt2)]
    >>> expr = delta10*tau*taurt4*(0.018 - 0.0034*delta2)*exp(-delta2)
    >>> singleton_variables_inline(assignments, expressions, expr)
    ([delta2, delta4, delta8, taurt2], [delta**2, delta2**2, delta4**2, sqrt(tau)], delta2*delta8*tau*sqrt(taurt2)*(0.018 - 0.0034*delta2)*exp(-delta2))
    '''
    pow_count = str(expr).count('**')
    assignment_use_in_expressions = []
    assignment_use_in_expr = []
    
    new_assignments = []
    new_expressions = []
    for v in assignments:
        assignment_use_in_expr.append(expr.count(v))
        assignment_use_in_expressions.append(sum(token.count(v) for token in expressions))
    
    for assignment, expression, count_expressions, count_expr in zip(assignments, expressions, assignment_use_in_expressions, assignment_use_in_expr):
        # This code won't work because sympy will consolidate terms
#         if count_expr + count_expressions > 1:
#             new_assignments.append(assignment)
#             new_expressions.append(expression)
#         elif count_expressions == 1:
#             for i in range(len(expressions)):
#                 expressions[i] = expressions[i].replace(assignment, expression)
#             for i in range(len(new_expressions)):
#                 new_expressions[i] = new_expressions[i].replace(assignment, expression)
#         elif count_expr == 1:
#             expr = expr.replace(assignment, expression)

        # This implementation only removes wasted things from the out expression
        if count_expr == 1 and count_expressions == 0:
            expr_tmp = expr.replace(assignment, expression)
            pow_count_tmp = str(expr_tmp).count('**')
            
            if pow_count_tmp > pow_count:
                # Abort! we accidentally caused a power
                new_assignments.append(assignment)
                new_expressions.append(expression)
            else:
                pow_count = pow_count_tmp
                expr = expr_tmp
        else:
            new_assignments.append(assignment)
            new_expressions.append(expression)
    return new_assignments, new_expressions, expr

def optimize_expression_for_var(expr, var, horner=True, intpows=True, fracpows=True):
    var_inv = symbols(var.name + '_inv') # Make it even if we don't need it

    expr, assignments, expressions = replace_inv(expr, var)
    
    expr, assign_tmp, expr_tmp = replace_power_sqrts(expr, var)
    assignments += assign_tmp
    expressions += expr_tmp
    expr, assign_tmp, expr_tmp = replace_power_sqrts(expr, var_inv)
    assignments += assign_tmp
    expressions += expr_tmp
    if horner:
        if var in expr.free_symbols:
            expr = horner_expr(expr, var)
        if var_inv in expr.free_symbols:
            expr = horner_expr(expr, var_inv)
    if intpows:
        expr, assign_tmp, expr_tmp = replace_intpowers(expr, var)
        assignments += assign_tmp
        expressions += expr_tmp
        expr, assign_tmp, expr_tmp = replace_intpowers(expr, var_inv)
        assignments += assign_tmp
        expressions += expr_tmp
    if fracpows:
        expr, assign_tmp, expr_tmp = replace_fracpowers(expr, var)
        assignments += assign_tmp
        expressions += expr_tmp
        expr, assign_tmp, expr_tmp = replace_fracpowers(expr, var_inv)
        assignments += assign_tmp
        expressions += expr_tmp
        
    return expr, assignments, expressions

def optimize_expression(expr, variables, horner=True,
                        intpows=True, fracpows=True):
    '''
    >>> tau, delta, tau_inv, delta_inv = symbols('tau, delta, tau_inv, delta_inv')
    >>> expr = 17.2752665749999998*tau - 0.000195363419999999995*tau**1.5 + log(delta) + 2.49088803199999997*log(tau) + 0.791309508999999967*log(1 - exp(-25.36365*tau)) + 0.212236767999999992*log(1 - exp(-16.90741*tau)) - 0.197938903999999999*log(exp(87.31279*tau) + 0.666666666666667) - 13.8419280760000003 - 0.000158860715999999992/tau - 0.0000210274769000000003/tau**2 + 6.05719400000000021e-8/tau**3
    >>> optimize_expression(expr, [tau,delta])[0]
    -0.00019536342*tau*taurt2 + 17.275266575*tau + tau_inv*(tau_inv*(6.057194e-8*tau_inv - 2.10274769e-5) - 0.000158860716) + log(delta) + 2.490888032*log(tau) + 0.791309509*log(1 - exp(-25.36365*tau)) + 0.212236768*log(1 - exp(-16.90741*tau)) - 0.197938904*log(exp(87.31279*tau) + 0.666666666666667) - 13.841928076
    
    >>> expr =  2.490888032/tau + 0.000158860716/tau**2 + 4.20549538e-5/tau**3 - 1.8171582e-7/tau**4
    >>> optimize_expression(expr, [tau,delta])[0]
    tau_inv*(tau_inv*(tau_inv*(4.20549538e-5 - 1.8171582e-7*tau_inv) + 0.000158860716) + 2.490888032)
    
    >>> tau, delta, tau_inv, delta_inv = symbols('tau, delta, tau_inv, delta_inv')
    >>> expr = delta*(0.1*delta**10*tau**(5/4)*exp(-delta**2) - 0.03*delta**5*tau_inv**(3/4)*exp(-delta))
    >>> optimize_expression(expr, [delta, tau], horner=False)
    (delta*(0.1*delta5*tau*sqrt(taurt2)*exp(-delta2)*delta5 - 0.03*delta5*tau_invrt2*tau_invrt4*exp(-delta)), [delta2, delta4, delta5, tau_inv, taurt2, tau_invrt2, tau_invrt4], [delta*delta, delta2*delta2, delta*delta4, 1.0/tau, sqrt(tau), sqrt(tau_inv), sqrt(tau_invrt2)])
    '''
    assignments = []
    expressions = []
    
    for var in variables:
        expr = simplify_powers_as_fractions(expr, var)
    
    for var in variables:
        expr, assign_tmp, expr_tmp = optimize_expression_for_var(expr, var, horner=horner, intpows=intpows, fracpows=fracpows)
        assignments += assign_tmp
        expressions += expr_tmp
        
    assignments, expressions, expr = singleton_variables_inline(assignments, expressions, expr)
    return expr, assignments, expressions
