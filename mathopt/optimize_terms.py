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
           'integer_chain_symbolic_path', 'simplify_powers_as_fractions']

def replace_inv(expr, var, var_inv):
    '''Accepts and expression, and replaces a specified variable and replaces
    it by its inverse where ever its inverse is used.
    
    Cases where replacement happens:
        
    >>> x, x_inv, y = symbols('x, x_inv, y')
    >>> replace_inv(x + 1/x, x, x_inv)
    x + x_inv
    
    >>> replace_inv(sin(x) + sin(1)/(5*x**1*y), x, x_inv)
    x_inv*sin(1)/(5*y) + sin(x)
    
    >>> tau, delta, tau_inv = symbols('tau, delta, tau_inv')
    >>> expr = 0.000233594806142*delta**11*tau**3.25*exp(-delta**2)
    >>> replace_inv(expr, tau, tau_inv)
    0.000233594806142*delta**11*tau**3.25*exp(-delta**2)
    
    Cases where replacement does not happen
    
    >>> replace_inv(sin(x) + 1/sin(x), x, x_inv)
    sin(x) + 1/sin(x)
    
    
    
    
    '''
    new = 0
    find_pow_inv = str(var)+'**-'
    find_inv = str(var)

    def change_term(arg):
        numer, denom = fraction(arg)
        str_denom = str(denom)
        if find_pow_inv in str_denom or find_inv in str_denom and not '(' in str_denom:
            arg = arg.subs(1/var, var_inv)
        return arg

    if isinstance(expr, Add):
        for arg in expr.args:
            new += change_term(arg)
    elif isinstance(expr, Mul):
        new = change_term(expr)
    return new


def replace_power_sqrts(expr, var):
    '''
    >>> x, y = symbols('x, y')
    >>> replace_power_sqrts(x**Rational(3,2)*y, x)
    x*y*sqrt(x)
    >>> replace_power_sqrts(x**Rational(7,2)*y, x)
    x**3*y*sqrt(x)
    >>> replace_power_sqrts(x**35.5*y, x)
    x**35*y*sqrt(x)
    >>> replace_power_sqrts(x**-35.5*y, x)
    y*sqrt(x)/x**36
    >>> replace_power_sqrts(x**-.5*y, x)
    y/sqrt(x)
    >>> replace_power_sqrts(x**35.25*y, x)
    x**35*y*sqrt(sqrt(x))
    >>> replace_power_sqrts(x**.25*y, x)
    y*sqrt(sqrt(x))
    >>> replace_power_sqrts(x**.75*y, x)
    y*sqrt(sqrt(x))*sqrt(x)
    
    Can't figure out how to make -1 be a division
    
    >>> replace_power_sqrts(x**-.25*y, x)
    y*(sqrt(sqrt(x)))**(-1)
    >>> replace_power_sqrts(x**-.75*y, x)
    y*(sqrt(sqrt(x))*sqrt(x))**(-1)
    '''
    new = 0
    def change_term(arg):
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
                
            if is05:
                if power == -0.5:
                    # Removing power completely
                    unev = 1/sqrt(var)
                else:
                    unev = UnevaluatedExpr(sqrt(var))
            elif is025:
                if power == -0.25:
                    unev = UnevaluatedExpr(1/UnevaluatedExpr(sqrt(UnevaluatedExpr(sqrt(var)))))
                else:
                    unev = UnevaluatedExpr(sqrt(UnevaluatedExpr(sqrt(var))))
            elif is075:
                if power == -0.75:
                    unev = UnevaluatedExpr(1/((UnevaluatedExpr(sqrt(UnevaluatedExpr(sqrt(var))))*UnevaluatedExpr(sqrt(var)))))
                else:
                    unev = UnevaluatedExpr(sqrt(UnevaluatedExpr(sqrt(var))))*UnevaluatedExpr(sqrt(var))
                
            if is05 or is025 or is075:
                arg = factor*unev*var**(new_power)
        return arg
    
    if isinstance(expr, Add):
        for arg in expr.args:
            new += change_term(arg)
    elif isinstance(expr, Mul):
        new = change_term(expr)
    return new

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
    
    >>> x = symbols('x')
    >>> expr = x**.15 + x**.2 + x**.33 + x**.35 + x**.8 + x**1.01 + x**1.6
    >>> simplify_powers_as_fractions(expr, x)
    x**(101/100) + x**(33/100) + x**(7/20) + x**(3/20) + x**(8/5) + x**(4/5) + x**(1/5)
    >>> expr = x**.15
    >>> simplify_powers_as_fractions(expr, x)
    x**(3/20)
    '''
    def change_term(arg):
        coeff, exponent = arg.as_coeff_exponent(var)
        exponent_simplified = nsimplify(exponent)
        if exponent_simplified.denominator() <= max_denom:
            return expr.replace(var**exponent, var**exponent_simplified)
        return expr
    
    if isinstance(expr, Add):
        for i in range(len(expr.args)):
            expr = change_term(expr.args[i])
    elif isinstance(expr, Mul) or isinstance(expr, Pow):
        expr = change_term(expr)
    return expr

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

def make_pow_sym(var, power):
    '''Create a new symbol for a specified symbol.
    
    >>> x = symbols('x')
    >>> make_pow_sym(x, 100)
    x100
    >>> make_pow_sym(x, 1)
    x
    '''
    if power == 1:
        return var
    name = var.name + str(power)
    sym = symbols(name)
    return sym

def integer_chain_symbolic_path(chain, var):
    '''Returns a tuple of assignments, expressions which can be joined together
    to calculate all of the necessary powers for an operation.
    
    >>> x = symbols('x')
    >>> chain = [[2], [2, 3], [2, 3, 5, 10, 13], [2, 3, 5, 10, 20]]
    >>> integer_chain_symbolic_path(chain, x)
    ([x2, x3, x5, x10, x13, x20], [x*x, x*x2, x2*x3, x5*x5, x3*x10, x10*x10])
    '''
    assignments = []
    expressions = []
    for l in chain:
        for i, v in enumerate(l):
            to_add_asign = make_pow_sym(var, v)
            if i == 0:
                assert v == 2
                to_add_expr = UnevaluatedExpr(var)*var
            else:
                prev = l[i-1]
                delta = v-l[i-1]
                to_add_expr = UnevaluatedExpr(make_pow_sym(var, prev))*make_pow_sym(var, delta)
            
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
    chain_length, chain = minimum_addition_chain_multi_heuristic(powers_int, small_chain_length=6)
    assignments, expressions = integer_chain_symbolic_path(chain, var)
    replacement_vars = [make_pow_sym(var, p) for p in powers]
    for power, replacement in zip(powers[::-1], replacement_vars[::-1]):
        # iterate from highest to low
        expr = expr.subs(var**power, replacement)
    return expr, assignments, expressions


def optimize_expression_for_var(expr, var, var_inv, horner=True, intpows=True):
    expr = replace_inv(expr, var, var_inv)
    if horner:
        expr = horner_expr(expr, var)
        expr = horner_expr(expr, var_inv)
    expr = replace_power_sqrts(expr, var)
    assignments = []
    expressions = []
    if intpows:
        expr, assign_tmp, expr_tmp = replace_intpowers(expr, var)
        assignments += assign_tmp
        expressions += expr_tmp
        expr, assign_tmp, expr_tmp = replace_intpowers(expr, var_inv)
        assignments += assign_tmp
        expressions += expr_tmp
    return expr, assignments, expressions

def optimize_expression(expr, variables, inverse_variables, horner=True,
                        intpows=True):
    '''
    >>> tau, delta, tau_inv, delta_inv = symbols('tau, delta, tau_inv, delta_inv')
    >>> expr = 17.2752665749999998*tau - 0.000195363419999999995*tau**1.5 + log(delta) + 2.49088803199999997*log(tau) + 0.791309508999999967*log(1 - exp(-25.36365*tau)) + 0.212236767999999992*log(1 - exp(-16.90741*tau)) - 0.197938903999999999*log(exp(87.31279*tau) + 0.666666666666667) - 13.8419280760000003 - 0.000158860715999999992/tau - 0.0000210274769000000003/tau**2 + 6.05719400000000021e-8/tau**3
    >>> optimize_expression(expr, [tau,delta], [tau_inv, delta_inv])[0]
    17.275266575*tau - 0.00019536342*tau*sqrt(tau) + tau_inv*(tau_inv*(6.057194e-8*tau_inv - 2.10274769e-5) - 0.000158860716) + log(delta) + 2.490888032*log(tau) + 0.791309509*log(1 - exp(-25.36365*tau)) + 0.212236768*log(1 - exp(-16.90741*tau)) - 0.197938904*log(exp(87.31279*tau) + 0.666666666666667) - 13.841928076
    
    '''
    assignments = []
    expressions = []
    
    for var in variables:
        expr = simplify_powers_as_fractions(expr, var)
    
    for var, var_inv in zip(variables, inverse_variables):
        expr, assign_tmp, expr_tmp = optimize_expression_for_var(expr, var, var_inv, horner=horner, intpows=intpows)
        assignments += assign_tmp
        expressions += expr_tmp
    return expr, assignments, expressions
