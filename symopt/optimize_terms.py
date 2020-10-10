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

def replace_inv(expr, var, var_inv):
    '''Accepts and expression, and replaces a specified variable and replaces
    it by its inverse where ever its inverse is used.
    
    Cases where replacement happens:
        
    >>> x, x_inv, y = symbols('x, x_inv, y')
    >>> replace_inv(x + 1/x, x, x_inv)
    x + x_inv
    
    >>> replace_inv(sin(x) + sin(1)/(5*x**1*y), x, x_inv)
    x_inv*sin(1)/(5*y) + sin(x)
    
    Cases where replacement does not happen
    
    >>> replace_inv(sin(x) + 1/sin(x), x, x_inv)
    sin(x) + 1/sin(x)
    
    
    '''
    new = 0
    find_pow_inv = str(var)+'**-'
    find_inv = str(var)
    
    for arg in expr.args:
        numer, denom = fraction(arg)
        str_denom = str(denom)
        if find_pow_inv in str_denom or find_inv in str_denom and not '(' in str_denom:
            arg = arg.subs(1/var, var_inv)
        new += arg
    return new


def replace_power_sqrts(expr, var):
    new = 0
    for arg in expr.args:
        factor, power = arg.as_coeff_exponent(var)
        if isinstance(power, Number):
            
            is05 = isclose(float(power %1), 0.5, rel_tol=1e-12)
            if isclose((power-.5) %1, 0, rel_tol=1e-12):
                new_power = int(power-.5)
            else:
                new_power = Rational("0.5")

    #         print(power, factor,  power %1 == 0.5)
            if is05:
                arg = factor*UnevaluatedExpr(sqrt(var))*var**(new_power)
        new += arg
    return new
