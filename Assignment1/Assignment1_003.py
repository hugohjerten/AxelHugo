# -*- coding: utf-8 -*-
"""
Created on Sun Apr 05 17:14:41 2015
@author: Hugo
"""
from pylab import *
from scipy import *
from scipy.misc import derivative
from scipy.optimize import fsolve

##----Assignment nr 1----##
   
## Task 1 a)##
def plotFunction(f):
    t = linspace(0, 2, 100)
    plot(t,f(t))
    grid(True)

## Task 1 b)
def bisec(f, a, b, tol):
    i = 0
    a = float(a)
    b = float(b)    
    if f(a)*f(b) > 0:
        raise Exception('This interval does not contain the root or this function is not only ascending or descending in this interval.')
    while True:
        i += 1        
        c = (a + b) / 2
        sizeNewInterval = abs(a - c)
        val = f(c)
        if val*f(a) > 0:
            a = c
        else:
            b = c
        if i > 1000:
            raise Exception('Loop too big... Something wrong!')
        if tol > sizeNewInterval or val == 0:
            break
    return c, i

## Task 1 c)
def func(x):
    return x**2 -2

def funcThreeRoots(x):
    return (x - 0.2)**3 + (x - 0.2)**2 - x + 0.2

def funcWithNoRoot(x):
    return x**2 + 1

## Task 3 c)
def iterationMethod(f,x0,tol):
    i = 0
    x0 = float(x0)
    while True:
        i += 1
        x1 = (x0 * (x0**2 + 6) )/ (3 * x0**2 + 2)
        changePrevValue = abs(x1 - x0)
#        print xnp1, f(xnp1), changePrevValue
        x0 = x1
        if tol > changePrevValue:
            if f(x0) < tol:
                break
        if i > 100:
            print 'With iterationMethod (which is the same as newton method??) a root was NOT found.'
            return
    print 'With iterationMethod (which is the same as newton method??) a root was found at {} with {} iterations.'.format(x0, i)
    return

## Task 4

def newton(f,xn,tol):
    i = 0    
    xn = float(xn)
    while True:
        i += 1
        val = f(xn)
        deriv = derivative(f,xn)
        xnp1 = xn - val/deriv
        changePrevValue = abs(xn-xnp1)
        xn = xnp1
        if tol > changePrevValue:
            if f(xn) < tol:
                break
        if i > 100:
            print 'With newton method a root was NOT found.'
            return
    print 'With newton method a root was found at {} with {} iterations.'.format(xn, i)
    return



#Main program/code

tolerance = 1.e-8
plotFunction(func)
show()
print 'For function func, x**2 -2:'
bisecReturn = bisec(func, 0, 2, tolerance)
print '--bisec method gives  {} in {} iterations.'.format(bisecReturn[0], bisecReturn[1])
print '--fsolve method gives {}.'.format(fsolve(func, 1.5)[0])
print 'For function funcThreeRoots, (x - 0.2)**3 + (x - 0.2)**2 - x + 0.2:'
bisecReturn = bisec(funcThreeRoots, -1, 0.5, tolerance)
print '--bisec method gives  {} in {} iterations.'.format(bisecReturn[0], bisecReturn[1])
print '--fsolve method gives {}.'.format(fsolve(funcThreeRoots, 0.3)[0])

iterationMethod(func, 50, tolerance)
iterationMethod(funcThreeRoots, 50, tolerance)
iterationMethod(funcWithNoRoot, 50, tolerance)

newton(func, 50, tolerance)
newton(funcThreeRoots, 50, tolerance)
newton(funcWithNoRoot, 50, tolerance)











