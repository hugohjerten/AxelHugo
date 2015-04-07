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
	iterations = []
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
		iterations.append(c)
        	if val*f(a) > 0:
            		a = c
        	else:
        		b = c
        	if i > 1000:
            		raise Exception('Loop too big... Something wrong!')
        	if tol > sizeNewInterval or val == 0:
            		break
    	return c, i, iterations

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
	iterations = []
	iterations.append(x0)
    	while True:
        	i += 1
        	x1 = (x0 * (x0**2 + 6) )/ (3 * x0**2 + 2)
        	changePrevValue = abs(x1 - x0)
        	x0 = x1
		iterations.append(x0)
        	if tol > changePrevValue:
			if f(x0) < tol:
				break
        		if i > 100:
            			print 'With iterationMethod (which is the same as newton method??) a root was NOT found.'
        			return
    	print 'With iterationMethod (which is the same as newton method??) a root was found at {} with {} iterations.'.format(x0, i)
    	return x0, i, iterations

## Task 4

def newton(f,xn,tol):
    	i = 0    
    	xn = float(xn)
	iterations = []
	iterations.append(xn)
    	while True:
        	i += 1
        	val = f(xn)
        	deriv = derivative(f,xn)
        	xnp1 = xn - val/deriv
        	changePrevValue = abs(xn-xnp1)
        	xn = xnp1
		iterations.append(xn)
        	if tol > changePrevValue:
            		if f(xn) < tol:
                		print 'With newton method a root was found at {} with {} iterations.'.format(xn, i)
				break
        	if i > 100:
            		print 'Session timed out with newton method, or root was NOT found.'
            		return
    	return xn, i, iterations
	
def secant(f, x0, x1, tol):
	i = 0
	x0 = float(x0)
	x1 = float(x1)
	iterations = []
	iterations.extend([x0,x1])
	while True:
		i += 1
		val0 = f(x0)
		val1 = f(x1)
		x2 = x1 - (val1 * (x1 - x0))/(val1 - val0)
		x0 = x1
		x1 = x2
		iterations.append(x1)
		interval = abs(x1-x0)
		if tol > interval or f(x1) < tol:
			print 'With secant method a root was found at {} with {} iterations.'.format(x1, i)
			break
		if i > 100:
			print 'Session timed out with secant method, or root was NOT found.'
			return
	return x1, i, iterations

## bisec function without Exception when root not found

def bisec2(f, a, b, tol):
        iterations = []
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
                iterations.append(c)
                if val*f(a) > 0:
                        a = c
                else:
                        b = c
                if i > 1000:
                        print 'Session timed out with the bisect method, or a root was NOT found.'
                        return
                if tol > sizeNewInterval or val == 0:
                        print 'With bisec method a root was found at {} with {} iterations'.format(c,i)
			break
        return c, i, iterations

#Task 5

def iterationPlot(name, iterations):
	plot(iterations)
	grid(True)
	title(name)
	xlabel('Iterations')
	ylabel('Value')
	show()
	return

def plotIncrementApproximation(name, iterations):
	tol = 1.e-14
	values = []
	i = 1
	val = 1
	while True:
		if i >= len(iterations) or val < tol:
			break
		val = abs(iterations[i] - iterations[i - 1])
		values.append(val)
		i += 1
	return

#Main program/code
tolerance = 1.e-8

mainmenu = str(input('For task 1-4 please press 1, for selection of specific iteration with details (task 5) press 2: '))
while True:
	if mainmenu == str(1):

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
		
		secant(func, 5, 10, tolerance)
		secant(func, -10, 0, tolerance)
		secant(funcThreeRoots, 5, 10, tolerance)
		secant(funcWithNoRoot, 5, 10, tolerance)
		break
	elif mainmenu == str(2):
		selection = str(input('For bisection press 1, for unknown iteration method press 2, for newton press 3, for secant press 4: '))
		while True:
			if selection == 1:
				bisecReturn = bisec2(func, 0, 2, tolerance)
				iterationPlot('Bisec method', bisecReturn[2])
				break
			elif selection == 2:
				break
			elif selection == 3:
				break
			elif selection == 4:
				break
			else:
				selection = input('Please type a value between 1 and 4: ')
		break
	else:
		mainmenu = str(input('Please type 1 or 2: '))
print 'Iterations complete.'








