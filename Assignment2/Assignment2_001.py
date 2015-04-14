# -*- coding: utf-8 -*-
"""
Created on Sun Apr 05 17:14:41 2015
@author: Hugo
"""
from pylab import *
from scipy import *
from scipy.misc import derivative
from scipy.optimize import fsolve
from scipy.linalg import solve

##----Assignment nr 2----##

##---Task 1--####

def newton():
	return

def Jacobian(f, x, h =1.e-8, fx = []):
	"""
	Computes the Jacobian of f
	eps is the increment and fx a reference value f(x)
	"""
	if fx == []:
		fx = f(x)
	jac = empty([len(fx),len(x)])
	for i in range(len(x)):
		x[i] += h
		jac[:,i] = (f(x)-fx)/h
		x[i] -= h
	return jac

def robotArm(x):
	alfa = x[0]
	beta = x[1]
	xf = 5 * cos(alfa) + 6 * cos(alfa + beta)
	yf = 5 * sin(alfa) + 6 * sin(alfa + beta)
	return array([xf, yf])

def LUFactorize(A):
	if A.shape[0] != A.shape[1]:
		print 'Matrix A is not of the size n*n.'
		return
	n =len(A)
	L = zeros((n,n))
	U = zeros((n,n))
	for i in range(n):
		for j in range(n):
			if i == j:
				L[i][j] = 1
	print L	
	return

x = [0.7, 0.7]
print Jacobian(robotArm, x)
LUFactorize(Jacobian(robotArm, x))
#print robotArm(x)
