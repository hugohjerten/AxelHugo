# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 10:09:23 2015
@author: sebastianelm
"""

#from __future__ import division # 2st. Understreck på varsin sida
from scipy import *
from pylab import *
from scipy.linalg import solve

def cubeSpline(xint, yint):
	h = abs(xint[1]-xint[0])
	d = yint
	m = len(xint)
	coeff = zeros((m-1, 4))
	Y = zeros((m, 1))
	allSigmas = zeros((m, m))
	allSigmas[0][0] = 1
	allSigmas[m-1][m-1] = 1
	col = 0
	for row in range(1,m-1):
		allSigmas[row][col] = 1
		col += 1
		allSigmas[row][col] = 4
		col += 1
		allSigmas[row][col] = 1
		col -= 1
		Y[row] = (yint[row + 1] - 2*yint[row] + yint[row - 1])*6/h**2
	sigmas = solve(allSigmas, Y)
	b = sigmas/2
	a = []
	c = []
	for i in range(m-1):
		a.append((sigmas[i+1]-sigmas[i])/6*h)
		c.append((yint[i+1]-yint[i])/h - (2*sigmas[i] + sigmas[i+1])*h/6)
	for i in range(m-1):
		coeff[i][0] = a[i]
		coeff[i][1] = b[i]
		coeff[i][2] = c[i]
		coeff[i][3] = d[i]
	print 'Size of coeff is: {}'.format(len(coeff))
	print 'Coeff-matrix\n {}'.format(coeff)
	return coeff

def cubeSplineVal(coeff, xint, xval):
	if xval <= xint[0]:
		row = 0
		xdiff = xval - xint[row]
		return coeff[row][0]*xdiff**3 + coeff[row][1]*xdiff**2 + coeff[row][2]*xdiff + coeff[row][3]
	for i in range(1, len(xint)):
		if xval <= xint[i]:
			row = i - 1
			xdiff = xval - xint[row]
			return coeff[row][0]*xdiff**3 + coeff[row][1]*xdiff**2 + coeff[row][2]*xdiff + coeff[row][3]
	row = len(xint)-1
	xdiff = xval - xint[row]
	return coeff[row][0]*xdiff**3 + coeff[row][1]*xdiff**2 + coeff[row][2]*xdiff + coeff[row][3]

x = linspace(0,9,10)
print x
y = [4, 8, 3, 6, 7, 9, 2, 5, 6, 11]
S = cubeSpline(x, y)
for i in range(10):
	print cubeSplineVal(S, x, i)

#print cubeSplineVal(cubeSpline(x, y), x, 6)









