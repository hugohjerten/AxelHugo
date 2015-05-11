
from scipy import *
from pylab import *
from scipy.linalg import solve
from numpy import *

def simpsons(f, a, b, n):
	h = abs(b-a)/n
	print 'Stepsize = {}'.format(h)
	tmp = 0
	limit = int(n)
	for j in range(1,limit):
		xj = a + j*h
		print xj
		tmp += 2.0*f(xj) + 4.0*f(xj + h/2.0)
		print tmp
	integral = (h/6.0)*(f(a) + 4.0*f(a + h/2.0) + f(b) + tmp)
	
	
	return integral 

def threePointGauss(f, a, b, n):
	h = abs(b-a)/n
	tmp = 0
	limit = int(n)
	for j in range(0, limit):
		xj = a + j*h
		tmp1 = f(xj + h*(5.0 + math.sqrt(15.0))/10.0)
		tmp2 = f(xj + 5.0*h/10.0)
		tmp3 = f(xj + h*(5.0 + math.sqrt(15.0))/10.0)
		tmp += 5.0*tmp1 + 8.0*tmp2 + 5.0*tmp3
	integral = h*tmp/18.0
	
	return

def xsquare(x):
	return x**2

def xcube(x):
	return x**3


print 'Simpons integral of x^2 between 0 and 3 = {}'.format(simpsons(xsquare, 0.0, 3.0, 1.0))

print '3PointGauss integral of x^2 between 0 and 3 = {}'.format(simpsons(xsquare, 0.0, 3.0, 1.0))

print 'Simpons integral of x^3 between 0 and 3 = {}'.format(simpsons(xcube, 0.0, 3.0, 1.0))

print '3PointGauss integral of x^3 between 0 and 3 = {}'.format(simpsons(xcube, 0.0, 3.0, 1.0))



