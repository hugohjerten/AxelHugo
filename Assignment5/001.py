
from scipy import *
from pylab import *
from scipy.linalg import solve
from numpy import *

def simpsons(f, a, b, n):
	h = abs(b-a)/n
	tmp = 0
	for j in range(1,n):
		xj = a + j*h
		tmp += 2*f(xj) + 4*f(xj + h/2)
	integral = (h/6)*(f(a) + 4*f(a + h/2) + f(b) + tmp)
	
	
	return integral 

def 3pointGauss(f, a, b, n):
	h = abs(b-a)/n
	tmp = 0
	for j in range(0,n):
		xj = a + j*h
		tmp1 = f(xj + h*(5 + math.sqrt(15))/10)
		
	return

