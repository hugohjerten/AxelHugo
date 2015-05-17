from scipy import *
from pylab import *
#from scipy.linalg import solve
from numpy import *




def explicitLinearTest(tol, eigenvalue):
	h = 5
	cnt = 0
	while True:
		h -= tol
		if abs(1 + h*eigenvalue) <= 1 :
			return h
		cnt += 1
		if cnt > 1000000:
			raise Exception('Stepsize not small enough for the excplicitLinearTest.')
		elif h <= 0:
			raise Exception('Stepsize cannot be smaller than zero.')

def implicitLinearTest(tol, eigenvalue):
	h = 0
	cnt = 0
	while True:
		h += tol
		if abs(1 - h*eigenvalue) >= 1 :
			return h
		cnt += 1
		if cnt > 1000000:
			raise Exception('Stepsize not small enough for the implicitLinearTest.')






## -- Task #1 -- ##

#dy/dt = lambda*y(t), with y(0) = 1, lambda = -2, -1 and 4. Find the maximum possible stepvalue (h) whilst still staying within stable region.

tol = 1e-5
eigenvalues = [-2, -1, 4]
for i in range(len(eigenvalues)):
	print 'Lambda = {}, Excplicit max h = {}, Implicit max h = {}'.format(eigenvalues[i], explicitLinearTest(tol, eigenvalues[i]), implicitLinearTest(tol, eigenvalues[i]))





