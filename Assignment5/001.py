
from scipy import *
from pylab import *
from scipy.linalg import solve
from numpy import *

def simpsons(f, a, b, n):
	h = abs(b-a)/n
	#print 'Stepsize = {}'.format(h)
	tmp = 0
	limit = int(n)
	for j in range(1,limit):
		xj = a + j*h
		#print xj
		tmp += 2.0*f(xj) + 4.0*f(xj + h/2.0)
		#print tmp
	integral = (h/6.0)*(f(a) + 4.0*f(a + h/2.0) + f(b) + tmp)
	
	
	return integral 

def threePointGauss(f, a, b, n):
	h = abs(b-a)/n
	tmp = 0
	limit = int(n)
	for j in range(0, limit):
		xj = a + j*h
		tmp1 = f(xj + h*(5.0 - math.sqrt(15.0))/10.0)
		tmp2 = f(xj + 5.0*h/10.0)
		tmp3 = f(xj + h*(5.0 + math.sqrt(15.0))/10.0)
		tmp += 5.0*tmp1 + 8.0*tmp2 + 5.0*tmp3
	integral = h*tmp/18.0
	return integral

def simpsonsError(a, b, n):
	h = abs(b-a)/n
	error = (abs(b-a)*96*h**4)/2880
	return error


def estimatedStep(a, b, error):
	n = 0
	est = 1
	while abs(est) > error:
		est = simpsonsError(0.0, 1.0, n+1)
		n += 1
	return n

def xsquare(x):
	return x**2

def xcube(x):
	return x**3

def xfourth(x):
	return x**4 + x

def taskThree(x):
	return 4/(1+x**2)

## -- Task 2 --##
a = 1.0
b = 3.0
n = 1.0

print 'Simpons integral of x^2 between {} and {} = {}'.format(a, b, simpsons(xsquare, a, b, n))
print '3PointGauss integral of x^2 between {} and {} = {}'.format(a, b, threePointGauss(xsquare, a, b, n))
print 'Simpons integral of x^3 between {} and {} = {}'.format(a, b, simpsons(xcube, a, b, n))
print '3PointGauss integral of x^3 between {} and {} = {}'.format(a, b, threePointGauss(xcube, a, b, n))
print 'Simpons integral of x^4 between {} and {} = {}'.format(a, b, simpsons(xfourth, a, b, n))
print '3PointGauss integral of x^4 between {} and {} = {}'.format(a, b, threePointGauss(xfourth, a, b, n))

simpsonsError = []
gaussError = []
nPlot = [2, 4, 8, 16, 32, 64]

for i in range(len(nPlot)):
	simpsonsError.append(abs(simpsons(taskThree, 0.0, 1.0, nPlot[i])-pi))
	gaussError.append(abs(threePointGauss(taskThree, 0.0, 1.0, nPlot[i])-pi))

figure(1)
grid(True)
loglog(nPlot, simpsonsError, label = "Simpsons")
loglog(nPlot, gaussError, label = "Gauss")
legend(loc ="upper right")
title('Order of convergence')
show()

## -- Task 3 --##

tol = 1e-8
n = 0

print 'Estimated number of steps for function task 3a = {}'.format(estimatedStep(0.0, 1.0, tol))

stepsSimpsons = 0
stepsGauss = 0
solSimpsons = 0
solGauss = 0

while abs(solSimpsons - pi) > tol:
	stepsSimpsons += 1
	solSimpsons = simpsons(taskThree, 0.0, 1.0, stepsSimpsons)

while abs(solGauss - pi) > tol:
	stepsGauss += 1
	solGauss = threePointGauss(taskThree, 0.0, 1.0, stepsGauss)

print 'Number of steps required for tolerance {} with Simpsons is n = {}'.format(tol, stepsSimpsons)
print 'Number of steps required for tolerance {} with Gauss is n = {}'.format(tol, stepsGauss)

print 'Task 3 simpsons integral of function integral = {}'.format(simpsons(taskThree, 0.0, 1.0, stepsSimpsons))
print 'Task 3 Gauss integral of function integral = {}'.format(threePointGauss(taskThree, 0.0, 1.0, stepsGauss))










