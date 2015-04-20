from scipy import *
from pylab import *


energy = array([27.93, 46.98, 31.95, 31.68, 21.00])
days = []
for i in range(len(energy)):
	days.append(i+1)
print energy
print days


## --Task 1-- ""

def polyfitter(x, y, name):
	n = len(x)-1
	p = polyfit(x, y, n)
	plotter(x1, polyval(p,x1), name)
	return p

def plotter(x, y, name = '', xname = '', yname = ''):
	plot(x, y)
	grid(True)
	if len(name) > 0:
		title(name)
	if len(xname) > 0:
		xlabel(xname)
	if len(yname) > 0:
		ylabel(yname)
	show()
	return

def lagrange(x, i, xm):
	"""
	Evaluates the i-th Lagrange polynomial at x
	based on grid data xm
	"""
	n = len(xm) - 1
	y = 1
	for j in range(n+1):
		if i != j:
			y *= (x - xm[j])/(xm[i]-xm[j])
	return y

def interpolationLagrange(x, xm, ym):
	n = len(xm) - 1
	p = array([lagrange(x, i, xm) for i in range(n+1)])
	y = dot(ym, p)
	print y
	plotter(x, polyval(y, x), 'Lagrange')
	return y





x1 = linspace(0.9,5.1,1000)
#print polyfitter(days, energy, 'Polyfitter')
print interpolationLagrange(x1, days, energy)
