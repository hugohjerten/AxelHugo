from scipy import *
from pylab import *


energy = array([27.93, 46.98, 31.95, 31.68, 21.00])
days = []
for i in range(len(energy)):
	days.append(i+1)
print(energy)
print(days)

## --Task 1-- ""

def polyfitter(x, y, name):
	n = len(x)-1
	p = polyfit(x, y, n)
	return p

def plotter(x, y, name = '', xname = '', yname = ''):
	subplot(211)
	plot(x,y)
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
	#print y
	return y

def interpolationVander(x, y):
	V = vander(x,len(x))
	a = solve(V,y)
	return a

x1 = linspace(0.9,5.1,1000)
polyfitterReturn =  polyfitter(days, energy, 'Polyfitter')
interpolationLagrangeReturn = interpolationLagrange(x1, days, energy)
interpolationVanderReturn = interpolationVander(days, energy)

polyFitEnergy0 = polyval(polyfitterReturn, 8.0)
lagrangeEnergy0 = interpolationLagrange(8.0, days, energy)
vanderEnergy0 = polyval(interpolationVanderReturn, 8.0)

print 'Polyfit value: {}'.format(polyFitEnergy0)
print 'Lagrange value: {}'.format(lagrangeEnergy0)
print 'Vandermonde value: {}'.format(vanderEnergy0)

figure(1)
subplot(311)
plot(x1, polyval(polyfitterReturn, x1))
grid(True)
title('Polyfitter')
subplot(312)
plot(x1, interpolationLagrangeReturn)
grid(True)
title('Lagrange')
subplot(313)
plot(x1, polyval(interpolationVanderReturn, x1))
grid(True)
title('Vandermonde')
show()
for i in range(len(energy)):
	dx = 0.01
	energy[i] += dx
	polyfitterReturn =  polyfitter(days, energy, 'Polyfitter')
	interpolationLagrangeReturn = interpolationLagrange(x1, days, energy)
	interpolationVanderReturn = interpolationVander(days, energy)
	
	polyFitEnergy = polyval(polyfitterReturn, 8.0)
	lagrangeEnergy = interpolationLagrange(8.0, days, energy)
	vanderEnergy = polyval(interpolationVanderReturn, 8.0)
	dYPoly = abs(polyFitEnergy0 - polyFitEnergy)
	dYLagrange = abs(lagrangeEnergy0 - lagrangeEnergy)
	dYVander = abs(vanderEnergy0 - vanderEnergy)
	relPol = (dYPoly/max(energy)) / (dx/max(days))
	relLag = (dYLagrange/max(energy)) / (dx/max(days))
	relVan = (dYVander/max(energy)) / (dx/max(days))
	print 'Change the data point #{}, the relative sensitivty of polyfitter: {}'.format(i+1, relPol)
	print 'Change the data point #{}, the relative sensitivty of lagrange: {}'.format(i+1, relLag)	
	print 'Change the data point #{}, the relative sensitivty of vander: {}'.format(i+1, relVan)



## Task 4 ##

def function(x):
	y = 1.0 / (1.0 + 25*x**2)
	return y


n = 3
x2 = linspace(-1,1,n)



