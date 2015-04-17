from scipy import *
from pylab import *
from scipy.interpolate import interp1d
from scipy.linalg import solve
from numpy import *

## --Task 1-- ""

def polyfitter(x, y, name):
	n = len(x)-1
	p = polyfit(x, y, n)
	return p

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

def interpolationVander(x, y):
	V = vander(x,len(x))
	a = solve(V,y)
	return a

def interpolationLagrange(x, xm, ym):
	n = len(xm) - 1
	p = array([lagrange(x, i, xm) for i in range(n+1)])
	y = dot(ym, p)
	#print y
	return y

energy = array([27.93, 46.98, 31.95, 31.68, 21.00])
days = []
for i in range(len(energy)):
	days.append(i+1)

x1 = linspace(0.9,5.1,1000)
polyfitterReturn =  polyfitter(days, energy, 'Polyfitter')
interpolationLagrangeReturn = interpolationLagrange(x1, days, energy)
interpolationVanderReturn = interpolationVander(days, energy)

polyFitEnergy0 = polyval(polyfitterReturn, 9.0)
lagrangeEnergy0 = interpolationLagrange(9.0, days, energy)
vanderEnergy0 = polyval(interpolationVanderReturn, 9.0)

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
	
	polyFitEnergy = polyval(polyfitterReturn, 9.0)
	lagrangeEnergy = interpolationLagrange(9.0, days, energy)
	vanderEnergy = polyval(interpolationVanderReturn, 9.0)
	dYPoly = abs(polyFitEnergy0 - polyFitEnergy)
	dYLagrange = abs(lagrangeEnergy0 - lagrangeEnergy)
	dYVander = abs(vanderEnergy0 - vanderEnergy)
	relPol = (dYPoly/max(energy)) / (dx/max(days))
	relLag = (dYLagrange/max(energy)) / (dx/max(days))
	relVan = (dYVander/max(energy)) / (dx/max(days))
	print 'Change the data point #{}, the relative sensitivty of polyfitter: {}'.format(i+1, relPol)
	print 'Change the data point #{}, the relative sensitivty of lagrange: {}'.format(i+1, relLag)	
	print 'Change the data point #{}, the relative sensitivty of vander: {}'.format(i+1, relVan)

## -- Task 2 -- ##

def w(x, xn):

	w = 1
	for j in range(len(xn)):
		w *= (x-xn[j])
	return w

#xn = [-0.25, -0.5, 0, 0.25, 0.5]
#figure(2)
xn = linspace(-1,1, 8)
x = linspace(-1, 1, 10000)
f = -x**3 + x 
f2 = interp1d(x, f) 
plot(x, f, x, w(x,xn), xn, f2(xn))
title('Interpolation Error')
plot(x, w(x, xn))
show()

## -- Task 3 -- ##

def chebyshev(n):
	
	xc = [] 
	for i in range(n):
		xc.append(cos(((2*(i+1) - 1)*pi/(2*n))))
	xc = xc[::-1]
	return xc

xc = chebyshev(8)
figure(3)
subplot(411)
plot(x, f)
#grid(True)
title('Function')
subplot(412)
plot(xc ,f2(xc))
#grid(True)
title('chebyshev')
subplot(413)
plot(xn, f2(xn))
#grid(True)
title('Regular interpolation')
subplot(414)
plot(x, f, xc, f2(xc), 'x', xn, f2(xn), 'o')
title('all in all')

show()
print xc
#plot(x, f, 


## -- Task 4 -- ##

def function(x):
	y = 1.0 / (1.0 + 25*x**2)
	return y
n = 3
figure(5)
for i in range(3):
	x2 = linspace(-1,1,n)
	xc = chebyshev(n)
	y2 = []
	yc = []
	for j in range(len(x2)):
		y2.append(function(x2[j]))
		yc.append(function(xc[j]))
	subplot(311+i)
	plot(x2, y2, xc, yc)
	n += 6
show()
