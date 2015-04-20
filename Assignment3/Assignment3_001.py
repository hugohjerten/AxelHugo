from scipy import *
from pylab import *


energy = array([27.93, 46.98, 31.95, 31.68, 21.00])
days = []
for i in range(len(energy)):
	days.append(i+1)
print energy
print days


## --Task 1-- ""

def polyfitter(x, y, name, xname = '', yname = ''):
	n = len(x)-1
	p = polyfit(x, y, n)
	x1 = linspace(1,5,1000)
	plotter(x1, polyval(p,x1), name, xname, yname)
	return p

def plotter(x, y, name, xname, yname):
	plot(x, y)
	grid(True)
	title(name)
	if len(xname) > 0:
		xlabel(xname)
	if len(yname) > 0:
		ylabel(yname)
	show()
	return

def 



print polyfitter(days,energy, 'Polyfitter')

