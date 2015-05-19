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

def explicitEuler(Un, Vn, En, Fn, h):
	n = linspace(0, 100, 100/h)
	x = []
	y = []
	G = 1.0
	m = 3.0
	x.append(Un)
	y.append(En)
	for i in range(len(n)):
		r = sqrt(Un**2+En**2)
		r3 = r**3
		Unp1 = Un + h*Vn
		Vnp1 = Vn - h*G*m*Un/r3
		Enp1 = En + h*Fn
		Fnp1 = Fn - h*G*m*En/r3
		Un = Unp1
		Vn = Vnp1
		En = Enp1
		Fn = Fnp1
		x.append(Un)
		y.append(En)
	return x,y


def implicitEuler(Un, Vn, En, Fn, h, tol):
	n = linspace(0, 100, 100/h)
	x = []
	y = []
	W = []
	G = 1.0
	m = 3.0
	x.append(Un)
	y.append(En)
	W.append([Un, Vn, En, Fn])
	for i in range(len(n)):
		print 'Loop i = {}'.format(i)
		Upred = Un
		Vpred = Vn
		Epred = En
		Fpred = Fn
		cnt = 0
		while True:
			print 'while loop n = {}'.format(cnt)
			print
			r = sqrt(Upred**2 + Epred**2)
			r3 = r**3
			Unp1 = Un + h*Vpred
			Vnp1 = Vn - h*G*m*Upred/r3
			Enp1 = En + h*Fpred
			Fnp1 = Fn - h*G*m*Epred/r3
			Upred = Unp1
			Vpred = Vnp1
			Epred = Enp1
			Fpred = Fnp1
			if abs(Upred - Un) < tol and abs(Epred - En) < tol:
				x.append(Upred)
				y.append(Epred)
				Un = Upred
				Vn = Vpred
				En = Epred
				Fn = Fpred
				break
			cnt += 1
	return x,y

## -- Task #1 -- ##

#dy/dt = lambda*y(t), with y(0) = 1, lambda = -2, -1 and 4. Find the maximum possible stepvalue (h) whilst still staying within stable region.

tol = 1e-5
eigenvalues = [-2, -1, 4]
for i in range(len(eigenvalues)):
	print 'Lambda = {}, Excplicit max h = {}, Implicit max h = {}'.format(eigenvalues[i], explicitLinearTest(tol, eigenvalues[i]), implicitLinearTest(tol, eigenvalues[i]))

## -- Task #2 -- ##
# (x0, y0) = (0,2)
# (x0', y0') = (1,0)

explicitReturn1 = explicitEuler(0.0, 1.0, 2.0, 0.0, 0.1)
explicitReturn2 = explicitEuler(0.0, 1.0, 2.0, 0.0, 0.01)
explicitReturn3 = explicitEuler(0.0, 1.0, 2.0, 0.0, 0.001)
figure(1)
plot(explicitReturn1[0], explicitReturn1[1], label = 'h = 0.1')
plot(explicitReturn2[0], explicitReturn2[1], label = 'h = 0.01')
plot(explicitReturn3[0], explicitReturn3[1], label = 'h = 0.001')
legend(loc = 'upper right')
title('Explicit Euler Method')
grid(True)
show()

tol = 1e-6

implicitReturn1 = implicitEuler(0.0, 1.0, 2.0, 0.0, 0.1, tol)
figure(2)
plot(implicitReturn1[0], implicitReturn[1], label = 'h = 0.1')
legend(loc = 'upper right')
title('Implicit Euler Method')
grid(True)
show()

















show()



















