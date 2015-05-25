from scipy import *
from pylab import *
from numpy import *

## -- Task 3 -- ##

def pendulum(duration, stepsize):
	h = stepsize
	time = linspace(0.0, duration, duration/h)
	l = 1
	g = 9.81
	Pn = pi/4.0
	Vn = 0.0
	P = [Pn] #position
	V = [Vn] #velocity
	Ptmp = Pn
	Vtmp = Vn
	for i in range(100):
		Pnp1 = Pn + h*Vtmp
		Vnp1 = Vn - h*g*math.sin(Ptmp)
		Ptmp = Pnp1
		Vtmp = Vnp1
	P.append(Pnp1)
	V.append(Vnp1)
	for mainIt in range(2, len(time)):
		Pn = P[mainIt - 1]
		Pnm1 = P[mainIt - 2]
		Vn = V[mainIt - 1]
		Vnm1 = V[mainIt - 2]
		Ptmp = Pn
		Vtmp = Vn
		for i in range(100):
			Pnp1 = (4.0/3.0)*Pn - (1.0/3.0)*Pnm1 + h*(2.0/3.0)*Vtmp
			Vnp1 = (4.0/3.0)*Vn - (1.0/3.0)*Vnm1 - h*(2.0/3.0)*(g/l)*math.sin(Ptmp)
			Ptmp = Pnp1
			Vtmp = Vnp1
		P.append(Pnp1)
		V.append(Vnp1)
	return time, P, V



"""
	a'' = -g/l*sin(a)   t >= 0
	a(0) = pi/4
	a'(0) = 0
	
"""
duration = 25.0
stepsize = 0.1

time, alpha, alphaSpeed = pendulum(duration, stepsize)
print len(time)
print len(alpha)
print len(alphaSpeed)
figure(1)
subplot(121)
plot(time, alpha, label = 'Angle')
plot(time, alphaSpeed, label = 'Angular Velocity')
legend(loc = 'upper right')
title('Pendulum (h = {})'.format(stepsize))
xlabel('Time [s]')
grid(True)
subplot(122)
plot(alpha, alphaSpeed)
title('Phase Plot (h = {})'.format(stepsize))
xlabel('Angle [rad]')
ylabel('Angular Velocity [rad/s]')
grid(True)
show()














