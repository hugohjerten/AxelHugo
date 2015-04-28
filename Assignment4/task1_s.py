# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 10:09:23 2015
@author: sebastianelm
"""
from __future__ import division 
from scipy import *
from matplotlib.pylab import *
from scipy.linalg import solve
from pylab import *
from numpy import *

def cubeSpline(xint,yint):
    h = abs(xint[1]-xint[0])
    di = yint 
    m = len(xint) # antal punkter
    coeff = zeros((m-1,4))
    allSigmas = zeros((m,m))
    allSigmas[0][0]=1    
    allSigmas[m-1][m-1]=1
    Y = zeros((m,1))
    col=0
    for row in range(1,m-1):
        allSigmas[row][col]=1
        col+=1
        allSigmas[row][col]=4
        col+=1
        allSigmas[row][col]
        col-=1
        Y[row][0] = (6/(h**2))*(yint[row+1]-2*yint[row]+yint[row-1])
    sigmas = solve(allSigmas,Y)
    bi = sigmas/2
    ai=[]
    ci=[]
    for i in range(m-1):
        ai.append((sigmas[i+1]-sigmas[i])/(6*h))
        ci.append((yint[i+1]-yint[i])/h - (2*sigmas[i]+sigmas[i+1])*h/6) 
    for i in range(m-1):
        coeff[i][0]=ai[i]
        coeff[i][1]=bi[i]
        coeff[i][2]=ci[i]
        coeff[i][3]=di[i]
    return(coeff)
   
def cubeSplineVal(coeff,xint,xval):
    for i in range(1,len(xint)):
        if xval<=xint[i]:
            xdiff = xval-xint[i-1] # ex. ai*(xi - x)^3, räknar diffen
            return coeff[i-1][0]*xdiff**3 + coeff[i-1][1]*xdiff**2 + coeff[i-1][2]*xdiff + coeff[i-1][3]
    xdiff = xval-xint[len(xint)-1] 
    return coeff[len(xint)-1][0]*xdiff**3 + coeff[len(xint)-1][1]*xdiff**2 + coeff[len(xint)-1][2]*xdiff + coeff[len(xint)-1][3]

def plotSpline(coeff, xint, yint, fignr = 1,name = 'UNNAMDED'):
    a = xint[0]
    b = xint[-1]
    xplot = linspace(a, b, 1000)
    yplot = []
    for i in range(len(xplot)):	
        yplot.append(cubeSplineVal(coeff, xint, xplot[i]))
    figure(fignr)
    title(name)
    plot(xplot, yplot)
    grid(True)
    plot(xint, yint, 'ro')
    grid(True)
    show()
    return

x = linspace(0,9,10)
y = [4, 8, 3, 6, 7, 9, 2, 5, 6, 11]
S = cubeSpline(x, y)
for i in range(10):
	print (cubeSplineVal(S, x, i))
plotSpline(S, x, y,1, 'Task #1b, random points')


# -- Task 2a -- #

x = [0, 1, 2, 3, 4, 5, 6]
y = [1., 3., -2., 0., 1., 0., 1.]
plotSpline(cubeSpline(x, y), x, y,2, 'Task #2a')

xi= [0,0,0,0, 1, 2, 3, 4, 5, 6,6,6,6] # The x-vector with 3 extra nodes on each sides, value of the endponts.
di= [1.,2.3,5.15,-4.35,0.8,1.4,-0.73,0.6,1.] # De Boor "control points" to the B-spline function two extra on each side.
dx= 0.02 # the interval between calculated values for the B-spline
BsplYval = Bsplbasis(xi,di,dx) # A vector with the B-splines value, points with distans dx between.

def Bsplbasis(xi,di,dx):
    """
    This code evaluates a cubic spline given by its knots and de Boor points
    at equidistant points
    On input:
    =========
    xi.... list or array of length m+7 where 
           xi is composed out of the m+1 knots x and 6 additional 
           points at both boundaries. Preferably one takes the first and last
           knot and repeats it so that it gets multiplicity 4 each.
    di.... list or an array of m+3 de Boor points
    dx.... float which is used to define the evaluation points, which form an
           equidistant grid starting at xi[3]=x[0] and having an interval length dx
    On return:
    ==========
    q..... list with the values of the spline evaluated at
           xi[3]+i*dt, i=0,... as long as xi[3]+i*dt <= xi[-4]  
    """
    
    eps = 1.e-14
    m = len(xi) #number knots
    i = 4      #index of first knot
    q=[]
    for u in arange(xi[3],xi[m-3]+dx,dx):
        
        # check if u value has moved to the next knot interval
        # include small tolerance on knots to avoid round-off error in comparisons.
        
        while (u>(xi[i]+eps)):
            i+=1
        # Now evaluate the spline at u using the deBoor algorithm.
        # Start with the relevant control points.
        # w used here to simplify indices.
        im4 = i-4
        qq = zeros(len(di))
        for j in range(1,5,1):
            qq[j-1]=di[im4+j-1]
        for j in range(1,4,1):
            for k in range(1,4-j+1,1):
                qq[k-1] = ((xi[im4 + k + 4-1] - u)/(xi[im4 + k + 4-1] - xi[im4 + k + j-1])) * qq[k-1] + \
                            ((u - xi[im4 + k + j-1])/(xi[im4 + k + 4-1] - xi[im4 + k + j-1])) * qq[k+1-1]
        #Create vector of points on the B-spline curve.
        q.append(qq[0])
    return q


figure(2)
BsplXval=linspace(min(xi),max(xi),len(BsplYval))
plot(BsplXval,BsplYval,'y')
grid(True)
show()

    

