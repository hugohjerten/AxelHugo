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
    return

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


x = linspace(0,9,10)
y = [4, 8, 3, 6, 7, 9, 2, 5, 6, 11]
S = cubeSpline(x, y)
for i in range(10):
	print (cubeSplineVal(S, x, i))
plotSpline(S, x, y,1, 'Task #1b, random points')
show()

# -- Task 2a -- #

x = [0, 1, 2, 3, 4, 5, 6]
y = [1., 3., -2., 0., 1., 0., 1.]
plotSpline(cubeSpline(x, y), x, y,2, 'Task #2a and b')

# -- Task 2b -- #

xi= [0,0,0,0, 1, 2, 3, 4, 5, 6,6,6,6] # The x-vector with 3 extra nodes on each sides, value of the endponts.
di= [1.,2.3,5.15,-4.35,0.8,1.4,-0.73,0.6,1.] # De Boor "control points" to the B-spline function two extra on each side.
dx= 0.02 # the interval between calculated values for the B-spline
BsplYval = Bsplbasis(xi,di,dx) # A vector with the B-splines value, points with distans dx between.

figure(2)
BsplXval=linspace(min(xi),max(xi),len(BsplYval))
plot(BsplXval,BsplYval,'y')
grid(True)
show()

def s1002(S):
     """
     this function describes the wheel profile s1002
     according to the standard. 
     S  independent variable in mm bewteen -69 and 60.
     wheel   wheel profile value
     (courtesy to Dr.H.Netter, DLR Oberpfaffenhofen)
                                               I
                                               I
                       IIIIIIIIIIIIIIIIIIIIIIII
                     II  D  C       B       A
                    I 
         I         I   E
          I       I
       H   I     I   F
            IIIII
 
              G
 
 
      FUNCTIONS:
      ---------- 
      Section A:   F(S) =   AA - BA * S                 
      Section B:   F(S) =   AB - BB * S    + CB * S^2 - DB * S^3
                               + EB * S^4 - FB * S^5 + GB * S^6
 

                              - HB * S^7 + IB * S^8
      Section C:   F(S) = - AC - BC * S    - CC * S^2 - DC * S^3
                               - EC * S^4 - FC * S^5 - GC * S^6
                               - HC * S^7
      Section D:   F(S) = + AD - SQRT( BD^2 - ( S + CD )^2 )
      Section E:   F(S) = - AE - BE * S
      Section F:   F(S) =   AF + SQRT( BF^2 - ( S + CF )^2 )
      Section G:   F(S) =   AG + SQRT( BG^2 - ( S + CG )^2 )
      Section H:   F(S) =   AH + SQRT( BH^2 - ( S + CH )^2 )
     """

     #Polynom coefficients:
     #-------------------- 
     #Section A:    
     AA =  1.364323640
     BA =  0.066666667

     #Section B:     
     AB =  0.000000000;
     BB =  3.358537058e-02;
     CB =  1.565681624e-03;
     DB =  2.810427944e-05;
     EB =  5.844240864e-08;
     FB =  1.562379023e-08;
     GB =  5.309217349e-15;
     HB =  5.957839843e-12;
     IB =  2.646656573e-13;
     #Section C:     
     AC =  4.320221063e+03;
     BC =  1.038384026e+03;
     CC =  1.065501873e+02;
     DC =  6.051367875;
     EC =  2.054332446e-01;
     FC =  4.169739389e-03;
     GC =  4.687195829e-05;
     HC =  2.252755540e-07;
     #Section D:     
     AD = 16.446;
     BD = 13.;
     CD = 26.210665;
     #Section E: 
     AE = 93.576667419;
     BE =  2.747477419;
     #Section F:     
     AF =  8.834924130;
     BF = 20.;
     CF = 58.558326413;
     #Section G:   
     AG = 16.;
     BG = 12.;
     CG = 55.;
     #Section H:   
     AH =  9.519259302;
     BH = 20.5;
     CH = 49.5;
     circarc = lambda S, A,B,C: A + sqrt(B**2 - ( S + C )**2)
     #
     #     Bounds
     #     ------- 
     #                        from                    to
     #     Section A:      Y = + 60               Y = + 32.15796
     #     Section B:      Y = + 32.15796         Y = - 26.
     #     Section C:      Y = - 26.              Y = - 35.
     #     Section D:      Y = - 35.              Y = - 38.426669071
     #     Section E:      Y = - 38.426669071     Y = - 39.764473993
     #     Section F:      Y = - 39.764473993     Y = - 49.662510381
     #     Section G:      Y = - 49.662510381     Y = - 62.764705882
     #     Section H:      Y = - 62.764705882     Y = - 70.

     YS=[-70.0,-62.764705882,-49.662510381,
             -39.764473993, -38.426669071, -35.0,
             -26.0,32.15796,60.0];

     if (S < YS[1]):                        #  Section H (Circular arc) 
        wheel = circarc(S,AH,BH,CH)
     elif (S  < YS[2]):                     #  Section G (Circular arc) 
        wheel = circarc(S,AG,BG,CG)
     elif (S < YS[3]):                      #  Section F (Circular arc) 
        wheel = circarc(S,AF,BF,CF)
     elif (S < YS[4]):                      #        Section E (LINEAR)
        wheel   = -BE*S-AE;
     elif (S < YS[5]):                      #  Section D (Circular arc) 
        wheel = -circarc(S,-AD,BD,CD)   
     elif (S < YS[6]):                      #                 Section C 
        wheel = - polyval([HC,GC,FC,EC,DC,CC,BC,AC],S)  
     elif (S < YS[7]):                      #                 Section B
        wheel = polyval([IB,-HB,GB,-FB,EB,-DB,CB,-BB,AB],S)
     else:                                  #        Section A (LINEAR)            
        wheel   =   -BA*S + AA;
     return wheel



  

