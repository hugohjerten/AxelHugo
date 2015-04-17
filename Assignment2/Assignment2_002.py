from pylab import*
from scipy import*
from scipy.linalg import solve



	##----------------Task 1-------------------##

def newton(f, x, eps = 1.e-8):
	xf, yf = f(x)
	n = 0
	while abs(xf) and abs(yf) > eps:
		deltax = solve(jacobian(f, x), -f(x))
		x += deltax
		xf, yf = f(x)
		n +=1
	return x


def jacobian (f , x , eps = 1.e-8 , fx = []):
	
	#Computes the Jacobian of f 
	#eps is the increment and fx a reference value f(x)
	
	if fx==[]: fx=f(x)
	jac = empty([len(fx) ,len(x)])
	for i in range(len(x)):
		x[i]+=eps
		jac[:,i]=(f(x)-fx)/eps
		x[i]-=eps
	return jac

def robotArm(x):
	alfa, beta = x
	xf = 5 * cos(alfa) + 6 * cos(alfa + beta) - 10 
	yf = 5 * sin(alfa) + 6 * sin(alfa + beta) - 4
	
	return array([xf, yf])

def LUFact(m):
	if not m.shape[0] == m.shape[1]:
		print "Matrix is not square"
		return	
	n = m.shape[0]
	L = zeros((n,n))
	U = zeros((n,n))
	for i in range(n):
		for j in range(n):
			if i==j:
				L[i,j] = 1

	values = []
        for col in range(len(m[0])+1):
                for row in range(col+1, len(m)):
                        r = []
                        for rowValue in m[col]: 
                                r.append((rowValue * (-(m[row][col] / m[col][col]))))
                #       r = [(rowValue * (-(m[row][col] / m[col][col]))) for rowValue in m[col]]
                        values.append((m[row][col] / m[col][col])) 
                        m[row] = [sum(pair) for pair in zip(m[row], r)]
	
	for col in range(0, n-1):
		for row in range(1, n):	
				L[row][col] = values[col]
				
	return m, L

def myGauss(m):
    # eliminate columns
	values = []
	for col in range(len(m[0])+1):
        	for row in range(col+1, len(m)):
			r = []
			for rowValue in m[col]:
				r.append((rowValue * (-(m[row][col] / m[col][col]))))
            	#	r = [(rowValue * (-(m[row][col] / m[col][col]))) for rowValue in m[col]]
			values.append((m[row][col] / m[col][col]))			
			m[row] = [sum(pair) for pair in zip(m[row], r)]		
	return m, values


#### -- Task 2 -- #####

def membrane(alpha = 1):
	#A = [[0]*30 for j in range(30)]
	A = zeros((30,30))
	A[0][0] = -(2 + alpha)
	A[0][1] = 1
	for i in range(1, 29):
		A[1][i - 1] = 1
		A[i][i] = -(2 + alpha)
		A[i][i + 1] = 1
	A[29][28] = 1
	A[29][29] = -(2 + alpha)
	#B = [[0]*30]
	B = zeros(30)
	B[5] = 2
	print B
	return A, B

def LUsolver(A, b):
	U, L = LUFact(A)
	n = len(A)	
	y = []
	for row in range(n):
		if row == 0:
			y.append(b[row])
		else:
			tmpSum = 0
			for k in range(row):
				tmpSum += L[row][k]*y[k]
			y.append(b[row] - tmpSum)
	print y
	x = zeros(n)
	for i in range(n):
		row = n - i - 1
		if row == n - 1:
			x[row] = y[row] / U[row][row]
		else:
			tmpSum = 0
			q = n - 1 - row
			for k in range(q):
				r = n - 1 - k
				tmpSum += U[row][r]*x[r]
			x[row] = (y[row] - tmpSum) / U[row][row]
	print x
	return

def test001():
	A = array([[1, 2, 3], [3, 4, 1], [1, 0, 1]])
	B = array([0, 6, 1])
	return A, B

def iterativeJacobi(A, b):
	Q = diag(A)
	return



x = array([0.7, 0.7])
m = array([[1,2,3], [3,4,1], [1,0,1]])
#print robotArm(x)
#print jacobian(robotArm, x)
#print LUFact(jacobian(robotArm, x))
#print LUFact(m)

print newton(robotArm, x)
#membrane()
membraneReturn = membrane()
LUsolver(membraneReturn[0], membraneReturn[1])
testReturn = test001()
LUsolver(testReturn[0], testReturn[1])


#print LUFact(jacobian(robotArm, x))
