import numpy as np
import matplotlib.pyplot as plt

Data_pts = np.array([[0,0,0], [1,1,0], [2,1,0],[3,2,0]])

print(Data_pts)

n = Data_pts.shape[0]

p = 3

w = np.ones(n)

def coxDeBoor(zeta, i, d, knot_vect):
    
    if d==0:
        if knot_vect[i] <= zeta and zeta <= knot_vect[i+1]:
            return 1
        return 0
    
    D1 = knot_vect[i+d] - knot_vect[i]
  
    D2 = knot_vect[i+d+1] - knot_vect[i+1]
    
    E1 = 0
    E2 = 0
    
    if D1 > 0:
        E1 = (zeta - knot_vect[i])/D1 * coxDeBoor(zeta, i, d-1, knot_vect)
    
    if D2 > 0:
        E2 = (knot_vect[i+d+1] - zeta)/D2 * coxDeBoor(zeta, i+1, d-1,knot_vect)
        
    return E1+E2


def NURBS_interpolation(Data_pts, p):

	'''
	the funnction takes in the input of data points and degree and outputs the control points
	and knot vector
	Data points should be in 3D
	'''
	m = Data_pts.shape[0] - 1
	
	U = np.zeros((m+p+2))
	
	t = np.zeros((m+1))
	
	dist = np.zeros((m))
	total_dist = 0
	
	for i in range(m):
		dist[i] = ((Data_pts[i+1,0]-Data_pts[i,0])**2
                           + (Data_pts[i+1,1]-Data_pts[i,1])**2
                           + (Data_pts[i+1,2]-Data_pts[i,2])**2)**0.5
		total_dist = total_dist + dist[i]	
	
	t[-1] = 1.
	
	for i in range(1,t.shape[0]-1):
		t[i] = t[i-1] + dist[i-1]/total_dist

	#Check there is a big here

	for j in range(1, m - p + 1):
		for i in range(j, j+p):
			U[j+p] = U[j+p] + (1/p)*t[i]

	for j in range(p+1):
		U[-j-1] = 1.
	
	N = np.zeros((m+1,m+1))
	for i in range(m+1):
		for j in range(m+1):
			N[i,j] = coxDeBoor(t[i], j, p, U)
	
	print(t)
	print(U)
	P = np.matmul(np.linalg.inv(N),Data_pts)

	return U,P	 
	
Xi, ctrlpts = NURBS_interpolation(Data_pts,3)

knot_vect = Xi

zeta = np.linspace(0,1,300)
N = np.zeros((300,n))
Wtd_N = np.zeros((300,n))
Curve = np.zeros((300,3))



for z in range(len(zeta)):
    for i in range(n):
        N[z,i] = coxDeBoor(zeta[z], i, p, Xi)
        Wtd_N[z,i] = N[z,i]*w[i]

		
	
	
Wtd_sum = np.sum(Wtd_N,1)
for z in range(len(zeta)):        
    for i in range(n):
        if Wtd_sum[z] != 0:
            Curve[z] = Curve[z] + Wtd_N[z,i]/Wtd_sum[z] * ctrlpts[i]


plt.plot(Curve[:,0],Curve[:,1],label='NURBS curve')
plt.scatter(Data_pts[:,0], Data_pts[:,1], label='Interpolating points', color = 'r')
plt.plot(ctrlpts[:,0],ctrlpts[:,1], label='Control points interpolation', marker='o')
plt.axis('equal')
plt.legend()
plt.show()



