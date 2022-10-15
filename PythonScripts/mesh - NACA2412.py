    
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

c = 0

with open('NACA-2412.csv') as csvfile:
	readCSV = csv.reader(csvfile, delimiter=',')
	for row in readCSV:
		c = c + 1


Data_pts = np.zeros((c,3))

c = 0 

with open('NACA-2412.csv') as csvfile:
	readCSV = csv.reader(csvfile, delimiter=',')
	for row in readCSV:
		Data_pts[c,0] = float(row[0])
		Data_pts[c,1] = float(row[1])
		c = c + 1

n1 = Data_pts.shape[0]
p1 = 3

p2 = 2 #order of the polynomial for second basis function

n2 = 4 #numer of basis vector for second basis function

w = np.ones((n1,n2))

knot_vect2 = np.array([0,0,0,0.5,1,1,1])

knot_vect1 = np.array([0,0,0,0.25,0.25,0.5,0.5,0.75,0.75,1,1,1])

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
        E2 = (knot_vect[i+d+1] - zeta)/D2 * coxDeBoor(zeta, i+1, d-1, knot_vect)
        
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
	
	P = np.matmul(np.linalg.inv(N),Data_pts)

	return U,P	 
	
knot_vect1, C = NURBS_interpolation(Data_pts,3)

ctrlpts = np.zeros((n1,n2,3))

for i in range(n1):    
    for j in range(n2):
        z = -3 + 2*j
        ctrlpts[i,j,0] = C[i,0]
        ctrlpts[i,j,2] = C[i,1]
        ctrlpts[i,j,1] = z

print(ctrlpts)


x = np.linspace(0,1,200)


def Surface(xi, neta):
    N1 = np.zeros((n1))
    N2 = np.zeros((n2))
    S = 0
    N_wtd = np.zeros((n1, n2)) 
    
    for i in range(n1):
        N1[i] = coxDeBoor(xi, i, p1, knot_vect1)

    for i in range(n2):
        N2[i] = coxDeBoor(neta, i, p2, knot_vect2)

    for i in range(n1):
        for j in range(n2):
            N_wtd[i,j] = N1[i]*N2[j]*w[i,j]

    Wtd_sum = np.sum(N_wtd)
    for i in range(n1):
        for j in range(n2):
            if Wtd_sum != 0:
                S = S + N_wtd[i,j]/Wtd_sum * ctrlpts[i,j]

    return S

fig = plt.figure()
ax = Axes3D(fig)

for i in range(knot_vect1.shape[0]-1):
    if knot_vect1[i] != knot_vect1[i-1]:
        Line = np.zeros((x.shape[0],3))
        for j in range(x.shape[0]):
            Line[j] = Surface(knot_vect1[i],x[j])
        ax.plot(Line[:,0], Line[:,1], Line[:,2], color='b')

for i in range(knot_vect2.shape[0]-1):
    if knot_vect2[i] != knot_vect2[i-1]:
        Line = np.zeros((x.shape[0],3))
        for j in range(x.shape[0]):
            Line[j] = Surface(x[j],knot_vect2[i])
        ax.plot(Line[:,0], Line[:,1], Line[:,2], color='b')


ax.set_zlim(-3,3)
ax.set_xlim(-3,3)

plt.axis('off')

plt.show()
        



