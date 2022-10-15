    
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

p1 = 2 #order of the polynomial for first basis function
p2 = 1 #order of the polynomial for second basis function

n1 = 3 #number of basis vector for first basis function
n2 = 2 #number of basis vector for second basis function

knot_vect2 = np.array([0,0,1,1])

knot_vect1 = np.array([0,0,0,1,1,1])

A = np.array([[1,0,1],[1,-1,1],[0,-1,1],[-1,-1,1],[-1,0,1],[-1,1,1],[0,1,1],[1,1,1],[1,0,1],
              [1,0,-1],[1,-1,-1],[0,-1,-1],[-1,-1,-1],[-1,0,-1],[-1,1,-1],[0,1,-1],[1,1,-1],[1,0,-1]])

B = np.array([[1,1/(2**(0.5)),1],
              [1,1/(2**(0.5)),1]])


ctrlpts = np.array([[[1,0,2],[1,0,-2]],
                    [[1,-1,2],[1,-1,-2]],
                    [[0,-1,2],[0,-1,-2]]])
                    
w = np.transpose(B)

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

x = np.linspace(0,1,200)
n = np.linspace(0,1,200)
xi, neta = np.meshgrid(x,n)

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


Surf = np.zeros((len(n),len(x),3))

for i in range(len(n)):
    for j in range(len(x)):
        Surf[i,j] = Surface(xi[i,j], neta[i,j])



fig = plt.figure()

ax = Axes3D(fig)

ax.plot_surface(Surf[:,:,0], Surf[:,:,1], Surf[:,:,2],color='r')
'''
for i in range(n1):
    if i == 0:
        ax.plot(ctrlpts[i,:,0],ctrlpts[i,:,1],ctrlpts[i,:,2], marker='o', color='orange', label='control net')
    ax.plot(ctrlpts[i,:,0],ctrlpts[i,:,1],ctrlpts[i,:,2], marker='o', color='orange')

for j in range(n2):
    ax.plot(ctrlpts[:,j,0],ctrlpts[:,j,1],ctrlpts[:,j,2], marker='o', color='orange')
'''
plt.show()


