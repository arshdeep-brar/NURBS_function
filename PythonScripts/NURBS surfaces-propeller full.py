    
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def Rotation(theta, X):
    rad = theta*math.pi/180
    R = np.array([[math.cos(rad),0,-math.sin(rad)],
                  [0,1,0],
                  [math.sin(rad), 0, math.cos(rad)]])
    x = np.array(X)
    x_new = R.dot(x)
    return x_new


p1 = 2 #order of the polynomial for first basis function
p2 = 4 #order of the polynomial for second basis function
p3 = 1

n1 = 3 #number of basis vector for first basis function
n2 = 9 #number of basis vector for second basis function
n3 = 2

knot_vect2 = np.array([0,0,0,0,0,0.2,0.4,0.6,0.8,1,1,1,1,1])

knot_vect1 = np.array([0,0,0,1,1,1])

knot_vect3 = np.array([0,0,1,1])
A = np.array([[0,0,0],[0,1,0],[1,2,0],[2,2,0],
              [2,1,0],[1.5,1,0],[1,0.5,0],[1,0,0]])

B = np.array([[1,1/(2**(0.5)),1],[1/(2**(0.5)),1/2.,1/(2**(0.5))],[1,1/(2**(0.5)),1]])

ctrlpts1 = np.array([[[-3,-20,0],[-3,-20,0.5],[0,-20,0.5],[3,-20,0.5],[3,-20,0],[3,-20,-0.5],[0,-20,0],[-3,-20,-0.5],[-3,-20,0]],
                    [Rotation(-5,[-4,0,0]),Rotation(-5,[-4,0,0.6]),Rotation(-5,[0,0,0.6]),Rotation(-5,[4,0,0.6]),Rotation(-5,[4,0,0]),Rotation(-5,[4,0,-0.6]),Rotation(-5,[0,0,0]),Rotation(-5,[-4,0,-0.6]),Rotation(-5,[-4,0,0])],
                    [Rotation(-10,[-1.2,20,0]),Rotation(-10,[-1.2,20,0.2]),Rotation(-10,[0,20,0.2]),Rotation(-10,[1.2,20,0.2]),Rotation(-10,[1.2,20,0]),Rotation(-10,[1.2,20,-0.2]),Rotation(-10,[0,20,0]),Rotation(-10,[-1.2,20,-0.2]),Rotation(-10,[-1.2,20,0])]])
ctrlpts2 = np.array([[[0,-20,0.25],[0,-20,0.25],[0,-20,0.25],[0,-20,0.25],[0,-20,0.25],[0,-20,0.25],[0,-20,0.25],[0,-20,0.25],[0,-20,0.25]],
                    [[-3,-20,0],[-3,-20,0.5],[0,-20,0.5],[3,-20,0.5],[3,-20,0],[3,-20,-0.5],[0,-20,0],[-3,-20,-0.5],[-3,-20,0]]])

ctrlpts3 = np.array([[[0,20,0.25],[0,20,0.25],[0,20,0.25],[0,20,0.25],[0,20,0.25],[0,20,0.25],[0,20,0.25],[0,20,0.25],[0,20,0.25]],
                     [Rotation(-10,[-1.2,20,0]),Rotation(-10,[-1.2,20,0.2]),Rotation(-10,[0,20,0.2]),Rotation(-10,[1.2,20,0.2]),Rotation(-10,[1.2,20,0]),Rotation(-10,[1.2,20,-0.2]),Rotation(-10,[0,20,0]),Rotation(-10,[-1.2,20,-0.2]),Rotation(-10,[-1.2,20,0])]])
                     

w = np.array([[1,0.75,1,0.75,1,1,5,1,1],
              [1,0.75,1,0.75,1,1,5,1,1],
              [1,0.75,1,0.75,1,1,5,1,1]])

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

x = np.linspace(0,1,100)
n = np.linspace(0,1,100)
xi, neta = np.meshgrid(x,n)



def Surface(xi, neta, l, m, k1, k2, p, q, CP):
    N1 = np.zeros((l))
    N2 = np.zeros((m))
    S = 0
    N_wtd = np.zeros((l, m)) 
    
    for i in range(l):
        N1[i] = coxDeBoor(xi, i, p, k1)

    for i in range(m):
        N2[i] = coxDeBoor(neta, i, q, k2)

    for i in range(l):
        for j in range(m):
            N_wtd[i,j] = N1[i]*N2[j]*w[i,j]

    Wtd_sum = np.sum(N_wtd)
    for i in range(l):
        for j in range(m):
            if Wtd_sum != 0:
                S = S + N_wtd[i,j]/Wtd_sum * CP[i,j]

    return S



Surf1 = np.zeros((len(n),len(x),3))

for i in range(len(n)):
    for j in range(len(x)):
        Surf1[i,j] = Surface(xi[i,j], neta[i,j], n1, n2, knot_vect1, knot_vect2 , p1, p2, ctrlpts1)

Surf2 = np.zeros((len(n),len(x),3))

for i in range(len(n)):
    for j in range(len(x)):
        Surf2[i,j] = Surface(xi[i,j], neta[i,j], n3, n2, knot_vect3, knot_vect2 , p3, p2, ctrlpts2)

Surf3 = np.zeros((len(n),len(x),3))

for i in range(len(n)):
    for j in range(len(x)):
        Surf3[i,j] = Surface(xi[i,j], neta[i,j], n3, n2, knot_vect3, knot_vect2 , p3, p2, ctrlpts3)

fig = plt.figure()

ax = Axes3D(fig)

ax.plot_surface(Surf1[:,:,0], Surf1[:,:,1], Surf1[:,:,2],color='r')
ax.plot_surface(Surf2[:,:,0], Surf2[:,:,1], Surf2[:,:,2],color='r')
ax.plot_surface(Surf3[:,:,0], Surf3[:,:,1], Surf3[:,:,2],color='r')

for i in range(n1):
    if i == 0:
        ax.plot(ctrlpts1[i,:,0],ctrlpts1[i,:,1],ctrlpts1[i,:,2], marker='o', color='orange', label='control net')
    ax.plot(ctrlpts1[i,:,0],ctrlpts1[i,:,1],ctrlpts1[i,:,2], marker='o', color='orange')
    ax.plot(ctrlpts1[:,i,0],ctrlpts1[:,i,1],ctrlpts1[:,i,2], marker='o', color='orange')

plt.axis('equal')
ax.set_zlim(-6,6)

plt.show()




