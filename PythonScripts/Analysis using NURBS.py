
#importing the required libraries for python

import numpy as np
import matplotlib.pyplot as plt


# In[41]:



def Rotation(theta, X):
    import math
    rad = theta*math.pi/180
    R = np.array([[math.cos(rad),0,-math.sin(rad)],
                  [0,1,0],
                  [math.sin(rad), 0, math.cos(rad)]])
    X_n = np.zeros(X.shape)
    for i in range(X.shape[0]):
        X_n[i] = R.dot(X[i])
    return X_n


#Settingg up parameters for NURBS surface

#Order of NURBS basis
p = 2
q = 1
r = 4

#Number of Basis function for NURBS
n = 3
m = 2
l = 9


# In[43]:


#Defining the knot Vectot for NURBS function

Xi = np.array([0,0,0,1,1,1])
Eta = np.array([0,0,1,1])
Zeta = np.array([0,0,0,0,0,0.2,0.4,0.6,0.8,1,1,1,1,1])


# In[44]:


#Control points array and weights 

Mod1 = np.array([[0,-20,0.5],[0,-20,0.5],[0,-20,0.5],[0,-20,0.5],[0,-20,0.5],[0,-20,0.5],[0,-20,0.5],[0,-20,0.5],[0,-20,0.5],
                    [-3,-20,0],[-3,-20,0.5],[0,-20,0.5],[3,-20,0.5],[3,-20,0],[3,-20,-0.5],[0,-20,0],[-3,-20,-0.5],[-3,-20,0.5]])

Mod2 = np.array([[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5],
                [-3.75,0,0],[-3.75,0,0.6],[0,0,0.6],[3.75,0,0.6],[3.75,0,0],[3.75,0,-0.6],[0,0,0],[-3.75,0,-0.6],[-3.75,0,0]])

Mod3 = np.array([[0,20,0.5],[0,20,0.5],[0,20,0.5],[0,20,0.5],[0,20,0.5],[0,20,0.5],[0,20,0.5],[0,20,0.5],[0,20,0.5],
             [-2.25,20,0],[-2.25,20,0.4],[0,20,0.4],[2.25,20,0.4],[2.25,20,0],[2.25,20,-0.4],[0,20,0],[-2.25,20,-0.4],[-2.25,20,0]])

ctrlpts = np.append(Mod1, Rotation(5,Mod2),axis=0)
ctrlpts = np.append(ctrlpts,Rotation(10,Mod3),axis=0)

w = np.array([1,0.75,1,0.75,1,1,5,1,1,
              1,0.75,1,0.75,1,1,5,1,1,
              1,0.75,1,0.75,1,1,5,1,1,
              1,0.75,1,0.75,1,1,5,1,1,
              1,0.75,1,0.75,1,1,5,1,1,
              1,0.75,1,0.75,1,1,5,1,1])



# In[45]:


#Converting the local indexing to the global indexing

def IEN2D(i,j,m):
    return (m*(i) + j)

def IEN3D(i,j,k,l,m):
    return ((l*m)*(i)+l*(j)+k)


# In[56]:


#Finding the basis functions for the NURBS

def coxDeBoor(zeta, i, d, knot_vect):
    '''
    zeta - value of the point we need to find
    i - i'th basis function
    d - order of basis function
    knot_vect - knot vector for NURBS    
    '''
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


# In[57]:



# In[58]:


def BasisFunction(xi, eta, zeta):
    N = np.zeros(l*m*n)
    Tot_weight = 0
    for i in range(n):
        for j in range(m):
            for k in range(l):
                N[IEN3D(i,j,k,l,m)] = coxDeBoor(xi, i, p, Xi)*coxDeBoor(eta, j, q, Eta)*coxDeBoor(zeta, k, r, Zeta)*w[IEN3D(i,j,k,l,m)]
                Tot_weight = Tot_weight + N[IEN3D(i,j,k,l,m)]
    return (N/Tot_weight)

print(IEN3D(1,1,1,l,m))

x = np.linspace(0,1,100)
e = np.linspace(0,1,100)
z = np.linspace(0,1,100)
xi, eta, zeta = np.meshgrid(x, e, z)

print(xi.shape)

P_t = np.transpose(ctrlpts)

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = Axes3D(fig)

print(xi.shape)
print(eta.shape)
print(zeta.shape)

for i in range(len(x)):
    S = np.zeros((len(e), len(z), 3))
##    for j in range(len(e)):
    for k in range(len(z)):
        S[i,k] = P_t.dot(BasisFunction(xi[99,i,k], eta[99,i,k], zeta[99,i,k]))
    ax.plot_surface(S[:,:,0], S[:,:,1], S[:,:,2], color='b')
    print(i)

ax.plot(Mod1[:,0], Mod1[:,1], Mod1[:,2], color='r')
ax.plot(Rotation(5,Mod2)[:,0], Rotation(5,Mod2)[:,1], Rotation(5,Mod2)[:,2], color='b')
ax.plot(Rotation(10,Mod3)[:,0], Rotation(10,Mod3)[:,1], Rotation(10,Mod3)[:,2], color='g')

##ax.plot(ctrlpts[:,0], ctrlpts[:,1], ctrlpts[:,2], color = 'r')

plt.axis('equal')  
ax.set_zlim(-6,6)

plt.show()

