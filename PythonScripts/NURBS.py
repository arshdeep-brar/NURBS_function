
# coding: utf-8

# In[17]:


import numpy as np
import matplotlib.pyplot as plt

p = 3 #order of the polynomial

n = 5 #number of basis vector

knot_vect = np.array([0,0,0,0,0.5,1,1,1,1])

ctrlpts = np.array([[0,0],[0.2,1],[0.6,0.5],[0.8,1],[1,0]])

s = 2**(0.5)/2
print(s)

w = np.array([1,1,2,1.5,1])

#Calculating the basis vector

def coxDeBoor(zeta, i, d):
    
    if d==0:
        if knot_vect[i] <= zeta and zeta <= knot_vect[i+1]:
            return 1
        return 0
    
    D1 = knot_vect[i+d] - knot_vect[i]
    D2 = knot_vect[i+d+1] - knot_vect[i+1]
    
    E1 = 0
    E2 = 0
    
    if D1 > 0:
        E1 = (zeta - knot_vect[i])/D1 * coxDeBoor(zeta, i, d-1)
    
    if D2 > 0:
        E2 = (knot_vect[i+d+1] - zeta)/D2 * coxDeBoor(zeta, i+1, d-1)
        
    return E1+E2

zeta = np.linspace(0,1,100)
N = np.zeros((100,n))
Wtd_N = np.zeros((100,n))
Curve = np.zeros((100,2))

for i in range(n):
    coxDeBoor(1,i,2)


for z in range(len(zeta)):
    for i in range(n):
        N[z,i] = coxDeBoor(zeta[z], i, p)
        Wtd_N[z,i] = N[z,i]*w[i]

Wtd_sum = np.sum(Wtd_N,1)
for z in range(len(zeta)):        
    for i in range(n):
        if Wtd_sum[z] != 0:
            Curve[z] = Curve[z] + Wtd_N[z,i]/Wtd_sum[z] * ctrlpts[i]


plt.plot(Curve[:,0],Curve[:,1],label='NURBS curve')
plt.plot(ctrlpts[:,0],ctrlpts[:,1], label='Control points interpolation', marker='o')
plt.axis('equal')
plt.legend()
plt.show()

