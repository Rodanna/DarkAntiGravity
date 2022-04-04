# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:35:03 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special

nr, root = np.loadtxt('roots.txt', unpack=True)
a = np.ones(len(nr),float)
b = np.ones(len(nr),float)

rmax = 1
mass = 100

u = np.linspace(-1,1,256)
X,Y = np.meshgrid(u,u)
r = np.sqrt(X**2+Y**2)
phi = np.arctan2(Y,X)

   
for i in range(0,len(nr)):
    m = nr[i]
    km = root[i]
    z = special.jv(m,km*r/rmax) * (a[i]*np.cos(m*phi)+b[i]*np.sin(m*phi))
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.contourf(X,Y,z,cmap='RdYlBu')
    plt.pause(0.5)


'''
rho = np.zeros((len(x),len(y)),float) # point mass distribution
rho[180][130] = mass
rho[130][120] = mass 
rho[150][100] = mass
rho[110][155] = mass
rho[120][90] = mass
rho[178][75] = mass
rho[120][100] = mass
rho[150][170] = mass
rho[89][160] = mass

plt.figure()
plt.imshow(rho)
plt.show()'''

'''def laplace(x,y):
    psi = 0
    r = np.sqrt(x**2+y**2)
    phi = np.arctan2(y,x)
    for M in range(0,len(number)):
        psi += (a[M]*np.sin(M*phi)+ b[M]*np.cos(M*phi))*(root[M]/rmax)**2*(special.jvp(number[M],r*root[M]/rmax,n=2)+(rmax/r*alpha[M]))*special.jvp(number[M],r*root[M]/rmax,n=1)
        psi += a[M]*M*np.cos(M*phi) - b[M]*M*np.sin(M*phi)*special.jv(number[M],r*root[M]/rmax)/r
    return psi

def gradient(x,y):
    xgrad = 0
    ygrad = 0
    r = np.sqrt(x**2+y**2)
    phi = np.arctan2(y,x)
    for M in range(0,len(number)):
        xgrad += (a[M]*np.sin(M*phi) + b[M]*np.cos(M*phi))*(r*root[M]/rmax)*special.jvp(number[M],r*root[M]/rmax,n=1)
        ygrad += a[M]*M*np.cos(M*phi) - b[M]*M*np.sin(M*phi)*(r*root[M]/rmax)*special.jv(number[M],r*root[M]/rmax)/r
    return xgrad,ygrad'''

