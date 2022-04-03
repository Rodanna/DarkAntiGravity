# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:35:03 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special

m = np.arange(0,10)
a = np.ones(len(m),float)
b = np.ones(len(m),float)
alpha = np.ones(len(m),float)
rmax = 80
mass = 100

v = np.arange(0,20)
x = np.linspace(-100,100,250)
y = np.linspace(-100,100,250)
X,Y = np.meshgrid(x,y)
Z = np.zeros((len(x),len(y)),float)
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
plt.show()

def fourierbessel(x,y,m):
    psi = 0
    for M in m:
        psi += a[M]*np.sin(M/np.tan(y/x)) + b[M]*np.cos(M/np.tan(y/x))*special.jv(v[M],np.sqrt(x**2+y**2)*alpha[M]/rmax)
    return psi

def laplace(x,y,m):
    psi = 0
    for M in m:
        psi += (a[M]*np.sin(M/np.tan(y/x))+ b[M]*np.cos(M/np.tan(y/x)))*(alpha[M]/rmax)**2*(special.jvp(v[M],np.sqrt(x**2+y**2)*alpha[M]/rmax,n=2)+(rmax/(np.sqrt(x**2+y**2)*alpha[M]))*special.jvp(v[M],np.sqrt(x**2+y**2)*alpha[M]/rmax,n=1))
        psi += a[M]*M*np.cos(M/np.tan(y/x)) - b[M]*M*np.sin(M/np.tan(y/x))*special.jv(v[M],np.sqrt(x**2+y**2)*alpha[M]/rmax)/np.sqrt(x**2+y**2)
    return psi

def gradient(x,y,m):
    xgrad = 0
    ygrad = 0
    for M in m:
        xgrad += (a[M]*np.sin(M/np.tan(y/x)) + b[M]*np.cos(M/np.tan(y/x)))*(np.sqrt(x**2+y**2)*alpha[M]/rmax)*special.jvp(v[M],np.sqrt(x**2+y**2)*alpha[M]/rmax,n=1)
        ygrad += a[M]*M*np.cos(M/np.tan(y/x)) - b[M]*M*np.sin(M/np.tan(y/x))*(np.sqrt(x**2+y**2)*alpha[M]/rmax)*special.jv(v[M],np.sqrt(x**2+y**2)*alpha[M]/rmax)/np.sqrt(x**2+y**2)
    return xgrad,ygrad


for i in range(0,len(x)):
    for j in range(0,len(y)):
        if np.sqrt(X[0,i]**2+Y[j,0]**2) <= rmax:
            Z[i,j] = fourierbessel(X[0,i],Y[j,0],m)
            
plt.figure()
plt.imshow(Z)
plt.show()

for i in range(0,len(x)):
    for j in range(0,len(y)):
        if np.sqrt(X[0,i]**2+Y[j,0]**2) <= rmax:
            Z[i,j] = laplace(X[0,i],Y[j,0],m)

plt.figure()
plt.imshow(Z)
plt.show()