# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:35:03 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special

rmax = 100

u = np.linspace(-rmax,rmax,256)
X,Y = np.meshgrid(u,u)
r = np.sqrt(X**2+Y**2)
phi = np.arctan2(Y,X)
z = 0

rho = np.zeros((len(X),len(Y)),float)  #point mass distribution
mass = np.array([100,90,50,34,57,31,64,87,11])*100
p1 = np.array([180,150,110,120,150,89,178,120,130])
p2 = np.array([130,120,120,-155,90,75,90,170,160])
for i in range(0,9):
    rho[p1[i]][p2[i]] = mass[i] 

plt.figure()
plt.imshow(rho)
plt.show()

nr, root = np.loadtxt('roots.txt', unpack=True) #import bessel roots
nr = nr[:]
root = root[:]
a = np.ones(len(nr),float)
b = np.ones(len(nr),float)
nr = [round(m) for m in nr]


Ba,Bb = 0,0
area = 4*rmax**2/(256*256)
for n in range(0,len(nr)): #coefficients of Fourier Bessel series
    m = nr[n]
    alpha = root[n]
    N = rmax**2*np.pi/2*(special.jv(m+1,alpha))**2
    if m == 0:
        N0 = N*2
    bess = special.jv(m,alpha/rmax*r)
    fnc = bess*np.cos(m*phi)
    fns = bess*np.sin(m*phi)
    Ba = area * np.sum(rho*fnc)
    Bb = area * np.sum(rho*fns)
    if m != 0:
        a[n] = Ba/N
        b[n] = Bb/N
    else: 
        a[n] = Ba/N0
        b[n] = Bb/N0


for i in range(0,len(nr)): #Fourier Bessel series with coefficients
    m = nr[i]
    alpha = root[i]
    angpart = a[i]*np.cos(m*phi)
    if m > 0:
        angpart += b[i]*np.sin(m*phi)
    z += special.jv(m,alpha*r/rmax)*angpart
    z[r > rmax] = 0

plt.clf()
plt.gca().set_aspect('equal')
plt.contourf(X,Y,z,cmap='RdYlBu')
plt.colorbar()


rgrad = 0
phigrad = 0

for i in range(0,len(nr)):
    m = nr[i]
    alpha = root[i]
    rgrad += (b[i]*np.sin(m*phi) + a[i]*np.cos(m*phi))*(alpha/rmax)*special.jvp(m,r*alpha/rmax,n=1)
    phigrad += (b[i]*m*np.cos(m*phi) - a[i]*m*np.sin(m*phi))*special.jv(m,r*alpha/rmax)/r

np.savetxt('rgrad.txt',rgrad)
np.savetxt('phigrad.txt',phigrad)

plt.figure()
plt.quiver(X,Y,rgrad*X/r-phigrad*Y/r,rgrad*Y/r+phigrad*X/r)
plt.show()