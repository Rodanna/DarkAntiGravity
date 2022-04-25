# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:35:03 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special

rmax = 150

u = np.linspace(-rmax,rmax,256)
X,Y = np.meshgrid(u,u)
r = np.sqrt(X**2+Y**2)
phi = np.arctan2(Y,X)
z = 0

rho = np.zeros((len(X),len(Y)),float)  


mass = np.array([100,90,50,34,57,31,64,87,11,54,78,43,57,99,93,65,45,33,22,79,95,93]) #point mass distribution
p1 = np.array([180,150,110,120,150,89,178,120,130,98,130,167,129,230,111,98,176,180,3,88,33,70])
p2 = np.array([130,120,120,155,90,75,90,170,160,78,129,230,180,67,130,10,120,170,189,200,210,167])
for i in range(0,22):
    rho[p1[i]][p2[i]] = mass[i] 
'''

for i in range (78,178): #square
    for j in range(78,178):
        rho[i][j] = 1
'''

plt.clf()
plt.gca().set_aspect('equal')
plt.contourf(X,Y,rho,cmap='RdYlBu')
plt.colorbar()
plt.show()

nr, root = np.loadtxt('roots.txt', unpack=True) #import bessel roots
nr = nr[:150]
root = root[:150]
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

np.savetxt('acoefficients.txt',a)
np.savetxt('bcoefficients.txt',b)

for i in range(0,len(nr)): #Fourier Bessel series with coefficients
    m = nr[i]
    alpha = root[i]
    angpart = a[i]*np.cos(m*phi)
    if m > 0:
        angpart += b[i]*np.sin(m*phi)
    z += -2*(rmax/alpha)**2*special.jv(m,alpha*r/rmax)*angpart
    z[r > rmax] = 0

plt.clf()
plt.gca().set_aspect('equal')
plt.contourf(X,Y,z,cmap='RdYlBu')
plt.colorbar()
plt.show()

rgrad = 0
phigrad = 0

for i in range(0,len(nr)):
    m = nr[i]
    alpha = root[i]
    rgrad += -2*(rmax/alpha)**2*(b[i]*np.sin(m*phi) + a[i]*np.cos(m*phi))*(alpha/rmax)*special.jvp(m,r*alpha/rmax,n=1)
    phigrad += -2*(rmax/alpha)**2*(b[i]*m*np.cos(m*phi) - a[i]*m*np.sin(m*phi))*special.jv(m,r*alpha/rmax)/r

np.savetxt('rgrad.txt',rgrad)
np.savetxt('phigrad.txt',phigrad)

Xgrad = rgrad*X/r-phigrad*Y/r
Ygrad = rgrad*Y/r+phigrad*X/r

plt.figure()
plt.gca().set_aspect('equal')
plt.quiver(X[::10,::10],Y[::10,::10],Xgrad[::10,::10],Ygrad[::10,::10])
plt.show()