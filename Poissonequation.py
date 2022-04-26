# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:35:03 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special
import Distances
from scipy.interpolate import RectBivariateSpline

rmax = 150

u = np.linspace(-rmax,rmax,256)
X,Y = np.meshgrid(u,u)
r = np.sqrt(X**2+Y**2)
phi = np.arctan2(Y,X)
z = 0

rho = np.zeros((len(X),len(Y)),float)  


mass = 10 #point mass distribution
p1 = np.array([180,150,110,120,150,89,3,120,130,98,130,167,29,130,111,98,176,1,150,88,133,70])
p2 = np.array([230,120,120,155,90,5,90,170,160,78,129,130,180,67,130,110,120,70,189,10,110,167])
for i in range(0,22):
    rho[p1[i]][p2[i]] = mass
    
'''
for i in range (78,177): #square
    for j in range(78,177):
        rho[i][j] = 1
'''

plt.clf()
plt.gca().set_aspect('equal')
plt.contourf(X,Y,rho,cmap='RdYlBu')
plt.colorbar()
plt.show()

nr, root = np.loadtxt('roots.txt', unpack=True) #import bessel roots
nr = nr[:50]
root = root[:50]
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
    z -= 2*(rmax/alpha)**2*special.jv(m,alpha*r/rmax)*angpart


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
    rgrad -= 2*(rmax/alpha)**2*(b[i]*np.sin(m*phi)+a[i]*np.cos(m*phi))*(alpha/rmax)*special.jvp(m,r*alpha/rmax,n=1)
    phigrad -= 2*(rmax/alpha)**2*(b[i]*m*np.cos(m*phi)-a[i]*m*np.sin(m*phi))*special.jv(m,r*alpha/rmax)/r

np.savetxt('rgrad.txt',rgrad)
np.savetxt('phigrad.txt',phigrad)

Xgrad = rgrad*X/r-phigrad*Y/r
Ygrad = rgrad*Y/r+phigrad*X/r

np.savetxt('Xgrad.txt',Xgrad)
np.savetxt('Ygrad.txt', Ygrad)

Xgrad = Xgrad*1000
Ygrad = Ygrad*1000

plt.figure()
plt.gca().set_aspect('equal')
plt.quiver(X[::5,::5],Y[::5,::5],Xgrad[::5,::5],Ygrad[::5,::5])
plt.show()