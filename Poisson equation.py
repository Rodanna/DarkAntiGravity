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
pr = np.sqrt(p1**2+p2**2)
pphi = np.arctan2(p2,p1)
for index in range(0,9):
    rho[p1[index]][p2[index]] = mass[index] 

plt.figure()
plt.imshow(rho)
plt.show()

nr, root = np.loadtxt('roots.txt', unpack=True) #import bessel roots
a = np.ones(len(nr),float)
b = np.ones(len(nr),float)

 
Ba,Bb = 0,0
area = 4*rmax**2/(256*256)
for m in range(0,len(nr)): #coefficients of Fourier Bessel series
    if m != 0:
        Na = rmax**2*np.pi/2*(special.jv(nr[m],root[m]))**2/m
    N0 = rmax**2*np.pi*special.jv(1,root[m])**2
    for k in range(0,9):
        Ba +=  area*pr[k]*mass[k]*special.jv(nr[m],root[m]*pr[k]/rmax)*np.cos(nr[m]*pphi[k])
        Bb +=  area*pr[k]*mass[k]*special.jv(nr[m],root[m]*pr[k]/rmax)*np.sin(nr[m]*pphi[k])
    if m != 0:
        a[m] = Ba/Na
        b[m] = Bb/Na
    else: 
        a[m] = Ba/N0
        b[m] = Bb/N0


for i in range(0,len(nr)): #Fourier Bessel series with coefficients
    z += special.jv(nr[i],root[i]*r/rmax)*(a[i]*np.cos(nr[i]*phi)+b[i]*np.sin(nr[i]*phi))

plt.clf()
plt.gca().set_aspect('equal')
plt.contourf(X,Y,z,cmap='RdYlBu')

def gradient(r,phi):
    rgrad = 0
    phigrad = 0
    for M in range(0,len(nr)):
        rgrad += (a[M]*np.sin(nr[M]*phi) + b[M]*np.cos(nr[M]*phi))*(r*root[M]/rmax)*special.jvp(nr[M],r*root[M]/rmax,n=1)
        phigrad += a[M]*nr[M]*np.cos(nr[M]*phi) - b[M]*nr[M]*np.sin(nr[M]*phi)*(r*root[M]/rmax)*special.jv(nr[M],r*root[M]/rmax)/r
    return rgrad,phigrad    

plt.figure()
plt.quiver(X,Y,gradient(r,phi)[0]*X/r-gradient(r,phi)[1]*Y/r,gradient(r,phi)[0]*Y/r+gradient(r,phi)[1]*X/r)
plt.show()