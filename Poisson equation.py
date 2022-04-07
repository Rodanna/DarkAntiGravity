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
    Na = rmax**2*np.pi/2*(special.jv(nr[m],root[m]))**2
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

    

#printing all Bessel functions individually
'''
for i in range(0,len(nr)):
    m = nr[i]
    km = root[i]
    z = special.jv(nr[i],root[i]*r/rmax) * (a[i]*np.cos(nr[i]*phi)+b[i]*np.sin(nr[i]*phi))
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.contourf(X,Y,z,cmap='RdYlBu')
    plt.pause(0.5)'''


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
    return xgrad,ygrad

Ca = integrate.quad(lambda x: np.cos(m*x),0,2*np.pi)[0]
Cb = integrate.quad(lambda x: np.sin(m*x),0,2*np.pi)[0]
    '''

