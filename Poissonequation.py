# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:35:03 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special

rmax = 150 #microradians
pixeltomicrorad = 0.11*4.848

u = np.linspace(-rmax,rmax,256)
X,Y = np.meshgrid(u,u) #microrad
r = np.sqrt(X**2+Y**2) #microrad
phi = np.arctan2(Y,X)
z = 0
dens = 2 #kg/m^2
k = np.zeros((len(X),len(Y)),float)

'''
#point mass distribution
p1 = np.array([94,96,99,212,180,150,110,120,150,89,3,120,130,98,130,167,29,130,111,98,176,1,150,88,133,70])
p2 = np.array([116,69,230,144,230,120,120,155,90,5,90,170,160,78,129,130,180,67,130,110,120,70,189,10,110,167])
for i in range(0,23):
    k[p1[i]][p2[i]] = dens
'''    
 
'''
#square
for i in range (78,177): 
    for j in range(78,177):
        k[i][j] = dens
'''

#circle
for i in range(0,len(X)):
    for j in range(0,len(Y)):
        mi = 128 - i
        mj = 128 - j
        if np.sqrt(mi**2+mj**2) <= 50:
            k[i][j] = dens

'''
#rectangle
for i in range (55,200): 
    for j in range(78,177):
        k[i][j] = dens
'''
'''
#diamond
for i in range (58,157): 
    for j in range(98,197):
        if i <j:
            k[i][j] = dens 
'''
'''
#triangle
for i in range (58,157): 
    for j in range(58,157):
        if i <j:
            k[i][j] = dens  
'''            
            
plt.clf()
plt.gca().set_aspect('equal')
plt.contourf(X,Y,k,cmap='RdYlBu')
plt.colorbar()
plt.show()

nr, root = np.loadtxt('roots.txt', unpack=True) #import bessel roots
nr = nr[:50]
root = root[:50]
a = np.ones(len(nr),float)
b = np.ones(len(nr),float)
nr = [int(m) for m in nr]

Ba,Bb = 0,0
area = 4*rmax**2/(256*256) # microrad^2/pixel^2
area = area / pixeltomicrorad**2 # microrad^2
for n in range(0,len(nr)): #coefficients of Fourier Bessel series
    m = nr[n]
    alpha = root[n]
    N = rmax**2*np.pi/2*(special.jv(m+1,alpha))**2 #microrads^2
    if m == 0:
        N0 = N*2 #microrads^2
    bess = special.jv(m,alpha/rmax*r)
    fnc = bess*np.cos(m*phi)
    fns = bess*np.sin(m*phi)
    Ba = area * np.sum(k*fnc)  #microrad^2
    Bb = area * np.sum(k*fns)
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
    z -= 2*(rmax/alpha)**2*special.jv(m,alpha*r/rmax)*angpart #microrad^2


plt.clf()
plt.gca().set_aspect('equal')
plt.contourf(X,Y,z,cmap='RdYlBu')
plt.colorbar()
plt.show()

rgrad = 0
phigrad = 0

for i in range(0,len(nr)):
    m = nr[i]
    alpha = root[i] #microrad
    rgrad -= 2*(rmax/alpha)**2*(b[i]*np.sin(m*phi)+a[i]*np.cos(m*phi))*(alpha/rmax)*special.jvp(m,r*alpha/rmax,n=1)
    phigrad -= 2*(rmax/alpha)**2*(b[i]*m*np.cos(m*phi)-a[i]*m*np.sin(m*phi))*special.jv(m,r*alpha/rmax)/r

Xgrad = rgrad*X/r-phigrad*Y/r #microrad
Ygrad = rgrad*Y/r+phigrad*X/r #microrad

np.savetxt('Xgrad.txt',Xgrad)
np.savetxt('Ygrad.txt', Ygrad)

Xgrad = Xgrad*1000
Ygrad = Ygrad*1000

plt.figure()
plt.gca().set_aspect('equal')
plt.quiver(X[::5,::5],Y[::5,::5],Xgrad[::5,::5],Ygrad[::5,::5])
plt.show()