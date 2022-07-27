# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:35:03 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special
from mpl_toolkits import mplot3d

rmax = 150 #microradians
res = 100 #1280

u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u) 
r = np.sqrt(X**2+Y**2) 
phi = np.arctan2(Y,X)
z = 0
w = 0
dens = 2 
k = np.zeros((len(X),len(Y)),float) #dimensionless


'''
#point mass distribution
p1 = np.array([94,96,99,212,180,150,110,120,150,89,3,120,130,98,130,167,29,130,111,98,176,1,150,88,133,70])
p2 = np.array([116,69,230,144,230,120,120,155,90,5,90,170,160,78,129,130,180,67,130,110,120,70,189,10,110,167])
for i in range(0,23):
    k[p1[i]][p2[i]] = dens
'''

'''
#square
for i in range (103,152): 
    for j in range(103,152):
        k[i][j] = dens
'''

k = np.loadtxt('galaxy2.txt')
'''

#circle
for i in range(0,len(X)):
    for j in range(0,len(Y)):
        mi = int(res/2) - i
        mj = int(res/2) - j
        if np.sqrt(mi**2+mj**2) <= 200:
            k[i][j] = dens
'''
 
            
plt.clf()
plt.title('mass distribution')
plt.gca().set_aspect('equal')
plt.contourf(X,Y,k,cmap='RdYlBu')
plt.colorbar()
plt.show()

nr, root = np.loadtxt('roots.txt', unpack=True) #import bessel roots
nr = nr[:]
root = root[:]
a = np.ones(len(nr),float)
b = np.ones(len(nr),float)
nr = [int(m) for m in nr]

Ba,Bb = 0,0
area = 4*rmax**2/res**2 #microrad**2 per pixel
for n in range(0,len(nr)): #coefficients of Fourier Bessel series
    m = nr[n]
    alpha = root[n]
    N = rmax**2*np.pi/2*(special.jv(m+1,alpha))**2 #microrad**2
    if m == 0:
        N0 = N*2 
    bess = special.jv(m,alpha/rmax*r)
    fnc = bess*np.cos(m*phi)
    fns = bess*np.sin(m*phi)
    Ba = area * np.sum(k*fnc) #microrad**2 (integral over pixel)
    Bb = area * np.sum(k*fns)
    if m != 0:
        a[n] = Ba/N #dimensionless
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
    z += 2*(rmax/alpha)**2*special.jv(m,alpha*r/rmax)*angpart #microrad**2
    w += special.jv(m,alpha*r/rmax)*angpart 

np.savetxt('potential2.txt',z)

plt.clf()
plt.title('potential')
plt.gca().set_aspect('equal')
plt.contourf(X,Y,z,cmap='RdYlBu')
plt.colorbar()
plt.show()

plt.clf()
plt.title('reconstructed density')
plt.gca().set_aspect('equal')
plt.contourf(X,Y,w,cmap='RdYlBu')
plt.colorbar()
plt.show()

rgrad = 0
phigrad = 0

for i in range(0,len(nr)):
    m = nr[i]
    alpha = root[i] #microrad
    rgrad += 2*(rmax/alpha)**2*(b[i]*np.sin(m*phi)+a[i]*np.cos(m*phi))*(alpha/rmax)*special.jvp(m,r*alpha/rmax,n=1)
    phigrad += 2*(rmax/alpha)**2*(b[i]*m*np.cos(m*phi)-a[i]*m*np.sin(m*phi))*special.jv(m,r*alpha/rmax)/r

Xgrad = rgrad*X/r-phigrad*Y/r #microradian
Ygrad = rgrad*Y/r+phigrad*X/r


np.savetxt('Xgrad2.txt',Xgrad)
np.savetxt('Ygrad2.txt', Ygrad)

u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

plt.contour(X,Y,Xgrad)
plt.title('Xgrad')
plt.show()
plt.contour(X,Y,Ygrad)
plt.title('Ygrad')
plt.show()

plt.figure()
plt.gca().set_aspect('equal')
plt.quiver(X[::5,::5],Y[::5,::5],Xgrad[::5,::5],Ygrad[::5,::5])
plt.show()


'''
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, w, rstride=1, cstride=1,cmap='viridis', edgecolor='none')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')'''