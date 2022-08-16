# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:35:03 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special

rmax = 150 #microradians
res = 120
critdens0 = 0.49823288405658067

u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u) 
r = np.sqrt(X**2+Y**2) 
phi = np.arctan2(Y,X)
pot = 0
fit = 0
d = 8 
Dd = 1.495717907193709e+17
dens = np.zeros((len(X),len(Y)),float) #kg/m^2

'''
#point mass distribution
p1 = np.array([94,96,99,212,180,150,110,120,150,89,3,120,130,98,130,167,29,130,111,98,176,1,150,88,133,70])
p2 = np.array([116,69,230,144,230,120,120,155,90,5,90,170,160,78,129,130,180,67,130,110,120,70,189,10,110,167])
for i in range(0,23):
    dens[p1[i]][p2[i]] = d
'''

dens = np.loadtxt('galaxy2.txt')

levs = np.linspace(1,8,2) #for kappa == 1 surface
lev = np.linspace(0,12,20)
           
plt.clf()
plt.title('galaxy cluster')
plt.gca().set_aspect('equal')
plt.contourf(X,Y,dens*critdens0,levels = lev, cmap='RdYlBu')
plt.colorbar()
plt.contour(X,Y,dens*critdens0,levels=levs,cmap='gist_gray',linewidths=0.75)
plt.xlabel('x in microradians')
plt.ylabel('y in microradians')
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, dens*critdens0, rstride=1, cstride=1,cmap='RdYlBu', edgecolor='none')
ax.set_xlabel('x in microradians')
ax.set_ylabel('y in microradians')
ax.set_zlabel('convergence')
plt.show()

nr, root = np.loadtxt('roots.txt', unpack=True) #import bessel roots
nr = nr[:500]
root = root[:500]
a = np.ones(len(nr),float)
b = np.ones(len(nr),float)
nr = [int(m) for m in nr]

Ba,Bb = 0,0
area = 4*rmax**2/res**2 #microrad^2 per pixel
for n in range(0,len(nr)): #coefficients of Fourier Bessel series
    m = nr[n]
    alpha = root[n]
    N = rmax**2*np.pi/2*(special.jv(m+1,alpha))**2 #microrad^2
    if m == 0:
        N0 = N*2 
    bess = special.jv(m,alpha/rmax*r)
    fnc = bess*np.cos(m*phi)
    fns = bess*np.sin(m*phi)
    Ba = area*np.sum(dens*fnc) #microrad^2*kg/m^2 (integral over pixel)
    Bb = area*np.sum(dens*fns)
    if m != 0:
        a[n] = Ba/N
        b[n] = Bb/N
    else: 
        a[n] = Ba/N0
        b[n] = Bb/N0

res2 = 500
u = np.linspace(-rmax,rmax,res2)
X,Y = np.meshgrid(u,u) 
r = np.sqrt(X**2+Y**2) 
phi = np.arctan2(Y,X)

for i in range(0,len(nr)): #Fourier Bessel series with coefficients
    m = nr[i]
    alpha = root[i]
    angpart = a[i]*np.cos(m*phi)
    if m > 0:
        angpart += b[i]*np.sin(m*phi)
    pot += 2*(rmax/alpha)**2*special.jv(m,alpha*r/rmax)*angpart #microrad**2*kg/m´2
    fit += special.jv(m,alpha*r/rmax)*angpart  #kg/m^2

np.savetxt('potential2_256.txt',pot)

plt.clf()
plt.title('potential')
plt.gca().set_aspect('equal')
plt.contourf(X,Y,pot,cmap='RdYlBu')
plt.colorbar()
plt.show()

plt.clf()
plt.title('galaxy cluster')
plt.gca().set_aspect('equal')
plt.contourf(X,Y,fit*critdens0,levels = lev, cmap='RdYlBu')
plt.colorbar()
plt.contour(X,Y,fit*critdens0,levels=levs,cmap='gist_gray',linewidths=0.75)
plt.xlabel('x in microradians')
plt.ylabel('y in microradians')
plt.show()

rgrad = 0
phigrad = 0

for i in range(0,len(nr)):
    m = nr[i]
    alpha = root[i] #microrad
    rgrad += 2*(rmax/alpha)**2*(b[i]*np.sin(m*phi)+a[i]*np.cos(m*phi))*(alpha/rmax)*special.jvp(m,r*alpha/rmax,n=1)
    phigrad += 2*(rmax/alpha)**2*(b[i]*m*np.cos(m*phi)-a[i]*m*np.sin(m*phi))*special.jv(m,r*alpha/rmax)/r

Xgrad = rgrad*X/r-phigrad*Y/r #microradian*kg/m^2
Ygrad = rgrad*Y/r+phigrad*X/r

np.savetxt('Xgrad2_256.txt',Xgrad)
np.savetxt('Ygrad2_256.txt', Ygrad)

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

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, fit*critdens0, rstride=1, cstride=1,cmap='RdYlBu', edgecolor='none')
ax.set_xlabel('x in microradians')
ax.set_ylabel('y in microradians')
ax.set_zlabel('convergence')
plt.show()