# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 18:48:05 2022
@author: Anna Rodriguez
"""

import numpy as np
from scipy import special
import matplotlib.pyplot as plt

def bisection(x,v,tol):
    a = x[0]
    b = x[-1]
    y = int(len(x)/2)
    if np.sign(special.jv(v,a)) == np.sign(special.jv(v,b)):
        if np.sign(special.jv(v,a)) != np.sign(special.jv(v,(a+b)/2)):
            return bisection(x[0:y],v,tol)
        elif np.sign(special.jv(v,b)) != np.sign(special.jv(v,(a+b)/2)):
            return bisection(x[y:-1],v,tol)
        else:
            return 999
    m = (a+b)/2
    if np.abs(special.jv(v,m)) < tol:
        return m
    elif np.sign(special.jv(v,a)) == np.sign(special.jv(v,m)):
        return bisection(x[y:-1],v,tol) # m closer than a
    elif np.sign(special.jv(v,b)) == np.sign(special.jv(v,m)):
        return bisection(x[0:y],v,tol) # m closer than b

rmax = 150 #microradians
res = 250
invcritdens0 = 0.49823288405658067 #(kg/m^2)^-1

u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u) 
radius = np.sqrt(X**2+Y**2) 
phi = np.arctan2(Y,X)

v = np.arange(0,70)
x = np.linspace(1,150,1000)
nrroots = 300

plt.figure()
plt.grid(True)
for i in range (0,len(v)):
    plt.plot(x,special.jv(v[i],x)) #bessel functions of first order
plt.show()

m = np.zeros(len(v))
xmin = np.empty((len(v),nrroots),float)
xmax = np.empty((len(v),nrroots),float)
roots0 = np.array([])

for n in range(0,nrroots):
    xmin[0][n] = np.pi*(n+0.75)
    xmax[0][n] = np.pi*(n+0.8)
    if xmin[0][n] >= x[0] and xmax[0][n] <= x[-1]:
        r = bisection(np.linspace(xmin[0][n],xmax[0][n],1000),v[0],0.001)
        m[0] += 1
        if r != 999:
            roots0 = np.append(roots0,r)
    elif xmin[0][n] <= x[0] and xmax[0][n] >= x[0]:
        r1 = bisection(np.linspace(x[0],xmax[0][n],1000),v[0],0.001)
        if r1 != 999:
            m[0] += 1
            roots0 = np.append(roots0,r1)
    elif xmin[0][n] <= x[-1] and xmax[0][n] >= x[-1]:
        r2 = bisection(np.linspace(xmin[0][n],x[-1],1000),v[0],0.001)
        if r2 != 999:
            m[0] += 1
            roots0 = np.append(roots0,r2)
#print('J 0 has',m[0],'roots')

lyst = []
root = np.zeros((len(v),nrroots),float)

for k in v:
    roots = np.array([[]])
    for n in range(0,nrroots):
        if k == 1 and n < len(roots0):
            xmin[k][n] = roots0[n] + 1
            xmax[k][n] = roots0[n] + np.pi/2
        else:
            xmin[k][n] = root[k-1][n] + 1
            xmax[k][n] = root[k-1][n]+ np.pi/2
        if xmin[k][n] >= x[0] and xmax[k][n] <= x[-1]:  
            m = (xmax[k][n] + xmin[k][n])/2
            rr = bisection(np.linspace(xmin[k][n],m,10000),v[k],0.0001)
            r = bisection(np.linspace(m,xmax[k][n],10000),v[k],0.001)
            if r != 999:
                roots = np.append(roots,r)
                root[k][n] = r
            elif rr != 999:
                roots = np.append(roots,rr)
                root[k][n] = rr
        elif xmin[k][n] <= x[0] and xmax[k][n] >= x[0]:
            r1 = bisection(np.linspace(x[0],xmax[k][n],10000),v[k],0.0001)
            if r1 != 999:
                roots = np.append(roots,r1)
                root[k][n] = r1
        elif xmin[k][n] <= x[-1] and xmax[k][n] >= x[-1]:
            r2 = bisection(np.linspace(xmin[k][n],x[-1],10000),v[k],0.0001)
            if r2 != 999:
                roots = np.append(roots,r2)
                root[k][n] = r2
    if k == 0:
        for xr in roots0:
            lyst.append((0,xr))
    elif k != 0:
        for xr in roots:
            lyst.append((k,xr))

        
lyst.sort(key = lambda v: v[1])            
#print('Here comes the long list!')
#print(lyst)
np.savetxt('roots.txt', lyst, fmt="%.3f", header="nr  root")

k = [v[0] for v in lyst]
xr = [v[1] for v in lyst]
plt.plot(k,xr,'.')
plt.show()


for n in range(0,16): #coefficients of Fourier Bessel series
    bess = special.jv(n,radius) #alpha/rmax*r
    fnc = bess*np.cos(m*phi)
    plt.figure()
    plt.plot(x,fnc)
    plt.show()