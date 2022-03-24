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
        return 999
    m = (a+b)/2
    if np.abs(special.jv(v,m)) < tol:
        return m
    elif np.sign(special.jv(v,a)) == np.sign(special.jv(v,m)):
        return bisection(x[y:-1],v,tol) # m closer than a
    elif np.sign(special.jv(v,b)) == np.sign(special.jv(v,m)):
        return bisection(x[0:y],v,tol) # m closer than b

v = np.arange(0,6)
x = np.linspace(0,65,1000)
nrroots = 100

plt.figure()
plt.grid(True)
for i in range (2,3):
    plt.plot(x,special.jv(v[i],x)) #bessel functions of first order
plt.show()

'''plt.figure() #first derivative
plt.grid(True)
for i in range (0,len(v)):
    plt.plot(x,special.jvp(v[i],x,n=1)) 
plt.show()'''

m = np.zeros(len(v))
xmin = np.empty((len(v),nrroots),float)
xmax = np.empty((len(v),nrroots),float)

for n in range(0,nrroots):
    xmin[0][n] = np.pi*(n+0.75)
    xmax[0][n] = np.pi*(n+0.8)
    if xmin[0][n] >= x[0] and xmax[0][n] <= x[-1]:
        m[0] += 1
    elif xmin[0][n] <= x[0] and xmax[0][n] >= x[0]:
        r1 = bisection(np.linspace(x[0],xmax[0][n],1000),v[0],0.001)
        if r1 != 999:
            m[0] += 1
    elif xmin[0][n] <= x[-1] and xmax[0][n] >= x[-1]:
        r1 = bisection(np.linspace(xmin[0][n],x[-1],1000),v[0],0.001)
        if r1 != 999:
            m[0] += 1
print(m[0])


for k in range(1,len(v)):
    for n in range(0,nrroots):
        xmin[k][n] = xmin[k-1][n]+ 1  # should it be pi?
        xmax[k][n] = xmax[k-1][n]+ np.pi/2
        if xmin[k][n] >= x[0] and xmax[k][n] <= x[-1]:
            m[k] += 1
        elif xmin[k][n] <= x[0] and xmax[k][n] >= x[0]:
            r1 = bisection(np.linspace(x[0],xmax[k][n],1000),v[k],0.01)
            if r1 != 999:
                m[k] += 1
        elif xmin[k][n] <= x[-1] and xmax[k][n] >= x[-1]:
            r1 = bisection(np.linspace(xmin[k][n],x[-1],1000),v[k],0.01)
            if r1 != 999:
                m[k] += 1
    print(m[k])
    