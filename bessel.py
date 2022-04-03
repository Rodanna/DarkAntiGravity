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

v = np.arange(0,60)
x = np.linspace(1,50,1000)
nrroots = 200

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
        r = bisection(np.linspace(xmin[0][n],xmax[0][n],1000),v[0],0.01)
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
print('J 0 has',m[0],'roots')


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
    print('J',k,'has',m[k],'roots')

rel = np.argpartition(m,-10)[-10:]   #gives the indeces corresponding to the largest elements

                
for k in rel:
    roots = np.array([])
    for n in range(0,nrroots):
        xmin[k][n] = xmin[k-1][n]+ 1  # should it be pi?
        xmax[k][n] = xmax[k-1][n]+ np.pi/2
        if xmin[k][n] >= x[0] and xmax[k][n] <= x[-1]:   
            r = bisection(np.linspace(xmin[k][n],xmax[k][n],1000),v[k],0.01)
            if r != 999:
                roots = np.append(roots,r)
        elif xmin[k][n] <= x[0] and xmax[k][n] >= x[0]:
            r1 = bisection(np.linspace(x[0],xmax[k][n],1000),v[k],0.01)
            if r1 != 999:
                roots = np.append(roots,r1)
                 #roots = np.append(roots[k],r1)
        elif xmin[k][n] <= x[-1] and xmax[k][n] >= x[-1]:
            r2 = bisection(np.linspace(xmin[k][n],x[-1],1000),v[k],0.01)
            if r2 != 999:
                roots = np.append(roots,r2)
    if k == 0:
        print('The zeros of the',0,'th. Bessel function are',roots0)
    elif k != 0:
        print('The zeros of the',k,'th. Bessel function are',roots)
