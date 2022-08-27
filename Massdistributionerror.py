# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 12:11:34 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt

rmax = 150 #microradians
res = 120
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)
r = np.sqrt(X**2+Y**2) 
phi = np.arctan2(Y,X)
lev = np.linspace(0,12,20)

galaxy = np.loadtxt('galaxy2.txt')
for k in range(0,len(X)):
    for l in range(0,len(Y)):
        mk = int(res/2) - k
        ml = int(res/2) - l
        if np.sqrt(mk**2+ml**2) > 60:
            galaxy[k,l] = 0
            
fit = np.zeros((42,120,120))
error = np.zeros((42,120,120))

for i in range (0,42):
    fit[i,:,:] = np.loadtxt(f'fit2_120(nr={(i+1)*50}).txt')
    for k in range(0,len(X)):
        for l in range(0,len(Y)):
            mk = int(res/2) - k
            ml = int(res/2) - l
            if np.sqrt(mk**2+ml**2) > 60:
                fit[i,k,l] = 0
    error[i,:,:] = galaxy - fit[i,:,:]

error1d  = np.zeros(42)

for i in range (0,42):
    plt.clf()
    plt.title('galaxy cluster')
    plt.gca().set_aspect('equal')
    plt.contourf(X,Y,error[i,:,:],levels = lev, cmap='RdYlBu') 
    plt.colorbar()
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.show()
    error1d[i] = np.mean(error[i])

x = np.arange(0,42)

plt.figure()
plt.plot(x,error1d,'o')
plt.show()