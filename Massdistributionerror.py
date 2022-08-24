# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 12:11:34 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special

rmax = 150 #microradians
res = 120
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u) 
lev = np.linspace(0,12,20)

galaxy = np.loadtxt('galaxy2.txt')
fit = np.zeros((42,120,120))
error = np.zeros((42,120,120))

for i in range (0,42):
    fit[i,:,:] = np.loadtxt(f'fit2_120(nr={(i+1)*50}).txt')
    error[i,:,:] = galaxy - fit[i,:,:]

for i in range (0,42):
    plt.clf()
    plt.title('galaxy cluster')
    plt.gca().set_aspect('equal')
    plt.contourf(X,Y,error[i,:,:],levels = lev, cmap='RdYlBu')
    plt.colorbar()
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.show()
    print(np.max(error[i,:,:]))