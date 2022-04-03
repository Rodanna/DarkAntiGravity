# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:35:03 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special
from scipy.interpolate import RectBivariateSpline

m = 5
a = 3
b = 5
rmax = 80

v = np.arange(0,20)
x = np.linspace(-100,100,500)
y = np.linspace(-100,100,500)
X,Y = np.meshgrid(x,y)
Z = np.zeros((len(x),len(y)),float)

def angle(x,y,m):
    psi = 0
    for M in range(0,m):
        psi += a*np.sin(M/np.tan(y/x)) + b*np.sin(M/np.tan(y/x))
    return psi

def bessel(x,y):
    psi = 0
    for i in range(0,len(v)):
        psi += special.jv(v[i],np.sqrt(x**2+y**2)/rmax) #root is still needed
    return psi


for i in range(0,len(x)):
    for j in range(0,len(y)):
        if np.sqrt(X[0,i]**2+Y[j,0]**2) <= rmax:
            Z[i,j] = angle(X[0,i],Y[j,0],m)*bessel(X[0,i],Y[i,0])
        

        
plt.figure()
plt.imshow(Z)
plt.show()