# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 11:30:27 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt

c = 300000 # km/s
H0 = 69.32 # km/(Mpc s)
Om = 0.3089
OL = 0.6911
Orad = 10**(-4)
steps = 1000
a0 = 1
a1 = 6
a = np.linspace(a0,a1,steps)

def redshift(a):
    z = (1-a)/a
    return z

def H(a):
    return H0*np.sqrt((Om/a**3 + Orad/a**4 + OL)) # km/(Mpc s)

def comoving(a):
    da = (a1-a0)/steps
    INT = 0
    for i in range (0,steps):
        INT += c*1/(a[i]**2*H(a[i]))*da # Mpc
    return INT

def angular(a):
    com = comoving(a)
    z = redshift(a)
    return com/(1+z)

plt.figure()
plt.grid(False)
plt.plot(a,H(a))
plt.xlabel('scale factor a')
plt.ylabel('Hubble parameter H')
plt.show()


print('The comoving distance is:',comoving(a),'Mpc, hence:', comoving(a)*1.03*10**14,'lightseconds')
print('The angular diameter distance is:',angular(a)[-1],'Mpc, hence:',angular(a)[-1]*1.03*10**14, 'lightseconds')