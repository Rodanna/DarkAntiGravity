# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 11:30:27 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt

c = 2.99792458e8 # m/s
pc = 1/np.tan(4.848*10**(-6))*1.495978707*10**11/c # cs
Mpc = pc*1e6 # cs
H0 = 69.32/Mpc*1000/c # 1/s 
Om = 0.3089
OL = 0.6911
Orad = 1e-4
steps = 1000
a0 = 0.0827129859387924
a1 = 1
a = np.linspace(a0,a1,steps)

def scalefactor(z):
    return 1/(1+z)

def redshift(a):
    z = (1-a)/a
    return z

def H(a):
    return H0*np.sqrt((Om/a**3 + Orad/a**4 + OL)) # 1/s

def comoving(a0,a1):
    da = (a0-a1)/steps
    INT = 0
    for i in range (0,steps):
        INT += 1/(a[i]**2*H(a[i]))*da # cs
    return INT

def lighttravel(a0,a1):
    da = (a0-a1)/steps
    INT = 0
    for i in range (0,steps):
        INT += 1/(a[i]*H(a[i]))*da # cs
    return INT

def angular(a0,a1):
    com = comoving(a0,a1)
    z = redshift(a1)
    return com/(1+z)

def luminosity(a):
    com = comoving(a)
    z = redshift(a)
    return com*(1+z)
    

'''plt.figure()
plt.grid(False)
plt.plot(a,H(a))
plt.xlabel('scale factor a')
plt.ylabel('Hubble parameter H')
plt.show()'''


#print('The comoving distance is:',comoving(a),'cs')
#print('The angular diameter distance is:',angular(a)[0],'cs')
#print('The luminosity distance is:',luminosity(a)[0],'cs')
#print('The light travel distance is:',lighttravel(a),'cs')
#print('The lighttravel distance is:',lighttravel(a)*3.171e-17,'Gly')