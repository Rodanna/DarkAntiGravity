# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 16:13:36 2022

@author: Anna Rodriguez
"""
import h5py
import numpy as np
import matplotlib.pyplot as plt

with h5py.File("data.h5","r") as f:
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]
    ds_arr = f[a_group_key][()]  # returns as a numpy array

rmax = 150 #microradians
res = 4000
Npix = res**2

Mpctom = 3.08567758128e+22 #m/Mpc
mdensity = 3.21956112e-27 #kg/m^3
L = Mpctom*187 #m
avden = mdensity*L #kg/m^2
mass = np.sum(10**(ds_arr)) #kg/m^2
density = 10**(ds_arr)*Npix/mass*avden #kg/m^2

u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

plt.clf()
plt.title('total mass distribution')
plt.gca().set_aspect('equal')
plt.contourf(X,Y,density,cmap='RdYlBu')
plt.colorbar()
plt.show()



for i in range(0,len(X)):
    for j in range(0,len(Y)):
        mi = int(res/2) - i + 1965
        mj = int(res/2) - j + 1353
        if np.sqrt(mi**2+mj**2) > 200:
            density[i,j] = 0
            

plt.clf()
plt.title('mass distribution')
plt.gca().set_aspect('equal')
plt.contourf(X[3640:4000,3173:3533],Y[3640:4000,3173:3533],density[3640:4000,3173:3533],cmap='RdYlBu')
plt.colorbar()
plt.show()

np.savetxt('totalgalaxy.txt',density)            
np.savetxt('galaxy4.txt',density[3640:4000,3173:3533])