# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 16:13:36 2022

@author: Anna Rodriguez
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt

filename = "data.h5"
Mpctom = 3.08567758128e22

with h5py.File(filename, "r") as f:
    # Print all root level object names (aka keys)
    print("Keys: %s" % f.keys())
    # get first object name/key; may or may NOT be a group
    a_group_key = list(f.keys())[0]
    ds_arr = f[a_group_key][()]  # returns as a numpy array

dens = np.zeros_like(ds_arr)
rmax = 150 #microradians
res = 4000
mass = 0

u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

for i in range(0,len(X)):
    for j in range(0,len(Y)):
        mi = int(res/2) - i + 400
        mj = int(res/2) - j - 1175
        mass += np.exp(ds_arr[i,j])
        if np.sqrt(mi**2+mj**2) <= 500:
            dens[i,j] = np.exp(ds_arr[i,j])

density = mass / (187*Mpctom)**3 # kg/m^3
criticaldensity = 1.028613775e-26 #kg/m^3
print(density)
N = criticaldensity/density
print(N)

plt.clf()
plt.title('total mass distribution')
plt.gca().set_aspect('equal')
plt.contourf(X,Y,dens,cmap='RdYlBu')
plt.colorbar()
plt.show()

np.savetxt('totalgalaxy.txt',dens)

plt.clf()
plt.title('mass distribution')
plt.gca().set_aspect('equal')
plt.contourf(X[1760:3040,200:1480],Y[1760:3040,200:1480],dens[1760:3040,200:1480],cmap='RdYlBu')
plt.colorbar()
plt.show()

np.savetxt('galaxy.txt',dens[1760:3040,200:1480])