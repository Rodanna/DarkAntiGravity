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
mass = np.sum(np.exp(ds_arr)) #kg/m^2
density = np.exp(ds_arr)*Npix/mass*avden #kg/m^2

u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

plt.clf()
plt.title('total mass distribution')
plt.gca().set_aspect('equal')
plt.contourf(X,Y,density,cmap='RdYlBu')
plt.colorbar()
plt.show()


'''
for i in range(0,len(X)):
    for j in range(0,len(Y)):
        mi = int(res/2) - i + 400
        mj = int(res/2) - j - 1175
        if np.sqrt(mi**2+mj**2) > 500:
            density[i,j] = 0

plt.clf()
plt.title('mass distribution')
plt.gca().set_aspect('equal')
plt.contourf(X[1760:3040,200:1480],Y[1760:3040,200:1480],density[1760:3040,200:1480],cmap='RdYlBu')
plt.colorbar()
plt.show()

np.savetxt('totalgalaxy.txt',density)            
np.savetxt('galaxy.txt',density[1760:3040,200:1480])'''