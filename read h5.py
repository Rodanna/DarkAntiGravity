# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 16:13:36 2022

@author: Anna Rodriguez
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt

filename = "data.h5"

with h5py.File(filename, "r") as f:
    # Print all root level object names (aka keys) 
    # these can be group or dataset names 
    print("Keys: %s" % f.keys())
    # get first object name/key; may or may NOT be a group
    a_group_key = list(f.keys())[0]
    ds_arr = f[a_group_key][()]  # returns as a numpy array

rmax = 150 #microradians
res = 4000

u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)
r = np.sqrt(X**2+Y**2) 
phi = np.arctan2(Y,X)

plt.clf()
plt.title('mass distribution')
plt.gca().set_aspect('equal')
plt.contourf(X,Y,ds_arr,cmap='RdYlBu')
plt.colorbar()
plt.show()