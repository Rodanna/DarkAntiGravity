# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 15:54:57 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt

rmax = 150
f = plt.imread('HUBBLE.jpg')
f = f[:,:,0]
ny,nx = f.shape

a = np.array([630,626,518,980,0])
b = np.array([680,651,530,1050,1024])
c = np.array([640,497,712,965,0])
d = np.array([690,521,724,1035,1024])

ff = np.empty((5,nx,ny),float)
ff[0][a[0]:b[0],c[0]:d[0]] = f[a[0]:b[0],c[0]:d[0]]
ff[1][a[1]:b[1],c[1]:d[1]] = f[a[1]:b[1],c[1]:d[1]]
ff[2][a[2]:b[2],c[2]:d[2]] = f[a[2]:b[2],c[2]:d[2]]
ff[3][a[3]:b[3],c[3]:d[3]] = f[a[3]:b[3],c[3]:d[3]] 
ff[4] = f

for i in range (0,len(ff)):
    plt.title('source image')
    plt.imshow(ff[i].T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.pause(1)
    plt.show()
    plt.clf()


for i in range (0,len(ff)):
    f = plt.imread('HUBBLE.jpg')[:,:,0]
    f = f[a[i]:b[i],c[i]:d[i]]
    plt.figure(frameon=False)
    plt.axis('off')
    plt.imshow(f.T)
    plt.savefig(f'cutout{i}.jpg',dpi = 4, bbox_inches='tight',pad_inches = 0)
    plt.pause(1)
    plt.clf()