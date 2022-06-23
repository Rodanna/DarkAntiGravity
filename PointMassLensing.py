# -*- coding: utf-8 -*-
"""
Spyder Editor

Dies ist eine temporäre Skriptdatei.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import Distances

a0 = 1
aL = 0.6
steps = 1000
rmax = 150
res = 1280

f = plt.imread('HUBBLE.jpg')/res
f = f[:,:,0]
ny,nx = f.shape

ff = np.empty((8,nx,ny),float)

ff[0][630:680,640:690] = f[630:680,640:690]
ff[1][600:650,500:520] = f[600:650,500:520]
ff[2][500:550,700:750] = f[500:550,700:750]
ff[3][530:580,700:750] = f[530:580,700:750]
ff[4][620:640,650:700] = f[620:640,650:700]
ff[5][600:610,600:620] = f[600:610,600:620]
ff[6][630:680,640:690] = f[630:680,640:690]
ff[6][600:650,500:520] = f[600:650,500:520]
ff[6][500:550,700:750] = f[500:550,700:750]
ff[6][530:580,700:750] = f[530:580,700:750]
ff[6][620:640,650:700] = f[620:640,650:700]
ff[6][600:610,600:620] = f[600:610,600:620]
ff[7] = f



z = np.array([9,1.8,2,5,2.1,7,6,3,8])
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

Xgrad = -np.loadtxt('Xgrad.txt')
Ygrad = -np.loadtxt('Ygrad.txt')

'''
cs = plt.contour(X,Y,Xgrad)
plt.clabel(cs)
plt.show()
plt.contour(X,Y,Ygrad)
plt.show()'''

for i in range (0,len(ff)-1):
    asrc = Distances.scalefactor(z[i])
    Ds = Distances.angular(1,asrc)
    Dds = Distances.angular(aL,asrc)
    x,y = X-(Dds/Ds)*Xgrad, Y-(Dds/Ds)*Ygrad
    spline = RectBivariateSpline(u,u,ff[i].T)
    g = spline.ev(x,y)
    plt.imshow(ff[i].T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.pause(1)
    plt.clf()
    plt.imshow(g.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.pause(1)
    plt.clf()