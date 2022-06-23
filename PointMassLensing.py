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

f = plt.imread('monsters.png')/256
f = f[-257:-1,:256,1]
ny,nx = f.shape
#plt.imshow(f)
#plt.show()

ff = np.empty((6,nx,ny),float)
ff[0][100:105,105:110] = f[100:105,105:110]
ff[1][125:131,125:131] = f[125:131,125:131]
ff[2][100:110,90:100] = f[100:110,90:100]
ff[3][100:110,140:145] = f[100:110,140:145]
ff[4][100:105,105:115], ff[4][130:135,120:125] = f[100:105,105:115], f[130:135,120:125]
ff[4][100:110,90:100], ff[4][100:110,140:145] = f[100:110,90:100], f[100:110,140:145]
ff[5] = f


z = np.array([9,1.8,2,5,2.1,7,6,3,8])
u = np.linspace(-rmax,rmax,256)
X,Y = np.meshgrid(u,u)

Xgrad = np.loadtxt('Xgrad.txt')
Ygrad = np.loadtxt('Ygrad.txt')
Xgrad = -Xgrad
Ygrad = -Ygrad
cs = plt.contour(X,Y,Xgrad)
plt.clabel(cs)
plt.show()
plt.contour(X,Y,Ygrad)
plt.show()

for i in range (0,len(ff)):
    asrc = Distances.scalefactor(z[i])
    Ds = Distances.angular(1,asrc)
    Dds = Distances.angular(aL,asrc)
    for t in range (0,250,20):
        ff[i][0,t:t+10] = ff[i][0,t:t+10]
    x,y = X-(Dds/Ds)*Xgrad, Y-(Dds/Ds)*Ygrad
    spline = RectBivariateSpline(u,u,ff[i].T)
    g = spline.ev(x,y)
    plt.imshow(ff[i].T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.pause(1)
    plt.clf()
    plt.imshow(g.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.pause(1)
    plt.clf()