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
plt.imshow(f)
plt.show()

ff = np.empty((6,nx,ny),float)
ff[0][100:105,105:110] = f[100:105,105:110]
ff[1][130:135,120:125] = f[130:135,120:125]
ff[2][100:110,90:100] = f[100:110,90:100]
ff[3][100:110,140:145] = f[100:110,140:145]
ff[4][100:105,105:115], ff[4][130:135,120:125] = f[100:105,105:115], f[130:135,120:125]
ff[4][100:110,90:100], ff[4][100:110,140:145] = f[100:110,90:100], f[100:110,140:145]
ff[5] = f

z = np.array([9,1.8,2,5,2.1,7,6,3,8])
u = np.linspace(-rmax,rmax,256)
X,Y = np.meshgrid(u,u)

Xgrad = np.loadtxt('Xgrad.txt', unpack=True)*2
Ygrad = np.loadtxt('Ygrad.txt', unpack = True)*2

for i in range (0,len(ff)):
    a = Distances.scalefactor(z[i])
    A = np.linspace(a,a0,steps)
    B = np.linspace(a,aL,steps)
    Ds = Distances.angular(A)[0]
    Dds = Distances.angular(B)[0]    
    x,y = X-(Ds/Dds)*Xgrad, Y-(Ds/Dds)*Ygrad
    spline = RectBivariateSpline(X[0,:],Y[:,0],ff[i].T)
    g = spline.ev(x,y)
    plt.imshow(ff[i])
    plt.show()
    plt.imshow(g)
    plt.show()