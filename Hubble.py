# -*- coding: utf-8 -*-
"""
Spyder Editor

Dies ist eine temporäre Skriptdatei.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import Distances


q = 0.8
C= 100000
a0 = 1
aL = 0.6
steps = 1000

f = plt.imread('HUBBLE.jpg')
f = f[:,:,0]
ny,nx = f.shape
plt.imshow(f)
plt.show()

ff = np.empty((8,nx,ny),float)
ff[0][980:1050,970:1050] = f[980:1050,970:1050]
ff[1][980:1070,210:300] = f[980:1070,210:300]
ff[2][470:510,500:550] = f[470:510,500:550]
ff[3][620:690,630:700] = f[620:690,630:700]
ff[4][620:700,470:550] = f[620:700,470:550]
ff[5][490:570,710:780] = f[490:570,710:780]
ff[6][730:800,400:450] = f[730:800,400:450]
ff[7][980:1050,970:1050] = f[980:1050,970:1050]
ff[7][470:510,500:550] = f[470:510,500:550]
ff[7][390:450,110:200] = f[390:450,110:200]
ff[7][620:690,630:700] = f[620:690,630:700]
ff[7][620:700,470:550] = f[620:700,470:550]
ff[7][490:570,710:780] = f[490:570,710:780]
ff[7][730:800,400:450] = f[730:800,400:450]

z = np.array([9,1.8,2,5,2.1,7,6,3,8])
u = np.linspace(-nx/2,nx/2,nx)
v = np.linspace(-ny/2,ny/2,ny)
X,Y = np.meshgrid(u,v)
r = np.sqrt(X**2+Y**2)
phi = np.arctan2(Y,X)

Xgrad = 2*X/(1+X**2+q*Y**2)*C #analytic functions
Ygrad = 2*q*Y/(1+X**2+q*Y**2)*C

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