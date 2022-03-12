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
C= 60000
R = 0.1
s = 50
a0 = 1
aL = 1.5
steps = 1000


f = plt.imread('HUBBLE.jpg')
f = f[:,:,0]
plt.imshow(f)
plt.show()


ff = np.array([f*0,f*0,f*0,f*0,f*0,f*0,f*0])
nr = len(ff)
ff[0][980:1050,970:1050] = f[980:1050,970:1050]
ff[1][980:1070,210:300] = f[980:1070,210:300]
ff[2][470:510,500:550] = f[470:510,500:550]
ff[3][620:690,630:700] = f[620:690,630:700]
ff[4][620:700,470:550] = f[620:700,470:550]
ff[5][490:570,710:780] = f[490:570,710:780]
ff[6][980:1050,970:1050] = f[980:1050,970:1050]
ff[6][470:510,500:550] = f[470:510,500:550]
ff[6][390:450,110:200] = f[390:450,110:200]
ff[6][620:690,630:700] = f[620:690,630:700]
ff[6][620:700,470:550] = f[620:700,470:550]
ff[6][490:570,710:780] = f[490:570,710:780]

z = np.array([9,1.8,2,5,2.1,7,6,3])
a = z*0
Ds = z*0
Dds = z*0
A = np.empty((nr,steps), float)
B = np.empty((nr,steps), float)

for i in range (0,nr):
    a[i] = Distances.scalefactor(z[i])
    A[i] = np.linspace(a[i],a0,steps)
    B[i] = np.linspace(a[i],aL,steps)
    Ds[i] = Distances.angular(A[i])[0]
    Dds[i] = Distances.angular(B[i])[0]


ny,nx = f.shape
u = np.linspace(-nx/2,nx/2,nx)
v = np.linspace(-ny/2,ny/2,ny)
X,Y = np.meshgrid(u,v)
x,y = np.meshgrid(u,v)
g = np.empty((nr,nx,ny),float)
x = np.empty((nr,nx,ny),float)
y = np.empty((nr,nx,ny),float)
spline = np.array([RectBivariateSpline(X[0,:],Y[:,0],ff[0].T),RectBivariateSpline(X[0,:],Y[:,0],ff[1].T),RectBivariateSpline(X[0,:],Y[:,0],ff[2].T),RectBivariateSpline(X[0,:],Y[:,0],ff[3].T),RectBivariateSpline(X[0,:],Y[:,0],ff[4].T),RectBivariateSpline(X[0,:],Y[:,0],ff[5].T),RectBivariateSpline(X[0,:],Y[:,0],ff[6].T)])


for i in range (0,nr):
    x[i],y[i] = X-(Ds[i]/Dds[i])*2*X/(1+X**2+q*Y**2)*C, Y-(Ds[i]/Dds[i])*2*q*Y/(1+X**2+q*Y**2)*C
    g[i] = spline[i].ev(x[i],y[i])
    plt.imshow(ff[i])
    plt.show()
    plt.imshow(g[i])
    plt.show()