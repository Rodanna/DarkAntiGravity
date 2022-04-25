# -*- coding: utf-8 -*-
"""
Spyder Editor

Dies ist eine temporäre Skriptdatei.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import Distances

def poten_x(x,y,a):
    xm, xp, ym, yp = x - a/2, x + a/2, y - a/2, y + a/2
    val = xm*np.arctan(ym/xm) + xp*np.arctan(yp/xp) \
        - xm*np.arctan(yp/xm) - xp*np.arctan(ym/xp) \
        + ym*np.log(xm*xm + ym*ym)/2 + yp*np.log(xp*xp + yp*yp)/2 \
        - ym*np.log(xp*xp + ym*ym)/2 - yp*np.log(xm*xm + yp*yp)/2
    return val/np.pi

def poten_y(x,y,a):
    xm, xp, ym, yp = x - a/2, x + a/2, y - a/2, y + a/2
    val = ym*np.arctan(xm/ym) + yp*np.arctan(xp/yp) \
        - ym*np.arctan(xp/ym) - yp*np.arctan(xm/yp) \
        + xm*np.log(xm*xm + ym*ym)/2 + xp*np.log(xp*xp + yp*yp)/2 \
        - xm*np.log(xm*xm + yp*yp)/2 - xp*np.log(xp*xp + ym*ym)/2
    return val/np.pi

a0 = 1
aL = 0.6
steps = 1000

f = plt.imread('monsters.png')/256
f = f[-257:-1,:256,0]
ny,nx = f.shape
plt.imshow(f)
plt.show()

ff = np.empty((4,nx,ny),float)
ff[0][100:105,30:35] = f[100:105,30:35]
ff[1][130:135,130:135] = f[130:135,130:135]
ff[2][100:110,90:100] = f[100:110,90:100]
ff[3][55:60,140:145] = f[55:60,140:145]

z = np.array([9,1.8,2,5,2.1,7,6,3,8])
u = np.linspace(-nx/2,nx/2,nx)
v = np.linspace(-ny/2,ny/2,ny)
X,Y = np.meshgrid(u,v)
Xgrad, Ygrad = np.zeros_like(X), np.zeros_like(Y)


r = np.sqrt(X**2+Y**2)
phi = np.arctan2(Y,X)
#r = np.sqrt(X[512:768,512:768]**2+Y[512:768,512:768]**2)
#phi = np.arctan2(Y[512:768,512:768],X[512:768,512:768])

rgrad = np.loadtxt('rgrad.txt', unpack=True)
phigrad = np.loadtxt('phigrad.txt', unpack = True)

Xgrad = (rgrad*X/r-phigrad*Y/r)*5 #256 pixel
Ygrad = (rgrad*Y/r+phigrad*X/r)*5

#Xgrad[512:768,512:768] = rgrad*X[512:768,512:768]/r-phigrad*Y[512:768,512:768]/r #1280 pixel image
#Ygrad[512:768,512:768] = rgrad*Y[512:768,512:768]/r+phigrad*X[512:768,512:768]/r
Xgrad = poten_x(X,Y,50)*5 #analytical potential
Ygrad = poten_y(X,Y,50)*5


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