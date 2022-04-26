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
rmax = 150

f = plt.imread('monsters.png')/256
f = f[-257:-1,:256,1]
#f = f[:,:,0]
ny,nx = f.shape
plt.imshow(f)
plt.show()


ff = np.empty((6,nx,ny),float)
ff[0][100:105,105:110] = f[100:105,105:110]
ff[1][130:135,120:125] = f[130:135,120:125]
ff[2][100:110,90:100] = f[100:110,90:100]
ff[3][100:110,140:145] = f[100:110,140:145]

ff[4][100:105,105:115] = f[100:105,105:115]
ff[4][130:135,120:125] = f[130:135,120:125]
ff[4][100:110,90:100] = f[100:110,90:100]
ff[4][100:110,140:145] = f[100:110,140:145]
ff[5] = f

z = np.array([9,1.8,2,5,2.1,7,6,3,8])
u = np.linspace(-rmax,rmax,256)
X,Y = np.meshgrid(u,u)

Ygrad = np.loadtxt('Xgrad.txt', unpack=True)
Xgrad = np.loadtxt('Ygrad.txt', unpack = True)

Xgrad = Xgrad*100
Ygrad = Ygrad*100

pi = 1
plt.figure()
plt.gca().set_aspect('equal')
plt.quiver(X[::pi,::pi],Y[::pi,::pi],Xgrad[::pi,::pi],Ygrad[::pi,::pi])
plt.show()

#Xgrad[512:768,512:768] = rgrad*X[512:768,512:768]/r-phigrad*Y[512:768,512:768]/r*100 #1280 pixel image
#Ygrad[512:768,512:768] = rgrad*Y[512:768,512:768]/r+phigrad*X[512:768,512:768]/r*100
#Xgrad = poten_x(X,Y,100)*5 #analytical potential
#Ygrad = poten_y(X,Y,100)*5


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