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

f = plt.imread('HUBBLE.jpg')/256
f = f[:,:,0]
ny,nx = f.shape
plt.imshow(f)
plt.show()

ff = np.empty((8,nx,ny),float)
''''ff[0][90:100,110:120] = f[90:100,110:120]'''

ff[0][512:562,720:760] = f[512:562,720:760]
ff[1][612:662,730:768] = f[612:662,730:768]
ff[2][690:720,512:530] = f[690:720,512:530]
ff[3][620:672,652:680] = f[620:672,652:680]
ff[4][720:750,610:640] = f[720:750,610:640]


z = np.array([9,1.8,2,5,2.1,7,6,3,8])
u = np.linspace(-nx/2,nx/2,nx)
v = np.linspace(-ny/2,ny/2,ny)
X,Y = np.meshgrid(u,v)
Xgrad, Ygrad = np.zeros_like(X), np.zeros_like(Y)
#r = np.sqrt(X**2+Y**2)
#phi = np.arctan2(Y,X)
r = np.sqrt(X[512:768,512:768]**2+Y[512:768,512:768]**2)
phi = np.arctan2(Y[512:768,512:768],X[512:768,512:768])

rgrad = np.loadtxt('rgrad.txt', unpack=True)*10
phigrad = np.loadtxt('phigrad.txt', unpack = True)*10

Xgrad[512:768,512:768] = rgrad*Y[512:768,512:768]/r+phigrad*X[512:768,512:768]/r #1280 pixel image
Ygrad[512:768,512:768] = rgrad*Y[512:768,512:768]/r+phigrad*X[512:768,512:768]/r

#Xgrad = rgrad*Y/r+phigrad*X/r #256 pixel
#Ygrad = rgrad*Y/r+phigrad*X/r

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