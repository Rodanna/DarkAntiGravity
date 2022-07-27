# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 17:30:02 2022

@author: Anna Rodriguez
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import Distances

a0 = 1
aL = np.array([0.75])
steps = 1000
rmax = 150
res = 100
t = 0
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
z = 2
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

Xgrad = -np.loadtxt('Xgrad2.txt')
Ygrad = -np.loadtxt('Ygrad2.txt')
potential = np.loadtxt('potential2.txt')

f = plt.imread('HUBBLE.jpg')/res
f = f[:100,:100,0]*0
f[50,50] = 1

for i in range(0,len(aL)):
    asrc = Distances.scalefactor(z)
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL[i],asrc)
    Dd = Distances.angular(a0,aL[i])
    critdens = 4*np.pi*G*Dd*Dds/(c*Ds) #(kg/m^2)^-1
    x,y = X-(Dds/Ds)*Xgrad*critdens, Y-(Dds/Ds)*Ygrad*critdens  
    spline = RectBivariateSpline(u,u,f.T)
    g = spline.ev(x,y)
    plt.title('source image')
    plt.imshow(f.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.pause(1)
    plt.clf()
    plt.title('lensed images')
    plt.imshow(g.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.pause(1)
    plt.clf()
