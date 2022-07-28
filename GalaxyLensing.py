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
aL = 1/3
steps = 1000
rmax = 150
res = 100
t = 0
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
z = np.linspace(3,13,9)
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

Xgrad = -np.loadtxt('Xgrad2.txt')
Ygrad = -np.loadtxt('Ygrad2.txt')
potential = np.loadtxt('potential2.txt')

f = plt.imread('HUBBLE.jpg')/res
f = f[:100,:100,0]

for i in range(0,len(z)):
    asrc = Distances.scalefactor(z[i])
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL,asrc)
    Dd = Distances.angular(a0,aL)
    critdens = 4*np.pi*G*Dd*Dds/(c*Ds) #(kg/m^2)^-1
    critdens0 = 4*np.pi*G*Dd/(c)
    x,y = X-(Dds/Ds)*Xgrad*critdens, Y-(Dds/Ds)*Ygrad*critdens  
    
    for k in range(0,3):
        for t in range(0,3):
            f = f*0
            f[49+k,49+t] = 1
            spline = RectBivariateSpline(u,u,f.T)
            g = spline.ev(x,y)
            plt.title('source image')
            plt.imshow(f.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
            plt.pause(0.1)
            plt.clf()
            plt.title('lensed images')
            plt.imshow(g.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
            plt.pause(0.1)
            plt.clf()
