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
aL = 0.6
zL = 1/aL-1
steps = 1000
rmax = 150
res = 120
t = 0
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
z = np.linspace(6,13,8)
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

Xgrad = -np.loadtxt('Xgrad1.txt')
Ygrad = -np.loadtxt('Ygrad1.txt')
potential = np.loadtxt('potential1.txt')

f = plt.imread('HUBBLE.jpg')/res
f = f[:120,:120,0]

for i in range(0,len(z)):
    asrc = Distances.scalefactor(z[i])
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL,asrc)
    Dd = Distances.angular(a0,aL)
    critdens = 4*np.pi*G*Dd*Dds/(c*Ds) #(kg/m^2)^-1
    x,y = X-(Dds/Ds)*Xgrad*critdens, Y-(Dds/Ds)*Ygrad*critdens  
    
    for k in range(0,10):
        for t in range(0,10):
            f = f*0
            f[55+k,55+t] = 1
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
            
            x0 = -5 + t
            y0 = -5 + k
            tnodim = (Ds/Dds)**2*((x-x0)**2 + (y-y0)**2)/2 + potential*critdens
            t = (1+zL)*Ds*Dd/Dds*tnodim #arrival time surface in s
            tmin = np.min(t)
            tmax = np.max(t)
            levs = np.linspace(tmin,tmin + (tmax-tmin)/5,40)
            plt.contour(Y,X,t,levels=levs)
            plt.title('arrival time surface')
            plt.gca().set_aspect('equal')
            plt.pause(0.1)
            plt.clf()
