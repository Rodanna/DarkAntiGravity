#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 12:37:31 2022

@author: annarodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import Distances

a0 = 1
aL = np.linspace(0.5,0.8,4)
zL = 1/aL-1
zL = np.around(zL,2)
zS = 6 #np.linspace(6,13,7)
rmax = 150
res = 250
res2 = int(res/2)
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
pc = 1/np.tan(4.848*10**(-6))*1.495978707*10**11/c # cs
Mpc = pc*1e6 # cs
invcritdens0 = 0.49823288405658067 #(kg/m^2)^-1
H0 = 67.6 # km/(s Mpc)
kmtoMpc = 3.24078e-20

u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)



Xgrad = -np.loadtxt('Xgrad_PM.txt')*50#4
Ygrad = -np.loadtxt('Ygrad_PM.txt')*50
pot = np.loadtxt('potential_PM.txt')*50


for i in range(0,len(zL)):
    aL[i] = round(aL[i],1)
    asrc = Distances.scalefactor(zS)
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL[i],asrc)
    Dd = Distances.angular(a0,aL[i])
    critdens = (c*Ds)/(4*np.pi*G*Dd*Dds) #kg/m^2
    x,y = X-(Dds/Ds)*Xgrad/critdens, Y-(Dds/Ds)*Ygrad/critdens
        

    x0 = 0
    y0 = 0
    rad2 = (X-x0)**2 + (Y-y0)**2
    f = np.exp(-rad2/50)
    spline = RectBivariateSpline(u,u,f.T)
    g = spline.ev(x,y)
    

    plt.imshow(f.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    #plt.savefig('sourceat6')
    plt.pause(0.1)
    plt.clf()
    plt.imshow(g,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.savefig('S1.pdf')
    plt.pause(0.1)
    plt.clf()
    
    tnodim = ((X-x0)**2 + (Y-y0)**2)/2 + pot/critdens*(Dds/Ds)
    tau = (1+zL[i])*Ds*Dd/Dds*tnodim #10^-9 s due to unit microrad
    tau /= 1e12  #arrival time surface in s
                 
    tmin = np.min(tau)
    tmax = np.max(tau)
    levs = np.linspace(tmin,tmin + (tmax-tmin)/5,60)
    plt.gca().set_aspect('equal')
    


    cs = plt.gca().set_aspect('equal')
    plt.contour(Y,X,tau.T,levels=levs,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.colorbar(cs)
    plt.grid()
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.savefig('S2.pdf')
    plt.pause(0.1)
    plt.clf()
