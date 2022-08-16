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
res = 500
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
z = np.linspace(9,13,8)
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

Xgrad = -np.loadtxt('Xgrad2_256.txt')
Ygrad = -np.loadtxt('Ygrad2_256.txt')
potential = np.loadtxt('potential2_256.txt')

hf = plt.imread('HUBBLE.jpg')
hf = hf[1038:1154,239:355,0]
f = 0*X
xshape = int(np.shape(hf)[0]/2) #only works if difference is an even number
yshape = int(np.shape(hf)[1]/2)
f[res//2-yshape:res//2+yshape,res//2-xshape:res//2+xshape] = hf
plt.imshow(f)
plt.show()

for i in range(0,len(z)):
    asrc = Distances.scalefactor(z[i])
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL,asrc)
    Dd = Distances.angular(a0,aL)
    critdens = (c*Ds)/(4*np.pi*G*Dd*Dds) #kg/m^2
    x,y = X-(Dds/Ds)*Xgrad/critdens, Y-(Dds/Ds)*Ygrad/critdens  
    
    for k in range(0,10):
        for t in range(0,10):
            if not (t==8 and k==7):
                continue
            #f = f*0
            #f[55+k,55+t] = 1
            x0 = -5 + t
            y0 = -5 + k
            rad2 = (X-x0)**2 + (Y-y0)**2
            #f = np.exp(-rad2/5)
            spline = RectBivariateSpline(u,u,f.T)
            g = spline.ev(x,y)
            plt.title('source image')
            plt.imshow(f.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
            plt.pause(1)
            plt.clf()
            plt.title('lensed images')
            plt.imshow(g.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
            plt.pause(2)
            plt.clf()
            
            tnodim = ((X-x0)**2 + (Y-y0)**2)/2 + potential/critdens*(Dds/Ds)
            tau = (1+zL)*Ds*Dd/Dd*tnodim #10^-9 s due to unit microrad
            tau /= 1e12  #arrival time surface in s
            tmin = np.min(tau)
            tmax = np.max(tau)
            levs = np.linspace(tmin,tmin + (tmax-tmin)/5,200)
            cs = plt.contour(Y,X,tau,levels=levs)
            #plt.colorbar(cs)
            plt.title('arrival time surface')
            plt.gca().set_aspect('equal')
            plt.title('%i %i' % (t,k))
            plt.pause(2)
            plt.clf()