# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 17:30:02 2022
@author: Anna Rodriguez
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import Distances
import math

a0 = 1
aL = 0.65
zL = 1/aL-1
rmax = 150
rmax2 = 150
res = 2048
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
z = np.linspace(9,13,8)
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

Xgrad = -np.loadtxt('Xgrad2_2048.txt')
Ygrad = -np.loadtxt('Ygrad2_2048.txt')
potential = np.loadtxt('potential2_2048.txt')

#hf = plt.imread('HUBBLE.jpg')
#hf = hf[:,:,0]
#hf = hf[1050:1154,251:355,0]
#hf = hf[626:651,497:521,0]
hf = plt.imread('cutout2.jpg')[:,:,1]
f = 0*X
yshape1 = int(np.shape(hf)[0]/2) #round down
yshape2 = int(math.ceil(np.shape(hf)[0]/2)) #round up
xshape1 = int(np.shape(hf)[1]/2)
xshape2 = int(math.ceil(np.shape(hf)[1]/2))
x0 = 0
y0 = 0
xdisp = res//2 + x0
ydisp = res//2 + y0

f[ydisp-yshape1:ydisp+yshape2,xdisp-xshape1:xdisp+xshape2] = hf

plt.title('Hubble Deep Field')
plt.imshow(hf,origin='lower',extent=[-rmax2,rmax2,-rmax2,rmax2])
plt.xlabel('x in microradians')
plt.ylabel('y in microradians')
plt.show()
plt.clf()

for i in range(len(z)-1,len(z)):
    asrc = Distances.scalefactor(z[i])
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL,asrc)
    Dd = Distances.angular(a0,aL)
    critdens = (c*Ds)/(4*np.pi*G*Dd*Dds) #kg/m^2
    x,y = X-(Dds/Ds)*Xgrad/critdens, Y-(Dds/Ds)*Ygrad/critdens  
    
    spline = RectBivariateSpline(u,u,f.T)
    g = spline.ev(x,y)
    plt.title('source image')
    plt.imshow(f.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.pause(1)
    plt.clf()
    plt.title('Lensed Hubble Deep Field')
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.imshow(g.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.pause(2)
    plt.clf()
    
    tnodim = ((X-x0)**2 + (Y-y0)**2)/2 + potential/critdens*(Dds/Ds)
    tau = (1+zL)*Ds*Dd/Dds*tnodim #10^-9 s due to unit microrad
    tau /= 1e12  #arrival time surface in s                           
    tmin = np.min(tau)
    tmax = np.max(tau)
    levs = np.linspace(tmin,tmin + (tmax-tmin)/5,100)
    plt.gca().set_aspect('equal')
    plt.contour(Y,X,tau,levels=levs)
    plt.grid()
    plt.title('arrival time surface')
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    #plt.savefig(f'timedelay{i}')
    plt.pause(0.1)
    plt.clf()
    