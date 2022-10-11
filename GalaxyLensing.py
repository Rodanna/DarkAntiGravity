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
res2 = int(res/2)
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
z = np.linspace(5,14,7)
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

Xgrad = -np.loadtxt('Xgrad2_2048.txt')
Ygrad = -np.loadtxt('Ygrad2_2048.txt')
potential = np.loadtxt('potential2_2048.txt')


hf = plt.imread('HUBBLE.jpg')[:2048,:2048,:]/256
#hf = plt.imread('cutout2.jpg')

f = np.zeros((res,res,3))
g = np.zeros((res,res,3))
yshape1 = int(np.shape(hf)[0]/2) #round down
yshape2 = int(math.ceil(np.shape(hf)[0]/2)) #round up
xshape1 = int(np.shape(hf)[1]/2)
xshape2 = int(math.ceil(np.shape(hf)[1]/2))
x0 = 0
y0 = 0
xdisp = res//2 + x0
ydisp = res//2 + y0

f[ydisp-yshape1:ydisp+yshape2,xdisp-xshape1:xdisp+xshape2,0] = hf[:,:,0]
f[ydisp-yshape1:ydisp+yshape2,xdisp-xshape1:xdisp+xshape2,1] = hf[:,:,1]
f[ydisp-yshape1:ydisp+yshape2,xdisp-xshape1:xdisp+xshape2,2] = hf[:,:,2]

ff = np.empty((7,res,res,3),float)
ff[0][924:1124,924:1124,:] = f[924:1124,924:1124,:]
ff[1][980:1070,210:300,:] = f[980:1070,210:300,:]
ff[2][470:510,500:550,:] = f[470:510,500:550,:]
ff[3][620:690,630:700,:] = f[620:690,630:700,:]
ff[4][620:700,470:550,:] = f[620:700,470:550,:]
ff[5][490:570,710:780,:] = f[490:570,710:780,:]
ff[6][730:800,400:450,:] = f[730:800,400:450,:]

plt.figure(frameon=False)
plt.axis('off')
plt.title('source image')
plt.imshow(f,origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
plt.pause(1)
plt.clf()

for i in range (0,len(ff)):
    plt.figure(frameon=False)
    plt.axis('off')
    plt.title('source image')
    plt.imshow(ff[i],origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
    plt.pause(1)
    plt.clf()

'''for i in range(0,len(z)):
    asrc = Distances.scalefactor(z[i])
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL,asrc)
    Dd = Distances.angular(a0,aL)
    critdens = (c*Ds)/(4*np.pi*G*Dd*Dds) #kg/m^2
    x,y = X-(Dds/Ds)*Xgrad/critdens, Y-(Dds/Ds)*Ygrad/critdens  
    
    spline0 = RectBivariateSpline(u,u,f[:,:,0].T)
    spline1 = RectBivariateSpline(u,u,f[:,:,1].T)
    spline2 = RectBivariateSpline(u,u,f[:,:,2].T)
    
    g[:,:,0] = spline0.ev(x,y)
    g[:,:,1] = spline1.ev(x,y)
    g[:,:,2] = spline2.ev(x,y)
    plt.figure(frameon=False)
    plt.axis('off')
    plt.title('source image')
    plt.imshow(f,origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
    plt.pause(1)
    plt.clf()
    plt.figure(frameon=False)
    plt.axis('off')
    #plt.title('Lensed Images')
    #plt.xlabel('x in microradians')
    #plt.ylabel('y in microradians')
    plt.imshow(g,origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
    plt.savefig('LensedHDF.jpg',dpi = 1500,bbox_inches='tight',pad_inches = 0)
    plt.pause(2)
    plt.clf()
        
    tnodim = ((X-x0)**2 + (Y-y0)**2)/2 + potential/critdens*(Dds/Ds)
    tau = (1+zL)*Ds*Dd/Dds*tnodim #10^-9 s due to unit microrad
    tau /= 1e12  #arrival time surface in s                           
    tmin = np.min(tau)
    tmax = np.max(tau)
    levs = np.linspace(tmin,tmin + (tmax-tmin)/5,50)
    plt.gca().set_aspect('equal')
    plt.contour(Y,X,tau,levels=levs)
    plt.grid()
    plt.title('Arrival Time Surface')
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.pause(0.1)
    plt.clf()'''