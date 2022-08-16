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
z = np.linspace(9,13,8)
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

Xgrad = -np.loadtxt('Xgrad2.txt')
Ygrad = -np.loadtxt('Ygrad2.txt')
potential = np.loadtxt('potential2.txt')

'''
dX = dY = u[1]-u[0]
Xgrad = 0*potential
Ygrad = 0*potential
Xgrad[1:,:] = -(potential[1:,:] - potential[:-1,:])/dX
Ygrad[:,1:] = -(potential[:,1:] - potential[:,:-1])/dY
'''
hf = plt.imread('HUBBLE.jpg')
#f = f[235:355,1034:1154,0]
f = 0*X
f[res//2-60:res//2+60,res//2-60:res//2+60] = hf[1034:1154,235:355,0]
plt.imshow(f)
plt.show()

'''
q = 0.5
C = 1e3
potential = C*np.log(1+X*X+q*Y*Y)
Xgrad = -2*C*X/(1+X*X+q*Y*Y)
Ygrad = -2*C*q*Y/(1+X*X+q*Y*Y)
'''

for i in range(0,len(z)):
    asrc = Distances.scalefactor(z[i])
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL,asrc)
    Dd = Distances.angular(a0,aL)
    critdens = 4*np.pi*G*Dd*Dds/(c*Ds) #(kg/m^2)^-1 
    #critdens *= 3
    x,y = X-(Dds/Ds)*Xgrad*critdens, Y-(Dds/Ds)*Ygrad*critdens  
    
    for k in range(0,10):
        
        for t in range(0,10):
            if not (t==8 and k==7):
                continue
            '''
            f = f*0
            f[55+k,55+t] = 1
            '''
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
            
            #x0 = -5 + t
            #y0 = -5 + k
            tnodim = ((X-x0)**2 + (Y-y0)**2)/2 + potential*critdens*Dds/Ds
            tau = (1+zL)*Ds*Dd/Dd*tnodim #arrival time surface in s
            tau /= 1e12  # because was in micro-radians
            tmin = np.min(tau)
            tmax = np.max(tau)
            levs = np.linspace(tmin,tmin + (tmax-tmin)/5,200)
            cs = plt.contour(Y,X,tau,levels=levs)
            plt.colorbar(cs)
            plt.title('arrival time surface')
            plt.gca().set_aspect('equal')
            plt.title('%i %i' % (t,k))
            plt.pause(2)
            plt.clf()
