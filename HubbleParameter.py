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
rmax = 150
res = 120
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
z = np.linspace(9,13,8)
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

Xgrad = -np.loadtxt('Xgrad2_120.txt')
Ygrad = -np.loadtxt('Ygrad2_120.txt')
potential = np.loadtxt('potential2_120.txt')

f = plt.imread('HUBBLE.jpg')
f = f[:120,:120,0]
index1 = np.array([0,0,7,7])
index2 = np.array([0,7,0,7])
timedelay = np.zeros(3)

for i in range(0,1):
    asrc = Distances.scalefactor(z[i])
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL,asrc)
    Dd = Distances.angular(a0,aL)
    critdens = (c*Ds)/(4*np.pi*G*Dd*Dds) #kg/m^2
    x,y = X-(Dds/Ds)*Xgrad/critdens, Y-(Dds/Ds)*Ygrad/critdens  
    
    for k in index1:
        for t in index2:
            f = f*0
            f[55+k,55+t] = 1
            plt.title('source image')
            plt.imshow(f.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
            plt.pause(0.1)
            plt.clf()
            x0 = -5 + t
            y0 = -5 + k
            rad2 = (X-x0)**2 + (Y-y0)**2
            f = np.exp(-rad2/5)
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
            
            tnodim = ((X-x0)**2 + (Y-y0)**2)/2 + potential/critdens*(Dds/Ds)
            tau = (1+zL)*Ds*Dd/Dd*tnodim #10^-9 s due to unit microrad
            tau /= 1e12  #arrival time surface in s
            tmin = np.min(tau)
            tmax = np.max(tau)
            levs = np.linspace(tmin,tmin + (tmax-tmin)/5,100)
            plt.gca().set_aspect('equal')
            if k == 0 and t == 0:
                plt.plot(-20,-15, 'b+')
                plt.plot(8,8,'b+')
                print('t1:',tau[-20+60,-15+60],'t2:',tau[8+60,8+60])
            if k == 7 and t == 7:
                plt.plot(-6,-2, 'b+') 
                plt.plot(12,20,'b+') 
                print('t1:',tau[-6+60,-2+60],'t2:',tau[12+60,20+60])
            if k == 0 and t == 7:
                plt.plot(-23,16,'b+')
                plt.plot(3,3,'b+')
            if t == 0 and k == 7:
                plt.plot(5,-15,'b+')
                plt.plot(-3,3,'b+')
            plt.contour(Y,X,tau,levels=levs)
            #plt.colorbar(cs)
            plt.title('arrival time surface')
            plt.title('%i %i' % (t,k))
            plt.pause(0.1)
            plt.clf()
            
            
t1 = np.array([530846213.63775915,556582512.2111806])
t2 = np.array([505811182.8377406,468241022.4749451])
z = np.array([9,9])

print(z[0]/(t1[0]-t2[0]))
print(z[1]/(t1[1]-t2[1]))