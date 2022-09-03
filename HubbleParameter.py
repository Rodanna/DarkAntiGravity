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
aL = np.linspace(0.5,0.8,4)
zL = 1/aL-1
rmax = 150
res = 500
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
kmtoMpc = 3.24078e-20
z = 6 #np.linspace(6,13,7)
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

Xgrad = -np.loadtxt('Xgrad2_500.txt')
Ygrad = -np.loadtxt('Ygrad2_500.txt')
potential = np.loadtxt('potential2_500.txt')

f = plt.imread('HUBBLE.jpg')
f = f[:res,:res,0]
index1 = np.array([3,7])
index2 = np.array([3,7])
timedelay = []


run1,run2,run3 = 0,0,0
for i in range(0,len(zL)):
    asrc = Distances.scalefactor(z)
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL[i],asrc)
    Dd = Distances.angular(a0,aL[i])
    critdens = (c*Ds)/(4*np.pi*G*Dd*Dds) #kg/m^2
    x,y = X-(Dds/Ds)*Xgrad/critdens, Y-(Dds/Ds)*Ygrad/critdens  
    for k in index1:
        for t in index2:
            f = f*0
            f[55+k,55+t] = 1
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
            tau = (1+zL[i])*Ds*Dd/Dd*tnodim #10^-9 s due to unit microrad
            tau /= 1e12  #arrival time surface in s
            tmin = np.min(tau)
            tmax = np.max(tau)
            levs = np.linspace(tmin,tmin + (tmax-tmin)/5,100)
            plt.gca().set_aspect('equal')
            if k == 3 and t == 3:
                if aL[i] == 0.5:
                    plt.plot(-20,-15, 'b+')
                    plt.plot(3,4,'b+')
                    plt.plot(15,-5,'b+')
                    print('t1:',tau[-20+60,-15+60],'t2:',tau[8+60,8+60])
                    timedelay.append((zL[run1],tau[-20+60,-15+60],tau[8+60,8+60]))
            if k == 7 and t == 7:
                if aL[i] == 0.5:
                    plt.plot(-13,-7, 'b+') 
                    plt.plot(17,0,'b+') 
                    plt.plot(-5,-10,'b+')
                    plt.plot(-5,15,'b+')
                    print('t1:',tau[-13+60,-7+60],'t2:',tau[17+60,5+60])
                    timedelay.append((zL[run1],tau[-6+60,-2+60],tau[12+60,20+60]))
            if k == 3 and t == 7:
                if aL[i] == 0.5:
                    plt.plot(-20,0,'b+')
                    plt.plot(3,3,'b+')
                    plt.plot(10,3,'b+')
                    print('t1:',tau[-23+60,16+60],'t2:',tau[3+60,3+60])
                    timedelay.append((zL[run1],tau[-23+60,16+60],tau[3+60,3+60]))
            if k == 7 and t == 3 :
                if aL[i] == 0.5:
                    plt.plot(5,-15,'b+')
                    plt.plot(0,6,'b+')
                    print('t1:',tau[5+60,-15+60],'t2:',tau[-3+60,3+60])
                    timedelay.append((zL[run1],tau[5+60,-15+60],tau[-3+60,3+60]))
                
            plt.contour(Y,X,tau,levels=levs)
            #plt.colorbar(cs)
            plt.title('arrival time surface')
            plt.title('%i %i' % (k,t))
            plt.pause(0.1)
            plt.clf()
            run3 += 1
        run2 += 1
    run1 +=1


for i in range(0,4):
    print(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][1]-timedelay[i][2]))