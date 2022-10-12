#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 13:27:27 2022

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
res = 120
res2 = int(res/2)
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
pc = 1/np.tan(4.848*10**(-6))*1.495978707*10**11/c # cs
Mpc = pc*1e6 # cs
invcritdens0 = 0.49823288405658067 #(kg/m^2)^-1
H0 = 67.6 # km/(s Mpc)
kmtoMpc = 3.24078e-20

u = np.linspace(-rmax,rmax,res)
ufine = np.linspace(-rmax,rmax,res*10)
Xfine,Yfine = np.meshgrid(ufine,ufine)
X,Y = np.meshgrid(u,u)
stoday = 3600*24 #s/day

pot = np.zeros((6,120,120))
Xgrad = np.zeros((6,120,120))
Ygrad = np.zeros((6,120,120))


'''Xgrad = -np.loadtxt('Xgrad2_120(nr=150).txt')#0
Ygrad= -np.loadtxt('Ygrad2_120(nr=150).txt')
pot = np.loadtxt('potential2_120(nr=150).txt')

Xgrad = -np.loadtxt('Xgrad2_120(nr=500).txt')#1
Ygrad = -np.loadtxt('Ygrad2_120(nr=500).txt')
pot = np.loadtxt('potential2_120(nr=500).txt')

Xgrad = -np.loadtxt('Xgrad2_120(nr=1000).txt')#2
Ygrad = -np.loadtxt('Ygrad2_120(nr=1000).txt')
pot = np.loadtxt('potential2_120(nr=1000).txt')

Xgrad = -np.loadtxt('Xgrad2_120(nr=1500).txt')#3
Ygrad= -np.loadtxt('Ygrad2_120(nr=1500).txt')
pot = np.loadtxt('potential2_120(nr=1500).txt')

Xgrad = -np.loadtxt('Xgrad2_120(nr=100).txt')#4
Ygrad = -np.loadtxt('Ygrad2_120(nr=100).txt')
pot = np.loadtxt('potential2_120(nr=100).txt')
'''
Xgrad= -np.loadtxt('Xgrad2_120(nr=2238).txt')#5
Ygrad = -np.loadtxt('Ygrad2_120(nr=2238).txt')
pot = np.loadtxt('potential2_120(nr=2238).txt')



f = plt.imread('HUBBLE.jpg')
f = f[:res,:res,0]
index = np.array([7]) #index = np.array([3,7]) =[]
Timedelay5 = []
Timedelay6 = []
Timedelay7 = []
Timedelay8 = []
k = 7
t = 7
x0 = -5 +t
y0 = -5 +k

splineX = RectBivariateSpline(u,u,Xgrad)
splineY = RectBivariateSpline(u,u,Ygrad)

gradX = splineX.ev(Xfine,Yfine)
gradY = splineY.ev(Xfine,Yfine)

plt.title('SplineX')
plt.imshow(gradX.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
plt.xlabel('x in microradians')
plt.pause(0.1)
plt.clf()

plt.title('SplineY')
plt.imshow(gradY.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
plt.xlabel('x in microradians')
plt.pause(0.1)
plt.clf()

tol = 0.001
indices = []
for i in range (0,len(gradX)):
    for j in range(0,len(gradY)):
        if np.abs(gradX[i,j]) < tol  :
            print('X:',gradX[i,j],gradY[i,j],i,j)
            indices.append((i/10-res2,j/10-res2))
        if np.abs(gradY[i,j]) < tol:
            print('Y:',gradX[i,j],gradY[i,j],i,j)
            indices.append((i/10-res2,j/10-res2))


for i in range(0,len(zL)):
    aL[i] = round(aL[i],1)
    asrc = Distances.scalefactor(zS)
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL[i],asrc)
    Dd = Distances.angular(a0,aL[i])
    critdens = (c*Ds)/(4*np.pi*G*Dd*Dds) #kg/m^2
    x,y = X-(Dds/Ds)*Xgrad/critdens, Y-(Dds/Ds)*Ygrad/critdens  

    f = f*0
    x0 = -5 + t
    y0 = -5 + k
    rad2 = (X-x0)**2 + (Y-y0)**2
    f = np.exp(-rad2/5)
    spline = RectBivariateSpline(u,u,f.T)
    g = spline.ev(x,y)
    
    
    plt.title('source image at z = 6')
    plt.imshow(f.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.savefig('sourceat6')
    plt.pause(0.1)
    plt.clf()
    plt.title(f'lensed image for lens at z = {zL[i]}')
    plt.imshow(g.T,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.savefig(f'lensedimages{i}')
    plt.pause(0.1)
    plt.clf()
    
    tnodim = ((X-x0)**2 + (Y-y0)**2)/2 + pot/critdens*(Dds/Ds)
    tau = (1+zL[i])*Ds*Dd/Dds*tnodim #10^-9 s due to unit microrad
    tau /= 1e12  #arrival time surface in s                           
    tmin = np.min(tau)
    tmax = np.max(tau)
    levs = np.linspace(tmin,tmin + (tmax-tmin)/5,500)
    plt.gca().set_aspect('equal')

    if aL[i] == 0.5:
        x1 = np.array([-10,14,-2,-3])#4
        y1 = np.array([-4,4,-7,9])
        #x1 = np.array([-10,14,-2,-3])#2 & 3
        #y1 = np.array([-4,4,-7,9])
        #x1 = np.array([-10,14,-2,-3])#1
        #y1 = np.array([-4,4,-7,9])
        #x1 = np.array([-10,14,-3,-4])#0
        #y1 = np.array([-4,4,-5,3])
        #x1 = np.array([-9,14,0,-2])# 5
        #y1 = np.array([-4,4,-7,7]) 
        for m in range (0,len(x1)):
            plt.plot(x1[m],y1[m],'r+',markersize = 3)   
            for n in range(m+1,len(x1)):
                timedelay = np.abs(tau[res2+y1[m],res2+x1[m]]-tau[res2+y1[n],res2+x1[n]])
                Timedelay5.append(timedelay)
            
            
    if aL[i] == 0.6:
        x1 = np.array([-13,17,0,-3])#4
        y1 = np.array([-5,4,-9,10])
        #x1 = np.array([-14,19,0,-3])#2 & 3
        #y1 = np.array([-5,4,-10,10])
        #x1 = np.array([-14,19,0,-3])#1
        #y1 = np.array([-5,4,-10,10])
        #x1 = np.array([-17,19,0,-3])#0
        #y1 = np.array([-5,4,-10,9])
        #x1 = np.array([-13,17,0,-3])#5
        #y1 = np.array([-5,4,-9,10])
        for m in range (0,len(x1)):
            plt.plot(x1[m],y1[m],'r+',markersize = 3)  
            for n in range(m+1,len(x1)):
                timedelay = np.abs(tau[res2+y1[m],res2+x1[m]]-tau[res2+y1[n],res2+x1[n]])
                Timedelay6.append(timedelay)
                
                
    if aL[i] == 0.7:
        x1 = np.array([-14,19,0,-3])#4
        y1 = np.array([-7,4,-11,11])
        #x1 = np.array([-16,19,0,-3])#2 & 3
        #y1 = np.array([-6,5,-11,10])
        #x1 = np.array([-16,19,0,-3])#1
        #y1 = np.array([-6,5,-11,10])
        #x1 = np.array([-16,19,0,-3])#0
        #y1 = np.array([-6,4,-11,10])
        #x1 = np.array([-14,19,0,-3])#5
        #y1 = np.array([-7,4,-11,11])
        for m in range (0,len(x1)):
            plt.plot(x1[m],y1[m],'r+',markersize = 3)  
            for n in range(m+1,len(x1)):
                timedelay = np.abs(tau[res2+y1[m],res2+x1[m]]-tau[res2+y1[n],res2+x1[n]])
                Timedelay7.append(timedelay)
                
                                                    
    if aL[i] == 0.8:
        x1 = np.array([-10,15,0,-3])#4
        y1 = np.array([-4,4,-8,8])
        #x1 = np.array([-10,16,-2,-3])#2 & 3
        #y1 = np.array([-4,4,-8,8])
        #x1 = np.array([-11,16,-2,-3])#1
        #y1 = np.array([-4,4,-8,8])
        #x1 = np.array([-11,16,-2,-3])#0
        #y1 = np.array([-4,4,-6,8])
        #x1 = np.array([-10,15,0,-3])#5
        #y1 = np.array([-4,4,-8,8])
        for m in range (0,len(x1)):
            plt.plot(x1[m],y1[m],'r+',markersize = 3)
            for n in range(m+1,len(x1)):
                timedelay = np.abs(tau[res2+y1[m],res2+x1[m]]-tau[res2+y1[n],res2+x1[n]])
                Timedelay8.append(timedelay)
                
                                    
    #cs = 
    plt.contour(Y,X,tau,levels=levs)
    #plt.colorbar(cs)
    plt.grid()
    plt.title(f'arrival time surface for lens at z = {zL[i]}')
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.savefig(f'timedelay{i}')
    plt.pause(0.1)
    plt.clf()

Timedelay5 = np.asarray(Timedelay5)/stoday
Timedelay6 = np.asarray(Timedelay6)/stoday
Timedelay7 = np.asarray(Timedelay7)/stoday
Timedelay8 = np.asarray(Timedelay8)/stoday

#print(Timedelay5,Timedelay6,Timedelay7,Timedelay8)