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
res = 500
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
stoday = 3600*24 #s/day


Xgrad = -np.loadtxt('Xgrad2_500.txt')
Ygrad = -np.loadtxt('Ygrad2_500.txt')
potential = np.loadtxt('potential2_500.txt')
noisypotential1 = np.loadtxt('noisypotential2_1000.txt')
noisypotential2 = np.loadtxt('noisypotential2_1000(2).txt')
noisypotential3 = np.loadtxt('noisypotential2_1000(3).txt')
noisypotential4 = np.loadtxt('noisypotential2_1000(4).txt')
noisypotential5 = np.loadtxt('noisypotential2_1000(5).txt')
noisypotential = np.array([noisypotential1,noisypotential2,noisypotential3,noisypotential4,noisypotential5])

Timedelay = np.zeros((4,10))
NoisyTimedelay = np.zeros((5,4,10))

p = 0
for i in range(5,9):
    Timedelay[p] = np.loadtxt(f'Timedelay0.{i}(nr=2238).txt')
    p += 1

factor = 3 #2 for 700 #29 for 1000 #2 for 2000 #2 for 1500 #13 for 2238
resfine = res*factor
resfine2 = int(resfine/2)
ufine = np.linspace(-rmax,rmax,resfine)[190*factor:310*factor] #cuts out central square and refines it
Xfine,Yfine = np.meshgrid(ufine,ufine)

f = np.zeros((resfine,resfine))
noisypot = np.zeros((5,360,360))

for p in range(0,5):
    noisypotspline = RectBivariateSpline(u, u, noisypotential[p])
    noisypot[p] = noisypotspline.ev(Xfine,Yfine)

potspline = RectBivariateSpline(u, u, potential)
potential = potspline.ev(Xfine,Yfine)
    


index = np.array([7]) #index = np.array([3,7]) =[]
H = []
timedelay = []
result = []

NoisyTimedelay5 = []
NoisyTimedelay6 = []
NoisyTimedelay7 = []
NoisyTimedelay8 = []

for p in range (0,5):
    for i in range(0,len(zL)):
        aL[i] = round(aL[i],1)
        asrc = Distances.scalefactor(zS)
        Ds = Distances.angular(a0,asrc)
        Dds = Distances.angular(aL[i],asrc)
        Dd = Distances.angular(a0,aL[i])
        critdens = (c*Ds)/(4*np.pi*G*Dd*Dds) #kg/m^2
        x,y = X-(Dds/Ds)*Xgrad/critdens, Y-(Dds/Ds)*Ygrad/critdens
            
        k = 7
        t = 7
        x0 = -5 + t
        y0 = -5 + k
        rad2 = (Xfine-x0)**2 + (Yfine-y0)**2
        f = np.exp(-rad2/5)
        spline = RectBivariateSpline(ufine,ufine,f.T)
        g = spline.ev(x,y)
        
        
        plt.title('source image at z = 6')
        plt.imshow(f.T,origin='lower',extent=[-rmax/6,rmax/6,-rmax/6,rmax/6])
        plt.xlabel('x in microradians')
        plt.ylabel('y in microradians')
        #plt.savefig('sourceat6')
        plt.pause(0.1)
        plt.clf()
        plt.title(f'lensed image for lens at z = {zL[i]}')
        plt.imshow(g.T[40:80,40:80],origin='lower',extent=[-rmax/6,rmax/6,-rmax/6,rmax/6])
        plt.xlabel('x in microradians')
        plt.ylabel('y in microradians')
        #plt.savefig(f'lensedimages{i}')
        plt.pause(0.1)
        plt.clf()
        
        tnodim = ((X-x0)**2 + (Y-y0)**2)/2 + noisypotential[p]/critdens*(Dds/Ds)
        tau = (1+zL[i])*Ds*Dd/Dds*tnodim #10^-9 s due to unit microrad
        tau /= 1e12  #arrival time surface in s     
        tspline = RectBivariateSpline(u,u,tau)
        tau = tspline.ev(Xfine,Yfine)
                     
        tmin = np.min(tau)
        tmax = np.max(tau)
        levs = np.linspace(tmin,tmin + (tmax-tmin)/5,40)
        plt.gca().set_aspect('equal')
        
        x1 = np.zeros(10)
        y1 = np.zeros(10)
        loc_min = []
        loc_max = []
        loc_sadd = []
        extremum = []
        for ix in range(1, 120*factor-1):
            for iy in range(1, 120*factor-1):
                if tau[ix, iy] < tau[ix, iy + 1] and tau[ix, iy] < tau[ix, iy - 1] and \
                   tau[ix, iy] < tau[ix + 1, iy] and tau[ix, iy] < tau[ix + 1, iy - 1] and \
                   tau[ix, iy] < tau[ix + 1, iy + 1] and tau[ix, iy] < tau[ix - 1, iy] and \
                   tau[ix, iy] < tau[ix - 1, iy - 1] and tau[ix, iy] < tau[ix - 1, iy + 1]:
                   loc_min.append((ix, iy))
                   extremum.append((ix,iy))
        for ix in range(1, 120*factor-1):
            for iy in range(1, 120*factor-1):
                if tau[ix, iy] > tau[ix, iy + 1] and tau[ix, iy] > tau[ix, iy - 1] and \
                   tau[ix, iy] > tau[ix + 1, iy] and tau[ix, iy] > tau[ix + 1, iy - 1] and \
                   tau[ix, iy] > tau[ix + 1, iy + 1] and tau[ix, iy] > tau[ix - 1, iy] and \
                   tau[ix, iy] > tau[ix - 1, iy - 1] and tau[ix, iy] > tau[ix - 1, iy + 1]:
                   loc_max.append((ix, iy))
                   extremum.append((ix ,iy))
        for ix in range(1, 120*factor-1):
            for iy in range(1, 120*factor-1):
                if tau[ix, iy] > tau[ix, iy + 1] and tau[ix, iy] > tau[ix, iy - 1] and \
                   tau[ix, iy] < tau[ix + 1, iy] and tau[ix, iy] < tau[ix - 1, iy]:
                   loc_sadd.append((ix, iy))
                   extremum.append((ix,iy))
        for ix in range(1, 120*factor-1):
            for iy in range(1, 120*factor-1):
                if tau[ix, iy] < tau[ix, iy + 1] and tau[ix, iy] < tau[ix, iy - 1] and \
                   tau[ix, iy] > tau[ix + 1, iy] and tau[ix, iy] > tau[ix - 1, iy]:
                   loc_sadd.append((ix, iy))
                   extremum.append((ix,iy))
        for ix in range(1, 120*factor-1):
            for iy in range(1, 120*factor-1):
                if tau[ix, iy] > tau[ix + 1, iy + 1] and tau[ix, iy] < tau[ix + 1, iy - 1] and \
                   tau[ix, iy] > tau[ix - 1, iy - 1] and tau[ix, iy] < tau[ix - 1, iy + 1]:
                   loc_sadd.append((ix, iy))
                   extremum.append((ix,iy))
        for ix in range(1, 120*factor-1):
            for iy in range(1, 120*factor-1):
                if tau[ix, iy] < tau[ix + 1, iy + 1] and tau[ix, iy] > tau[ix + 1, iy - 1] and \
                   tau[ix, iy] < tau[ix - 1, iy - 1] and tau[ix, iy] > tau[ix - 1, iy + 1]:
                   loc_sadd.append((ix, iy))
                   extremum.append((ix,iy))    
        extremum = list(set(extremum))
        loc_sadd = list(set(loc_sadd))
        print('minima:',loc_min,'maxima:',loc_max,'saddle:',loc_sadd)  
        
        for kik in range (0,len(extremum)):
            y1[kik] = int(extremum[kik][0])
            x1[kik] = int(extremum[kik][1])
        for m in range (0,len(extremum)):
            plt.plot(Xfine[0,int(x1[m])],Yfine[int(y1[m]),0],'r+',markersize = 3)
            #print(tau[int(y1[m])-1:int(y1[m]+2),int(x1[m])-1:int(x1[m])+2])
            for n in range(m+1,len(extremum)):
                timedelay = np.abs(tau[int(y1[m]),int(x1[m])]-tau[int(y1[n]),int(x1[n])])    
                if aL[i] == 0.5:         
                        NoisyTimedelay5.append(timedelay/stoday)                
                if aL[i] == 0.6:
                        NoisyTimedelay6.append(timedelay/stoday)                      
                if aL[i] == 0.7:
                        NoisyTimedelay7.append(timedelay/stoday)            
                if aL[i] == 0.8:
                        NoisyTimedelay8.append(timedelay/stoday)
        
    
        cs = plt.gca().set_aspect('equal')
        plt.contour(Yfine,Xfine,tau.T,levels=levs,origin='lower',extent=[-rmax,rmax,-rmax,rmax])
        plt.colorbar(cs)
        plt.grid()
        plt.title(r'lens at $z_L$'f' = {zL[i]}')
        plt.xlabel('x in microradians')
        plt.ylabel('y in microradians')
        plt.pause(0.1)
        plt.clf()

for i in range (1,6):
    NoisyTimedelay[i-1][0] = sorted(np.asarray(NoisyTimedelay5[(i-1)*10:i*10]))
    NoisyTimedelay[i-1][1] = sorted(np.asarray(NoisyTimedelay6[(i-1)*10:i*10]))
    NoisyTimedelay[i-1][2] = sorted(np.asarray(NoisyTimedelay7[(i-1)*10:i*10]))
    NoisyTimedelay[i-1][3] = sorted(np.asarray(NoisyTimedelay8[(i-1)*10:i*10]))

Hubble =[]


for i in range(0,5): 
    for j in range (0,4):
        for k in range (0,10):
            timedelay1 = NoisyTimedelay[i][j][k]
            timedelay2 = Timedelay[j][k]
            Hubble.append(timedelay1/timedelay2) 

Hubble = np.asarray(Hubble)

np.savetxt('Hubble1000.txt',Hubble)
'''
plt.figure()
plt.title('Hubble Parameter in ')
plt.grid()
plt.hist(Hubble, bins = 200,rwidth = 0.9)
plt.show()'''
