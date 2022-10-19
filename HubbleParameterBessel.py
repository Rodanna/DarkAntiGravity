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
stoday = 3600*24

x1 = np.zeros((5,5))  
y1 = np.zeros((5,5))
Timedelay = np.zeros((5,4,10))
extremum = np.zeros((5,5,2))

p = 0
for i in range (5,9):
    Timedelay[0][p] = sorted(np.loadtxt(f'Timedelay0.{i}(nr=700).txt')*stoday)
    Timedelay[1][p] = sorted(np.loadtxt(f'Timedelay0.{i}(nr=1000).txt')*stoday)
    Timedelay[2][p] = sorted(np.loadtxt(f'Timedelay0.{i}(nr=1500).txt')*stoday)
    Timedelay[3][p] = sorted(np.loadtxt(f'Timedelay0.{i}(nr=2000).txt')*stoday)
    Timedelay[4][p] = sorted(np.loadtxt(f'Timedelay0.{i}(nr=2238).txt')*stoday)
    p +=1        

    
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)
stoday = 3600*24 #s/day

pot = np.zeros((6,120,120))
Xgrad = np.zeros((6,120,120))
Ygrad = np.zeros((6,120,120))


Xgrad = -np.loadtxt('Xgrad2_120(nr=2238).txt')#4
Ygrad = -np.loadtxt('Ygrad2_120(nr=2238).txt')
pot = np.loadtxt('potential2_120(nr=2238).txt')


index = np.array([7]) #index = np.array([3,7]) =[]
Hubble =[]

factor = 2 #2 for 700 #29 for 1000 #2 for 2000 #2 for 1500 #13 for 2238
resfine = res*factor
resfine2 = int(resfine/2)
ufine = np.linspace(-rmax,rmax,resfine)[40*factor:80*factor] #cuts out central square and refines it
Xfine,Yfine = np.meshgrid(ufine,ufine)

f = np.zeros((resfine,resfine))

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
    
    tnodim = ((X-x0)**2 + (Y-y0)**2)/2 + pot/critdens*(Dds/Ds)
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
    timedelay = []
    for ix in range(1, 40*factor-1):
        for iy in range(1, 40*factor-1):
            if tau[ix, iy] < tau[ix, iy + 1] and tau[ix, iy] < tau[ix, iy - 1] and \
               tau[ix, iy] < tau[ix + 1, iy] and tau[ix, iy] < tau[ix + 1, iy - 1] and \
               tau[ix, iy] < tau[ix + 1, iy + 1] and tau[ix, iy] < tau[ix - 1, iy] and \
               tau[ix, iy] < tau[ix - 1, iy - 1] and tau[ix, iy] < tau[ix - 1, iy + 1]:
               loc_min.append((ix, iy))
               extremum.append((ix,iy))
    for ix in range(1, 40*factor-1):
        for iy in range(1, 40*factor-1):
            if tau[ix, iy] > tau[ix, iy + 1] and tau[ix, iy] > tau[ix, iy - 1] and \
               tau[ix, iy] > tau[ix + 1, iy] and tau[ix, iy] > tau[ix + 1, iy - 1] and \
               tau[ix, iy] > tau[ix + 1, iy + 1] and tau[ix, iy] > tau[ix - 1, iy] and \
               tau[ix, iy] > tau[ix - 1, iy - 1] and tau[ix, iy] > tau[ix - 1, iy + 1]:
               loc_max.append((ix, iy))
               extremum.append((ix ,iy))
    for ix in range(1, 40*factor-1):
        for iy in range(1, 40*factor-1):
            if tau[ix, iy] > tau[ix, iy + 1] and tau[ix, iy] > tau[ix, iy - 1] and \
               tau[ix, iy] < tau[ix + 1, iy] and tau[ix, iy] < tau[ix - 1, iy]:
               loc_sadd.append((ix, iy))
               extremum.append((ix,iy))
    for ix in range(1, 40*factor-1):
        for iy in range(1, 40*factor-1):
            if tau[ix, iy] < tau[ix, iy + 1] and tau[ix, iy] < tau[ix, iy - 1] and \
               tau[ix, iy] > tau[ix + 1, iy] and tau[ix, iy] > tau[ix - 1, iy]:
               loc_sadd.append((ix, iy))
               extremum.append((ix,iy))
    for ix in range(1, 40*factor-1):
        for iy in range(1, 40*factor-1):
            if tau[ix, iy] > tau[ix + 1, iy + 1] and tau[ix, iy] < tau[ix + 1, iy - 1] and \
               tau[ix, iy] > tau[ix - 1, iy - 1] and tau[ix, iy] < tau[ix - 1, iy + 1]:
               loc_sadd.append((ix, iy))
               extremum.append((ix,iy))
    for ix in range(1, 40*factor-1):
        for iy in range(1, 40*factor-1):
            if tau[ix, iy] < tau[ix + 1, iy + 1] and tau[ix, iy] > tau[ix + 1, iy - 1] and \
               tau[ix, iy] < tau[ix - 1, iy - 1] and tau[ix, iy] > tau[ix - 1, iy + 1]:
               loc_sadd.append((ix, iy))
               extremum.append((ix,iy))    
    extremum = sorted(list(set(extremum)))
    loc_sadd = list(set(loc_sadd))
    print('minima:',loc_min,'maxima:',loc_max,'saddle:',loc_sadd)  
    
    for kik in range (0,len(extremum)):
        y1[kik] = int(extremum[kik][0])
        x1[kik] = int(extremum[kik][1])
    hah = 0
    for m in range (0,len(extremum)):
        plt.plot(Xfine[0,int(x1[m])],Yfine[int(y1[m]),0],'r+',markersize = 3)
        for n in range(m+1,len(extremum)):
            timedelay.append(p.abs(tau[int(y1[m]),int(x1[m])]-tau[int(y1[n]),int(x1[n])]))
    for m in range (0,len(extremum)):
        plt.plot(Xfine[0,int(x1[m])],Yfine[int(y1[m]),0],'r+',markersize = 3)
        for n in range(m+1,len(extremum)):
            timedelay = np.abs(tau[int(y1[m]),int(x1[m])]-tau[int(y1[n]),int(x1[n])])
            timedelay2 = Timedelay[4][i][hah]
            print(timedelay2,timedelay)
            H = timedelay/timedelay2
            print(H)
            Hubble.append(H)
            hah += 1


    cs = plt.gca().set_aspect('equal')
    plt.contour(Yfine,Xfine,tau.T,levels=levs,origin='lower',extent=[-rmax/6,rmax/6,-rmax/6,rmax/6])
    plt.colorbar(cs)
    plt.grid()
    plt.title(r'lens at $z_L$'f' = {zL[i]}')
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    #plt.savefig(f'timedelay{i}(nr=2238)')
    plt.pause(0.1)
    plt.clf()


HubbleA = np.asarray(Hubble)
plt.figure()
plt.title('Hubble Parameter in ')
plt.grid()
plt.hist(HubbleA,bins=len(HubbleA),rwidth = 0.6)
#plt.savefig('HubbleParameter.png')
plt.show()
#np.savetxt('HubbleBessel(nr=2238).txt', HubbleA)

