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
X,Y = np.meshgrid(u,u)
stoday = 3600*24 #s/day

pot = np.zeros((6,120,120))
Xgrad = np.zeros((6,120,120))
Ygrad = np.zeros((6,120,120))


Xgrad = -np.loadtxt('Xgrad2_120(nr=500).txt')#4
Ygrad = -np.loadtxt('Ygrad2_120(nr=500).txt')

pot = np.loadtxt('potential2_120(nr=500).txt')

'''
Xgrad= -np.loadtxt('Xgrad2_120(nr=2238).txt')#5
Ygrad = -np.loadtxt('Ygrad2_120(nr=2238).txt')
pot = np.loadtxt('potential2_120(nr=2238).txt')
'''


index = np.array([7]) #index = np.array([3,7]) =[]
Timedelay5 = []
Timedelay6 = []
Timedelay7 = []
Timedelay8 = []


factor = 31 #29 for 1000 #10 for 2000 #20 for 1500 #13 for 2238
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
    plt.savefig('sourceat6')
    plt.pause(0.1)
    plt.clf()
    plt.title(f'lensed image for lens at z = {zL[i]}')
    plt.imshow(g.T[40:80,40:80],origin='lower',extent=[-rmax/6,rmax/6,-rmax/6,rmax/6])
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.savefig(f'lensedimages{i}')
    plt.pause(0.1)
    plt.clf()
    
    tnodim = ((X-x0)**2 + (Y-y0)**2)/2 + pot/critdens*(Dds/Ds)
    tau = (1+zL[i])*Ds*Dd/Dds*tnodim #10^-9 s due to unit microrad
    tau /= 1e12  #arrival time surface in s     
    tspline = RectBivariateSpline(u, u, tau)
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
    for ix in range(1, 40*factor-1):
        for iy in range(1, 40*factor-1):
            if tau[ix, iy] <= tau[ix, iy + 1] and tau[ix, iy] <= tau[ix, iy - 1] and \
               tau[ix, iy] <= tau[ix + 1, iy] and tau[ix, iy] <= tau[ix + 1, iy - 1] and \
               tau[ix, iy] <= tau[ix + 1, iy + 1] and tau[ix, iy] <= tau[ix - 1, iy] and \
               tau[ix, iy] <= tau[ix - 1, iy - 1] and tau[ix, iy] <= tau[ix - 1, iy + 1]:
               loc_min.append((int(ix/factor)+50, int(iy/factor)+50))
               extremum.append((ix,iy))
    for ix in range(1, 40*factor-1):
        for iy in range(1, 40*factor-1):
            if tau[ix, iy] >= tau[ix, iy + 1] and tau[ix, iy] >= tau[ix, iy - 1] and \
               tau[ix, iy] >= tau[ix + 1, iy] and tau[ix, iy] >= tau[ix + 1, iy - 1] and \
               tau[ix, iy] >= tau[ix + 1, iy + 1] and tau[ix, iy] >= tau[ix - 1, iy] and \
               tau[ix, iy] >= tau[ix - 1, iy - 1] and tau[ix, iy] >= tau[ix - 1, iy + 1]:
               loc_max.append((int(ix/factor)+50, int(iy/factor)+50))
               extremum.append((ix,iy))
    for ix in range(1, 40*factor-1):
        for iy in range(1, 40*factor-1):
            if tau[ix, iy] >= tau[ix, iy + 1] and tau[ix, iy] >= tau[ix, iy - 1] and \
               tau[ix, iy] <= tau[ix + 1, iy] and tau[ix, iy] <= tau[ix - 1, iy]:
               loc_sadd.append((int(ix/factor)+50, int(iy/factor)+50))
               extremum.append((ix,iy))
    for ix in range(1, 40*factor-1):
        for iy in range(1, 40*factor-1):
            if tau[ix, iy] <= tau[ix, iy + 1] and tau[ix, iy] <= tau[ix, iy - 1] and \
               tau[ix, iy] >= tau[ix + 1, iy] and tau[ix, iy] >= tau[ix - 1, iy]:
               loc_sadd.append((int(ix/factor)+50, int(iy/factor)+50))
               extremum.append((ix,iy))
    print('minima:',loc_min,'maxima:',loc_max,'saddle:',loc_sadd)  
    
    for kik in range (0,len(extremum)):
        y1[kik] = int(extremum[kik][0])
        x1[kik] = int(extremum[kik][1])
    for m in range (0,len(extremum)):
        plt.plot(Xfine[0,int(x1[m])],Yfine[int(y1[m]),0],'r+',markersize = 3)
        for n in range(m+1,len(extremum)):
            timedelay = np.abs(tau[int(x1[m]),int(y1[m])]-tau[int(y1[n]),int(x1[n])])    
            if aL[i] == 0.5:         
                    Timedelay5.append(timedelay)                
            if aL[i] == 0.6:
                    Timedelay6.append(timedelay)                      
            if aL[i] == 0.7:
                    Timedelay7.append(timedelay)            
            if aL[i] == 0.8:
                    Timedelay8.append(timedelay)

            
    #cs = 
    plt.gca().set_aspect('equal')
    plt.contour(Yfine,Xfine,tau.T,levels=levs,origin='lower',extent=[-rmax/6,rmax/6,-rmax/6,rmax/6])
    #plt.colorbar(cs)
    plt.grid()
    plt.title(f'arrival time surface for lens at z = {zL[i]}')
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.savefig(f'timedelay{i}(nr=500)')
    plt.pause(0.1)
    plt.clf()

    
Timedelay5 = np.asarray(Timedelay5)/stoday
Timedelay6 = np.asarray(Timedelay6)/stoday
Timedelay7 = np.asarray(Timedelay7)/stoday
Timedelay8 = np.asarray(Timedelay8)/stoday
np.savetxt('Timedelay0.5(nr=500).txt',Timedelay5)
np.savetxt('Timedelay0.6(nr=500).txt',Timedelay6)
np.savetxt('Timedelay0.7(nr=500).txt',Timedelay7)
np.savetxt('Timedelay0.8(nr=500).txt',Timedelay8)
#print(Timedelay5,Timedelay6,Timedelay7,Timedelay8)