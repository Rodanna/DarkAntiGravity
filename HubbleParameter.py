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
index = np.array([7]) #index = np.array([3,7])
timedelay = []


run1,run2,run3 = 0,0,0
for i in range(0,len(zL)):
    aL[i] = round(aL[i],1)
    asrc = Distances.scalefactor(z)
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL[i],asrc)
    Dd = Distances.angular(a0,aL[i])
    critdens = (c*Ds)/(4*np.pi*G*Dd*Dds) #kg/m^2
    x,y = X-(Dds/Ds)*Xgrad/critdens, Y-(Dds/Ds)*Ygrad/critdens  
    for k in index:
        for t in index:
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
                    plt.plot(-15,-10, 'b+')
                    plt.plot(3,5,'b+')
                    print('t1:',tau[-20+60,-15+60],'t2:',tau[8+60,8+60])
                    timedelay.append((zL[run1],tau[-20+60,-15+60],tau[8+60,8+60]))
                if aL[i] == 0.6:
                    plt.plot(-19,-12, 'b+')
                    plt.plot(2,6,'b+')
                    print('t1:',tau[-20+60,-15+60],'t2:',tau[8+60,8+60])
                    timedelay.append((zL[run1],tau[-20+60,-15+60],tau[8+60,8+60]))
                if aL[i] == 0.7:
                    plt.plot(-20,-12, 'b+')
                    plt.plot(1,7,'b+')
                    print('t1:',tau[-20+60,-15+60],'t2:',tau[8+60,8+60])
                    timedelay.append((zL[run1],tau[-20+60,-15+60],tau[8+60,8+60]))
                if aL[i] == 0.8:
                    plt.plot(-16,-11, 'b+')
                    plt.plot(2,6,'b+')
                    print('t1:',tau[-20+60,-15+60],'t2:',tau[8+60,8+60])
                    timedelay.append((zL[run1],tau[-20+60,-15+60],tau[8+60,8+60]))
            if k == 7 and t == 7:
                if aL[i] == 0.5:
                    plt.plot(-8,-5,'b+') 
                    plt.plot(14,5,'r+') 
                    plt.plot(1,-5,'r+')
                    plt.plot(-3,7,'g+')
                    print('t1:',tau[60-8,60-5],'t2:',tau[60+14,60+5],'t3:',tau[60+1,60-5],'t4:',tau[60-3,60+7])
                    timedelay.append((zL[run1],tau[60-8,60-5],tau[60+14,60+5],tau[60+1,60-5],tau[60-3,60+7]))
                if aL[i] == 0.6:
                    plt.plot(-13,-5,'b+') 
                    plt.plot(17,6,'r+') 
                    plt.plot(0,-9,'r+')
                    plt.plot(-4,10,'g+')
                    print('t1:',tau[60-13,60-5],'t2:',tau[60+17,60+6],'t3:',tau[60,60-9],'t4:',tau[60-4,60+10])
                    timedelay.append((zL[run1],tau[60-13,60-5],tau[60+17,60+6],tau[60,60-9],tau[60-4,60+10]))
                if aL[i] == 0.7:
                    plt.plot(-14,-5, 'b+') 
                    plt.plot(18,7,'r+') 
                    plt.plot(0,-11,'r+')
                    plt.plot(-4,11,'g+')
                    print('t1:',tau[60-14,60-5],'t2:',tau[60+18,60+7],'t3:',tau[60,60-11],'t4:',tau[60-4,60+11])
                    timedelay.append((zL[run1],tau[60-14,60-5],tau[60+18,60+7],tau[60,60-11],tau[60-4,60+11]))
                if aL[i] == 0.8:
                    plt.plot(-10,-4,'b+') 
                    plt.plot(15,6,'r+') 
                    plt.plot(0,-6,'r+')
                    plt.plot(-3,8,'g+')
                    print('t1:',tau[60-10,60-4],'t2:',tau[60+15,60+6],'t3:',tau[60,60-6],'t4:',tau[60-3,60+8])
                    timedelay.append((zL[run1],tau[60-10,60-4],tau[60+15,60+6],tau[60,60-6],tau[60-3,60+8]))
            if k == 3 and t == 7:
                if aL[i] == 0.5:
                    plt.plot(-15,-3,'b+')
                    plt.plot(3,3,'b+')
                    plt.plot(8,5,'b+')
                    print('t1:',tau[-23+60,16+60],'t2:',tau[3+60,3+60])
                    timedelay.append((zL[run1],tau[-23+60,16+60],tau[3+60,3+60]))
                if aL[i] == 0.6:
                    plt.plot(-20,-4, 'b+') 
                    plt.plot(15,6,'b+') 
                    plt.plot(2,3,'b+')
                    plt.plot(1,14,'b+')
                    print('t1:',tau[-13+60,-7+60],'t2:',tau[17+60,5+60])
                    timedelay.append((zL[run1],tau[-6+60,-2+60],tau[12+60,20+60]))
                if aL[i] == 0.7:
                    plt.plot(-26,-5, 'b+') 
                    plt.plot(15,6,'b+') 
                    plt.plot(1,-5,'b+')
                    plt.plot(1,10,'b+')
                    print('t1:',tau[-13+60,-7+60],'t2:',tau[17+60,5+60])
                    timedelay.append((zL[run1],tau[-6+60,-2+60],tau[12+60,20+60]))
                if aL[i] == 0.8:
                    plt.plot(-9,-5, 'b+') 
                    plt.plot(15,6,'b+') 
                    plt.plot(1,-5,'b+')
                    plt.plot(-2,5,'b+')
                    print('t1:',tau[-13+60,-7+60],'t2:',tau[17+60,5+60])
                    timedelay.append((zL[run1],tau[-6+60,-2+60],tau[12+60,20+60]))
            if k == 7 and t == 3 :
                if aL[i] == 0.5:
                    plt.plot(2,-13,'b+')
                    plt.plot(0,6,'b+')
                    print('t1:',tau[5+60,-15+60],'t2:',tau[-3+60,3+60])
                    timedelay.append((zL[run1],tau[5+60,-15+60],tau[-3+60,3+60]))
                if aL[i] == 0.6:
                    plt.plot(-9,-15, 'b+') 
                    plt.plot(15,-10,'b+') 
                    plt.plot(4,-20,'b+')
                    plt.plot(-2,7,'b+')
                    print('t1:',tau[-13+60,-7+60],'t2:',tau[17+60,5+60])
                    timedelay.append((zL[run1],tau[-6+60,-2+60],tau[12+60,20+60]))
                if aL[i] == 0.7:
                    plt.plot(-9,-10, 'b+') 
                    plt.plot(15,6,'b+') 
                    plt.plot(1,-15,'b+')
                    plt.plot(3,5,'b+')
                    print('t1:',tau[-13+60,-7+60],'t2:',tau[17+60,5+60])
                    timedelay.append((zL[run1],tau[-6+60,-2+60],tau[12+60,20+60]))
                if aL[i] == 0.8:
                    plt.plot(-9,-5, 'b+') 
                    plt.plot(15,6,'b+') 
                    plt.plot(1,-5,'b+')
                    plt.plot(-2,5,'b+')
                    print('t1:',tau[-13+60,-7+60],'t2:',tau[17+60,5+60])
                    timedelay.append((zL[run1],tau[-6+60,-2+60],tau[12+60,20+60]))
                
            plt.contour(Y,X,tau,levels=levs)
            #plt.colorbar(cs)
            plt.grid()
            plt.title('arrival time surface')
            plt.title('%f %i %i' % (aL[i],k,t))
            plt.pause(0.1)
            plt.clf()
            run3 += 1
        run2 += 1
    run1 +=1

Hubble =[]
for i in range(0,4):
    #print(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][1]-timedelay[i][2]))
    Hubble.append(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][1]-timedelay[i][2]))
    #print(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][1]-timedelay[i][3]))
    Hubble.append(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][1]-timedelay[i][3]))
    #print(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][1]-timedelay[i][4]))
    Hubble.append(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][1]-timedelay[i][4]))
    #print(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][2]-timedelay[i][3]))
    Hubble.append(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][2]-timedelay[i][3]))
    #print(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][2]-timedelay[i][4]))
    Hubble.append(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][2]-timedelay[i][4]))
    #print(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][3]-timedelay[i][4]))
    Hubble.append(kmtoMpc*timedelay[i][0]/np.abs(timedelay[i][3]-timedelay[i][4]))
    
HubbleA = np.asarray(Hubble)

plt.figure()
plt.title('Hubble Parameter in ')
plt.grid()
plt.hist(HubbleA, bins = 20,rwidth = 0.1)
plt.show()