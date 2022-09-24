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
res2 = int(res/2)
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
pc = 1/np.tan(4.848*10**(-6))*1.495978707*10**11/c # cs
Mpc = pc*1e6 # cs
H0 = 67.6 # km/(s Mpc)
kmtoMpc = 3.24078e-20
z = 6 #np.linspace(6,13,7)
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

Xgrad = -np.loadtxt('Xgrad2_500.txt')
Ygrad = -np.loadtxt('Ygrad2_500.txt')
potential = np.loadtxt('potential2_500.txt')
noisypotential1 = np.loadtxt('noisypotential2_10.txt')
noisypotential2 = np.loadtxt('noisypotential2_10(2).txt')
noisypotential3 = np.loadtxt('noisypotential2_10(3).txt')
noisypotential4 = np.loadtxt('noisypotential2_10(4).txt')
noisypotential5 = np.loadtxt('noisypotential2_10(5).txt')
noisypotential = np.array([noisypotential1,noisypotential2,noisypotential3,noisypotential4,noisypotential5])

f = plt.imread('HUBBLE.jpg')
f = f[:res,:res,0]
index = np.array([7]) #index = np.array([3,7])
Hubble =[]

run1,run2,run3 = 0,0,0
for i in range(0,len(zL)):
    aL[i] = round(aL[i],1)
    asrc = Distances.scalefactor(z)
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL[i],asrc)
    Dd = Distances.angular(a0,aL[i])
    critdens = (c*Ds)/(4*np.pi*G*Dd*Dds) #kg/m^2
    x,y = X-(Dds/Ds)*Xgrad/critdens, Y-(Dds/Ds)*Ygrad/critdens  
    for p in range (0,len(noisypotential)):
        for k in index:
            for t in index:
                f = f*0
                #f[245+k,245+t] = 1
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
                tau = (1+zL[i])*Ds*Dd/Dds*tnodim #10^-9 s due to unit microrad
                tau /= 1e12  #arrival time surface in s                           
                tmin = np.min(tau)
                tmax = np.max(tau)
                levs = np.linspace(tmin,tmin + (tmax-tmin)/5,500)
                plt.gca().set_aspect('equal')
                if k == 7 and t == 7:
                    if aL[i] == 0.5:
                        x1 = np.array([-9,14,0,-2])
                        y1 = np.array([-4,4,-7,7]) 
                        for m in range (0,len(x1)):
                            plt.plot(x1[m],y1[m],'r.',markersize = 1)                     
                            for n in range(m+1,len(x1)):
                                Hnodim = ((X[0,res2+x1[m]]-x0)**2 + (Y[res2+y1[m],0]-y0)**2 -(X[0,res2+x1[n]]-x0)**2 -(Y[res2+y1[n],0]-y0)**2)/2 - np.abs(noisypotential[p][res2+y1[m],res2+x1[m]]-noisypotential[p][res2+y1[n],res2+x1[n]])/critdens*(Dds/Ds) #x and y are flipped for potential
                                Hnodim /= 1e12
                                timedelay = tau[res2+y1[m],res2+x1[m]]-tau[res2+y1[n],res2+x1[n]]
                                H = H0*(1+zL[i])*Ds*Dd/Dds*Hnodim/timedelay
                                H = round(H,4)
                                Hubble.append(H)
            
            
                    if aL[i] == 0.6:
                        x1 = np.array([-13,17,0,-3])
                        y1 = np.array([-5,4,-9,10])
                        for m in range (0,len(x1)):
                            plt.plot(x1[m],y1[m],'r.',markersize = 1)                     
                            for n in range(m+1,len(x1)):
                                Hnodim = ((X[0,res2+x1[m]]-x0)**2 + (Y[res2+y1[m],0]-y0)**2 -(X[0,res2+x1[n]]-x0)**2 -(Y[res2+y1[n],0]-y0)**2)/2 - np.abs(noisypotential[p][res2+y1[m],res2+x1[m]]-noisypotential[p][res2+y1[n],res2+x1[n]])/critdens*(Dds/Ds) #x and y are flipped for potential
                                Hnodim /= 1e12
                                timedelay = tau[res2+y1[m],res2+x1[m]]-tau[res2+y1[n],res2+x1[n]]
                                H = H0*(1+zL[i])*Ds*Dd/Dds*Hnodim/timedelay
                                H = round(H,4)
                                Hubble.append(H)
                                                
                    if aL[i] == 0.7:
                        x1 = np.array([-14,19,0,-3])
                        y1 = np.array([-7,4,-11,11])
                        for m in range (0,len(x1)):
                            plt.plot(x1[m],y1[m],'r.',markersize = 1)                     
                            for n in range(m+1,len(x1)):
                                Hnodim = ((X[0,res2+x1[m]]-x0)**2 + (Y[res2+y1[m],0]-y0)**2 -(X[0,res2+x1[n]]-x0)**2 -(Y[res2+y1[n],0]-y0)**2)/2 - np.abs(noisypotential[p][res2+y1[m],res2+x1[m]]-noisypotential[p][res2+y1[n],res2+x1[n]])/critdens*(Dds/Ds) #x and y are flipped for potential
                                Hnodim /= 1e12
                                timedelay = tau[res2+y1[m],res2+x1[m]]-tau[res2+y1[n],res2+x1[n]]
                                H = H0*(1+zL[i])*Ds*Dd/Dds*Hnodim/timedelay
                                H = round(H,4)
                                Hubble.append(H)
                                                
                    if aL[i] == 0.8:
                        x1 = np.array([-10,15,0,-3])
                        y1 = np.array([-4,4,-8,8])
                        for m in range (0,len(x1)):
                            plt.plot(x1[m],y1[m],'r.',markersize = 1)                     
                            for n in range(m+1,len(x1)):
                                Hnodim = ((X[0,res2+x1[m]]-x0)**2 + (Y[res2+y1[m],0]-y0)**2 -(X[0,res2+x1[n]]-x0)**2 -(Y[res2+y1[n],0]-y0)**2)/2 - np.abs(noisypotential[p][res2+y1[m],res2+x1[m]]-noisypotential[p][res2+y1[n],res2+x1[n]])/critdens*(Dds/Ds) #x and y are flipped for potential
                                Hnodim /= 1e12
                                timedelay = tau[res2+y1[m],res2+x1[m]]-tau[res2+y1[n],res2+x1[n]]
                                H = H0*(1+zL[i])*Ds*Dd/Dds*Hnodim/timedelay
                                H = round(H,4)
                                Hubble.append(H)
                                                        
            plt.contour(Y,X,tau,levels=levs)
            #plt.colorbar(cs)
            plt.grid()
            plt.title('arrival time surface')
            plt.title('%f %i %i' % (aL[i],k,t))
            plt.pause(0.1)
            plt.clf()

    
HubbleA = np.asarray(Hubble)

plt.figure()
plt.title('Hubble Parameter in ')
plt.grid()
plt.hist(HubbleA, bins = 200,rwidth = 0.6)
plt.savefig('HubbleParameter.png')
plt.show()