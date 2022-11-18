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
aL = 0.65
zL = 1/aL-1
rmax = 150
rmax2 = 150
res = 2048
res2 = int(res/2)
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
zS = np.linspace(9,14,7)
rmaxS = 10
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

resfine = res*5
resfine2 = int(resfine/2)

ufine = np.linspace(-rmaxS,rmaxS,resfine)

Xgrad = -np.loadtxt('Xgrad2_2048.txt')
Ygrad = -np.loadtxt('Ygrad2_2048.txt')
g = np.zeros((7,res,res,3))
 
a = 2
len0y = int(2048/a)
len0x = int(2048/a)
len1y = int(2048/a)
len1x = int(1580/a)
len2y = int(1638/a)
len2x = int(2048/a)
len3y = int(2048/a) 
len3x = int(1168/a)
len4y = int(2048/a) 
len4x = int(1598/a)
len5y = int(2032/a)
len5x = int(2048/a)
len6y = int(2048/a)
len6x = int(1526/a)-1


ff = np.empty((7,resfine,resfine,3),float)
#ff[0][resfine2-len0x:resfine2+len0x,resfine2-len0y:resfine2+len0y] = plt.imread('HUBBLE.jpg')/256
ff[1][resfine2-len1x-3000:resfine2+len1x-3000,resfine2-len1y:resfine2+len1y] = plt.imread('Galaxy1.jpg')/256
ff[2][resfine2-len2x:resfine2+len2x,resfine2-len2y+3000:resfine2+len2y+3000] = plt.imread('Galaxy2.jpg')/256
ff[3][resfine2-len3x+20:resfine2+len3x+21,resfine2-len3y:resfine2+len3y] = plt.imread('Galaxy3.jpg')/256 #fix
ff[4][resfine2-len4x-2000:resfine2+len4x+1-2000,resfine2-len4y-3800:resfine2+len4y-3800] = plt.imread('Galaxy4.jpg')/256
ff[5][resfine2-len5x+3000:resfine2+len5x+3000,resfine2-len5y:resfine2+len5y] = plt.imread('Galaxy5.jpg')/256
ff[6][resfine2-len6x:resfine2+len6x+1,resfine2-len6y-3500:resfine2+len6y-3500] = plt.imread('Galaxy6.jpg')[:-1,:,:]/256


for i in range(1,len(zS)):
    asrc = Distances.scalefactor(zS[i])
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL,asrc)
    Dd = Distances.angular(a0,aL)
    critdens = (c*Ds)/(4*np.pi*G*Dd*Dds) #kg/m^2
    
    #x,y = X-(Dds/Ds)*Xgrad/critdens, Y-(Dds/Ds)*Ygrad/critdens
    x,y = X,Y
    spline0 = RectBivariateSpline(ufine,ufine,ff[i][:,:,0].T) #physicak length has to match if pieces are cut out
    spline1 = RectBivariateSpline(ufine,ufine,ff[i][:,:,1].T)
    spline2 = RectBivariateSpline(ufine,ufine,ff[i][:,:,2].T)
    #x,y = X,Y
    
    g[i][:,:,0] = np.clip(spline0.ev(x,y),0,1)
    g[i][:,:,1] = np.clip(spline1.ev(x,y),0,1)
    g[i][:,:,2] = np.clip(spline2.ev(x,y),0,1)
    
       
    plt.figure(frameon=False)
    ax = plt.gca()
    ax.set_xlabel([-150,150])
    ax.set_xlabel([-150,150])
    #plt.imsave(f'Lensedimage{i}nolens.png',g[i])
    plt.imsave(f'Lensedimage{i}.png',ff[i])
    plt.title('source image')
    plt.xlabel('x in microradians')
    plt.ylabel('x in microradians')
    plt.imshow(ff[i],origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
    plt.pause(1)
    plt.clf()
    plt.figure(frameon=False)
    ax = plt.gca()
    ax.set_xlabel([-150,150])
    ax.set_xlabel([-150,150])
    plt.xlabel('x in microradians')
    plt.ylabel('x in microradians')
    plt.imshow(g[i],origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
    plt.pause(2)
    plt.clf()


totalf = np.zeros((2048,2048,3)) 
    
plt.figure()
for i in range(0,7):
    totalf += plt.imread(f'Lensedimage{i}.png')[:,:,:3]
ax = plt.gca()
ax.set_xlabel([-150,150])
ax.set_xlabel([-150,150])
plt.xlabel('x in microradians')
plt.ylabel('x in microradians')
plt.imshow(totalf,origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
plt.savefig('LensedHDF.pdf')
#plt.imsave('LensedHDF.jpg',totalf)
plt.show()