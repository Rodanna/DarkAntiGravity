# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 17:30:02 2022
@author: Anna Rodriguez
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import Distances
import math

a0 = 1
aL = 0.65
zL = 1/aL-1
rmax = 150
rmax2 = 150
res = 2048
res2 = int(res/2)
G = 6.67408e-11 #m^3/kgs^2
c = 299792458 #m/s
z = np.linspace(5,14,7)
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)

resfine = res*15
resfine2 = int(resfine/2)

ufine = np.linspace(-rmax,rmax,resfine)

Xgrad = -np.loadtxt('Xgrad2_2048.txt')
Ygrad = -np.loadtxt('Ygrad2_2048.txt')
potential = np.loadtxt('potential2_2048.txt')
g = np.zeros((7,res,res,3))

a = 10
len0y = int(2048*resfine/res/2)
len0x = int(2048*resfine/res/2)
len1y = int(1024/a)
len1x = int(770/a)
len2y = int(1638/a)
len2x = int(2048/a)
len3y = int(2048/a) 
len3x = int(2268/a)
len4y = int(2048/a) 
len4x = int(1598/a)
len5y = int(2032/a)
len5x = int(2048/a)
len6y = int(2048/a)
len6x = int(2003/a)


ff = np.empty((7,resfine,resfine,3),float)
#ff[0][resfine2-len0x:resfine2+len0,resfine2-len0:resfine2+len0] = plt.imread('Cutout0.png')[:,:,:3]
f = plt.imread('Galaxy1.jpg')/256
ff[1][resfine2-len1x:resfine2+len1x,resfine2-len1y:resfine2+len1y] = f

'''
ff[2][resfine2-len2x:resfine2+len2x,resfine2-len2y+3800:resfine2+len2y+3800] = plt.imread('Galaxy2.jpg')[:,:,:3] 
ff[3][resfine2-len3x+5800:resfine2+len3x+5800,resfine2-len3y:resfine2+len3y] = plt.imread('Galaxy3.jpg')[:,:,:3] 
ff[4][resfine2-len4x:resfine2+len4x,resfine2-len4y-800:resfine2+len4y-800] = plt.imread('Galaxy4.jpg')[:,:,:3] 
ff[5][resfine2-len5x+4500:resfine2+len5x+4500,resfine2-len5y:resfine2+len5y] = plt.imread('Galaxy5.jpg')[:,:,:3] 
ff[6][resfine2-len6x:resfine2+len6x,resfine2-len6y-4600:resfine2+len6y-4600] = plt.imread('Galaxy6.jpg')[:,:,:3]'''


for i in range(1,2):
    asrc = Distances.scalefactor(z[i])
    Ds = Distances.angular(a0,asrc)
    Dds = Distances.angular(aL,asrc)
    Dd = Distances.angular(a0,aL)
    critdens = (c*Ds)/(4*np.pi*G*Dd*Dds) #kg/m^2
    x,y = X-(Dds/Ds)*Xgrad/critdens, Y-(Dds/Ds)*Ygrad/critdens  

    spline0 = RectBivariateSpline(ufine,ufine,ff[i][:,:,0].T) #physicak length has to match if pieces are cut out
    spline1 = RectBivariateSpline(ufine,ufine,ff[i][:,:,1].T)
    spline2 = RectBivariateSpline(ufine,ufine,ff[i][:,:,2].T)
    
    
    g[i][:,:,0] = np.clip(spline0.ev(x,y),0,1)
    g[i][:,:,1] = np.clip(spline1.ev(x,y),0,1)
    g[i][:,:,2] = np.clip(spline2.ev(x,y),0,1)
    
       
    plt.figure(frameon=False)
    plt.axis('off')
    plt.title('source image')
    plt.imshow(ff[i],origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
    plt.pause(1)
    plt.clf()
    plt.figure(frameon=False)
    plt.axis('off')
    #plt.title(Lensed Images)
    #plt.xlabel('x in microradians')
    #plt.ylabel('y in microradians')
    plt.imshow(g[i],origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
    plt.imsave(f'Lensedimage{i}.png',g[i])
    plt.pause(2)
    plt.clf()
    

'''
total = np.zeros((res,res,3))
totalf = np.zeros((res,res,3))
    
    
plt.figure(frameon=False)
plt.axis('off')
for i in range(0,7):
    totalf += plt.imread(f'Lensedimage{i}.png')[:,:,:3]
plt.imshow(totalf**(4/5),origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
plt.savefig('HDF.jpg',dpi = 1500,bbox_inches='tight',pad_inches = 0)
plt.show()
   
plt.figure(frameon=False)
plt.axis('off')
for i in range(0,7):
    g = plt.imread(f'Lensedimage{i}.png')[:,:,:3]
    if i < 1:
        total += g
    else:
        total += g**(7/4)
plt.imshow(total**(4/5),origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
plt.imsave('LensedHDF.jpg',np.clip(total,0,1))
plt.show()'''