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

resfine = res*10
resfine2 = int(resfine/2)

ufine = np.linspace(-rmax,rmax,resfine)

Xgrad = -np.loadtxt('Xgrad2_2048.txt')
Ygrad = -np.loadtxt('Ygrad2_2048.txt')
potential = np.loadtxt('potential2_2048.txt')
g = np.zeros((7,res,res,3))

hf = plt.imread('HUBBLE.jpg')[:2048,:2048,:]/256

len0 = int(2048*resfine/res/2)
len1 = int(22*resfine/res/2) #pixel: 44x44 before refining pixels
len2 = int(40*resfine/res/2)
len3 = int(50*resfine/res/2)
len4 = int(50*resfine/res/2)
len5 = int(20*resfine/res/2)
len6 = int(30*resfine/res/2)


ff = np.empty((7,resfine,resfine,3),float)
#ff[0][resfine2-len0:resfine2+len0,resfine2-len0:resfine2+len0] = plt.imread('Cutout0.png')[:,:,:3]
ff[1][resfine2-len1-320:resfine2+len1-320,resfine2-len1:resfine2+len1] = plt.imread('Cutout1.png')[:,:,:3] 
ff[2][resfine2-len2:resfine2+len2,resfine2-len2+380:resfine2+len2+380] = plt.imread('Cutout2.png')[:,:,:3] 
'''ff[3][resfine2-len3+580:resfine2+len3+580,resfine2-len3:resfine2+len3] = plt.imread('Cutout3.png')[:,:,:3] 
ff[4][resfine2-len4:resfine2+len4,resfine2-len4-80:resfine2+len4-80] = plt.imread('Cutout4.png')[:,:,:3] 
ff[5][resfine2-len5+450:resfine2+len5+450,resfine2-len5:resfine2+len5] = plt.imread('Cutout5.png')[:,:,:3] 
ff[6][resfine2-len6:resfine2+len6,resfine2-len6-460:resfine2+len6-460] = plt.imread('Cutout6.png')[:,:,:3]'''

'''
for i in range (1,len(ff)):
    plt.figure(frameon=False)
    plt.axis('off')
    plt.title('source image')
    plt.imshow(ff[i],origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
    plt.pause(1)
    plt.clf()
'''
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
    #plt.title('Lensed Images')
    #plt.xlabel('x in microradians')
    #plt.ylabel('y in microradians')
    plt.imshow(g[i],origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
    plt.imsave(f'Lensedimage{i}.png',g[i])
    plt.pause(2)
    plt.clf()


total = np.zeros((res,res,3))
totalf = np.zeros((res,res,3))

'''
for i in range (0,len(ff)):
    plt.figure(frameon=False)
    plt.axis('off')
    plt.title('source image')
    plt.imshow(ff[i],origin='lower',extent=[-rmax,rmax,-rmax,rmax],interpolation='none')
    plt.pause(1)
    plt.clf()
    
    
plt.figure(frameon=False)
plt.axis('off')
for i in range(0,7):
    totalf += np.imread(f'Lensedimage{i}.png')
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