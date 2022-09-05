# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 12:11:34 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def foo(x,y):
    try:
        return x/y
    except ZeroDivisionError:
        return 0
    
rmax = 150 #microradians
res = 120
invcritdens0 = 0.49823288405658067 #(kg/m^2)^-1
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)
r = np.sqrt(X**2+Y**2) 
phi = np.arctan2(Y,X)
lev = np.linspace(0,1,20)
levs = np.linspace(1,8,2)

galaxy = np.loadtxt('galaxy2.txt')
for k in range(0,len(X)):
    for l in range(0,len(Y)):
        mk = int(res/2) - k
        ml = int(res/2) - l
        if np.sqrt(mk**2+ml**2) > 60:
            galaxy[k,l] = 0
            
fit = np.zeros((42,120,120))
error = np.zeros((42,120,120))
massdispersion = np.zeros((42,120,120))


for i in range (0,42):
    fit[i,:,:] = np.loadtxt(f'fit2_120(nr={(i+1)*50}).txt')
    for k in range(0,len(X)):
        for l in range(0,len(Y)):
            mk = int(res/2) - k
            ml = int(res/2) - l
            if np.sqrt(mk**2+ml**2) > 60:
                fit[i,k,l] = 0
    error[i,:,:] = galaxy - fit[i,:,:]
'''
error1d = np.zeros(42)
for i in range (0,42):
    plt.clf()
    plt.title('total mass difference')
    plt.gca().set_aspect('equal')
    plt.contourf(X,Y,error[i,:,:]*invcritdens0,levels = lev, cmap='Reds') 
    plt.colorbar()
    plt.contour(X,Y,galaxy*invcritdens0,levels=levs,cmap='gist_gray',linewidths=0.5)
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.show()
    error1d[i] = np.mean(error[i,:,:])

x = np.arange(0,42)

plt.figure()
plt.grid()
plt.plot(x,error1d,'.')
plt.xlabel('number of bessel roots')
plt.ylabel('normalized mass error')
plt.show()

array = np.linspace(0,1,10)

for i in range (0,42):
    plt.clf()
    plt.title('fractional mass difference')
    plt.gca().set_aspect('equal')
    plt.contourf(X,Y,foo(error[i,:,:]*invcritdens0,galaxy),levels = lev, cmap='Reds') 
    plt.colorbar()
    plt.contour(X,Y,galaxy*invcritdens0,levels=levs,cmap='gist_gray',linewidths=0.5)
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.show()
'''

    
for i in range (0,42):
    err = (foo(error[i,:,:],galaxy)).flatten()
    plt.figure()
    plt.title('mass error count')
    plt.hist(err,bins = 50)
    #sns.distplot(err,hist_kws=dict(edgecolor="w", linewidth=6))
    plt.show()