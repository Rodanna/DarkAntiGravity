# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 12:11:34 2022

@author: Anna Rodriguez
"""

import numpy as np
import matplotlib.pyplot as plt

def foo(x,y):
    error = np.zeros_like(galaxy)
    for i in range(0,120):
        for j in range(0,120):
            if y[i,j] != 0:
                error[i,j] = x[i,j]/y[i,j]
            else:
                error[i,j] = 0
    return error
    
rmax = 150 #microradians
res = 120
invcritdens0 = 0.49823288405658067 #(kg/m^2)^-1
u = np.linspace(-rmax,rmax,res)
X,Y = np.meshgrid(u,u)
r = np.sqrt(X**2+Y**2) 
phi = np.arctan2(Y,X)

lev = np.linspace(0,0.3,20)
levs = np.linspace(1,8,2)

galaxy = np.loadtxt('galaxy2.txt')*invcritdens0 
potential = np.loadtxt('potential2_120(nr=2238).txt')*invcritdens0         
kappa = np.zeros((43,120,120))
pot = np.zeros((43,120,120))

errorkappa = np.zeros((43,120,120))
errorpot = np.zeros((43,120,120))


for i in range (0,43):
    kappa[i,:,:] = np.loadtxt(f'fit2_120(nr={(i+1)*50}).txt')*invcritdens0
    pot[i,:,:] = np.loadtxt(f'potential2_120(nr={(i+1)*50}).txt')*invcritdens0
    errorkappa[i,:,:] = kappa[i,:,:] - galaxy
    errorpot[i,:,:] = pot[i,:,:] - potential
    #errorpot[i,:,:] = pot[i,:,:] - np.mean(pot[i,:,:]) 

meankappa = np.zeros(43)
errorkappa1d = np.zeros(43)
rmserrorkappa = np.zeros(43)
rmskappa = np.zeros(43)
errorpot1d = np.zeros(43)
rmspot = np.zeros(43)
rmserrorpot = np.zeros(43)
meanpot = np.zeros(43)


for i in range (0,43):
    plt.clf()
    plt.title(r'Total error of $\kappa$ ' f' (nr = {(i+1)*50})')
    plt.gca().set_aspect('equal')
    plt.contourf(X,Y,errorkappa[i,:,:],levels = lev, cmap='Reds') 
    plt.colorbar()
    plt.contour(X,Y,galaxy,levels=levs,cmap='gist_gray',linewidths=0.5)
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.savefig(f'2total{(i+1)*50}.png')
    plt.show()
    meankappa[i] = np.mean(np.abs(kappa[i,:,:]))
    errorkappa1d[i] = np.mean(np.abs(errorkappa[i,:,:]))
    rmskappa[i] = np.sqrt(np.mean(kappa[i,:,:]**2))
    rmserrorkappa[i] = np.sqrt(np.mean(errorkappa[i,:,:]**2))

  
for i in range (0,43):
    plt.clf()
    plt.title(f'Potential error (nr = {(i+1)*50})')
    plt.gca().set_aspect('equal')
    plt.contourf(X,Y,errorpot[i,:,:],levels = lev, cmap='Reds') 
    plt.colorbar()
    plt.contour(X,Y,galaxy,levels=levs,cmap='gist_gray',linewidths=0.5)
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.savefig(f'2potential{(i+1)*50}.png')
    plt.show()
    meanpot[i] = np.mean(np.abs(pot[i,:,:]))
    errorpot1d[i] = np.mean(np.abs(errorpot[i,:,:]))
    rmserrorpot[i] = np.sqrt(np.mean(errorpot[i,:,:]**2))
    rmspot[i] = np.sqrt(np.mean(pot[i,:,:]**2))


x = np.arange(0,43)

plt.figure()
plt.title('Mean mass error in dependence of number of Bessel roots')
plt.grid()
plt.plot((x+1)*50,errorkappa1d,'b+',label='Mean error of $\kappa$')
plt.plot((x+1)*50,meankappa,'b+',label=r'Mean $\kappa$')
plt.xlabel('number of bessel roots')
plt.ylabel(r'Error of $\kappa$')
plt.legend()
plt.savefig('Meanmasserror.png')
plt.show()

plt.figure()
plt.title('RMS of mass in dependence of number of Bessel roots')
plt.grid()
plt.plot((x+1)*50,rmserrorkappa,'b+',label ='RMS of error of $\kappa$')
plt.plot((x+1)*50,rmskappa,'b+',label ='RMS $\kappa$')
plt.xlabel('number of bessel roots')
plt.ylabel('RMS of $\kappa$')
plt.legend()
plt.savefig('RMSmasserror.png')
plt.show()

plt.figure()
plt.title('Mean potential error in dependence of number of Bessel roots')
plt.grid()
plt.plot((x+1)*50,errorpot1d,'b+',label='Mean error of potential')
plt.plot((x+1)*50,meanpot,'b+',label='Mean potential')
plt.xlabel('number of bessel roots')
plt.ylabel('Potential Error')
plt.legend()
plt.savefig('Meanpotentialerror.png')
plt.show()

plt.figure()
plt.title('RMS of potential in dependence of number of Bessel roots')
plt.grid()
plt.plot((x+1)*50,rmserrorpot,'b+',label='RMS of error of potential')
plt.plot((x+1)*50,rmspot,'b+',label='RMS of potential')
plt.xlabel('number of bessel roots')
plt.ylabel('RMS of Potential')
plt.legend()
plt.savefig('RMSpotentialerror.png')
plt.show()


'''
for i in range (0,43):
    plt.clf()
    plt.title(r'Fractional error of $\kappa$ 'f' (nr = {(i+1)*50})')
    plt.gca().set_aspect('equal')
    plt.contourf(X,Y,foo(error[i,:,:],galaxy),levels = lev, cmap='Reds') 
    plt.colorbar()
    plt.contour(X,Y,galaxy,levels=levs,cmap='gist_gray',linewidths=0.5)
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.savefig(f'2fractional{(i+1)*50}.png')
    plt.show()


numbers = np.array([0,4,10,30,41])
plt.figure() 
plt.grid()   
for i in numbers:
    err = (foo(error[i,:,:],galaxy)).flatten()
    #plt.title('Fractional mass difference in dependence of number of Bessel roots')
    plt.hist(err,bins = 100,histtype='step',label=f'nr = {(i+1)*50}')
plt.xlabel('Fractional Mass Error')
plt.ylabel('Count')
plt.legend()
plt.savefig('2Fractionalmasserrorcount.png')
plt.show()'''