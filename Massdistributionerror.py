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

galaxy = np.loadtxt('galaxy1.txt')*invcritdens0 
potential = np.loadtxt('potential1_120(nr=2238).txt')*invcritdens0         
kappa = np.zeros((43,120,120))
pot = np.zeros((43,120,120))

errorkappa = np.zeros((43,120,120))
errorpot = np.zeros((43,120,120))


for i in range (0,43):
    kappa[i,:,:] = np.loadtxt(f'fit1_120(nr={(i+1)*50}).txt')*invcritdens0
    pot[i,:,:] = np.loadtxt(f'potential1_120(nr={(i+1)*50}).txt')*invcritdens0
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
    plt.savefig(f'1total{(i+1)*50}.pdf')
    plt.show()
    meankappa[i] = np.mean(np.abs(kappa[i,:,:]))
    errorkappa1d[i] = np.mean(np.abs(errorkappa[i,:,:]))
    rmskappa[i] = np.sqrt(np.mean((kappa[i,:,:]-meankappa[i])**2))
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
    plt.savefig(f'1potential{(i+1)*50}.pdf')
    plt.show()
    meanpot[i] = np.mean(np.abs(pot[i,:,:]))
    errorpot1d[i] = np.mean(np.abs(errorpot[i,:,:]))
    rmserrorpot[i] = np.sqrt(np.mean(errorpot[i,:,:]**2))
    rmspot[i] = np.std(pot[i,:,:])


x = np.arange(0,43)

fig, axs = plt.subplots(2)
fig.suptitle('Mean $\kappa$ and $\kappa$ error')
axs[1].plot((x+1)*50,errorkappa1d,'b+',label='Mean error of $\kappa$')
axs[1].label_outer()
axs[1].grid(True)
axs[1].legend()
axs[0].plot((x+1)*50,meankappa,'r+',label=r'Mean of $\kappa$')
axs[0].label_outer()
axs[0].grid(True)
axs[0].legend()
plt.xlabel('number of Bessel roots')
plt.savefig('1Meanmasserror.pdf')
plt.show()

fig, axs = plt.subplots(2)
fig.suptitle('RMS of $\kappa$ and $\kappa$ error')
axs[1].plot((x+1)*50,rmserrorkappa,'b+',label ='RMS of error of $\kappa$')
axs[1].label_outer()
axs[1].grid(True)
axs[1].legend()
axs[0].plot((x+1)*50,rmskappa,'r+',label ='RMS of $\kappa$')
axs[0].label_outer()
axs[0].grid(True)
axs[0].legend()
plt.xlabel('number of Bessel roots')
plt.savefig('1RMSmasserror.pdf')
plt.show()

fig, axs = plt.subplots(2)
fig.suptitle('Mean potential and potential error in pikoradian')
axs[1].semilogy((x+1)*50, errorpot1d,'r+',label='Mean error of potential')
axs[1].label_outer()
axs[1].grid(True)
axs[1].legend()
axs[0].plot((x+1)*50,meanpot,'b+',label='Mean of potentiall')
axs[0].label_outer()
axs[0].grid(True)
axs[0].legend()
plt.xlabel('number of Bessel roots')
plt.savefig('1Meanpotentialerror.pdf')
plt.show()

fig, axs = plt.subplots(2)
fig.suptitle('RMS of potential and potential error in pikoradian')
axs[0].plot((x+1)*50, rmspot,'r+',label='RMS of potential')
axs[0].label_outer()
axs[0].grid(True)
axs[0].legend()
axs[1].semilogy((x+1)*50,rmserrorpot,'b+',label='RMS of error of potential')
axs[1].label_outer()
axs[1].grid(True)
axs[1].legend()
plt.xlabel('number of Bessel roots')
plt.savefig('1RMSpotentialerror.pdf')
plt.show()


for i in range (0,43):
    plt.clf()
    plt.title(r'Fractional error of $\kappa$ 'f' (nr = {(i+1)*50})')
    plt.gca().set_aspect('equal')
    plt.contourf(X,Y,foo(errorkappa[i,:,:],galaxy),levels = lev, cmap='Reds') 
    plt.colorbar()
    plt.contour(X,Y,galaxy,levels=levs,cmap='gist_gray',linewidths=0.5)
    plt.xlabel('x in microradians')
    plt.ylabel('y in microradians')
    plt.savefig(f'1fractional{(i+1)*50}.pdf')
    plt.show()


numbers = np.array([0,4,10,30,41])
plt.figure() 
plt.grid()   
for i in numbers:
    err = (foo(errorkappa[i,:,:],galaxy)).flatten()
    #plt.title('Fractional mass difference in dependence of number of Bessel roots')
    plt.hist(err,bins = 100,histtype='step',label=f'nr = {(i+1)*50}')
plt.xlabel('$\Delta$m/m')
plt.legend()
plt.savefig('Fractionalmasserrorcount1.pdf')
plt.show()