#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 19:48:59 2022

@author: annarodriguez
"""

import numpy as np
import matplotlib.pyplot as plt


Hubble10 = np.loadtxt('Hubble10.txt')
Hubble10 = np.sort(Hubble10)
median10 = Hubble10[100]
Hubble100 = np.loadtxt('Hubble100.txt')
Hubble100 = np.sort(Hubble100)
median100 = Hubble100[100]
Hubble1000 = np.loadtxt('Hubble1000.txt')
Hubble1000 = np.sort(Hubble1000)
median1000 = Hubble1000[100]


x10 = np.linspace(40,140,1000)
plt.figure()
plt.title('Hubble Parameter for $\sigma$ = 10 %')
plt.grid()
plt.hist(Hubble10, bins = 100,rwidth = 0.6)
plt.vlines(median10,0,60,'r',linestyles='dotted')
plt.vlines(Hubble10[int(len(Hubble10)*0.16)],0,60,'r',linestyles='dotted')
plt.vlines(Hubble10[int(len(Hubble10)*0.84)],0,60,'r',linestyles='dotted')
plt.xlabel('Normalized estimation for Hubble parameter')
plt.ylabel('number count')
plt.savefig('HubbleParameter10.png')
plt.show()



x100 = np.linspace(65,72,1000)
plt.figure()
plt.title('Hubble Parameter for $\sigma$ = 1 % ')
plt.grid()
plt.hist(Hubble100, bins = 100,rwidth = 0.6)
plt.vlines(median100,0,30,'r',linestyles='dotted')
plt.vlines(Hubble100[int(len(Hubble10)*0.16)],0,30,'r',linestyles='dotted')
plt.vlines(Hubble100[int(len(Hubble10)*0.84)],0,30,'r',linestyles='dotted')
plt.xlabel('Normalized estimation for Hubble parameter')
plt.ylabel('number count')
plt.savefig('HubbleParameter100.png')
plt.show()



x1000 = np.linspace(67,68,1000)
plt.figure()
plt.title('Hubble Parameter for $\sigma$ = 0.1 % ')
plt.grid()
plt.hist(Hubble1000, bins = 100,rwidth = 0.6)
plt.vlines(median1000,0,40,'r',linestyles='dotted')
plt.vlines(Hubble1000[int(len(Hubble10)*0.16)],0,40,'r',linestyles='dotted')
plt.vlines(Hubble1000[int(len(Hubble10)*0.84)],0,40,'r',linestyles='dotted')
plt.xlabel('Normalized estimation for Hubble parameter')
plt.ylabel('number count')
plt.savefig('HubbleParameter1000.png')
plt.show()