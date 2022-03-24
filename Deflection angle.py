# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 15:35:24 2022

@author: Anna Rodriguez
"""

import numpy as np
import Distances

steps = 1000
G = 6.67e-11 # m**3/(kg s**2)
c = 300000000 # m/s
M = 5e14*2e30 # kg
beta = 5

zS = 1.41 #redshift of Q0957+561
aS = Distances.scalefactor(zS)
zL = 0.36 # redshift of lensing galaxy
aL = Distances.scalefactor(zL)
z0 = 0
a0 = 1

Dls = Distances.angular(np.linspace(aS,aL,steps))[0] #cs
Ds = Distances.angular(np.linspace(aS,aL,steps))[0] #cs
DL = Distances.angular(np.linspace(aS,aL,steps))[0] #cs


def alpha(x):
    return 4*G*M/c**2 *Dls/(Ds*DL)*x/np.abs(x)**2

print('An estimate for the deflection angle: Beta =', beta, 'degrees -> Theta =', alpha(beta),'degrees')