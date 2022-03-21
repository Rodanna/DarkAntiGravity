# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 18:48:05 2022

@author: Anna Rodriguez
"""

import numpy as np
from scipy import special
import matplotlib.pyplot as plt

def bisection(x,v,tol):
    a = x[0]
    b = x[-1]
    y = int(len(x)/2)
    if np.sign(special.jv(v,a)) == np.sign(special.jv(v,b)):
        return 999
    m = (a+b)/2
    if np.abs(special.jv(v,m)) < tol:
        return m
    elif np.sign(special.jv(v,a)) == np.sign(special.jv(v,m)):
        return bisection(x[y:-1],v,tol) # m closer than a
    elif np.sign(special.jv(v,b)) == np.sign(special.jv(v,m)):
        return bisection(x[0:y],v,tol) # m closer than b

v = np.array([0,1,2,3,4])
x = np.linspace(4,6,1000)

plt.figure()
plt.grid(True)
for i in range (0,len(v)):
    plt.plot(x,special.jv(v[i],x)) #bessel functions of first order
plt.show()

'''plt.figure()
plt.grid(True)
for i in range (0,len(v)):
    plt.plot(x,special.jvp(v[i],x,n=1)) #first derivative
plt.show()'''

for i in range (0,len(v)):
    r1 = bisection(x,v[i],0.01)
    if r1 != 999:
        print("FOR J",i,'THE ROOT IS', r1)

    #print("f(r1) =", f(r1))