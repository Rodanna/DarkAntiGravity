# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 21:02:36 2022

@author: Anna Rodriguez
"""

import matplotlib.pyplot as plt
import numpy as np

norm = 5000000000
def trfun(x,y):
    q = 0.8
    C= 20000
    R = 0.01
    s = 5
    #return cs*x-sn*y, sn*x+cs*y #ROTATION
    #return x+(2*x+y)*(x**2+y**2+x*y)**2/norm, 2*y+(2*y+x)*(x**2/10+y**2/5+x*y)**2/norm #drop
    #return x*np.sqrt(np.abs(1-y**2/2))/100, y*np.sqrt(np.abs(1-x**2/2))/100
    #return x*y/50,y*1.2
    #return 1/(np.sqrt(1-q**2))*np.arctan(np.sqrt(1-q**2)*x/(R+s)), 1/(np.sqrt(1-q**2))*np.arctan(np.sqrt(1-q**2)*y/(R+q**2*s))
    return 2*x/(1+x**2+q*y**2)*C, 2*q*y/(1+x**2+q*y**2)*C+10


def display(fname):
    src = plt.imread(fname)
    src = 1/255 * src
    src = src[:,:,:3]
    plt.style.use('dark_background')
    img = 0*src
    isize = img.shape[0]
    jsize = img.shape[1]
    imid = isize//2
    jmid = jsize//2
    for i in range(isize):
        for j in range(jsize):
            x = i - imid
            y = j - jmid
            X,Y = trfun(x,y)
            si = imid + int(X)
            sj = jmid + int(Y)
            if 0 <= si < isize and 0 <= sj < jsize:
                img[i,j,:] = src[si,sj,:]
    plt.grid(False)
    plt.imshow(img)
    plt.show()

print()
display('TREE.jpg')