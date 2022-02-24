# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 21:02:36 2022

@author: Anna Rodriguez
"""

import matplotlib.pyplot as plt
import numpy as np

norm = 10000000000
def trfun(x,y):
    ang = np.pi
    cs = np.cos(ang)
    sn = np.sin(ang)
    #return cs*x-sn*y, sn*x+cs*y #ROTATION
    return x+(2*x+y)*((x/2)**2+y**2+x*y)**2/norm, y+(2*y+x)*(x**2/50+y**2/5+30*x*y)**2/norm
    #return x*np.sqrt(np.abs(1-y**2/2))/100, y*np.sqrt(np.abs(1-x**2/2))/100

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

display('TREE.jpg')