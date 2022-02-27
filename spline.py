import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline

rmax = 2*np.pi
u = np.linspace(-rmax,rmax,32)
X,Y = np.meshgrid(u,u)
f = np.sin(X)*np.sin(2*Y)

spline = RectBivariateSpline(X[0,:],Y[:,0],f.T)
u = np.linspace(-rmax,rmax,256)
x,y = np.meshgrid(u,u)
g = spline.ev(x,y)

plt.imshow(f)
plt.show()

plt.imshow(g)
plt.show()