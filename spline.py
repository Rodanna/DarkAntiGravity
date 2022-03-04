import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline

f = plt.imread('monsters.png')
f = f[:,:,0]
plt.imshow(f)
plt.show()

ny,nx = f.shape
u = np.linspace(-nx/2,nx/2,nx)
v = np.linspace(-ny/2,ny/2,ny)
X,Y = np.meshgrid(u,v)

spline = RectBivariateSpline(X[0,:],Y[:,0],f.T)
x,y = 0.8*X,0.9*Y
g = spline.ev(x,y)

plt.imshow(g)
plt.show()