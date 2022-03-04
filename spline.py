import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline

q = 0.8
C= 500
R = 100
s = 500


f = plt.imread('monsters.png')
f = f[:,:,0]
plt.imshow(f)
plt.show()

ny,nx = f.shape
u = np.linspace(-nx/2,nx/2,nx)
v = np.linspace(-ny/2,ny/2,ny)
X,Y = np.meshgrid(u,v)

spline = RectBivariateSpline(X[0,:],Y[:,0],f.T)
x1,y1 = 2*X/(1+X**2+q*Y**2)*C, 2*q*Y/(1+X**2+q*Y**2)*C
x2,y2 = -2*X/(1+X**2+q*Y**2)*C, 2*q*Y/(1+X**2+q*Y**2)*C
g = spline.ev(x1,y1)
h = spline.ev(x2,y2)

plt.imshow(h)
plt.show()
plt.imshow(g)
plt.show()