import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline

q = 0.8
C= 6000
R = 0.01
s = 500
cs = np.cos(np.pi/4)
sn = np.sin(np.pi/4)

f = plt.imread('monsters.png')
f = f[:,:,0]
plt.imshow(f)
plt.show()

ny,nx = f.shape
u = np.linspace(-nx/2,nx/2,nx)
v = np.linspace(-ny/2,ny/2,ny)
X,Y = np.meshgrid(u,v)


spline = RectBivariateSpline(X[0,:],Y[:,0],f.T)
x1,y1 = X-2*X/(1+X**2+q*Y**2)*C, Y-2*q*Y/(1+X**2+q*Y**2)*C
x2,y2 = X+100*np.exp(X*1j)/(1+X**2+q*Y**2)*C, Y+100*np.exp(Y*1j)/(1+X**2+q*Y**2)*C
g = spline.ev(x1,y1)
h = spline.ev(x2,y2)


plt.imshow(g)
plt.show()
plt.imshow(h)
plt.show()