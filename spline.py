import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline

q = 0.5
C= 6000
R = 0.01
s = 500
cs = np.cos(np.pi/4)
sn = np.sin(np.pi/4)

f = plt.imread('monsters.png')
f = f[:,:,0]
t = 0*f
t[80:100,380:400] = f[80:100,380:400]
t[145:155,245:255] = 1
f = t
plt.imshow(f)
plt.show()

ny,nx = f.shape
u = np.linspace(-nx/2,nx/2,nx)
v = np.linspace(-ny/2,ny/2,ny)
X,Y = np.meshgrid(u,v)


spline = RectBivariateSpline(X[0,:],Y[:,0],f.T)
eps = 1
x1,y1 = X-2*X/(eps+X**2+q*Y**2)*C, Y-2*q*Y/(eps+X**2+q*Y**2)*C
x2,y2 = X+100*np.exp(X*1j)/(1+X**2+q*Y**2)*C, Y+100*np.exp(Y*1j)/(1+X**2+q*Y**2)*C
g = spline.ev(x1,y1)
print(np.max(f),np.max(g))
h = spline.ev(x2,y2)


plt.imshow(g)
plt.show()
plt.imshow(h)
plt.show()