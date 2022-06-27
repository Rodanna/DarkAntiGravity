import numpy as np
import matplotlib.pyplot as plt

def poten(x,y,a):
    xm, xp, ym, yp = x - a/2, x + a/2, y - a/2, y + a/2
    val = -3*a*a \
    + xm*xm*np.arctan(ym/xm) + ym*ym*np.arctan(xm/ym) \
    + xp*xp*np.arctan(yp/xp) + yp*yp*np.arctan(xp/yp) \
    - xm*xm*np.arctan(yp/xm) - yp*yp*np.arctan(xm/yp) \
    - xp*xp*np.arctan(ym/xp) - ym*ym*np.arctan(xp/ym) \
    + xm*ym*np.log(xm*xm + ym*ym) \
    + xp*yp*np.log(xp*xp + yp*yp) \
    - xp*ym*np.log(xp*xp + ym*ym) \
    - xm*yp*np.log(xm*xm + yp*yp)
    return val/(2*np.pi)

def poten_x(x,y,a):
    xm, xp, ym, yp = x - a/2, x + a/2, y - a/2, y + a/2
    val = xm*np.arctan(ym/xm) + xp*np.arctan(yp/xp) \
        - xm*np.arctan(yp/xm) - xp*np.arctan(ym/xp) \
        + ym*np.log(xm*xm + ym*ym)/2 + yp*np.log(xp*xp + yp*yp)/2 \
        - ym*np.log(xp*xp + ym*ym)/2 - yp*np.log(xm*xm + yp*yp)/2
    return val/np.pi;
       
def poten_y(x,y,a):
    xm, xp, ym, yp = x - a/2, x + a/2, y - a/2, y + a/2
    val = ym*np.arctan(xm/ym) + yp*np.arctan(xp/yp) \
        - ym*np.arctan(xp/ym) - yp*np.arctan(xm/yp) \
        + xm*np.log(xm*xm + ym*ym)/2 + xp*np.log(xp*xp + yp*yp)/2 \
        - xm*np.log(xm*xm + yp*yp)/2 - xp*np.log(xp*xp + ym*ym)/2
    return val/np.pi;

   
rmax = 150 
u = np.linspace(-rmax,rmax,256)
X,Y = np.meshgrid(u,u)
r = np.sqrt(X**2+Y**2)
phi = np.arctan2(Y,X)


pi = 5

Xgrad = -np.loadtxt('Xgrad.txt')
Ygrad = -np.loadtxt('Ygrad.txt')

plt.figure()
plt.gca().set_aspect('equal')
plt.quiver(X[::pi,::pi],Y[::pi,::pi],Xgrad[::pi,::pi],Ygrad[::pi,::pi])
plt.show()



Xgrad = poten_x(X,Y,50)
Ygrad = poten_y(X,Y,50)

plt.figure()
plt.gca().set_aspect('equal')
plt.quiver(X[::pi,::pi],Y[::pi,::pi],Xgrad[::pi,::pi],Ygrad[::pi,::pi])
plt.show()

