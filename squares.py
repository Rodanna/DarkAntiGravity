import numpy as np

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
    
def poten_xx(x,y,a):
    xm, xp, ym, yp = x - a/2, x + a/2, y - a/2, y + a/2
    val = np.arctan(yp/xp) + np.arctan(ym/xm) - np.arctan(yp/xm) - np.arctan(ym/xp)
    return val/np.pi
    
def poten_yy(x,y,a):
    xm, xp, ym, yp = x - a/2, x + a/2, y - a/2, y + a/2
    val = np.arctan(xp/yp) + np.arctan(xm/ym) - np.arctan(xp/ym) - np.arctan(xm/yp)
    return val/np.pi
    
def poten_xy(x,y,a):
    xm, xp, ym, yp = x - a/2, x + a/2, y - a/2, y + a/2
    val = np.log(xp*xp+yp*yp) + np.log(xm*xm+ym*ym) - np.log(xp*xp+ym*ym) - np.log(xm*xm+yp*yp)
    return val/(2*np.pi)

def kappa(x,y,a):
    xm, xp, ym, yp = x - a/2, x + a/2, y - a/2, y + a/2
    val = np.arctan(yp/xp) + np.arctan(ym/xm) \
        - np.arctan(yp/xm) - np.arctan(ym/xp) \
        + np.arctan(xp/yp) + np.arctan(xm/ym) \
        - np.arctan(xp/ym) - np.arctan(xm/yp)
    return val/(2*np.pi)

def sumoverpix(f,x,y):
    a = 1 - 2**.5/1e10
    g = f(x-1,y,a) + f(x,y,a) + f(x+1,y,a)
    return 4*g
