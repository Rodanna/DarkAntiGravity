import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import Distances

q = 0.8
C= 60000
R = 0.1
s = 50
a0 = 1
aL = 1.5
steps = 1000


f = plt.imread('HUBBLE.jpg')
f = f[:,:,0]
plt.imshow(f)
plt.show()

f1 = f*0
f1[980:1050,970:1050] = f[980:1050,970:1050]
z1 = 9
a1 = Distances.scalefactor(z1)
A1 = np.linspace(a1,a0,steps)
B1 = np.linspace(a1,aL,steps)
Ds1 = Distances.angular(A1)[0]
Dds1 = Distances.angular(B1)[0]

f2 = f*0
f2[980:1070,210:300] = f[980:1070,210:300]
z2 = 1.8
a2 = Distances.scalefactor(z2)
A2 = np.linspace(a2,a0,steps)
B2 = np.linspace(a2,aL,steps)
Ds2 = Distances.angular(A2)[0]
Dds2 = Distances.angular(B2)[0]


f3 = f*0
f3[390:450,110:200] = f[390:450,110:200]
z3 = 2
a3 = Distances.scalefactor(z3)
A3 = np.linspace(a3,a0,steps)
B3 = np.linspace(a3,aL,steps)
Ds3 = Distances.angular(A3)[0]
Dds3 = Distances.angular(B3)[0]

f4 = f*0
f4[620:690,630:700] = f[620:690,630:700]
z4 = 5
a4 = Distances.scalefactor(z4)
A4 = np.linspace(a4,a0,steps)
B4 = np.linspace(a4,aL,steps)
Ds4 = Distances.angular(A4)[0]
Dds4 = Distances.angular(B4)[0]

f5 = f*0
f5[620:700,470:550] = f[620:700,470:550]
z5 = 2.1
a5 = Distances.scalefactor(z5)
A5 = np.linspace(a5,a0,steps)
B5 = np.linspace(a5,aL,steps)
Ds5 = Distances.angular(A5)[0]
Dds5 = Distances.angular(B5)[0]

f6 = f*0
f6[490:570,710:780] = f[490:570,710:780]
z6 = 7.6
a6 = Distances.scalefactor(z6)
A6 = np.linspace(a6,a0,steps)
B6 = np.linspace(a6,aL,steps)
Ds6 = Distances.angular(A6)[0]
Dds6 = Distances.angular(B6)[0]


ny,nx = f.shape
u = np.linspace(-nx/2,nx/2,nx)
v = np.linspace(-ny/2,ny/2,ny)
X,Y = np.meshgrid(u,v)

spline1 = RectBivariateSpline(X[0,:],Y[:,0],f1.T)
x1,y1 = X-(Ds1/Dds1)*2*X/(1+X**2+q*Y**2)*C, Y-(Ds1/Dds1)*2*q*Y/(1+X**2+q*Y**2)*C
g1 = spline1.ev(x1,y1)
plt.imshow(f1)
plt.show()
plt.imshow(g1)
plt.show()


spline2 = RectBivariateSpline(X[0,:],Y[:,0],f2.T)
x2,y2 = X-(Ds2/Dds2)*2*X/(1+X**2+q*Y**2)*C, Y-(Ds2/Dds2)*2*q*Y/(1+X**2+q*Y**2)*C
g2 = spline2.ev(x2,y2)
plt.imshow(f2)
plt.show()
plt.imshow(g2)
plt.show()

spline3 = RectBivariateSpline(X[0,:],Y[:,0],f3.T)
x3,y3 = X-(Ds3/Dds3)*2*X/(1+X**2+q*Y**2)*C, Y-(Ds3/Dds3)*2*q*Y/(1+X**2+q*Y**2)*C
g3 = spline3.ev(x3,y3)
plt.imshow(f3)
plt.show()
plt.imshow(g3)
plt.show()

spline4 = RectBivariateSpline(X[0,:],Y[:,0],f4.T)
x4,y4 = X-(Ds4/Dds4)*2*X/(1+X**2+q*Y**2)*C, Y-(Ds4/Dds4)*2*q*Y/(1+X**2+q*Y**2)*C
g4 = spline4.ev(x4,y4)
plt.imshow(f4)
plt.show()
plt.imshow(g4)
plt.show()

spline5 = RectBivariateSpline(X[0,:],Y[:,0],f5.T)
x5,y5 = X-(Ds5/Dds5)*2*X/(1+X**2+q*Y**2)*C, Y-(Ds5/Dds5)*2*q*Y/(1+X**2+q*Y**2)*C
g5 = spline5.ev(x5,y5)
plt.imshow(f5)
plt.show()
plt.imshow(g5)
plt.show()

spline6 = RectBivariateSpline(X[0,:],Y[:,0],f6.T)
x6,y6 = X-(Ds6/Dds6)*2*X/(1+X**2+q*Y**2)*C, Y-(Ds6/Dds6)*2*q*Y/(1+X**2+q*Y**2)*C
g6 = spline6.ev(x6,y6)
plt.imshow(f6)
plt.show()
plt.imshow(g6)
plt.show()