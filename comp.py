#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

f = open("radmc3d/image.out")
AU = 1.496e+13
theta0 = 45./180*np.pi

iformat = int(f.readline())
im_nx, im_ny = np.array([int(i) for i in f.readline().split()])
nlam = int(f.readline())
pixsize_x, pixsize_y = np.array([float(i) for i in f.readline().split()]) / AU

if nlam==1:
    lam = float(f.readline())
else:
    print "You gotta be kidding me."
    exit(1)
f.readline()

data = np.loadtxt(f)
I1 = data[:,0]
Q1 = data[:,1]
U1 = data[:,2]
I1.shape = im_nx, im_ny
Q1.shape = im_nx, im_ny
U1.shape = im_nx, im_ny
I1 = I1.T
Q1 = Q1.T
U1 = U1.T

x, y = np.meshgrid(np.arange(-im_nx/2, (im_nx+1)/2)*pixsize_x,
        np.arange(-im_ny/2, (im_ny+1)/2)*pixsize_y)
extent = [x[0,0], x[-1,-1], y[0,0], y[-1,-1]]

p1 = np.sqrt(Q1**2+U1**2)/I1
#angle = 0.5*np.arctan2(U,Q)+np.pi/2
angle1 = 0.5*np.arctan2(U1,-Q1)

Ip1 = I1*p1

R = np.sqrt(x*x/(np.cos(theta0)**2)+y*y)
flag = R<150

plt.figure()
plt.quiver(x[flag], y[flag], p1[flag]*np.cos(angle1[flag]), p1[flag]*np.sin(angle1[flag]))
plt.savefig("pvec_radmc.pdf")
#plt.colorbar()

plt.figure()
plt.tricontourf(x[flag], y[flag], p1[flag], 20)
plt.colorbar()
plt.title("Polarization Degree")
plt.savefig("p_radmc.pdf")

"""
plt.figure()
plt.imshow(Ip1, extent=extent)
plt.colorbar()

plt.figure()
plt.imshow(I1, extent=extent)
plt.colorbar()
"""

plt.figure()
plt.tricontourf(x[flag], y[flag], Ip1[flag], 20)
plt.colorbar()
plt.title("Polarized Intensity")
plt.savefig("Ip_radmc.pdf")

"""
plt.figure()
plt.imshow(angle1)
plt.colorbar()
"""

"""
plt.figure()
plt.hist(p1[flag], bins=50)
"""

#####################################################
# Below starts the formal solution data.

data2 = np.loadtxt("fmlsol/result.dat")

R2 = data2[:,0]
phi2 = data2[:,1]

I2 = data2[:,2]
Q2 = data2[:,3]
U2 = data2[:,4]
V2 = data2[:,5]

p2 = np.sqrt(Q2**2+U2**2)/I2
theta2 = 0.5*np.arctan2(U2, Q2)

kwargs = {"headlength":0, "headaxislength":0, "headwidth":1, "pivot":'middle'}#, "scale":scale, "width":width}

plt.figure()
plt.quiver(R2*np.cos(phi2)*np.cos(theta0), R2*np.sin(phi2), p2*np.cos(theta2), p2*np.sin(theta2), angles='uv', **kwargs)
ax = plt.gca()

plt.savefig("pvec_fml.pdf")

"""
plt.figure()
plt.tricontourf(R2*np.cos(phi2)*np.cos(theta0), R2*np.sin(phi2), theta2, 20)
plt.colorbar()
plt.show()
"""

#ax.set_aspect(1/np.cos(theta0))
"""
"""

R2.shape = 15, -1
phi2.shape = 15, -1
p2.shape = 15, -1
I2.shape = 15, -1
Q2.shape = 15, -1
theta2.shape = 15, -1

phi2 = np.append(phi2.T, [[2*np.pi for _ in phi2[:,0]]],0).T
R2 = np.append(R2.T, [[i for i in R2[:,0]]],0).T
p2 = np.append(p2.T, [[i for i in p2[:,0]]],0).T
I2 = np.append(I2.T, [[i for i in I2[:,0]]],0).T
Q2 = np.append(Q2.T, [[i for i in Q2[:,0]]],0).T
theta2 = np.append(theta2.T, [[i for i in theta2[:,0]]],0).T

"""
plt.figure()
plt.contourf(R2*np.cos(phi2), R2*np.sin(phi2), p2*I2, 20) 
ax = plt.gca()
ax.set_aspect(1/np.cos(theta0))

plt.figure()
plt.hist(p2.ravel(), bins=50)
"""

#####################################################
# Let's plot them together

# x, y, p1, I1, Q1, U1, V1 are data from RADMC-3D
# R2, phi2, I2, Q2, U2, V2 are data from formal solution

# Firstly, p along major axis
plt.figure()

f1 = x==0
tny = y[f1]
tnp = p1[f1]
f2 = tny>0
nyy = tny[f2]
np1 = tnp[f2]
plt.plot(nyy, np1)

f1 = phi2==1.5708
nRR = R2[f1]
np2 = p2[f1]
plt.plot(nRR, np2)

plt.savefig("p_major.pdf")

# Let's study angles along major axis
plt.figure()

f1 = x==0
tny = y[f1]
tnp = angle1[f1]
f2 = tny>0
nyy = tny[f2]
np1 = tnp[f2]
plt.plot(nyy, -np1)

f1 = phi2==1.5708
nRR = R2[f1]
np2 = theta2[f1]
plt.plot(nRR, np2)

plt.savefig("pvecAngle_major.pdf")

# Secondly, pQ along minor axis
pQ1 = -Q1/I1
pQ2 = Q2/I2

plt.figure()

f1 = y==0
tnx = x[f1]
tnp = pQ1[f1]
f2 = tnx>0
nxx = tnx[f2]
npQ1 = tnp[f2]
plt.plot(nxx, npQ1)

#f1 = phi2==0.0
f1 = phi2==3.14159
nRR = R2[f1] * np.cos(theta0)
npQ2 = pQ2[f1]
plt.plot(nRR, npQ2)

plt.savefig("p_minor.pdf")

# Let's do something with I's
plt.figure()

f1 = x==0
tny = y[f1]
tnI = I1[f1]
f2 = tny>0
nyy = tny[f2]
nI1 = tnI[f2]
plt.plot(nyy, nI1, label='RADMC-3D')
#plt.plot(nyy, nI1*1.5)
#plt.plot(nyy, nI1*0.5)

f1 = phi2==1.5708
nRR = R2[f1]
nI2 = I2[f1]

#ratio = nI2[0]/nI1[0]
ratio = 1
plt.plot(nRR, nI2/ratio, label='formal_sol')
plt.yscale('log')
plt.legend(loc='best')
plt.title("Total intensity")

plt.savefig("Itot.pdf")

#plt.show()
