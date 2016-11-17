#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

f = open("image.out")
AU = 1.496e+13

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
I = data[:,0]
Q = data[:,1]
U = data[:,2]
I.shape = im_nx, im_ny
Q.shape = im_nx, im_ny
U.shape = im_nx, im_ny
I = I.T
Q = Q.T
U = U.T

x, y = np.meshgrid(np.arange(-im_nx/2, (im_nx+1)/2)*pixsize_x,
        np.arange(-im_ny/2, (im_ny+1)/2)*pixsize_y)
extent = [x[0,0], x[-1,-1], y[0,0], y[-1,-1]]

p = np.sqrt(Q**2+U**2)/I
#angle = 0.5*np.arctan2(U,Q)+np.pi/2
angle = 0.5*np.arctan2(U,-Q)

It = I.copy()
cap = 1e-15
It[It<cap]=cap
im = plt.imshow(It, norm=LogNorm(), origin='lower', extent=extent)
plt.colorbar(im)
plt.title('Total intensity')

plt.figure()
im = plt.imshow(p, origin='lower', extent=extent)
plt.colorbar(im)
plt.title('Polarization degree')

plt.figure()
im = plt.imshow(p*I, origin='lower', extent=extent)
plt.colorbar(im)
plt.title('Polarized intensity')

r = np.sqrt(x*x+y*y)
flag = r<120
x = x[flag]
y = y[flag]
angle = angle[flag]
p = p[flag]

plt.figure()
plt.quiver(x, y, p*np.cos(angle), p*np.sin(angle))

plt.show()
