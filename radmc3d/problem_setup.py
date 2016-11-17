#!/usr/bin/env python

import numpy as np
from scipy.integrate import quad

# Some natural constants
AU  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
MS  = 1.98892e33     # Solar mass              [g]
TS  = 5.78e3         # Solar temperature       [K]
LS  = 3.8525e33      # Solar luminosity        [erg/s]
RS  = 6.96e10        # Solar radius            [cm]

# Monte Carlo parameters
nphot = 1000000

# Grid parameters
nx       = 128
ny       = 128
nz       = 128
sizex    = 200*AU
sizey    = 200*AU
sizez    = 200*AU

# Model parameters
radius   = 5*AU
rho0     = 1e-16

# Star paramters
mstar    = 0.55*MS
rstar    = RS
tstar    = TS
pstar    = [0.,0.,0.]

# Make the coordinates
xi       = np.linspace(-sizex, sizex, nx+1)
yi       = np.linspace(-sizey, sizey, ny+1)
zi       = np.linspace(-sizez, sizez, nz+1)
xc       = 0.5 * ( xi[0:-1]+ xi[1:] )
yc       = 0.5 * ( yi[0:-1]+ yi[1:] )
zc       = 0.5 * ( zi[0:-1]+ zi[1:] )

# Make the dust density model
p = 1.064
q = 0.43
Rc = 79.*AU
T0 = 70.     # [K]
R0 = 10.*AU
mdisk = 0.1349*MS

Ts0 = 400.
rs0 = 3.*AU

gas_to_dust = 100.

gamma = -1.5+p+q/2
alpha = 7./2-p-q/2
def func(x):
  return x**(1-gamma)*np.exp(-x**alpha)
I = quad(func, 0, np.Inf)[0]

k = 1.38e-16   # erg/K
G = 6.673e-8   # cgs
mbar = 2.36 * 1.67e-24 # g

rho0 = mdisk/(I*np.sqrt(2*np.pi*k * T0/G/mstar/mbar)*(Rc/R0)**(-q/2.)*Rc**3.5)/gas_to_dust
print rho0

zz, yy, xx = np.meshgrid(zc,yc,xc, indexing = 'ij') # radmc is Fortran code!
rr       = np.sqrt(xx**2+yy**2+zz**2)
RR       = np.sqrt(xx**2+yy**2)

TR = T0*(R0/RR)**q
#HR = np.sqrt(2*k*TR*RR**3/G/mstar/mbar)
HR = 16.8*AU*(RR/Rc)**(1.5-q/2.)

rhod = rho0*(RR/Rc)**(-p)*np.exp(-(RR/Rc)**alpha)*np.exp(-(zz/HR)**2)
rhod[zz>3*HR]=0
rhod[zz<-3*HR]=0
rhod[RR>200*AU]=0

W = np.exp(-(zz/3/HR)**2)
Td = W*TR + (1-W) * Ts0*(rs0/rr)**q

# Make the wavelength array
lambda1 = 0.1
lambda2 = 7.0
lambda3 = 25.
lambda4 = 1.0e4
n12     = 20
n23     = 100
n34     = 30
lam12   = lambda1 * (lambda2/lambda1)**(np.arange(n12)/(1.*n12))
lam23   = lambda2 * (lambda3/lambda2)**(np.arange(n23)/(1.*n23))
lam34   = lambda3 * (lambda4/lambda3)**(np.arange(n34)/(1.*(n34-1.)))
lam12.shape = 1,-1
lam23.shape = 1,-1
lam34.shape = 1,-1
lamb  = np.column_stack((lam12,lam23,lam34))
lamb.shape = -1
"""
lamb    = np.logspace(1, 4, 200)
"""
nlam    = len(lamb)

# Write the wavelength_micron.inp file
f = open('wavelength_micron.inp', 'w')
#f.write()
print >> f, nlam
np.savetxt(f, lamb)
f.close()

# Write the stars.inp file
f = open('stars.inp', 'w')
print >> f, 2
print >> f, 1, nlam
print >> f 
print >> f, rstar, mstar, pstar[0], pstar[1], pstar[2]
print >> f 
np.savetxt(f, lamb)
print >> f 
print >> f, -tstar
f.close()

# Write the grid file
f = open('amr_grid.inp', 'w')
print >> f,1                      # iformat
print >> f,0                      # AMR grid style  (0=regular grid, no AMR)
print >> f,0                      # Coordinate system
print >> f,0                      # gridinfo
print >> f,1,1,1                  # Include x,y,z coordinate
print >> f,nx,ny,nz               # Size of grid
np.savetxt(f, xi)
np.savetxt(f, yi)
np.savetxt(f, zi)
f.close()

# Write the density file
f = open('dust_density.inp', 'w')
print >> f, 1
print >> f, nx*ny*nz
print >> f, 1
rhod.shape = -1
np.savetxt(f, rhod)
f.close()
"""
# Write the density file
fb = open('dust_density.binp', 'w')
print >> fb, bytearray(1)
print >> fb, bytearray(8)
print >> fb, bytearray(nx*ny*nz)
print >> fb, bytearray(1)
rhod.shape = -1
rhod.tofile(fb)
f.close()
"""

# Write the temperature file
f = open('dust_temperature.dat', 'w')
print >> f, 1
print >> f, nx*ny*nz
print >> f, 1
Td.shape = -1
np.savetxt(f, Td)
f.close()

# Dust opacity control file
f = open('dustopac.inp', 'w')
print >> f,'2               Format number of this file'
print >> f,'1               Nr of dust species'
print >> f,'============================================================================'
print >> f,'10              Way in which this dust species is read'
print >> f,'0               0=Thermal grain'
print >> f,'kataoka         Extension of name of dustkappa_***.inp file'
print >> f,'----------------------------------------------------------------------------'
f.close()

# Write the radmc3d.inp control file
f = open('radmc3d.inp', 'w')
print >> f,'nphot = ',nphot
print >> f,'nphot_scat = ',nphot
print >> f,'scattering_mode_max = 4'
#print >> f,'iranfreqmode = 1'
f.close()
