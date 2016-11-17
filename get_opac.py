#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

# Parameters: wavelength and file to open
lam = 1e3
try:
  f = open("opac/dustkapscatmat_kataoka.inp")
except:
  print "Cannot open opacity file."
  print "You need to run opac/makeopac first."
  exit(1)

# Move opacity file to radmc3d/
print "Copying dust data file from opac/ to radmc3d/..."
os.system("cp -i opac/dustkapscatmat_kataoka.inp radmc3d/")

# Read through basic parameters
for i in range(7): f.readline() 
nlams = int(f.readline())
ntheta = int(f.readline())
f.readline()

# Get Kappa's for different wavelengths
lines1 = []
for i in range(nlams): lines1.append(f.readline())
Kave = np.array([[float(t) for t in line.split()] for line in lines1])
lams = Kave[:,0]
Kabs = Kave[:,1]
Ksca = Kave[:,2]
g    = Kave[:,3]

# Interpolate to get the value we want at lam
fKabs = interp1d(lams, Kabs)
fKsca = interp1d(lams, Ksca)
kabs, ksca = float(fKabs(lam)), float(fKsca(lam))
print "Kabs=", kabs
print "Ksca=", ksca

# Get the angle grids in the file
f.readline()
thetas = []
for i in range(ntheta): thetas.append(f.readline())
thetas = np.array([float(i) for i in thetas])
f.readline()

def read_one_block():
  """Helper function to read through one block of data
  corresponding to one wavelength.
  """
  for i in range(ntheta+1):
    f.readline()

# Find the indices corresponding to the wavelengths closest to lam
index = np.where((lams>lam) == 1)[0][0]
ind1 = index-1
ind2 = index

for i in range(ind1): read_one_block() # skip uninteresting data

def plot_one_lam(tlam):
  """
  Helper function to investigate Muller Matrix at one wavelength tlam.
  Results are saved as pdf figures in working directory
  """
  # Get the Muller Matrix elements first
  lines = []
  for i in range(ntheta): lines.append(f.readline())
  f.readline()
  data = np.array([[float(t) for t in line.split()] for line in lines])
  Z11 = data[:,0]
  Z12 = data[:,1]
  Z22 = data[:,2]
  Z33 = data[:,3]
  Z34 = data[:,4]
  Z44 = data[:,5]

  # Plot Muller Matrix elements against angles
  plt.figure()
  plt.plot(thetas, Z11, label="Z11")
  plt.plot(thetas, -Z12, label="-Z12")
  plt.xlabel(r"Scattering angle \theta")
  plt.ylabel(r"Elements of Muller Matrix")
  plt.legend(loc='best')
  plt.title(r"Muller Matrix Elements for \lambda = {0:.3e}".format(tlam))
  plt.savefig("muller{0:.2e}.pdf".format(tlam))

  # Plot Polarization degree against angles
  plt.figure()
  plt.plot(thetas, -Z12/Z11)
  plt.xlabel(r"Scattering angle \theta")
  plt.ylabel(r"Polarization degree")
  plt.title(r"Polarization degree for \lambda = {0:.3e}".format(tlam))
  plt.savefig("p{0:.2e}.pdf".format(tlam))

# Call helper function to plot wavelength of interest. (not really)
# Note: You must call ind1 first before calling ind2. I didn't rewind the
# file so it goes in only one direction. You need to modify the code to 
# be able to do so.
plot_one_lam(lams[ind1])
#plot_one_lam(lams[ind2])

f.close()

# Let's generate param.inp for fml_sol
text = """hltau\t\t; Problem type
-1\t\t; rho0 at Rc
{0:e}\t; Cabs
{1:e}\t; Csca
45\t\t; viewing angle
result.dat\t; output data file name
"""

f = open("fmlsol/param.inp","w")
f.write(text.format(kabs, ksca))
