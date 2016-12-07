"""
PYTHON module for RADMC3D 
(c) Attila Juhasz 2011,2012,2013,2014

This sub-module contains natural constants in CGS units 
(Translated from RADMC's IDL function problem_natconst.pro)
"""
gg = 6.672e-8   # Gravitational constant
mp = 1.6726e-24 # Mass of proton [g]
me = 9.1095e-28 # Mass of electron [g]
kk = 1.3807e-16 # Bolzmann's constant [erg/K]
hh = 6.6262e-27 # Planck's constant [erg.s]
ee = 4.8032e-10 # Unit charge 
cc  = 2.99792458e10  # Light speed [cm/s]
st = 6.6524e-25 # Thmpson cross-section [cm^2]
ss = 5.6703e-5  # Stefan-Boltzmann const [erg/cm^2/K^4/s]
aa = 7.5657e-15 # 4 * self.ss / self.cc
#
# Gas constants
#
muh2 = 2.3000e0   # Mean molecular weight (H2 + He + metals)
#
# Alternative units
#
ev   = 1.6022e-12  # Electronvolt [erg]
kev  = 1.6022e-9   # Kilo electronvolt [erg]
micr = 1e-4        # Micron [cm]
km   = 1e5         # Kilometer [cm]
angs = 1e-8        # Angstrom [cm]
#
# Astronomical constants
#

ls   = 3.8525e33   # Solar luminosity [erg/s]
rs   = 6.96e10     # Solar radius [cm]
ms   = 1.99e33     # Solar mass [g]
ts   = 5.780e3     # Solar effective temperature [K]
au   = 1.496e13    # Astronomical unit [cm]
pc   = 3.08572e18  # Parsec [cm]
mea  = 5.9736e27   # Mass of Earth [g]
rea  = 6.375e8     # Equatorila radius of Earth [cm]
mmo  = 7.347e25    # Mass of Moon [g]
rmo  = 1.738e8     # Radius of Moon [cm]
dmo  = 3.844e10    # Distance earth-moon (center-to-center) [cm]
mju  = 1.899e30    # Mass of Jupiter
rju  = 7.1492e9    # Equatorial radius of Jupiter [cm]
dju  = 7.78412e13  # Distance Jupiter-Sun [cm]
#
# Time units
#

year = 3.1536e7    # Year [s]
hour = 3.6000e3    # Hour [s]
day  = 8.6400e4    # Day  [s]

