import sys
sys.path.insert(1, '../src/magpies/')
from magpies import *
from atmos import *
import numpy as np
from math import *
from scipy.optimize import minimize
import matplotlib as mpl
import matplotlib.pyplot as plt

## Radius and mass of neutron star
Rns = 12  ## km
Mns = 1.4 ## M_solar
Tb = pow(10, 8.3648)  ## K
Bp = 1e14 ## G

## Instrument and exposure
eph = np.linspace (0.20, 3.00, 32) ## Nicer soft X-ray range
nphot = 1e5

g14c = g14 (Rns, Mns) ## computing the free fall acceleration


atm_iron_2003 = NS_atmosphere ('Potekhin_2003_iron', g14c, Tb, Bp)

theta = np.linspace (0, pi, 100)  ## theta coordinates
phi   = np.linspace (0, 2*pi, 99) ## phi coordinates

theta1, phi1 = np.meshgrid (theta, phi)

Ts = atm_iron_2003.Ts (theta1) ## Surface temperatures

L = compute_L(theta, phi, Rns, Ts)
Teff = compute_Teff(theta, phi, Rns, Ts)

Ts_uniform = np.full((Ts.shape[0], Ts.shape[1]), Teff)
Teff1 = compute_Teff(theta, phi, Rns, Ts_uniform)
L1 = compute_L(theta, phi, Rns, Ts_uniform)

#print (abs(L - L1) / abs(L))
#print (abs(Teff - Teff1) / abs(Teff))

assert (abs(L - L1) / abs(L) < 0.03)
assert (abs(Teff - Teff1) / abs(Teff) < 0.03)

print ('Tests are successful')
