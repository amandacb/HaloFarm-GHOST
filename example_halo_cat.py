'''
Description: Provides a sample halo catalog of dark matter halos
at redshift zero, and calls functions from ghost_model.py to generate a 
halo galaxy catalog. Requires NumPy, SciPy, CosmoloPy, and the functions
from ghost_model.py

'''

import matplotlib.pyplot as plt
import numpy as np
from numpy import exp
from numpy import log
import cosmolopy.distance as cd
import scipy.interpolate
import cosmolopy.magnitudes as mg
from ghost_model import stellar_mass, galaxy_type_and_SFR
from ghost_model import comoving_volume_of_catalog, lum_to_flux
from ghost_model import radio_lum, h_alpha_lum

halos = np.genfromtxt("halo_cat_z_0.csv", delimiter=',').T

mass = halos[0]
z = halos[1]

print "catalog size", len(mass)

s_mass = stellar_mass(z, mass)

type, SFR = galaxy_type_and_SFR(z, s_mass)

radio_lum = radio_lum(SFR)
h_alpha_lum = h_alpha_lum(SFR)

h_alpha_flux = lum_to_flux(z, h_alpha_lum)
radio_flux = lum_to_flux(z, radio_lum)

halo_galaxy_catalog = np.column_stack((z, mass, 
							s_mass, type, SFR,
							h_alpha_lum, radio_lum, 
							h_alpha_flux, radio_flux
							))
							
np.save("save_halo_galaxy_cat", halo_galaxy_catalog)

'''
Code for plotting different relationships of the halo-galaxy catalog.

'''

fig1 = plt.figure(1)
plt.scatter(mass, s_mass, s=5, c='k')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Dark Matter Halo Mass $M_h$($M_\odot$$h^{-1}$)", fontsize = 12)
plt.ylabel("Stellar Mass $M_*$ ($M_\odot$/h)", fontsize = 12)
plt.title("Halo Mass $M_h$ vs. Stellar Mass $M_*$", fontsize = 12)

fig2 = plt.figure(2)
plt.scatter(s_mass, SFR, s=5, c='k')
plt.xscale("log")
plt.yscale("log")
plt.ylabel("$\psi$ SFR ($M_\odot$$yr^{-1}$)")
plt.xlabel("Stellar Mass $M_*$ ($M_\odot$$h^{-1}$)")
plt.title("$\psi$ SFR vs. $M_*$")
plt.ylim(10**-6, 10**2)

plt.show()