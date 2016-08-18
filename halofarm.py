"""
Description: Generates a dark matter halo catalog according to the 
Press-Schechter halo mass function described in Murray et al. 2013
and Tinker et al. 2008. Requires NumPy, a function from ghost_model.py
and a function from halofarm_mass_function.py.

"""
import matplotlib.pyplot as plt
import numpy as np
from numpy import log
from ghost_model import mass_function
from ghost_model import comoving_volume_of_catalog

# select a name for your halo catalog
with open("halo_cat_filename.csv", 'w') as file_mass:

	N = 40	# number of number density bins
	
	Z = 10  # redshift partitions
	
	# choose min/max mass of halos in catalog
	min = 10**10
	max = 10**15

	# desired range of redshifts
	z = np.linspace(0.0, 0.01, Z)
	
	for i in range(Z):
		if i < Z-1:
			volume_element = comoving_volume_of_catalog(z[i], z[i+1], 400)
			mass, number_dens = mass_function(z[i], min, max, N)
			how_many = number_dens*volume_element
			
			# deal with the fractional halos using a discrete uniform distribution
			decimal = how_many%1.0
			rands = np.random.uniform(size=len(how_many))
			#returns True if random # is less than decimal; if True, add a halo
			halo_or_not = np.less(rands, decimal)
			indx = np.where(halo_or_not == True)
			how_many[indx] = np.ceil(how_many[indx])
			how_many = np.floor(how_many)
	
			# conversion for h
			mass = mass*0.72
			for j in range(len(mass)-1):
				y = np.random.uniform(low=mass[j], high=mass[j+1], size=how_many[j])
				
				for item in y:
					file_mass.write(str(item)+ "," + str(z[i])+ "\n")
