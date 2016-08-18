'''
Description: Implementation of GHOST model. Five functions provided which 
can be used to create a halo-galaxy catalog. Sample code provided at the
end for getting started up with a halo catalog containing over 2,500 dark 
matter halos at redshift zero. Requires SciPy, NumPy, and CosmoloPy. 

'''
import numpy as np
from numpy import exp
from numpy import log
import cosmolopy.distance as cd
import scipy.interpolate
import cosmolopy.magnitudes as mg
import cosmolopy.perturbation as cp
from scipy.integrate import simps

############################################################################

# input parameter values

# parameters for computing mean stellar mass as a function of halo mass
ms_beta0 = 1.4
ms_beta1 = 0.17
ms_gamma0 = 0.65
ms_gamma1 = -0.26
M1 = 10**12.1
M2 = 10**11.8
mu0 = 0.019
nu = -0.72
sig_inf = 0.16
sig_1 = 0.05
xi = 4.25

sfms_sigma = 0.4
pass_sigma = 0.4
alpha_fpass0 = 10.8
alpha_fpass1 = 0.5
beta_fpass = -1.1
zeta_fpass = -0.3
c = 0.5*(1+np.tanh(zeta_fpass))

sfr_sfms_alpha0 = 0.005
sfr_sfms_beta = 0.9
sfr_sfms_alpha1 = 1.07

m_shift = 0.03

h_alpha_const = 1.266e41
radio_const = 1.812e28

cosmo = {
	'omega_M_0' : 0.27, 
	'omega_lambda_0' : .73, 
	'omega_b_0' : 0.044,
	'omega_n_0' : 0.0, 
	'N_nu' : 0, 
	'h' : 0.7, 
	'n' : 0.95, 
	'sigma_8' : 0.8
	}
	
cd.set_omega_k_0(cosmo)

############################################################################

def stellar_mass(z, mvir):
	"""
	Function returns an array of predicted stellar masses given an array 
	of halo masses and redshifts.
	
	Parameters
	----------
	z : numpy array of floats
		redshifts of the dark matter halos
	mvir: numpy array of floats
		masses of dark matter halos in M_sun/h
	
	Returns
	-------
	numpy array of stellar masses in M_sun/h
	"""
	
	# redshift-dependent parameters
	gamma = ms_gamma0*(1. + z)**ms_gamma1 
	beta = ms_beta0 + ms_beta1*z 
	mu = mu0*(1. + z)**nu
	
	# mean stellar mass given by broken power law
	stellar_mass_mean = 2.*mvir*mu/(((mvir/M1)**-beta)+((mvir/M1)**gamma))
	dev = sig_inf+sig_1*(1-((2./np.pi)*np.arctan(xi*np.log10(mvir/M2))))
	
	# generate array of stellar masses according to a lognormal distribution
	return np.random.lognormal(np.log(stellar_mass_mean), dev)


def galaxy_type_and_SFR(z, s_mass):
	"""
	Function takes redshift and predicted stellar mass arrays of a list of galaxies
	and returns each galaxy's predicted type (active or passive) and star 
	formation rate. An active type galaxy is denoted by the boolean value True, 
	while a passive type is denoted by the boolean value False. 
	
	Parameters
	----------
	z : numpy array of floats
		redshifts of the dark matter halos
	s_mass: numpy array of floats
		masses of dark matter halos in M_sun
	
	Returns
	-------
	boolean array of values describing galaxy type, 
	float array of star formation rates (M_sun/yr)
	"""
	
	# redshift-dependent parameters
	alpha = alpha_fpass0 + alpha_fpass1*z
	sfr_alpha = sfr_sfms_alpha0+ z*sfr_sfms_alpha1
	
	# fraction of galaxies at a certain stellar mass that are passive
	f_passive = c + (1. - c)/(1. + (s_mass/(10.**alpha))**beta_fpass)
	
	# select type according to a uniform distribution
	uni_random = np.random.uniform(size = len(f_passive)) 
	type =  np.greater(uni_random, f_passive)
	
	# create array of predicted mean star formation rates
	mean_sfr = (10.**sfr_alpha)*((s_mass/1e10)**sfr_sfms_beta)
	
	# for passive galaxies, multiple mean_sfr by some factor m_shift
	pass_indx = np.where((type == False))
	mean_sfr[pass_indx] = m_shift * mean_sfr[pass_indx]
	
	# assign appropriate sigma to each galaxy depending on type
	sigma = np.full(len(mean_sfr), sfms_sigma)
	sigma[pass_indx] = pass_sigma
	
	# generate array of sfr according to a lognormal distribution
	SFR_random = np.random.lognormal(np.log(mean_sfr), sigma)
	
	return type, SFR_random


def lum_to_flux(z, luminosity):
	"""
	Function takes an array of redshifts and intrinsic luminosities and
	computes the luminous flux (erg/s/cm^2) of each galaxy in the 
	specified frequency of emission.
	
	Parameters
	----------
	z : numpy array of floats
		redshifts of the dark matter halos
	luminosity: numpy array of floats
		intrinsic luminosity value of each dark matter in erg/s
	
	Returns
	-------
	float array of luminous fluxes in erg/s/cm^2
	"""
	
	#Import file containing table of redshifts and corresponding
	#comoving distances in Mpc.
	
	table = np.genfromtxt("dist_z_vals.txt").T
	redshifts = table[0]
	co_dist = table[1]

	d_interp = scipy.interpolate.interp1d(redshifts, co_dist, kind='linear')

	selected_co_dist = d_interp(z)
	
	# calculate luminosity distance (d_L) in Mpc
	lum_dist = (1 + z) * selected_co_dist

	#convert from Mpc to cm
	mpc_to_cm = 3.0856776e24
	lum_dist_cm = lum_dist * mpc_to_cm
	
	# consider cases where redshift is zero
	divisor = (4.*np.pi*lum_dist_cm**2.)
	zero_divide = np.where(lum_dist_cm == 0)
	
	divisor[zero_divide] = 1.

	return luminosity/divisor

	
def comoving_volume_of_catalog(z_min, z_max, sq_deg_of_sky):
	"""
	Computes the co-moving volume of the redshift space with
	a specified field of view.
	
	Parameters
	----------
	z_min, z_max : float
		Min/max redshift of dark matter halos
	sq_deg_: float
		Field of view in square degrees
	
	Returns
	-------
	float number describing the volume
	"""
	
	v_min = cd.comoving_volume(z_min, **cosmo)
	v_max = cd.comoving_volume(z_max, **cosmo)
	
	shell_volume = (v_max - v_min)/(.72**3)

	deg_to_rad = (np.pi*np.pi)/(180.*180.)

	#fraction of sky
	f_sky = (sq_deg_of_sky * deg_to_rad)/(4. * np.pi)

	return shell_volume * f_sky

def radio_lum(SFR):
	"""
	Returns the predicted intrinsic luminosity values in the radio frequency
	
	Parameters
	----------
	SFR: numpy array of floats
		star formation rate value for each galaxy
	
	Returns
	-------
	numpy array of floats describing the intrinsic luminosity in erg/s
	
	"""
	
	return radio_const*SFR
	
def h_alpha_lum(SFR):
	"""
	Returns the predicted intrinsic luminosity values in the H-alpha spectral band
	
	Parameters
	----------
	SFR: numpy array of floats
		star formation rate value for each galaxy
	
	Returns
	-------
	numpy array of floats describing the intrinsic luminosity in erg/s
	
	"""
	
	return h_alpha_const*SFR
	
def mass_function(z, min, max, n_bins):
	"""
	Returns a redshift-dependent halo mass function number density
	
	Parameters
	----------
	z: numpy array of floats
		array of redshift values
	min, max: float
		min/max halo masses
	n_bins: integer
		number of bins for halo mass function calculation
		
	Returns
	-------
	float numpy array of masses
	float numpy array of number densities
	
	"""
	
	# choose a delta vir
	delta_vir = 200.

	# define variables used in the fitting function (A, b, a, c)
	big_a = 0.186*((1 + z)**-0.14)
	little_a = 1.47*((1 + z)**-0.06)
	c = 1.19
	alpha = exp(-(0.75/np.log(delta_vir/75.))**1.2)
	b = 2.57*((1 + z)**-alpha)

	# high mass halo range 10^10 - 10^15 m_sun/h
	min_mass = min
	max_mass = max
    	
	# parameters used in HMF calc
	cosmo = {
		'omega_M_0' : 0.27, 
		'omega_lambda_0' : 1-0.27, 
		'omega_b_0' : 0.0455,
		'omega_n_0' : 0.0, 
		'N_nu' : 0, 
		'h' : 0.704, 
		'n' : 0.967, 
		'sigma_8' : 0.81
		}

	r_little = cp.mass_to_radius(min_mass, **cosmo)
	r_big = cp.mass_to_radius(max_mass, **cosmo)

	k = np.logspace(-4., 2., 1e4)

	radii = np.logspace(np.log10(r_little), np.log10(r_big), n_bins)

	mass = cp.radius_to_mass(radii, **cosmo)

	
	# mean density of the universe in solar mass/ Mpc units
	mean_dens = (3*min_mass)/(4*np.pi*(r_little**3))

	list_of_root_vari = []
	for item in radii:
	
		#top hat window function used in HMFcalc
		r = item
		window_sq = lambda k: ((3.*np.sin(r*k) - 3.*r*k*np.cos(k*r)) / (r*k)**3.)**2.

		# cosmolopy power spectrum function
		pow = cp.power_spectrum(k, z, **cosmo)

		integrand =(1./(2.*np.pi*np.pi))*(k**2.)*pow*window_sq(k)

		vari = simps(integrand, k)
	
		list_of_root_vari.append(vari**0.5)

	# convert list of deviances to an np array
	dev = np.array(list_of_root_vari)

	# convert to solar mass/h
	mass = mass
	
	dev_inverse = 1./dev
	ln_dev = np.log(dev_inverse)

	# acquiring d(ln(1/dev)/d(mass) term
	slope = []

	# numerical differentiation
	ln_m = np.log(mass)
	for i in range(len(mass)):
		if (i != len(mass)-1):
			m = (ln_dev[i+1] - ln_dev[i])/(ln_m[i+1] - ln_m[i])
			#append the absolute value of the slope
			slope.append(abs(m))
		else: 
			m = (ln_dev[i] - ln_dev[i-1])/(ln_m[i] - ln_m[i-1])
			# append the absolute value of the slope
			slope.append(abs(m))

	# insert root variance values into the function of collapsed densities
	fcoll = []
	for item in dev:
		f_coll = lambda item: big_a*((b/item)**(little_a) + 1)*exp(-c/(item**2))
		fcoll.append(f_coll(item))
	
	# calculate number density of halos as a function of mass
	number_dens = (mean_dens / mass) * fcoll * slope
	return mass, number_dens*np.log(10.)
