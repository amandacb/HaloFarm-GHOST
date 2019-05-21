# HaloFarm-GHOST
Python implementation of a galaxy-halo model developed by Phil Bull.

===========================================================
GHOST Model Implementation/HaloFarm Code
              ---
         Amanda Brown
   aclairebrown@hotmail.com
         August 2016
===========================================================

Statistically modeling the relationship between dark matter halos and the
galaxies that they contain.

Released under the MIT License.

ABOUT GHOST AND HALOFARM
-----------------

GHOST is an analytic model developed by Philip Bull that that forecasts the 
observable properties of galaxies given the mass and redshift of the dark 
matter halos that contain them. I have written a Python script to
implement this model. It takes a text file containing the mass and 
redshift of a number of dark matter halos and tabulates forecasted properties
of the galaxy contained in each. I also include a script, called HaloFarm, that 
generates a dark matter halo catalog that can be used as input an input file
when testing the GHOST script. 

The halo catalog has been generated using the Press-Schechter halo mass 
function described in Murray et al. 2013 and Tinker et al. 2008. This formalism
provides the redshift-dependent number density across a number of mass bins. 
After multiplying this value by the co-moving volume each the redshift partition,
we ascertain the number of halos expected to be found in each given 
redshift slice as a function of mass. According to this number, I draw from 
a uniform random distribution to select halos from each mass bin.  

The code can be used to forecast observations in the H-alpha and radio 
frequencies that will be made with upcoming sky-mapping surveys like SKA
and WFIRST.

REQUIREMENTS
------------

 - Python 2.7
 - NumPy, SciPy, CosmoloPy
 - matplotlib


INSTALLATION
------------

No installations required. The Python scripts can be copy-and-pasted into 
a text editor and used directly.

GETTING STARTED
---------------

To begin, install Python 2.7 and the latest versions of SciPy, NumPy, matplotlib,
and CosmoloPy. Next, select a halo catalog text file to be read in as input 
(or create a catalog directly using the HaloFarm script). If using HaloFarm, 
specify your desired redshift range and field of view (in degrees). If using 
your own halo catalog, make sure that your file contains a row for each halo 
with its mass (in solar mass/h) and redshift. 


WORKING WITH THE CODE
---------------------
I have provided a sample piece of code that describes how to read-in your halo 
catalog, use functions to make predictions about each halo's galaxy, and then 
output your halo-galaxy catalog to a file. 

 * halofarm.py: You can create your own halo catalog that distributes
 		halo masses according to the redshift-dependent
 		Press-Schechter Halo Mass Function described in Tinker
 		et al. 2008 and Murray et al. 2013.
 
 * ghost_model.py: Contains functions for implementing Phil Bull's GHOST model.

 * example_halo_cat.py: Get started right away using the GHOST model. In 
 			this script, a halo catalog has been provided containing
 			over 8,000 dark matter halos at redshift 0.


CITING THIS CODE
----------------

If you wish to use this code in a scientific paper, please cite Philip Bull, 
who developed the GHOST model, and Amanda Brown, who has helped to implement
it.

BUGS
----

Email me at aclairebrown@hotmail.com to let me know if there are problems 
with the code, or if you would like more information about the formalism
or implementation of the model.
