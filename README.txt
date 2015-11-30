!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 11/2015
! Objctv: README
! note  : Description of B0D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

This code solve the 0D Boltzmann equation, depending on time for the Helium gas.
It includes electron-electron collisions, electron-neutral collision, heating by electric field.
Furthermore, it takes into account inelastic and super-elastic collisions, ionization 
(direct, non-direct, Penning and Associative), dissociative recombination, 
l-exchange and s-exchange processes, ionic conversion, radiative transfert and diffusion.
The geometry concidered here is a cylinder (used for the absorbed power and the diffusion coefficient).

We invite you to look at these papers where most of the rates and cross-sections 
used in the code come from (our source of inspiration! :D ):

----------------------------------------------------------------!
   -o Luis Alves et al (doi:10.1088/0022-3727/25/12/007)	!
   -o M Santos et al (doi:10.1088/0022-3727/47/26/265201)	!
   -o G Hagelaar et al (doi:10.1088/0963-0252/14/4/011)		!
   -o T Belmonte et al (doi:10.1088/0022-3727/40/23/015)	!
----------------------------------------------------------------!

MAIN ROUTINES :
---------------

-o evolution.f90 :
	Contains the main loop. This is the loop in time where all reactions concidered in the code 
	are called.
	In this loop, time is updated but also some densities (charges, excited states). 
	The time step is adapted to the strongest collision rate between inelastic-superelastic collisions
	and ionization processes.
	
-o excit.f90 :
	Contains the excitation and de-excitation processes between all excited states.
	The subroutine "equil" check if the collision rate is not to high. If yes, it calculates an 
	equilibrium solution between the both excited states. This permit to keep a time step not too small.
	The subroutine "implic" calculates in a semi-implicit form, the densities at time k+1. 
	It allows to have a greater time-step without calculating an equilibrium solution between 
	some states. We advise to take care to have a "reasonable" time-step (i mean not too high!).

-o ioniz.f90 :
	Contains ionization (direct and undirect) processes. In the case of high ground state density, 
	in order to have a greater time-step, we perform subcycles (only for the fundamental state).

-o penn-asso.f90 :
	Contains Penning ionization and associative ionization.
	
-o l-xchange.f90 :
	Contains l-exchange and s-exchange processes for He.

-o recomb.f90 :
	Contains the dissociative recombination processes. Two routines are available here. 
	The first one depends on a fixed collision rate and the EEDF (Energy Electron Distribution Function)
	is renormalized. The second one consider a given cross-section and the EEDF is calculated 
	in a kinetic way.
	Contains also the 3 body ionic conversion process.
	
-o radiff.f90 :
	Contains the radiative transfert process and the loss of particles due to ambipolar diffusion.
	Two methods are concidered for the diffusion process.
	
-o heat.f90 :
	Contains heat by electric field and depends of the input power defined (== absorbed power).
	Contains electron-electron collisions by solving the Fokker-Planck equation 
	(full conservative scheme). 
	Contains electron-neutral collisions.
	
-o read-input.f90 :
	Contains initialization of all parameters, cross-sections, EEDF, charge and excited state densities,
	read input files, ..., etc.
	
-o param.f90 :
	Contains the definition of structures and types used in the code. 
	Contains all constants and physical parameters used (i.e qe, me, kB, ...,etc) in the code.
	Contains routines for dynamic allocation and deallocation.
	
-o main.f90 : the SOURCE! :)

------------------------------------------------------------------------------------------------------------

/!\ Inside of the folder "datFile", make sure you have a folder named "Rstart", else create it!

/!\ In the same folder ("datFile"), you should have a file called "input_he". 
    This file is needed to define (and read) the simulation parameters. If you don't have it, 
    create it by using the example below :

**** EXAMPLE OF INPUT FILE : (name_of_file = input_he)

500               Number of grid points 
0                 Restart simulation ( 0  - no restart | 1 - restart )
2.13000E+01       E (V/cm)      Electric field
3.00000E+18   0   Ng (cm-3) ( 1 - input parameter N | 0 - input parameter Prs)
2.45000E+09       f (Hz)        Heating Frequency
7.60000E+02       prs (Torr)    Pressure (760 Torr = Atmospheric pressure)
2000.000000       Tp (K)        Temperature
2.10000E+00       Tpe (eV)      Initial electron temperature
1.00000E-01       Ra (cm)       Cylinder Radius
2.30000E+00       L  (cm)       Tube Length
2.34000E+12       ne (cm-3)     Initial electron density
8.30000E+01   0   P  (W/cm3) ( 1 - input parameter E | 0 - input parameter P)
3.000000000       Simulation Time (micro-secnd)  
1.00000E-12       Dt ( Time-Step (secnd) )
1.00000E+02       Emax (Energy grid max (eV) ) (Maximum allowed = 1000 eV)
1.0000            save data's every XXX micro-sec

**** END EXAMPLE

