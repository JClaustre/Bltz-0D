!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 11/2015
! Objctv: README
! note  : Description of B0D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!**** Modification from previous version ****!
___________________________________ 24-11-2015

o--> In Excit.f90
     o Correction of the loop by removing "GoTo 1000"
     o Add rate calculation to compute the adaptative time-step
     o Two routine added : 
       	   o- Implicit calculation of densities and correction of EEDF
	   o- Equilibrium state for reactions with high rates 
	      (most probably for high energy level)

o--> In Ioniz.f90
     o Add subcycles for the ground state atom because of the high rate
     o Add rate calculation to compute the adaptative time-step

o--> In Penn-asso.f90
     o Add rate calculation to compute the adaptative time-step

o--> In evolution.f90
     o Add computation of new time-step.
     o Add MaxDt to define a max time step
_______________________________________________
