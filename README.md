0D Time Dependent Bltz Solver
===============
*Author: J. Claustre*

*Date  : 12/2015*

Description
-----------

This code solve the 0D Bltz equation, depending on time for the
electrons in the _pure_ helium gas. It includes *electron-electron*
collisions, *electron-neutral* collision, *heating* by electric field.

Furthermore, it takes into account *inelastic* and *superelastic*
collisions, *ionization* (direct, non-direct, Penning and
Associative), *dissociative recombination*, *l-exchange* and
*s-exchange* processes, *ionic conversion*, *radiative transfert*,
*diffusion* and so on.

The geometry concidered here is a cylinder (the geometry consideration
is used for the absorbed power calculation and the ambipolar diffusion).

Have a look at our papers, to have a detailled description of the model:

* J Claustre *et al.* (doi:10.1088/1361-6595/aa8a16)
  [(link article)](https://doi.org/10.1088/1361-6595/aa8a16)
* C Boukandou *et al.* (doi:10.1016/j.cpc.2017.07.004)
  [(link paper)](https://doi.org/10.1016/j.cpc.2017.07.004)

We invite you to look at these papers too (our source of inspiration!:D ):

* L Alves *et al.* (doi:10.1088/0022-3727/25/12/007)
  [(link article)](http://m.iopscience.iop.org/article/10.1088/0022-3727/25/12/007/meta;jsessionid=AE4353A7414EB307AA0214AD6A4BA223.c3.iopscience.cld.iop.org)
* M Santos *et al.* (doi:10.1088/0022-3727/47/26/265201)
  [(link article)](http://iopscience.iop.org/article/10.1088/0022-3727/47/26/265201#)
* G Hagelaar *et al.* (doi:10.1088/0963-0252/14/4/011)
  [(link article)](http://m.iopscience.iop.org/article/10.1088/0963-0252/14/4/011/meta)
* T Belmonte *et al.* (doi:10.1088/0022-3727/40/23/015)
  [(link article)](http://iopscience.iop.org/article/10.1088/0022-3727/40/23/015/meta)

----------------------------------------------------------------

Install it
-------

1. Download the zip file on your machine. (button on the right of the
   web page)
2. unzip the file

### On Linux
* To unzip with the terminal mode :

		unzip file.zip -d destination_folder
		
* In the folder, a makefile is included. Modify it with your own
        compiler (gcc by default) and compile it by using the terminal
        consol :

		make makedirectories
		make

* /!\ In the same folder ("datFile"), you should have a file called
"input_he". This file is needed to define (and read) the simulation
parameters. If you don't have it, create it by using the example below
:

----------------------------------------

*EXAMPLE OF INPUT FILE : (name_of_file =* input_he*)*

----------------------------------------
	
	2000              Number of grid points 
	0                 Restart simulation ( 0  - no restart | 1 - restart )
	2.320E+01         E (V/cm)      Electric field (see also the absorbed power parameter)
	4.000E+18   0     Ng (cm-3) ( 1 - input parameter neutral gas (Ng) | 0 - input parameter Pressure (prs))
	3.000E+06   0     f (Hz)        Heating Frequency (If 1 ==> "RF mode" Else "HF mode" )
	1.0000000         prs (Torr)    Pressure
	3.000E+02  600    Tp (K)        Temperature [+ Tp at tube boundary (if needed)]
	2.000E+00         Tpe (eV)      Electronic Temperature
	2.500E+00         Ra (cm)       Cylinder Radius
	5.000E+00         L  (cm)       Tube Length
	1.000E+03         ne (cm-3)     Initial electron density
	1.000E+00   0     P  (W/cm3)    Absorbed power ( 1 - input parameter E | 0 - input parameter P)
	1000.0001         Simulation Time (micro-sec)
	1.000E-09         Dt ( Time-Step (sec) ) ==> Max time-step allowed
	3.000E+02         Emax (Max Energy [for the grid] (eV) ) 
	5.000E+01         save data's every XXX micro-sec (includes Restart files + EEDF)
	*************************************************************************
	0                 Activate the laser (1) or don't (0) (if 0, don't care about the following lines!)
	1083              lenght wave of the laser (nm)
	1                 Polarization of the Laser (0=pi ; 1=sigma+ ; 2=sigma-)
	1.000E+00         Laser intensity (W)
	1.963E+01         Surface or laser section (cm2)
	5.000E+02         Start time of the Laser (micro-sec)
	1                 Numbers of transitions to use with the laser
	9                 Write the transitions to use (from 1 to 9) (write on the same line, ex: 1 5 7 9)

----------------------------------------
### On Windows
* Sorry, i really don't know! :X

Run it
-------

After the code compiled perfectly, you should have an executale
*run_BOD*. In a Linux terminal, wrote the command

	./run_BOD

Output Files and results
-------------------------

During the calculation and after the code end, you can check results
and outputs in the file directory you choose at the beginning of the
simulation (see *param.f90* file).

Want to Contribute (in GitHub)?
--------------------

1. Fork it
2. Create your own feature branch (`git checkout -b my-own-new-feature`)
3. Commit your changes (`git commit -am "Add some changes"`)
4. Push to the branch (`git push origin my-own-new-feature`)
5. Create new Pull Request

Routine descriptions
---------------

* Evolution.f90
	* Contains the main loop (time evolution). This is the loop in
	time where all reactions concidered and heating, in the code, are called.
	* In this loop, time is updated but also some densities (charges,
      excited states).
	* The time step is adapted to (mainly) the strongest collision
     rate between inelastic-superelastic collisions and ionization
     processes.
	
* Excit.f90
	* Contains the excitation and de-excitation processes between all
      excited states.
	* The subroutine "equil" check if the collision rate is not to
	high. If yes, it calculates an equilibrium solution between the
	both excited states. This permit to keep a time step not too
	small.
	* The subroutine "implic" calculates in a semi-implicit form, the
	densities at time k+1.  It allows to have a greater time-step
	without calculating an equilibrium solution between some states.
	* We advise to take care to have a "reasonable" time-step (*i.e*
      not too high!).

* Ioniz.f90
	* Contains ionization (direct and undirect) processes. In the case
	of high ground state density.
	* In order to have a greater time-step, we perform subcycles (only
	for the fundamental state).

* Penn-asso.f90
	* Contains Penning ionization and associative ionization.
	* Contains also the dimer Penning processes.
	* In the init subroutine, you have several rates. You can choose
      which one to use in the *read_input.f90* (look for "penning").
	
* L-xchange.f90
	* Contains l-exchange and s-exchange processes for He.

* Recomb.f90
	* Contains the dissociative recombination processes.
	* Two routines are available here :
		* The first one depends on a fixed collision rate and the EEDF
		(Energy Electron Distribution Function) is renormalized.
		* The second one consider a given cross-section and the EEDF
		is calculated in a kinetic way.
	* Contains also the 3 body ionic conversion process.
	
* Radiff.f90
	* Contains the radiative transfert process and the loss of
     particles due to ambipolar diffusion.
	* Two methods are concidered for the diffusion process.
	
* Heat.f90
	* Contains heat by electric field and depends of the input power
	defined (= absorbed power).
	* Contains electron-electron collisions by solving the
	Fokker-Planck equation (full conservative scheme).
	* Contains electron-neutral collisions.
	
* Pumping.f90
	* Contains the optical pumping part and the distribution of
      populations on the sublevels of the excited states 2^3^P and
      2^3^S (*e.g.* radiative and metastable states) without magnetic
      field.
	  * The polarization of Helium is also considered but it is done
        is two steps.
		  * The calculation of a~i~ populations for several
            polarization values (-> a~i~(P) ).
		  * The calculation of the polarization with the a~i~ values in
            an other files (*i.e.* python program).

* Gaz_Tp.f90
	* Contains the routine to calculate the 1D gas temperature
	evolution inside of the cylinder, considering:
		* The electron density [Ne] with a Bessel shape in the transverse direction.
		* The gas temperature at the boundary of the tube.
		* The conductivity
		
* Read-input.f90
	* Contains initialization of all parameters, cross-sections, EEDF,
	charge and excited state densities, read input files, ..., etc.
	
* Param.f90
	* Contains the definition of structures and types used in the
      code.
	* Contains all constants and physical parameters used (i.e qe, me,
	  kB, ...,etc) in the code.
	* Contains routines for dynamic allocation and deallocation.
	
* Main.f90 :
	* The top of the pyramid :)

Tips
---------------

* Outputs
	* In the file param.f90, you can choose the directories for your
	  outputs. If the directory doesn't exist, it will create it.
	* If you want to create manually your directory, _you_ _have_ _to_
      _create_ the folder *Rstart* inside of your new directory!
	* Be sure to *make clean* (and *make*) before to run the code to take into
      account of any changes in *param.f90*.

* Cross Sections
	* Several cross sections are available. We found some differences
      between the cross section given by the team of Alves (Lisboa)
      and Biagi. We advise to use the Biagi cross section if you take
      into account the _Associative_ _Ionization_ processes, else the
      Alves cross section.
	* You can change from one to the other in the file *excit.f90*, in
      the subroutine *Init_ExcDxc*. you'll find the variable *Switch_CS*.

			if (Switch_CS == 1) Then
				Alves Cross section
			else
				Biagi cross section
			end if

* Gas temperature
	* If you want to calculate the gas temperature, just uncomment the
	  call in the main loop (in *evolution.f90*).

------------------------------------------------------------------------------------------------------------
