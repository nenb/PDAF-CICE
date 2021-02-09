
User Documentation - Coming Soon!
=================================

---------------------------------

</br>
</br>

### CICE-PDAF - Development Notes

### 1. PDAF manual for ‘online mode’ implementation:  

http://pdaf.awi.de/trac/wiki/ImplementationGuide 

### 2. Building CICE-PDAF 

The coupled CICE-PDAF model can be compiled and run using the original CICE build and run scripts. Some extra environment variables have been added to the ‘comp_ice’ and ‘clean_ice’ scripts. These variables are used in building the PDAF library, selecting the appropriate makefile macro and in building the CICE-PDAF model.  

A new makefile macro has been created (‘... .pdaf’), which is very similar to the previous makefile macro (‘... .racc’). The main differences are that (i) the mathematical libraries (lapack,…) and header files required for PDAF have been included, as well as the PDAF library itself*, (ii) all real variables are now treated with double precision accuracy by default, using the –r8 compiler flag and (iii) a new preprocessor directive has been included (USE_PDAF). By removing this directive prior to the build, the original CICE model can be recovered.  

The CICE makefile has been modified. The major change is that a rule has been introduced for building/cleaning the PDAF library. Building the PDAF library occurs before building CICE-PDAF. 

PDAF also includes its own makefile and makefile macro. These are located in the subdirectories pdaf/src and pdaf/make.arch. These largely follow the guidelines from the PDAF website. Probably the most relevant sections for the CICE project are (i) the preprocessor directives and (ii) how the compiler treats ‘real’ variables. Regarding treatment of real variables, PDAF does not declare the ‘KIND’ specification in the source code, for flexibility purposes. Thus, when compiling, one has to tell the compiler how to treat real variables. At the moment, the makefile has been written to treat all real variables with double precision accuracy.  

In the run directory (rundir), the only change is that a namelist has been added (namelist.pdaf). More details on this namelist are discussed below. The most important thing to note at this stage is that PDAF has been **attached to CICE assuming that CICE will be run in serial mode. This means that PDAF-CICE should not be run using the parallel (MPI) CICE option as the PDAF assimilation step will almost certainly be performed incorrectly** when CICE is run this way.  

*The PDAF library contains all the necessary routines for performing the assimilation step.

### 3.a ‘Fully-parallel’ implementation of PDAF 

PDAF has been attached to CICE using the fully-parallel implementation. See the following schematic from the PDAF website: 

![PDAF_Schematic](/PDAF_Schematic.png)

### 3.b Overview of calls from CICE to PDAF (PDAF amendments denoted by (P) and italics) 

Changes have been made to the ‘driver’ program CICE.F90 and subroutine CICE_Run_Mod.F90, both located in the ‘drivers’ subdirectory. The changes can be summarised schematically as follows (these changes should be directly compared against the schematic from the previous section on the fully-parallel implementation of PDAF):  

<ins>**CICE.F90**</ins>:

(P) _**CALL** init_parallel()_ 

(P) _**CALL** init_parallel_pdaf()_                                    

**CALL** CICE_Initialize()

(P) _**CALL** init_pdaf()_  

...

**CALL** CICE_Finalize()

(P) _**CALL** finalize_pdaf()_

<br>

<ins>**CICE_Run_Mod.F90**</ins>: 
<br> ...

**CALL** ice_step 

(P) _**CALL** assimilate_pdaf()_ 

istep = istep + 1 

...

### 4 Description of calls from CICE to PDAF 

PDAF is designed in such a way that (in theory) the only required changes to the original model source code are those described in section 3b. The remaining work then takes place inside separate user-supplied routines. All such user-supplied routines live in the subdirectory pdaf/modelbindings. This subdirectory is where the bulk of the coding work for the CICE-PDAF project takes place.  The following is a brief description of most routines in this subdirectory, *excluding those required for the analysis step*.

<ins> **1. init_parallel() and init_parallel_pdaf()**</ins> 

The PDAF documentation for these routines can be found here: http://pdaf.awi.de/trac/wiki/AdaptParallelization . 

A crucial assumption for the current implementation is that PDAF is being attached to a single program (i.e. we are using the SPMD paradigm). If additional programs are coupled at a later stage (e.g., NEMO, CESM coupler) then this step of the attachment process will likely be invalid – **do not use CICE-PDAF in its current form with any external program.**  

The parallel step is necessary to create the N parallel ensemble members. Each ensemble member is referred to as a model task in PDAF jargon.   

The subroutine init_parallel() initialises the MPI environment.  Any attempt to initialise the MPI environment before this step (e.g. if CICE is called from an external program) will likely result in an error (see above warning about using CICE-PDAF with an external program). This subroutine is called outside the init_parallel_pdaf() subroutine so as to act as an explicit reminder about where the MPI environment should be initialised for correct use with PDAF. 

The subroutine init_parallel_pdaf() is responsible for creating the different MPI communicators. The best reference for the purpose of each of these communicators is the PDAF website documentation outlined in Section 1. The modifications to this routine that are specific to CICE-PDAF are as follows: i) The number of ensemble members (‘tasks’) is read from a namelist file (namelist.pdaf) and ii) CICE is assumed to operate only in serial mode (see Section 2 for further details). Thus, there will be N MPI communicators, and each communicator will represent a single copy of CICE.  

<ins> **2. init_pdaf()** </ins> 

The documentation for this step can be found here: http://pdaf.awi.de/trac/wiki/InitPdaf 

After the adaption of the parallelization, the initialization of PDAF has to be implemented. The various parameter values are read from a namelist file (namelist.pdaf) which lives in the ‘rundir’ subdirectory. The routine read_config_pdaf  is used to read the namelist file. The routine init_pdaf_info is used to print relevant information from the initialization process to STDOUT. A slight subtlety when running N ensemble members in parallel is that when each member attempts to print information, the information becomes unreadable. Hence it is common to only have one member print output (e.g. model member 1). 

The initial conditions for the N ensemble members are read by the routine init_ens_pdaf() called from within init_pdaf. The size of the state vector must be computed prior to calling ‘init_ens_pdaf’. This computation is performed by the routine calc_statevector_dim(). 

As CICE is a forced, dissipative model (i.e. it is not a chaotic model) it is necessary to supply each ensemble member with different boundary conditions, so as to preserve sufficient spread between the different members. To do this, the approach adopted here is to overwrite the CICE variable 'fyear' with a different value for each ensemble member. 

A routine next_observation_pdaf() is necessary to tell PDAF how often to perform the update step. For CICE-PDAF, next_observation_pdaf() uses the CICE variable ‘npt’ to determine how many timesteps have been performed, and how many timesteps remain. 

A routine prepoststep_pdaf() is necessary for output of data, both pre and post-analysis step. It also includes an option to output intial conditions.

A module mod_statevector has been written that is entirely specific to CICE. It contains all necessary information and routines for producing the statevector for CICE-PDAF. The module was originally written following an extremely crude format. Over time, it should be possible to considerably improve the readability of this module. It is important to note that this module contains a preprocessor directive ‘USE_STRESS’. If this directive is included, then the state variables ‘a11_1, …, a12_4’ will be included in the statevector. Otherwise, these state variables will not be included in the statevector. 

<ins> **3. assimilate_pdaf()** </ins>

PDAF-OMI v1.2 is used for performing the analysis step calculations in PDAF-CICE. Please see the git log and/or contact a member of the sea-ice team at Reading University for further details.

--------
