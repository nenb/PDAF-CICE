!$Id: obs_TYPE_pdafomi_TEMPLATE.F90 464 2020-05-28 10:19:36Z lnerger $
!> PDAF-OMI template observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
!! 
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF.
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs_f)
!! and for the observation operator (obs_op_f).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! The module uses two derived data types (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!! **Using this template:**
!! To be able to distinguish the observation type and the routines in this module,
!! we recommend to rename the module according to the observation module-type.
!! Further,we recommend to replace 'TYPE' in the routine names according to the
!! type of the observation so that they can be identified when calling them from 
!! the call-back routines.
!!
!! These 2 routines usually need to be adapted for the particular observation type:
!! * init_dim_obs_f_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_f_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE obs_ice_thickness_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_filter, abort_parallel    ! Rank of filter process
  USE PDAFomi, &
       ONLY: obs_f, obs_l   ! Declaration of observation data types
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_ice_thickness=.TRUE.        !< Whether to assimilate this data type
  LOGICAL :: twin_experiment=.FALSE.           ! Whether to perform an identical twin experiment
  REAL    :: rms_ice_thickness=0.1      !< Observation error standard deviation (for constant errors)
  REAL    :: noise_amp = 0.1  ! Standard deviation for Gaussian noise in twin experiment

  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.
  LOGICAL            :: obs_file=.FALSE. ! Are observations read from file or manually created
  CHARACTER(len=200) :: file_ice_thickness= &
       '/storage/silver/cpom/fm828007/CICE/cice_r1155_pondsnow/rundir_test/restart/iced.2008-01-01-00000.nc'  ! netcdf file holding observations

! ***********************************************************************
! *** The following two data types are used in PDAFomi                ***
! *** They are declared in PDAFomi and only listed here for reference ***
! ***********************************************************************

! Data type to define the full observations by internally shared variables of the module
!   TYPE obs_f
!           Mandatory variables to be set in init_dim_obs_f
!      INTEGER :: doassim                   ! Whether to assimilate this observation type
!      INTEGER :: disttype                  ! Type of distance computation to use for localization
!      INTEGER :: ncoord                    ! Number of coordinates use for distance computation
!      LOGICAL :: use_global_obs=.true.     ! Whether to use (T) global full obs. 
!                                           ! or (F) obs. restricted to those relevant for a process domain
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! indices of observed field in state vector
!           
!           Optional variables - they can be set in init_dim_obs_f
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!
!           The following variables are set in the routine PDAFomi_gather_obs_f
!      INTEGER :: dim_obs_p                 ! number of PE-local observations
!      INTEGER :: dim_obs_f                 ! number of full observations
!      REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
!      REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
!      REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
!      INTEGER :: dim_obs_g                 ! global number of observations 
!                                           ! (only if full obs. are restricted to process domain))
!      INTEGER, ALLOCATABLE :: id_obs_f_lim(:) ! Indices of domain-relevant full obs. in global vector of obs.
!                                           ! (only if full obs. are restricted to process domain))
!
!           Mandatory variable to be set in obs_op_f
!      INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
!   END TYPE obs_f

! Data type to define the local observations by internally shared variables of the module
!   TYPE obs_l
!      INTEGER :: dim_obs_l                 ! number of local observations
!      INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
!      INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector 
!      REAL, ALLOCATABLE :: distance_l(:)   ! Distances of local observations
!      REAL, ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
!   END TYPE obs_l
! ***********************************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could renamed the variables
  type(obs_f), public :: thisobs      ! full observation
  type(obs_l), public :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

CONTAINS

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations 
!! required for the analyses in the loop over all local 
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * thisobs\%doassim     - Whether to assimilate this type of observations
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!! * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: .true.: use global full observations)
!!
!! The following variables are set in the routine gather_obs_f
!! * thisobs\%dim_obs_p   - PE-local number of module-type observations
!! * thisobs\%dim_obs_f   - full number of module-type observations
!! * thisobs\%obs_f       - full vector of module-type observations
!! * thisobs\%ocoord_f    - coordinates of observations in OBS_MOD_F
!! * thisobs\%ivar_obs_f  - full vector of inverse obs. error variances of module-type
!! * thisobs\%dim_obs_g   - Number of global observations (only if if use_global_obs=.false)
!! * thisobs\%id_obs_f_lim - Ids of full observations in global observations (if use_global_obs=.false)
!!
!! **Adapting the template**
!! In this routine the variables listed above have to be initialized. One
!! can include modules from the model with 'use', e.g. for mesh information.
!! Alternatively one could include these as subroutine arguments
!!
  SUBROUTINE init_dim_obs_f_ice_thickness(step, dim_obs_f)

    USE PDAFomi_obs_f, &
         ONLY: PDAFomi_gather_obs_f
    USE mod_assimilation, &
         ONLY: filtertype, local_range
    USE netcdf
    USE ice_domain_size, &
         ONLY: nx_global, ny_global, ncat
    USE ice_grid, &
         ONLY: tlon, tlat
    USE ice_constants, &
         ONLY: pi, puny
    USE mod_statevector, &
         ONLY: aicen_offset, vicen_offset

    IMPLICIT NONE

    ! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs_f  !< Dimension of full observation vector

    ! *** Local variables ***
    INTEGER :: s, i, j, k                  ! Counters
    INTEGER :: cnt, cnt0                   ! Counters
    INTEGER :: stat(20000)                 ! Status flag for NetCDF commands
    INTEGER :: ncid_in                     ! ID for NetCDF file
    INTEGER :: id_3dvar                    ! IDs for fields
    INTEGER :: pos_nc(3),cnt_nc(3)         ! Vectors for reading 2D fields
    INTEGER :: dim_obs_p                   ! Number of process-local observations
    INTEGER :: status                      ! Status flag
    REAL, ALLOCATABLE :: obs_field1(:,:,:)    ! Observation field read from file
    REAL, ALLOCATABLE :: obs_field2(:,:,:)    ! Observation field read from file
    REAL, ALLOCATABLE :: ice_thick_field(:,:) ! Combination of aicen and vicen
    REAL, ALLOCATABLE :: obs_p(:)             ! PE-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)        ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)        ! PE-local observation coordinates


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - OBS_ice_thickness'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_ice_thickness) thisobs%doassim = 1

    ! Specify type of distance computation
    !thisobs%disttype = 0   ! 0=Cartesian
    thisobs%disttype = 3   ! 3=Haversine ??? TO BE CONFIRMED

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2

    ! Arrays for aicen and vicen
    ALLOCATE(obs_field1(nx_global, ny_global, ncat))
    ALLOCATE(obs_field2(nx_global, ny_global, ncat))
    
 !  ! We read in fields from a file
    IF (obs_file) THEN
       ! **********************************
       ! *** Read PE-local observations ***
       ! **********************************
       
       ! Read observation values and their coordinates,
       ! also read observation error information if available
       
       
       !1. First read fractional ice area 'aicen'
       s = 1
       stat(s) = NF90_OPEN(trim(file_ice_thickness), NF90_NOWRITE, ncid_in)

       DO i = 1, s
          IF (stat(i) .NE. NF90_NOERR) THEN
             WRITE(*,'(/9x, a, 3x, a)') &
                  'NetCDF error in opening file:', &
                  trim(file_ice_thickness)
             CALL abort_parallel()
          END IF
       END DO

       ! ******************************************
       ! *** Read file containing initial state ***
       ! ******************************************

       s=1
       stat(s) = NF90_INQ_VARID(ncid_in, 'aicen', id_3dvar)

       ! Read state variable data from file
       pos_nc = (/ 1, 1, 1 /)
       cnt_nc = (/ nx_global , ny_global, ncat /)
       s = s + 1
       stat(s) = NF90_GET_VAR(ncid_in, id_3dvar, obs_field1, start=pos_nc, count=cnt_nc)

       DO i = 1, s
          IF (stat(i) .NE. NF90_NOERR) THEN
             WRITE(*,'(/9x, a, 3x, a)') &
                  'NetCDF error in reading observations for aicen'
             CALL abort_parallel()
          END IF
       END DO

       ! *******************************************
       ! *** Close file containing initial state ***
       ! *******************************************

       s = 1
       stat(s) = NF90_CLOSE(ncid_in)

       DO i = 1, s
          IF (stat(i) .NE. NF90_NOERR) THEN
             WRITE(*,'(/9x, a, 3x, a)') &
                  'NetCDF error in closing file:', &
                  trim(file_ice_thickness)
             CALL abort_parallel()
          END IF
       END DO


       !2. Now read volume per unit ice area 'vicen'
       s = 1
       stat(s) = NF90_OPEN(trim(file_ice_thickness), NF90_NOWRITE, ncid_in)

       DO i = 1, s
          IF (stat(i) .NE. NF90_NOERR) THEN
             WRITE(*,'(/9x, a, 3x, a)') &
                  'NetCDF error in opening file:', &
                  trim(file_ice_thickness)
             CALL abort_parallel()
          END IF
       END DO

       ! ******************************************
       ! *** Read file containing initial state ***
       ! ******************************************

       s=1
       stat(s) = NF90_INQ_VARID(ncid_in, 'vicen', id_3dvar)

       ! Read state variable data from file
       pos_nc = (/ 1, 1, 1 /)
       cnt_nc = (/ nx_global , ny_global, ncat /)
       s = s + 1
       stat(s) = NF90_GET_VAR(ncid_in, id_3dvar, obs_field2, start=pos_nc, count=cnt_nc)

       DO i = 1, s
          IF (stat(i) .NE. NF90_NOERR) THEN
             WRITE(*,'(/9x, a, 3x, a)') &
                  'NetCDF error in reading observations for vicen'
             CALL abort_parallel()
          END IF
       END DO

       ! *******************************************
       ! *** Close file containing initial state ***
       ! *******************************************

       s = 1
       stat(s) = NF90_CLOSE(ncid_in)

       DO i = 1, s
          IF (stat(i) .NE. NF90_NOERR) THEN
             WRITE(*,'(/9x, a, 3x, a)') &
                  'NetCDF error in closing file:', &
                  trim(file_ice_thickness)
             CALL abort_parallel()
          END IF
       END DO

!   ! User manually specifies observation field.
    ELSE
       ! obs_field1 is aicen
       DO k = 1, ncat
          DO j = 60,60!1, ny_global
             DO i= 48,48!1, nx_global
                obs_field1(i,j,k) = 0.1
             END DO
          END DO
       END DO
       ! obs_field2 is vicen
       DO k = 1, ncat
          DO j = 60,60!1, ny_global
             DO i= 48,48!1, nx_global
                IF (k == 1) THEN
                   obs_field2(i,j,k) = 0.1!0.3
                ELSE
                   obs_field2(i,j,k) = 0.1!REAL(k) - 1.0
                END IF
             END DO
          END DO
       END DO
    END IF

! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! Compute ice thickness observation field by combining aicen
    ! and vicen fields and summing over categories.
    ALLOCATE(ice_thick_field(nx_global,ny_global))
    ice_thick_field=0.0
    DO j = 1, ny_global
       DO i= 1, nx_global
          DO k = 1, ncat
             ice_thick_field(i,j) = ice_thick_field(i,j) + &
                  obs_field2(i,j,k) * obs_field1(i,j,k)
          END DO
          IF( SUM(obs_field1(i,j,:)) > puny ) THEN
             ice_thick_field(i,j) = ice_thick_field(i,j) / &
                  SUM(obs_field1(i,j,:))
          END IF
       END DO
    END DO

    ! Count physically meaningful observations
    cnt = 0
    DO j = 1, ny_global
       DO i = 1, nx_global
          IF (ice_thick_field(i,j) > puny) cnt = cnt + 1
       END DO
    END DO
    dim_obs_p = cnt
    dim_obs_f = cnt

    IF (mype_filter==0) &
         WRITE (*,'(8x, a, i6)') '--- number of full observations', dim_obs_f


    ! *** Initialize vector of observations on the process sub-domain ***
    ! *** Initialize coordinate array of observations on the process sub-domain ***

    ! Allocate process-local observation arrays
    ALLOCATE(obs_p(dim_obs_p))
    ALLOCATE(ocoord_p(2, dim_obs_p))

    ! *** Initialize process local index array                         ***
    ! *** This array holds the information which elements of the state ***
    ! *** vector are used in the observation operator.                 ***
    ! *** It has a many rows as required for the observation operator, ***
    ! *** i.e. 1 if observations are at grid points; >1 if             ***
    ! *** interpolation is required                                    ***

    ! The initialization is done locally for each process sub-domain and later
    ! used in the observation operator. 
    ! Examples:
    ! 1. If the observations are model fields located at grid points, one should
    !   initialize the index array thisobs%id_obs_p with one row so that it contains 
    !   the indices of the observed field values in the process-local state vector
    !   (state_p). Then one can use the observation operator OBS_OP_F_GRIDPOINT 
    !   provided by MOD_OBS_OP_GENERAL_PDAF.
    ! 2. If the observations are the average of model fields located at grid points,
    !   one should initialize the index array thisobs%id_obs_p with as many rows as 
    !   values to be averaged. Each column of the arrays then contains the indices of
    !   the elements of the process-local state vector that have to be averaged. With
    !   this index array one can then use the observation operator OBS_OP_F_GRIDAVG
    !   provided by MOD_OBS_OP_GENERAL_PDAF.
    ! 3. If model values need to be interpolated to the observation location
    !   one should initialize the index array thisobs%id_obs_p with as many rows as 
    !   values are required in the interpolationto be averaged. Each column of the 
    !   array then contains the indices of elements of the process-local state vector 
    !   that are used in the interpolation.
    ! Below, you need to replace NROWS by the number of required rows

    ! Initialize with 10 rows = 2 variables * 5 categories.Observation operator then
    ! forms quotient of variables and sums over categories.
    ALLOCATE(thisobs%id_obs_p(2*ncat, dim_obs_p))

    cnt = 0
    cnt0 = 0
    DO j = 1, ny_global
       DO i = 1, nx_global
          cnt0 = cnt0 + 1
          IF (ice_thick_field(i,j) > puny) THEN
             cnt = cnt + 1
             DO k = 1, ncat
                ! Compute locations in state vector of (i,j,:) indices
                thisobs%id_obs_p(k, cnt)  = aicen_offset + cnt0 + &
                     ((k-1)*nx_global*ny_global)
                thisobs%id_obs_p(k+ncat, cnt)  = vicen_offset + cnt0 + &
                     ((k-1)*nx_global*ny_global)
             END DO
             obs_p(cnt) = ice_thick_field(i,j)
             ! Use (i+1, j+1) due to ghost cells
             ocoord_p(1, cnt) = tlon(i+1,j+1,1)*180.0/pi
             ocoord_p(2, cnt) = tlat(i+1,j+1,1)*180.0/pi
          END IF
       END DO
    END DO

! **********************************************************************
! *** Initialize interpolation coefficients for observation operator ***
! **********************************************************************

  ! This initialization is only required if an observation operator
  ! with interpolation is used. The coefficients should be determined
  ! here instead of the observation operator, because the operator is
  ! called for each ensemble member while init_dim_obs_f is only called
  ! once.

  ! Allocate array of interpolation coefficients. As thisobs%id_obs_p, the number
  ! of rows corresponds to the number of grid points using the the interpolation

!    ALLOCATE(thisobs%icoeff_p( NROWS , dim_obs_p))

  ! Ensure that the order of the coefficients is consistent with the
  ! indexing in thisobs%id_obs_p. Further ensure that the order is consistent
  ! with the assumptions used in the observation operator.

!    thisobs%icoeff_p = ...

! *************************************************************
! *** For twin experiment: Generate synthetic observations  ***
! *************************************************************

    IF (twin_experiment) THEN
       CALL add_noise(cnt, obs_p)
    END IF


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ALLOCATE(ivar_obs_p(dim_obs_p))

    ivar_obs_p = rms_ice_thickness


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    ! This routine is generic for the case that only the observations, 
    ! inverse variances and observation coordinates are gathered

    CALL PDAFomi_gather_obs_f(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, local_range, dim_obs_f)


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    DEALLOCATE(obs_field1, obs_field2)
    DEALLOCATE(ice_thick_field)
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE init_dim_obs_f_ice_thickness



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module
!! It has to append the observations to ostate_f from
!! position OFFSET_OBS+1. For the return value OFFSET_OBS
!! has to be incremented by the number of added observations.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The order of the calls to this routine for different modules
!! is important because it influences the offset of the 
!! module-type observation in the overall full observation vector.
!!
!! In general one could directly call the PDAFomi observation
!! operators in the routine obs_op_f_pdafomi in callback_obs_pdafomi.F90.
!! We leave this call here because it shows the place where one
!! would implement an own observation operator.
!!
!! Outputs for within the module are:
!! * thisobs\%off_obs_f - Offset of full module-type observation in overall full obs. vector
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE obs_op_f_ice_thickness(dim_p, dim_obs_f, state_p, ostate_f, offset_obs)

    !USE PDAFomi, &
    !     ONLY: PDAFomi_obs_op_f_gridpoint
    USE PDAFomi_obs_op_extra, &
         ONLY: PDAFomi_obs_op_f_ice_thickness
    USE ice_domain_size, &
         ONLY: ncat

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs_f             !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: ostate_f(dim_obs_f)   !< Full observed state
    INTEGER, INTENT(inout) :: offset_obs         !< input: offset of module-type observations in ostate_f
                                                 !< output: input + number of added observations


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    IF (thisobs%doassim==1) THEN
    
       !+++  Choose suitable observation operator from the
       !+++  module PDAFomi_obs_op or implement your own

       ! observation operator for observed grid point values
       !CALL PDAFomi_obs_op_f_gridpoint(thisobs, state_p, ostate_f, offset_obs)

       ! observation operator for observed grid point values
       CALL PDAFomi_obs_op_f_ice_thickness(thisobs, ncat, state_p, ostate_f, offset_obs)

    END IF

  END SUBROUTINE obs_op_f_ice_thickness
  
  SUBROUTINE add_noise(dim_state, x)

! !DESCRIPTION:
! Routine to add model error. The variance of the model
! error is  given by noise_amp times the time step size.
!
! !REVISION HISTORY:
! 2016-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: dim_state    ! State dimension
    REAL, INTENT(inout) :: x(dim_state) ! Observations
!EOP

! local variables
    REAL, ALLOCATABLE :: noise(:)      ! Random noise
    INTEGER, SAVE :: iseed(4)          ! Seed for random number generator
    INTEGER, SAVE :: firststep=1       ! Flag for first call


    ! ***********************
    ! *** Add model error ***
    ! ***********************

    ! Seeds taken from PDAF Lorenz96 routine
    IF (firststep==1) THEN
       WRITE (*,'(9x, a)') '--- Initialize seed for ice thickness noise'
       iseed(1)=2*220+1
       iseed(2)=2*100+5
       iseed(3)=2*10+7
       iseed(4)=2*30+9
       firststep=0
    ENDIF

    ! Generate random Gaussian noise
    ALLOCATE(noise(dim_state))
    CALL dlarnv(3, iseed, dim_state, noise)

    ! Add noise to generate observations
    x = x + (noise_amp * noise)

    DEALLOCATE(noise)

  END SUBROUTINE add_noise

END MODULE obs_ice_thickness_pdafomi
