!$Id: output_netcdf_asml.F90 61 2019-02-01 08:49:36Z lnerger $
!BOP
!
! !MODULE:
MODULE output_netcdf_asml

! !DESCRIPTION:
! This module provides routines to initialize
! NetCDF output files for the assimilation
! and to write output into the files.
!
! !REVISION HISTORY:
! 2010-01 - Lars Nerger - Initial code
! Later revisions - see SVN log
!
! !USES:
  USE mod_parallel_pdaf, &
       ONLY: abort_parallel

  IMPLICIT NONE
  SAVE
  PUBLIC

! !PUBLIC DATA MEMBERS:
! Name of the NetCDF output file
  CHARACTER(len=100) :: file_asml            ! Name of output file; set via namelist.pdaf
  INTEGER :: delt_write_asml = 1             ! Output interval in assimilation intervals
  LOGICAL :: write_states    = .TRUE.        ! Whether to write estimated states
  LOGICAL :: write_ens       = .TRUE.        ! Whether to write full ensemble

!EOP

! Private variables
  INTEGER, PRIVATE :: file_pos   ! File position to write to
  INTEGER, PRIVATE :: fileid     ! Id of netcdf file
  INTEGER, PRIVATE :: cnt_steps  ! Count time step for delt_write_asml

! Array of 2d/3d state variables that are to be output to file
  CHARACTER(len=20), DIMENSION(1) :: output_2dvar
  CHARACTER(len=20), DIMENSION(3) :: output_3dvar

  DATA output_2dvar / 'hi_m' /
  DATA output_3dvar / 'aicen', 'vicen', 'hi_dist' /

! Array of 2d/3d state variable offsets (actually *indices* for var2d_offset and
! var3d_offset) that are to be output to file.
  INTEGER, DIMENSION(1) :: index2d_offset
  INTEGER, DIMENSION(3) :: index3d_offset

! These values can be found in mod_statevector (search for var2d_offset/var3d_offset).
  DATA index2d_offset / 16 /
  DATA index3d_offset / 1, 2, 27 /

CONTAINS
!BOP
!
! !ROUTINE: init_netcdf_asml  --- initialize netcdf output
!
! !INTERFACE:
  SUBROUTINE init_netcdf_asml(step, dt, dimx, dimy, dimcat, &
       filtertype, subtype, dim_ens, forget, local_range, &
       locweight, srange, rms_obs, delt_obs, total_steps)

! !DESCRIPTION:
! This routine initializes the netcdf file

! !USES:
    USE netcdf

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(IN) :: step       ! Initial time step
    REAL, INTENT(IN)    :: dt         ! Size of time step
    INTEGER, INTENT(IN) :: dimx       ! Dimension: longitude
    INTEGER, INTENT(IN) :: dimy       ! Dimension: latitude
    INTEGER, INTENT(IN) :: dimcat     ! Dimension: categories
    INTEGER, INTENT(IN) :: filtertype ! Type of filter
    INTEGER, INTENT(IN) :: subtype    ! Sub-type of filter
    INTEGER, INTENT(IN) :: dim_ens    ! ensemble_size
    REAL, INTENT(IN)    :: forget     ! forgetting factor
    REAL, INTENT(IN)    :: local_range  ! Localization radius
    INTEGER, INTENT(IN) :: locweight    ! Type of localization
    REAL, INTENT(IN)    :: rms_obs      ! RMS error of observations
    REAL, INTENT(IN)    :: srange       ! Support range for 5th order polynomial
                                        !   and range for 1/e for exponential weighting
    INTEGER, INTENT(IN) :: delt_obs     ! Number of time steps between two analysis steps
    INTEGER, INTENT(IN) :: total_steps  ! Total number of time steps in experiment

! Local variables
    INTEGER :: i, s, idx              ! Counters
    INTEGER :: dimid_1                ! Dimension IDs
    INTEGER :: dimid_x, dimid_y       ! Dimension IDs
    INTEGER :: dimid_cat              ! Dimension IDs
    INTEGER :: dimid_step             ! Dimension ID
    INTEGER :: dimid_ens, dimid_ensp1 ! Dimension IDs
    INTEGER :: id_tmp                 ! Variable IDs
    INTEGER :: stat(250)              ! Array for status flag
    INTEGER :: dimarray(2)            ! Array for dimensions
    CHARACTER(len=100) :: attstr      ! String to write attributes

! *** Initialize file ***

! Print screen information
    WRITE (*, '(/1x, a)') 'Initialize NetCDF output'
    WRITE (*,'(5x,a,i6)') 'Writing interval (steps) ', delt_write_asml
    IF (write_states) WRITE (*,'(5x,a)') '--> Write model states'
    IF (write_ens) WRITE (*,'(5x,a)') '--> Write ensemble states'

! Initialize file position
    file_pos = 1

! Initialize counter for output interval
    cnt_steps = 1

! Initialize file and write global atributes
    s = 1

! Overwrite existing files if present (CLOBBER)
    stat(s) = NF90_CREATE(TRIM(file_asml), NF90_CLOBBER, fileid)
    s = s + 1

    attstr  = 'Assimilation performed using PDAF-CICE model'
    stat(s) = NF90_PUT_ATT(fileid, NF90_GLOBAL, 'title', &
         TRIM(attstr))
    s = s + 1

! Define Dimensions
    stat(s) = NF90_DEF_DIM(fileid, 'ni', dimx, dimid_x)
    s = s + 1
    stat(s) = NF90_DEF_DIM(fileid, 'nj', dimy, dimid_y)
    s = s + 1
    stat(s) = NF90_DEF_DIM(fileid, 'nc', dimcat, dimid_cat)
    s = s + 1
    stat(s) = NF90_DEF_DIM(fileid, 'mem', dim_ens, dimid_ens)
    s = s + 1
    stat(s) = NF90_DEF_DIM(fileid, 'iteration', NF90_UNLIMITED, dimid_step)
    s = s + 1
    stat(s) = NF90_DEF_DIM(fileid, 'one', 1, dimid_1)
    s = s + 1
    stat(s) = NF90_DEF_DIM(fileid, 'dim_ensp1', dim_ens+1, dimid_ensp1)
    s = s + 1

! Define variables characterizing the experiment
    stat(s) = NF90_DEF_VAR(fileid, 'filtertype', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'subtype', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'dim_ens', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'forget', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'step_null', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'total_steps', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'local_range', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'locweight', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'srange', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'rms_obs', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'delt_obs', NF90_INT, DimId_1, Id_tmp)
    s = s + 1

    stat(s) = NF90_DEF_VAR(fileid, 'step', NF90_INT, DimId_step, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'step_ini', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'time', NF90_DOUBLE, DimId_step, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'time_ini', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1

    writestates: IF (write_states) THEN

       DO idx = 1, size(output_2dvar)
          ! Initialise the different CICE 2D state variables
          CALL init_2dstatevar(dimid_x, dimid_y, dimid_1, dimid_step, &
               output_2dvar(idx))
       END DO

       DO idx = 1, size(output_3dvar)
       ! Initialise the different CICE 3D state variables
          CALL init_3dstatevar(dimid_x, dimid_y, dimid_cat, dimid_1, dimid_step, &
               output_3dvar(idx))
       END DO

    END IF writestates

    writeens: IF (write_ens) THEN

       DO idx = 1, size(output_2dvar)
          ! Initialise the ensemble for different CICE 2D state variables
          CALL init_2dens(dimid_x, dimid_y, dimid_ens, dimid_1, dimid_step, &
               output_2dvar(idx))
       END DO

       DO idx = 1, size(output_3dvar)
          ! Initialise the ensemble for different CICE 3D state variables
          CALL init_3dens(dimid_x, dimid_y, dimid_cat, dimid_ens, dimid_1, &
               dimid_step, output_3dvar(idx))
       END DO

    END IF writeens

    stat(s) = NF90_ENDDEF(fileid)
    s = s + 1

! Write variables characterizing the experiment
    stat(s) = NF90_INQ_VARID(fileid, 'filtertype', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, filtertype)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'subtype', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, subtype)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'dim_ens', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, dim_ens)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'forget', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, forget)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'step_null', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, step)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'total_steps', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, total_steps)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'local_range', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, local_range)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'locweight', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, locweight)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'srange', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, srange)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'rms_obs', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, rms_obs)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'delt_obs', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, delt_obs)
    s = s + 1

    DO i = 1,  s - 1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in file initialization, no.', i
          CALL abort_parallel()
       END IF
    END DO

  END SUBROUTINE init_netcdf_asml

  SUBROUTINE init_2dstatevar(id_x, id_y, id_1, id_step, statevar)

! !DESCRIPTION:
! ! This routine initializes the different 2D CICE state variables

! !USES:
    USE netcdf

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(IN)    :: id_x        ! ID for x dimension
    INTEGER, INTENT(IN)    :: id_y        ! ID for y dimension
    INTEGER, INTENT(IN)    :: id_1        ! ID for initial state (dim=1)
    INTEGER, INTENT(IN)    :: id_step     ! ID for timestep
    CHARACTER(len=*), INTENT(IN) :: statevar     ! Name of state variable

! Local variables
    INTEGER :: i, s                    ! Counters
    INTEGER :: dimarray(3)             ! Array for dimensions
    INTEGER :: stat(250)               ! Array for status flag
    INTEGER :: id_statevar             ! ID for state variable
    CHARACTER(len=100) :: statevar_ini ! Label for initial state
    CHARACTER(len=100) :: statevar_for ! Label for forecast
    CHARACTER(len=100) :: statevar_ana ! Label for analysis
    CHARACTER(len=4)   :: anastr       ! Dummy string for label subscript


    anastr='_ini'
    WRITE(statevar_ini,'(a)') TRIM(statevar)//TRIM(anastr)
    anastr='_for'
    WRITE(statevar_for,'(a)') TRIM(statevar)//TRIM(anastr)
    anastr='_ana'
    WRITE(statevar_ana,'(a)') TRIM(statevar)//TRIM(anastr)

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_1
    s=1
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_ini), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_step
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_for), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_step
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_ana), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    DO i = 1,  s-1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in CICE variable initialization:', statevar
          CALL abort_parallel()
       END IF
    END DO

  END SUBROUTINE init_2dstatevar

  SUBROUTINE init_3dstatevar(id_x, id_y, id_cat, id_1, id_step, statevar)

! !DESCRIPTION:
! ! This routine initializes the different 3D CICE state variables

! !USES:
    USE netcdf

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(IN)    :: id_x        ! ID for x dimension
    INTEGER, INTENT(IN)    :: id_y        ! ID for y dimension
    INTEGER, INTENT(IN)    :: id_cat      ! ID for category dimension
    INTEGER, INTENT(IN)    :: id_1        ! ID for initial state (dim=1)
    INTEGER, INTENT(IN)    :: id_step     ! ID for timestep
    CHARACTER(len=*), INTENT(IN) :: statevar     ! Name of state variable

! Local variables
    INTEGER :: i, s                    ! Counters
    INTEGER :: dimarray(4)             ! Array for dimensions
    INTEGER :: stat(250)               ! Array for status flag
    INTEGER :: id_statevar             ! ID for state variable
    CHARACTER(len=100) :: statevar_ini ! Label for initial state
    CHARACTER(len=100) :: statevar_for ! Label for forecast
    CHARACTER(len=100) :: statevar_ana ! Label for analysis
    CHARACTER(len=4)   :: anastr       ! Dummy string for label subscript


    anastr='_ini'
    WRITE(statevar_ini,'(a)') TRIM(statevar)//TRIM(anastr)
    anastr='_for'
    WRITE(statevar_for,'(a)') TRIM(statevar)//TRIM(anastr)
    anastr='_ana'
    WRITE(statevar_ana,'(a)') TRIM(statevar)//TRIM(anastr)

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_cat
    dimarray(4) = id_1
    s=1
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_ini), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_cat
    dimarray(4) = id_step
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_for), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_cat
    dimarray(4) = id_step
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_ana), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    DO i = 1,  s - 1
       IF (stat(i) /= NF90_NOERR) THEN
            WRITE(*, *) 'NetCDF error in CICE variable initialization:', statevar
            CALL abort_parallel()
         END IF
      END DO

  END SUBROUTINE init_3dstatevar

  SUBROUTINE init_2dens(id_x, id_y,id_ens, id_1, id_step, statevar)

! !DESCRIPTION:
! ! This routine initializes the ensemble for different 2D CICE state variables

! !USES:
    USE netcdf

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(IN)    :: id_x        ! ID for x dimension
    INTEGER, INTENT(IN)    :: id_y        ! ID for y dimension
    INTEGER, INTENT(IN)    :: id_ens      ! ID for ensemble
    INTEGER, INTENT(IN)    :: id_1        ! ID for initial state (dim=1)
    INTEGER, INTENT(IN)    :: id_step     ! ID for timestep
    CHARACTER(len=*), INTENT(IN) :: statevar     ! Name of state variable

! Local variables
    INTEGER :: i, s                    ! Counters
    INTEGER :: dimarray(4)             ! Array for dimensions
    INTEGER :: stat(250)               ! Array for status flag
    INTEGER :: id_statevar             ! ID for state variable
    CHARACTER(len=100) :: statevar_ini ! Label for initial state
    CHARACTER(len=100) :: statevar_for ! Label for forecast
    CHARACTER(len=100) :: statevar_ana ! Label for analysis
    CHARACTER(len=8)   :: anastr       ! Dummy string for label subscript


    anastr='_ens_ini'
    WRITE(statevar_ini,'(a)') TRIM(statevar)//TRIM(anastr)
    anastr='_ens_for'
    WRITE(statevar_for,'(a)') TRIM(statevar)//TRIM(anastr)
    anastr='_ens_ana'
    WRITE(statevar_ana,'(a)') TRIM(statevar)//TRIM(anastr)

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_ens
    dimarray(4) = id_1
    s=1
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_ini), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_ens
    dimarray(4) = id_step
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_for), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_ens
    dimarray(4) = id_step
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_ana), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    DO i = 1,  s - 1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in CICE ensemble variable initialization:', statevar
          CALL abort_parallel()
       END IF
    END DO

  END SUBROUTINE init_2dens

  SUBROUTINE init_3dens(id_x, id_y, id_cat, id_ens, id_1, id_step, statevar)

! !DESCRIPTION:
! ! This routine initializes the ensemble for different 3D CICE state variables

! !USES:
    USE netcdf

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(IN)    :: id_x        ! ID for x dimension
    INTEGER, INTENT(IN)    :: id_y        ! ID for y dimension
    INTEGER, INTENT(IN)    :: id_cat      ! ID for category dimension
    INTEGER, INTENT(IN)    :: id_ens      ! ID for ensemble
    INTEGER, INTENT(IN)    :: id_1        ! ID for initial state (dim=1)
    INTEGER, INTENT(IN)    :: id_step     ! ID for timestep
    CHARACTER(len=*), INTENT(IN) :: statevar     ! Name of state variable

! Local variables
    INTEGER :: i, s                    ! Counters
    INTEGER :: dimarray(5)             ! Array for dimensions
    INTEGER :: stat(250)               ! Array for status flag
    INTEGER :: id_statevar             ! ID for state variable
    CHARACTER(len=100) :: statevar_ini ! Label for initial state
    CHARACTER(len=100) :: statevar_for ! Label for forecast
    CHARACTER(len=100) :: statevar_ana ! Label for analysis
    CHARACTER(len=8)   :: anastr       ! Dummy string for label subscript


    anastr='_ens_ini'
    WRITE(statevar_ini,'(a)') TRIM(statevar)//TRIM(anastr)
    anastr='_ens_for'
    WRITE(statevar_for,'(a)') TRIM(statevar)//TRIM(anastr)
    anastr='_ens_ana'
    WRITE(statevar_ana,'(a)') TRIM(statevar)//TRIM(anastr)

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_cat
    dimarray(4) = id_ens
    dimarray(5) = id_1
    s=1
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_ini), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_cat
    dimarray(4) = id_ens
    dimarray(5) = id_step
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_for), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_cat
    dimarray(4) = id_ens
    dimarray(5) = id_step
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_ana), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    DO i = 1,  s - 1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in CICE ensemble variable initialization:', statevar
          CALL abort_parallel()
       END IF
    END DO

  END SUBROUTINE init_3dens

!BOP
!
! !ROUTINE: write_netcdf_asml  --- write netcdf output during assimilation
!
! !INTERFACE:
  SUBROUTINE write_netcdf_asml(calltype, step, time, dim, dimx, dimy, &
       dimcat, state, dim_ens, ens)

! !DESCRIPTION:
! This routine writes the netcdf file.

! !USES:
    USE netcdf
    USE mod_statevector, &
         ONLY: var2d_offset, var3d_offset

    IMPLICIT NONE

! !ARGUMENTS:
    CHARACTER(len=3)    :: calltype    ! Type of output call
    INTEGER, INTENT(IN) :: step        ! Current time step
    REAL, INTENT(IN)    :: time        ! Current model time
    INTEGER, INTENT(IN) :: dim         ! Dimension of model state
    INTEGER, INTENT(IN) :: dimx       ! Dimension: longitude
    INTEGER, INTENT(IN) :: dimy       ! Dimension: latitude
    INTEGER, INTENT(IN) :: dimcat     ! Dimension: categories
    REAL, INTENT(IN)    :: state(dim)  ! Model state
    INTEGER, INTENT(IN) :: dim_ens     ! Ensemble size
    REAL, INTENT(IN)    :: ens(dim, dim_ens)       ! Ensemble
!EOP

! Local variables
    INTEGER :: i, s, idx, idx_off      ! Counters
    INTEGER :: ID_time, ID_step        ! Variable IDs
    INTEGER :: stat(250)               ! Array for status flag
    INTEGER :: pos(1)                  ! Position index for writing
    LOGICAL :: dowrite                 ! Flag whether to write at the current call
    INTEGER, ALLOCATABLE :: id_2dstate(:)   ! Array of 2D variable IDs
    INTEGER, ALLOCATABLE :: id_3dstate(:)   ! Array of 3D variable IDs
    INTEGER, ALLOCATABLE :: id_2dens(:)     ! Array of ensemble 2D variable IDs
    INTEGER, ALLOCATABLE :: id_3dens(:)     ! Array of ensemble 3D variable IDs
    CHARACTER(len=100)   :: varstr          ! State variable name
    CHARACTER(len=8)     :: anastr          ! Dummy string for variable name subscript


     s = 1
     stat(s) = NF90_OPEN(file_asml , NF90_WRITE, fileid)

     WRITE (*,'(/9x, a, 3x, a)') 'Writing to file:', file_asml

     IF (stat(1) .NE. NF90_NOERR) THEN
        WRITE(*,'(/9x, a, 3x, a)') &
             'NetCDF error in opening assimilation file:', file_asml
        CALL abort_parallel()
     END IF

! Check if we have to write states at this time step
    IF (cnt_steps==delt_write_asml) THEN
       dowrite = .TRUE.
       IF (calltype == 'ana') THEN
          cnt_steps = 1
       END IF
    ELSE
       dowrite = .FALSE.
       IF (calltype == 'ana') THEN
          cnt_steps = cnt_steps + 1
       END IF
    END IF

! Inquire variable Ids
    s = 1
    IF (calltype == 'ini') THEN

       IF (write_states) THEN
          ! Inquire 2D variable IDs
          ALLOCATE( id_2dstate( size(output_2dvar) ) )
          DO idx = 1, size(output_2dvar)

             anastr='_ini'
             WRITE(varstr,'(a)') TRIM(output_2dvar(idx))//TRIM(anastr)

             stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_2dstate(idx))
             s = s + 1
          END DO
          ! Inquire 3D variable IDs
          ALLOCATE( id_3dstate( size(output_3dvar) ) )
          DO idx = 1, size(output_3dvar)

             anastr='_ini'
             WRITE(varstr,'(a)') TRIM(output_3dvar(idx))//TRIM(anastr)

             stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_3dstate(idx))
             s = s + 1
          END DO
       END IF

       IF (write_ens) THEN
          ! Inquire 2D variable IDs
          ALLOCATE( id_2dens( size(output_2dvar) ) )
          DO idx = 1, size(output_2dvar)

             anastr='_ens_ini'
             WRITE(varstr,'(a)') TRIM(output_2dvar(idx))//TRIM(anastr)

             stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_2dens(idx))
             s = s + 1
          END DO
          ! Inquire 3D variable IDs
          ALLOCATE( id_3dens( size(output_3dvar) ) )
          DO idx = 1, size(output_3dvar)

             anastr='_ens_ini'
             WRITE(varstr,'(a)') TRIM(output_3dvar(idx))//TRIM(anastr)

             stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_3dens(idx))
             s = s + 1
          END DO
       END IF

       stat(s) = NF90_INQ_VARID(fileid, "time_ini", Id_time)
       s = s + 1
       stat(s) = NF90_INQ_VARID(fileid, "step_ini", Id_step)
       s = s + 1

    ELSE IF (calltype == 'for') THEN
       IF (write_states) THEN
          ! Inquire 2D variable IDs
          ALLOCATE( id_2dstate( size(output_2dvar) ) )
          DO idx = 1, size(output_2dvar)

             anastr='_for'
             WRITE(varstr,'(a)') TRIM(output_2dvar(idx))//TRIM(anastr)
             stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_2dstate(idx))
             s = s + 1
          END DO
          ALLOCATE( id_3dstate( size(output_3dvar) ) )
          ! Inquire 3D variable IDs
          DO idx = 1, size(output_3dvar)

             anastr='_for'
             WRITE(varstr,'(a)') TRIM(output_3dvar(idx))//TRIM(anastr)

             stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_3dstate(idx))
             s = s + 1
          END DO
       END IF

       IF (write_ens) THEN
          ALLOCATE( id_2dens( size(output_2dvar) ) )
          ! Inquire 2D variable IDs
          DO idx = 1, size(output_2dvar)

             anastr='_ens_for'
             WRITE(varstr,'(a)') TRIM(output_2dvar(idx))//TRIM(anastr)

             stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_2dens(idx))
             s = s + 1
          END DO
          ALLOCATE( id_3dens( size(output_3dvar) ) )
          ! Inquire 3D variable IDs
          DO idx = 1, size(output_3dvar)

             anastr='_ens_for'
             WRITE(varstr,'(a)') TRIM(output_3dvar(idx))//TRIM(anastr)

             stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_3dens(idx))
             s = s + 1
          END DO
       END IF

       stat(s) = NF90_INQ_VARID(fileid, "time", Id_time)
       s = s + 1
       stat(s) = NF90_INQ_VARID(fileid, "step", Id_step)
       s = s + 1

    ELSE

       IF (write_states) THEN
          ! Inquire 2D variable IDs
          ALLOCATE( id_2dstate( size(output_2dvar) ) )
          DO idx = 1, size(output_2dvar)
             anastr='_ana'
             WRITE(varstr,'(a)') TRIM(output_2dvar(idx))//TRIM(anastr)

             stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_2dstate(idx))
             s = s + 1
          END DO
          ALLOCATE( id_3dstate( size(output_3dvar) ) )
          ! Inquire 3D variable IDs
          DO idx = 1, size(output_3dvar)
             anastr='_ana'
             WRITE(varstr,'(a)') TRIM(output_3dvar(idx))//TRIM(anastr)

             stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_3dstate(idx))
             s = s + 1
          END DO
       END IF

       IF (write_ens) THEN
          ALLOCATE( id_2dens( size(output_2dvar) ) )
          ! Inquire 2D variable IDs
          DO idx = 1, size(output_2dvar)
             anastr='_ens_ana'
             WRITE(varstr,'(a)') TRIM(output_2dvar(idx))//TRIM(anastr)

             stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_2dens(idx))
             s = s + 1
          END DO
          ALLOCATE( id_3dens( size(output_3dvar) ) )
          ! Inquire 3D variable IDs
          DO idx = 1, size(output_3dvar)
             anastr='_ens_ana'
             WRITE(varstr,'(a)') TRIM(output_3dvar(idx))//TRIM(anastr)

             stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_3dens(idx))
             s = s + 1
          END DO
       END IF

       stat(s) = NF90_INQ_VARID(fileid, "time", Id_time)
       s = s + 1
       stat(s) = NF90_INQ_VARID(fileid, "step", Id_step)
       s = s + 1

    END IF

    DO i = 1,  s - 1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in preparing output, no.', i
          CALL abort_parallel()
       END IF
    END DO
    s = 1

! Write variables
    IF (calltype == 'ini') THEN
       pos(1) = file_pos
       stat(s) = NF90_PUT_VAR(fileid, Id_time, time, pos)
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_step, step, pos)
       s = s + 1
    END IF

    IF (calltype == 'for' .AND. dowrite) THEN
       pos(1) = file_pos
       stat(s) = NF90_PUT_VAR(fileid, Id_time, time, pos)
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_step, -step, pos)
       s = s + 1
    END IF

    DO i = 1,  s - 1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in writing output, no.', i
          CALL abort_parallel()
       END IF
    END DO

    ! Write state information and RMS errors only at specified intervals
    ! and at initialization
    writetimedep: IF (dowrite .OR. calltype=='ini') THEN
       IF (write_states) THEN
          ! Write 2D states
          DO idx = 1, size(id_2dstate)
             idx_off=index2d_offset(idx)
             CALL write_2dstate(Id_2dstate(idx), var2d_offset(idx_off), dim, state, &
                  dimx, dimy, file_pos)
          END DO

          ! Write 3D states
          DO idx = 1, size(id_3dstate)
             idx_off=index3d_offset(idx)
             CALL write_3dstate(Id_3dstate(idx), var3d_offset(idx_off), dim, state, &
                  dimx, dimy, dimcat, file_pos)
          END DO
       END IF

       IF (write_ens) THEN
          ! Write 2D ensemble
          DO idx = 1, size(id_2dens)
             idx_off=index2d_offset(idx)
             CALL write_2dens(Id_2dens(idx), var2d_offset(idx_off), dim, dim_ens, &
                  ens, dimx, dimy, file_pos)
          END DO

          ! Write 3D ensemble
          DO idx = 1, size(id_3dens)
             idx_off=index3d_offset(idx)
             CALL write_3dens(Id_3dens(idx), var3d_offset(idx_off), dim, dim_ens, &
                  ens, dimx, dimy, dimcat, file_pos)
          END DO
       END IF

    END IF writetimedep

    ! Increment file position counter
    IF (dowrite .AND. calltype == 'ana') THEN
       file_pos = file_pos + 1
    END IF

    ! Tidy up arrays of state variable IDs
    DEALLOCATE(Id_2dstate, Id_3dstate, Id_2dens, Id_3dens)

  END SUBROUTINE write_netcdf_asml

  SUBROUTINE write_2dstate(idvar, offset, dim, state, dimx, dimy, iter)

! !DESCRIPTION:
! ! This routine writes the different 2D CICE state variables to file

! !USES:
    USE netcdf

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(IN) :: idvar       ! 2D variable ID
    INTEGER, INTENT(IN) :: offset      ! 2D variable offset in state vector
    INTEGER, INTENT(IN) :: dim         ! Dimension of model state
    REAL, INTENT(IN)    :: state(dim)  ! Model state
    INTEGER, INTENT(IN) :: dimx        ! Dimension: longitude
    INTEGER, INTENT(IN) :: dimy        ! Dimension: latitude
    INTEGER, INTENT(IN) :: iter        ! Iteration number

    ! Local variables
    INTEGER :: i, j, s                      ! Counters
    INTEGER :: stat(20)                     ! Array for status flag
    INTEGER :: pos3(3)                      ! Position index for writing
    INTEGER :: cnt3(3)                      ! Count index for writing
    REAL, ALLOCATABLE :: tmp_array(:,:)     ! Temporary array holding 2D data


    ! Reshape statevector into 2D format
    ALLOCATE(tmp_array(dimx, dimy))
    DO j = 1, dimy
       DO i = 1, dimx
          tmp_array(i,j) = state(offset + i + (j-1)*dimx)
       END DO
    END DO

    s=1
    pos3(1) = 1
    pos3(2) = 1
    pos3(3) = iter
    cnt3(1) = dimx
    cnt3(2) = dimy
    cnt3(3) = 1
    stat(s) = NF90_PUT_VAR(fileid, idvar, tmp_array, pos3, cnt3)
    s = s + 1

    DO i = 1,  s-1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in CICE 2D state write. Offset:', offset
          CALL abort_parallel()
       END IF
    END DO

    DEALLOCATE(tmp_array)

  END SUBROUTINE write_2dstate

  SUBROUTINE write_3dstate(idvar, offset, dim, state, dimx, dimy, dimcat, iter)

! !DESCRIPTION:
! ! This routine writes the different 2D CICE state variables to file

! !USES:
    USE netcdf

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(IN) :: idvar       ! 3D variable ID
    INTEGER, INTENT(IN) :: offset      ! 3D variable offset in state vector
    INTEGER, INTENT(IN) :: dim         ! Dimension of model state
    REAL, INTENT(IN)    :: state(dim)  ! Model state
    INTEGER, INTENT(IN) :: dimx        ! Dimension: longitude
    INTEGER, INTENT(IN) :: dimy        ! Dimension: latitude
    INTEGER, INTENT(IN) :: dimcat      ! Dimension: Number of categories
    INTEGER, INTENT(IN) :: iter        ! Iteration number

    ! Local variables
    INTEGER :: i, j, k, s                     ! Counters
    INTEGER :: stat(20)                       ! Array for status flag
    INTEGER :: pos4(4)                        ! Position index for writing
    INTEGER :: cnt4(4)                        ! Count index for writing
    REAL, ALLOCATABLE :: tmp_array(:,:,:)     ! Temporary array holding 3D data


    ! Reshape statevector into 3D format
    ALLOCATE(tmp_array(dimx, dimy, dimcat))
    DO k = 1 ,dimcat
       DO j = 1, dimy
          DO i = 1, dimx
             tmp_array(i,j,k) = state(offset + i + (j-1)*dimx + (k-1)*dimx*dimy)
          END DO
       END DO
    END DO

    s=1
    pos4(1) = 1
    pos4(2) = 1
    pos4(3) = 1
    pos4(4) = iter
    cnt4(1) = dimx
    cnt4(2) = dimy
    cnt4(3) = dimcat
    cnt4(4) = 1
    stat(s) = NF90_PUT_VAR(fileid, idvar, tmp_array,  pos4, cnt4)
    s = s + 1

    DO i = 1,  s-1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in CICE 3D state write. Offset:', offset
          CALL abort_parallel()
       END IF
    END DO

    DEALLOCATE(tmp_array)

  END SUBROUTINE write_3dstate

  SUBROUTINE write_2dens(idvar, offset, dim, dim_ens, ens, dimx, dimy, iter)

! !DESCRIPTION:
! ! This routine writes the different 2D CICE ensemble variables to file

! !USES:
    USE netcdf

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(IN) :: idvar       ! 2D variable ID
    INTEGER, INTENT(IN) :: offset      ! 2D variable offset in state vector
    INTEGER, INTENT(IN) :: dim         ! Dimension of model state
    INTEGER, INTENT(IN) :: dim_ens     ! Ensemble size
    REAL, INTENT(IN)    :: ens(dim, dim_ens)       ! Ensemble
    INTEGER, INTENT(IN) :: dimx        ! Dimension: longitude
    INTEGER, INTENT(IN) :: dimy        ! Dimension: latitude
    INTEGER, INTENT(IN) :: iter        ! Iteration number

    ! Local variables
    INTEGER :: i, j, k, s                     ! Counters
    INTEGER :: stat(20)                       ! Array for status flag
    INTEGER :: pos4(4)                        ! Position index for writing
    INTEGER :: cnt4(4)                        ! Count index for writing
    REAL, ALLOCATABLE :: tmp_array(:,:,:)     ! Temporary array holding ens data


    ! Reshape statevector into ensemble format
    ALLOCATE(tmp_array(dimx, dimy, dim_ens))
    DO k = 1 ,dim_ens
       DO j = 1, dimy
          DO i = 1, dimx
             tmp_array(i,j,k) = ens(offset + i + (j-1)*dimx, k)
          END DO
       END DO
    END DO

    s=1
    pos4(1) = 1
    pos4(2) = 1
    pos4(3) = 1
    pos4(4) = iter
    cnt4(1) = dimx
    cnt4(2) = dimy
    cnt4(3) = dim_ens
    cnt4(4) = 1
    stat(s) = NF90_PUT_VAR(fileid, idvar, tmp_array, pos4, cnt4)
    s = s + 1

    DO i = 1,  s-1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in CICE 2D ensemble write. Offset:', offset
          CALL abort_parallel()
       END IF
    END DO

    DEALLOCATE(tmp_array)

  END SUBROUTINE write_2dens

  SUBROUTINE write_3dens(idvar, offset, dim, dim_ens, ens, dimx, dimy, dimcat, iter)

! !DESCRIPTION:
! ! This routine writes the different 3D CICE ensemble variables to file

! !USES:
    USE netcdf

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(IN) :: idvar       ! 2D variable ID
    INTEGER, INTENT(IN) :: offset      ! 2D variable offset in state vector
    INTEGER, INTENT(IN) :: dim         ! Dimension of model state
    INTEGER, INTENT(IN) :: dim_ens     ! Ensemble size
    REAL, INTENT(IN)    :: ens(dim, dim_ens)       ! Ensemble
    INTEGER, INTENT(IN) :: dimx        ! Dimension: longitude
    INTEGER, INTENT(IN) :: dimy        ! Dimension: latitude
    INTEGER, INTENT(IN) :: dimcat      ! Dimension: Number of categories
    INTEGER, INTENT(IN) :: iter        ! Iteration number

    ! Local variables
    INTEGER :: i, j, k, l, s                    ! Counters
    INTEGER :: stat(20)                         ! Array for status flag
    INTEGER :: pos5(5)                          ! Position index for writing
    INTEGER :: cnt5(5)                          ! Count index for writing
    REAL, ALLOCATABLE :: tmp_array(:,:,:,:)     ! Temporary array holding ens data


    ! Reshape statevector into ensemble format
    ALLOCATE(tmp_array(dimx, dimy, dimcat, dim_ens))
    DO l = 1 , dim_ens
       DO k = 1 ,dimcat
          DO j = 1, dimy
             DO i = 1, dimx
                tmp_array(i,j,k,l) = ens(offset + i + (j-1)*dimx + dimx*dimy*(k-1), l)
             END DO
          END DO
       END DO
    END DO

    s=1
    pos5(1) = 1
    pos5(2) = 1
    pos5(3) = 1
    pos5(4) = 1
    pos5(5) = iter
    cnt5(1) = dimx
    cnt5(2) = dimy
    cnt5(3) = dimcat
    cnt5(4) = dim_ens
    cnt5(5) = 1
    stat(s) = NF90_PUT_VAR(fileid, idvar, tmp_array,  pos5, cnt5)
    s = s + 1

    DO i = 1,  s-1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in CICE 3D ensemble write. Offset:', offset
          CALL abort_parallel()
       END IF
    END DO

    DEALLOCATE(tmp_array)

  END SUBROUTINE write_3dens
!BOP
!
! !ROUTINE: close_netcdf_asml  --- close netcdf file
!
! !INTERFACE:
  SUBROUTINE close_netcdf_asml()

! !DESCRIPTION:
! This routine closes the netcdf file.

! !USES:
    USE netcdf

    IMPLICIT NONE

!EOP

! Local variables
    INTEGER :: stat(50)                ! Array for status flag

! Close file

    stat(1) = NF90_CLOSE(fileid)
    IF (stat(1) /= NF90_NOERR) THEN
       WRITE(*, *) 'NetCDF error in closing output file, no. 1'
       CALL abort_parallel()
    END IF

  END SUBROUTINE close_netcdf_asml

END MODULE output_netcdf_asml
