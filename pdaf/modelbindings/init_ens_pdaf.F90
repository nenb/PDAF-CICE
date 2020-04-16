!$Id: init_ens.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: init_ens --- Initialize ensemble
!
! !INTERFACE:
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize an ensemble of dim\_ens states.
! Typically, the ensemble will be directly read from files.
!
! The routine is called by all filter processes and
! initializes the ensemble for the PE-local domain.
!
! Implementation for the 2D online example
! without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:

  USE netcdf
  USE mod_assimilation, &
       ONLY: screen, istate_fname
  USE ice_domain_size, &
       ONLY: nx_global, ny_global, ncat
  USE mod_statevector, &
       ONLY: fill2d_ensarray, fill3d_ensarray

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  ! It is not necessary to initialize the array 'state_p' for SEIK.
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! *** local variables ***
  INTEGER :: s, i                         ! Counters
  INTEGER :: stat(20000)                  ! Status flag for NetCDF commands
  INTEGER :: ncid_in                      ! ID for NetCDF file
  INTEGER :: id_dimx,id_dimy,id_dim_cat   ! IDs for dimensions
  INTEGER :: dim_lon_ic, dim_lat_ic, dim_cat_ic ! Initial state dimensions
  CHARACTER(len=100) :: istate_ncfile     ! File holding initial state estimate


  ! ************************************
  ! *** Routine to fill state vector ***
  ! ************************************
  !
  ! Initial state file is first checked
  ! to confirm that file dimensions con-
  ! sistent with dimensions for model run.
  ! Then ensemble array of state vectors
  ! is filled via call to routines in
  ! state vector module.

  ! ***********************************
  ! ***********************************


  ! ******************************************


  ! ******************************************
  ! *** Open file containing initial state ***
  ! ******************************************

  istate_ncfile= trim(istate_fname)//'.nc'

  WRITE (*,'(/9x, a, 3x, a)') 'Initial state estimate file:', istate_ncfile

  s = 1
  stat(s) = NF90_OPEN(istate_ncfile , NF90_NOWRITE, ncid_in)

  DO i = 1, s
     IF (stat(i) .NE. NF90_NOERR) THEN
        WRITE(*,'(/9x, a, 3x, a)') &
             'NetCDF error in opening initial state file:', istate_ncfile
        STOP
     END IF
  END DO

  ! ************************
  ! *** Check dimensions ***
  ! ************************

  s=1
  stat(s) = NF90_INQ_DIMID(ncid_in, 'ni', id_dimx)
  s = s + 1
  stat(s) = NF90_INQUIRE_DIMENSION(ncid_in, id_dimx, len=dim_lon_ic)
  s = s + 1
  stat(s) = NF90_INQ_DIMID(ncid_in, 'nj', id_dimy)
  s = s + 1
  stat(s) = NF90_INQUIRE_DIMENSION(ncid_in, id_dimy, len=dim_lat_ic)
  s = s + 1
  stat(s) = NF90_INQ_DIMID(ncid_in, 'ncat', id_dim_cat)
  s = s + 1
  stat(s) = NF90_INQUIRE_DIMENSION(ncid_in, id_dim_cat, len=dim_cat_ic)

  DO i = 1, s
     IF (stat(i) .NE. NF90_NOERR) THEN
        WRITE(*, '(/9x, a, 3x,i1 )') &
             'NetCDF error in reading dimensions from initial state file, s=', i
        STOP
     END IF
  END DO

  ! Compare global dimensions from model to dimensions in file.
  IF (dim_lon_ic .NE. nx_global .OR. dim_lat_ic .NE. ny_global .OR. &
       dim_cat_ic .NE. ncat) THEN
     WRITE (*,'(/3x, a)') 'ERROR: Dimensions in initial state file not valid.'
     STOP
  END IF

  ! *******************************************
  ! *** Close file containing initial state ***
  ! *******************************************

  s = 1
  stat(s) = NF90_CLOSE(ncid_in)

  DO i = 1, s
     IF (stat(i) .NE. NF90_NOERR) THEN
        WRITE(*,'(/9x, a, 3x, a)') &
             'NetCDF error in closing initial state file:', istate_ncfile
        STOP
     END IF
  END DO

  ! ***************************************************
  ! *** Initialize ensemble array of state vectors  ***
  ! ***************************************************

  WRITE (*,'(/1x,a)') '------- Reading Initial State -------------'
  WRITE (*,'(/1x,a)') 'Calling 2dfill_ensarray'
  WRITE (*,'(/1x,a)') 'Calling 3dfill_ensarray'

  CALL fill2d_ensarray(dim_p, dim_ens, ens_p)
  CALL fill3d_ensarray(dim_p, dim_ens, ens_p)

END SUBROUTINE init_ens_pdaf
