MODULE mod_statevector

! !DESCRIPTION:
! This module provides variables & routines for
! manipulating the state vector.

! !USES:
  IMPLICIT NONE
  SAVE

! 2d state vector variables - start index
  INTEGER :: uvel_offset
  INTEGER :: vvel_offset
  INTEGER :: stressp_1_offset
  INTEGER :: stressp_2_offset
  INTEGER :: stressp_3_offset
  INTEGER :: stressp_4_offset
  INTEGER :: stressm_1_offset
  INTEGER :: stressm_2_offset
  INTEGER :: stressm_3_offset
  INTEGER :: stressm_4_offset
  INTEGER :: stress12_1_offset
  INTEGER :: stress12_2_offset
  INTEGER :: stress12_3_offset
  INTEGER :: stress12_4_offset
  INTEGER :: sst_offset
#ifdef USE_STRESS
  INTEGER :: a11_1_offset
  INTEGER :: a11_2_offset
  INTEGER :: a11_3_offset
  INTEGER :: a11_4_offset
  INTEGER :: a12_1_offset
  INTEGER :: a12_2_offset
  INTEGER :: a12_3_offset
  INTEGER :: a12_4_offset
#endif

! Array holding 2d state variable offsets
#ifdef USE_STRESS
  INTEGER :: var2d_offset(23)
#else
  INTEGER :: var2d_offset(15)
#endif

! 2d state vector variables - dimension size
  INTEGER :: uvel_dim_state
  INTEGER :: vvel_dim_state
  INTEGER :: stressp_1_dim_state
  INTEGER :: stressp_2_dim_state
  INTEGER :: stressp_3_dim_state
  INTEGER :: stressp_4_dim_state
  INTEGER :: stressm_1_dim_state
  INTEGER :: stressm_2_dim_state
  INTEGER :: stressm_3_dim_state
  INTEGER :: stressm_4_dim_state
  INTEGER :: stress12_1_dim_state
  INTEGER :: stress12_2_dim_state
  INTEGER :: stress12_3_dim_state
  INTEGER :: stress12_4_dim_state
  INTEGER :: sst_dim_state
#ifdef USE_STRESS
  INTEGER :: a11_1_dim_state
  INTEGER :: a11_2_dim_state
  INTEGER :: a11_3_dim_state
  INTEGER :: a11_4_dim_state
  INTEGER :: a12_1_dim_state
  INTEGER :: a12_2_dim_state
  INTEGER :: a12_3_dim_state
  INTEGER :: a12_4_dim_state
#endif

! Array holding 2d state variable dimensions
#ifdef USE_STRESS
  INTEGER :: var2d_dim_state(23)
#else
  INTEGER :: var2d_dim_state(15)
#endif

! 3d state vector variables - start index
  INTEGER :: aicen_offset
  INTEGER :: vicen_offset
  INTEGER :: vsnon_offset
  INTEGER :: Tsfcn_offset
  INTEGER :: iage_offset
  INTEGER :: FY_offset
  INTEGER :: alvl_offset
  INTEGER :: vlvl_offset
  INTEGER :: apnd_offset
  INTEGER :: hpnd_offset
  INTEGER :: ipnd_offset
  INTEGER :: sice001_offset
  INTEGER :: qice001_offset
  INTEGER :: sice002_offset
  INTEGER :: qice002_offset
  INTEGER :: sice003_offset
  INTEGER :: qice003_offset
  INTEGER :: sice004_offset
  INTEGER :: qice004_offset
  INTEGER :: sice005_offset
  INTEGER :: qice005_offset
  INTEGER :: sice006_offset
  INTEGER :: qice006_offset
  INTEGER :: sice007_offset
  INTEGER :: qice007_offset
  INTEGER :: qsno001_offset

! Array holding 3d state variable offsets
  INTEGER :: var3d_offset(26)

! 3d state vector variables - dimension size
  INTEGER :: aicen_dim_state
  INTEGER :: vicen_dim_state
  INTEGER :: vsnon_dim_state
  INTEGER :: Tsfcn_dim_state
  INTEGER :: iage_dim_state
  INTEGER :: FY_dim_state
  INTEGER :: alvl_dim_state
  INTEGER :: vlvl_dim_state
  INTEGER :: apnd_dim_state
  INTEGER :: hpnd_dim_state
  INTEGER :: ipnd_dim_state
  INTEGER :: sice001_dim_state
  INTEGER :: qice001_dim_state
  INTEGER :: sice002_dim_state
  INTEGER :: qice002_dim_state
  INTEGER :: sice003_dim_state
  INTEGER :: qice003_dim_state
  INTEGER :: sice004_dim_state
  INTEGER :: qice004_dim_state
  INTEGER :: sice005_dim_state
  INTEGER :: qice005_dim_state
  INTEGER :: sice006_dim_state
  INTEGER :: qice006_dim_state
  INTEGER :: sice007_dim_state
  INTEGER :: qice007_dim_state
  INTEGER :: qsno001_dim_state

! Array holding 3d state variable dimensions
  INTEGER :: var3d_dim_state(26)

! Array for 2d/3d state variable .NC ids
#ifdef USE_STRESS
  CHARACTER(len=20), DIMENSION(23) :: id2d_list
#else
  CHARACTER(len=20), DIMENSION(15) :: id2d_list
#endif
  CHARACTER(len=20), DIMENSION(26) :: id3d_list

! Fill array of 2d state variable .NC ids
#ifdef USE_STRESS
  DATA id2d_list / 'uvel', 'vvel', 'stressp_1', 'stressp_2',&
       'stressp_3', 'stressp_4', 'stressm_1', 'stressm_2',&
       'stressm_3', 'stressm_4', 'stress12_1', 'stress12_2',&
       'stress12_3', 'stress12_4', 'sst', 'a11_1', 'a11_2',&
       'a11_3', 'a11_4',  'a12_1', 'a12_2', 'a12_3',&
       'a11_2' /
#else
  DATA id2d_list / 'uvel', 'vvel', 'stressp_1', 'stressp_2',&
       'stressp_3', 'stressp_4', 'stressm_1', 'stressm_2',&
       'stressm_3', 'stressm_4', 'stress12_1', 'stress12_2',&
       'stress12_3', 'stress12_4', 'sst' /
#endif

! Fill array of 3d state variable .NC ids
  DATA id3d_list / 'aicen', 'vicen', 'vsnon', 'Tsfcn',&
       'iage', 'FY', 'alvl', 'vlvl', 'apnd',&
       'hpnd', 'ipnd', 'sice001', 'qice001', 'sice002',&
       'qice002', 'sice003', 'qice003', 'sice004', 'qice004',&
       'sice005', 'qice005', 'sice006', 'qice006', 'qice007',&
       'sice007', 'qsno001' /


CONTAINS


SUBROUTINE calc_2d_offset()

! !DESCRIPTION:
! This routine calculates the offset values for the 2d state variables.
! It then stores these offset values in a 1D array.

! !USES:
  USE ice_domain_size, &
       ONLY: nx_global, ny_global

  IMPLICIT NONE


  uvel_offset = 0
  vvel_offset = uvel_offset + nx_global*ny_global
  stressp_1_offset = vvel_offset + nx_global*ny_global
  stressp_2_offset = stressp_1_offset + nx_global*ny_global
  stressp_3_offset = stressp_2_offset + nx_global*ny_global
  stressp_4_offset = stressp_3_offset + nx_global*ny_global
  stressm_1_offset = stressp_4_offset + nx_global*ny_global
  stressm_2_offset = stressm_1_offset + nx_global*ny_global
  stressm_3_offset = stressm_2_offset + nx_global*ny_global
  stressm_4_offset = stressm_3_offset + nx_global*ny_global
  stress12_1_offset = stressm_4_offset + nx_global*ny_global
  stress12_2_offset = stress12_1_offset + nx_global*ny_global
  stress12_3_offset = stress12_2_offset + nx_global*ny_global
  stress12_4_offset = stress12_3_offset + nx_global*ny_global
  sst_offset = stress12_4_offset + nx_global*ny_global
#ifdef USE_STRESS
  a11_1_offset = sst_offset + nx_global*ny_global
  a11_2_offset = a11_1_offset + nx_global*ny_global
  a11_3_offset = a11_2_offset + nx_global*ny_global
  a11_4_offset = a11_3_offset + nx_global*ny_global
  a12_1_offset = a11_4_offset + nx_global*ny_global
  a12_2_offset = a12_1_offset + nx_global*ny_global
  a12_3_offset = a12_2_offset + nx_global*ny_global
  a12_4_offset = a12_3_offset + nx_global*ny_global
#endif

  ! Fill array of state variable offsets
  var2d_offset(1) = uvel_offset
  var2d_offset(2) = vvel_offset
  var2d_offset(3) = stressp_1_offset
  var2d_offset(4) = stressp_2_offset
  var2d_offset(5) = stressp_3_offset
  var2d_offset(6) = stressp_4_offset
  var2d_offset(7) = stressm_1_offset
  var2d_offset(8) = stressm_2_offset
  var2d_offset(9) = stressm_3_offset
  var2d_offset(10) = stressm_4_offset
  var2d_offset(11) = stress12_1_offset
  var2d_offset(12) = stress12_2_offset
  var2d_offset(13) = stress12_3_offset
  var2d_offset(14) = stress12_4_offset
  var2d_offset(15) = sst_offset
#ifdef USE_STRESS
  var2d_offset(16) = a11_1_offset
  var2d_offset(17) = a11_2_offset
  var2d_offset(18) = a11_3_offset
  var2d_offset(19) = a11_4_offset
  var2d_offset(20) = a12_1_offset
  var2d_offset(21) = a12_2_offset
  var2d_offset(22) = a12_3_offset
  var2d_offset(23) = a12_4_offset
#endif

END SUBROUTINE calc_2d_offset


SUBROUTINE calc_3d_offset()

! !DESCRIPTION:
! This routine calculates the offset values for the 3d state variables.
! It then stores these offset values in a 1D array.

! !USES:
  USE ice_domain_size, &
       ONLY: nx_global, ny_global, ncat

  IMPLICIT NONE


  ! Compute 2d offset values in case not already done
  CALL calc_2d_offset()

#ifdef USE_STRESS
  aicen_offset = a12_4_offset + nx_global*ny_global ! Continue
  ! from 2d state variable a12_4 offset
#else
  aicen_offset = sst_offset + nx_global*ny_global ! Continue
  ! from 2d state variable sst offset
#endif
  vicen_offset = aicen_offset + nx_global*ny_global*ncat
  vsnon_offset = vicen_offset + nx_global*ny_global*ncat
  Tsfcn_offset = vsnon_offset + nx_global*ny_global*ncat
  iage_offset = Tsfcn_offset + nx_global*ny_global*ncat
  FY_offset = iage_offset + nx_global*ny_global*ncat
  alvl_offset = FY_offset + nx_global*ny_global*ncat
  vlvl_offset = alvl_offset + nx_global*ny_global*ncat
  apnd_offset = vlvl_offset + nx_global*ny_global*ncat
  hpnd_offset = apnd_offset + nx_global*ny_global*ncat
  ipnd_offset = hpnd_offset + nx_global*ny_global*ncat
  sice001_offset = ipnd_offset + nx_global*ny_global*ncat
  qice001_offset = sice001_offset + nx_global*ny_global*ncat
  sice002_offset = qice001_offset + nx_global*ny_global*ncat
  qice002_offset = sice002_offset + nx_global*ny_global*ncat
  sice003_offset = qice002_offset + nx_global*ny_global*ncat
  qice003_offset = sice003_offset + nx_global*ny_global*ncat
  sice004_offset = qice003_offset + nx_global*ny_global*ncat
  qice004_offset = sice004_offset + nx_global*ny_global*ncat
  sice005_offset = qice004_offset + nx_global*ny_global*ncat
  qice005_offset = sice005_offset + nx_global*ny_global*ncat
  sice006_offset = qice005_offset + nx_global*ny_global*ncat
  qice006_offset = sice006_offset + nx_global*ny_global*ncat
  sice007_offset = qice006_offset + nx_global*ny_global*ncat
  qice007_offset = sice007_offset + nx_global*ny_global*ncat
  qsno001_offset = qice007_offset + nx_global*ny_global*ncat

  ! Fill  array of state variable offsets
  var3d_offset(1) = aicen_offset
  var3d_offset(2) = vicen_offset
  var3d_offset(3) = vsnon_offset
  var3d_offset(4) = Tsfcn_offset
  var3d_offset(5) = iage_offset
  var3d_offset(6) = FY_offset
  var3d_offset(7) = alvl_offset
  var3d_offset(8) = vlvl_offset
  var3d_offset(9) = apnd_offset
  var3d_offset(10) = hpnd_offset
  var3d_offset(11) = ipnd_offset
  var3d_offset(12) = sice001_offset
  var3d_offset(13) = qice001_offset
  var3d_offset(14) = sice002_offset
  var3d_offset(15) = qice002_offset
  var3d_offset(16) = sice003_offset
  var3d_offset(17) = qice003_offset
  var3d_offset(18) = sice004_offset
  var3d_offset(19) = qice004_offset
  var3d_offset(20) = sice005_offset
  var3d_offset(21) = qice005_offset
  var3d_offset(22) = sice006_offset
  var3d_offset(23) = qice006_offset
  var3d_offset(24) = qice007_offset
  var3d_offset(25) = sice007_offset
  var3d_offset(26) = qsno001_offset

END SUBROUTINE calc_3d_offset

SUBROUTINE calc_2d_dim()

! !DESCRIPTION:
! This routine calculates the dimension of the 2d state variables.

! !USES:
  USE ice_domain_size, &
       ONLY: nx_global, ny_global

  IMPLICIT NONE


  uvel_dim_state = nx_global*ny_global
  vvel_dim_state = nx_global*ny_global
  stressp_1_dim_state = nx_global*ny_global
  stressp_2_dim_state = nx_global*ny_global
  stressp_3_dim_state = nx_global*ny_global
  stressp_4_dim_state = nx_global*ny_global
  stressm_1_dim_state = nx_global*ny_global
  stressm_2_dim_state = nx_global*ny_global
  stressm_3_dim_state = nx_global*ny_global
  stressm_4_dim_state = nx_global*ny_global
  stress12_1_dim_state = nx_global*ny_global
  stress12_2_dim_state = nx_global*ny_global
  stress12_3_dim_state = nx_global*ny_global
  stress12_4_dim_state = nx_global*ny_global
  sst_dim_state = nx_global*ny_global
#ifdef USE_STRESS
  a11_1_dim_state = nx_global*ny_global
  a11_2_dim_state = nx_global*ny_global
  a11_3_dim_state = nx_global*ny_global
  a11_4_dim_state = nx_global*ny_global
  a12_1_dim_state = nx_global*ny_global
  a12_2_dim_state = nx_global*ny_global
  a12_3_dim_state = nx_global*ny_global
  a12_4_dim_state = nx_global*ny_global
#endif

  ! Fill array of state variable dimensions
  var2d_dim_state(1) = uvel_dim_state
  var2d_dim_state(2) = vvel_dim_state
  var2d_dim_state(3) = stressp_1_dim_state
  var2d_dim_state(4) = stressp_2_dim_state
  var2d_dim_state(5) = stressp_3_dim_state
  var2d_dim_state(6) = stressp_4_dim_state
  var2d_dim_state(7) = stressm_1_dim_state
  var2d_dim_state(8) = stressm_2_dim_state
  var2d_dim_state(9) = stressm_3_dim_state
  var2d_dim_state(10) = stressm_4_dim_state
  var2d_dim_state(11) = stress12_1_dim_state
  var2d_dim_state(12) = stress12_2_dim_state
  var2d_dim_state(13) = stress12_3_dim_state
  var2d_dim_state(14) = stress12_4_dim_state
  var2d_dim_state(15) = sst_dim_state
#ifdef USE_STRESS
  var2d_dim_state(16) = a11_1_dim_state
  var2d_dim_state(17) = a11_2_dim_state
  var2d_dim_state(18) = a11_3_dim_state
  var2d_dim_state(19) = a11_4_dim_state
  var2d_dim_state(20) = a12_1_dim_state
  var2d_dim_state(21) = a12_2_dim_state
  var2d_dim_state(22) = a12_3_dim_state
  var2d_dim_state(23) = a12_4_dim_state
#endif

END SUBROUTINE calc_2d_dim

SUBROUTINE calc_3d_dim()

! !DESCRIPTION:
! This routine calculates the dimension of the 3d state variables.

! !USES:
  USE ice_domain_size, &
       ONLY: nx_global, ny_global, ncat

  IMPLICIT NONE


  aicen_dim_state = nx_global*ny_global*ncat
  vicen_dim_state = nx_global*ny_global*ncat
  vsnon_dim_state = nx_global*ny_global*ncat
  Tsfcn_dim_state = nx_global*ny_global*ncat
  iage_dim_state = nx_global*ny_global*ncat
  FY_dim_state = nx_global*ny_global*ncat
  alvl_dim_state = nx_global*ny_global*ncat
  vlvl_dim_state = nx_global*ny_global*ncat
  apnd_dim_state = nx_global*ny_global*ncat
  hpnd_dim_state = nx_global*ny_global*ncat
  ipnd_dim_state = nx_global*ny_global*ncat
  sice001_dim_state = nx_global*ny_global*ncat
  qice001_dim_state = nx_global*ny_global*ncat
  sice002_dim_state = nx_global*ny_global*ncat
  qice002_dim_state = nx_global*ny_global*ncat
  sice003_dim_state = nx_global*ny_global*ncat
  qice003_dim_state = nx_global*ny_global*ncat
  sice004_dim_state = nx_global*ny_global*ncat
  qice004_dim_state = nx_global*ny_global*ncat
  sice005_dim_state = nx_global*ny_global*ncat
  qice005_dim_state = nx_global*ny_global*ncat
  sice006_dim_state = nx_global*ny_global*ncat
  qice006_dim_state = nx_global*ny_global*ncat
  sice007_dim_state = nx_global*ny_global*ncat
  qice007_dim_state = nx_global*ny_global*ncat
  qsno001_dim_state = nx_global*ny_global*ncat

  ! Fill  array of state variable dimensions
  var3d_dim_state(1) = aicen_dim_state
  var3d_dim_state(2) = vicen_dim_state
  var3d_dim_state(3) = vsnon_dim_state
  var3d_dim_state(4) = Tsfcn_dim_state
  var3d_dim_state(5) = iage_dim_state
  var3d_dim_state(6) = FY_dim_state
  var3d_dim_state(7) = alvl_dim_state
  var3d_dim_state(8) = vlvl_dim_state
  var3d_dim_state(9) = apnd_dim_state
  var3d_dim_state(10) = hpnd_dim_state
  var3d_dim_state(11) = ipnd_dim_state
  var3d_dim_state(12) = sice001_dim_state
  var3d_dim_state(13) = qice001_dim_state
  var3d_dim_state(14) = sice002_dim_state
  var3d_dim_state(15) = qice002_dim_state
  var3d_dim_state(16) = sice003_dim_state
  var3d_dim_state(17) = qice003_dim_state
  var3d_dim_state(18) = sice004_dim_state
  var3d_dim_state(19) = qice004_dim_state
  var3d_dim_state(20) = sice005_dim_state
  var3d_dim_state(21) = qice005_dim_state
  var3d_dim_state(22) = sice006_dim_state
  var3d_dim_state(23) = qice006_dim_state
  var3d_dim_state(24) = qice007_dim_state
  var3d_dim_state(25) = sice007_dim_state
  var3d_dim_state(26) = qsno001_dim_state

END SUBROUTINE calc_3d_dim

SUBROUTINE calc_statevector_dim(dim_state_p)

! !DESCRIPTION:
! This routine calculates the total state vector dimension.

! !USES:
  IMPLICIT NONE

! !ARGUMENTS
  INTEGER, INTENT(inout) :: dim_state_p   ! Dimension of state vector

  ! Local variables
  INTEGER :: i  ! Counter


  ! Calculate state variable dimensions in case not already calculated
  CALL calc_2d_dim()
  CALL calc_3d_dim()

  dim_state_p = 0

  ! Add 2d state variable dimensions
  DO i = 1, size(var2d_dim_state)
     dim_state_p = dim_state_p + var2d_dim_state(i)
  END DO

  ! Add 3d state variable dimensions
  DO i = 1, size(var3d_dim_state)
     dim_state_p = dim_state_p + var3d_dim_state(i)
  END DO

END SUBROUTINE calc_statevector_dim

INTEGER FUNCTION calc_local_dim()

! !DESCRIPTION:
! This function calculates the local domain dimension.

  ! !USES:
  USE ice_domain_size, &
       ONLY: ncat
  IMPLICIT NONE

  ! Dimension size is total number of state variables.
  ! NOTE: different categories are considered as different
  ! 'state variables'.
  calc_local_dim = size(id2d_list) + ( ncat*size(id3d_list) )

  RETURN

END FUNCTION

SUBROUTINE fill2d_ensarray(dim_p, dim_ens, ens_p)

! !DESCRIPTION:
! Fill ensemble array with 2d state variables from initial state file.

! !USES:
  USE netcdf
  USE ice_calendar, &
       ONLY: year_init
  USE mod_assimilation, &
       ONLY: istate_dir
  USE ice_domain_size, &
       ONLY: nx_global, ny_global
  USE mod_parallel_pdaf, &
       ONLY: abort_parallel

  IMPLICIT NONE

! !ARGUMENTS
  INTEGER, INTENT(in) :: dim_p                     ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                   ! Size of ensemble
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)     ! PE-local state ensemble

! *** local variables ***
  INTEGER :: s, i, j, idx, member,yr ! Counters
  INTEGER :: stat(20000)             ! Status flag for NetCDF commands
  INTEGER :: ncid_in                 ! ID for NetCDF file
  INTEGER :: id_2dvar                ! IDs for fields
  INTEGER :: pos(2),cnt(2)           ! Vectors for 2D reading fields
  CHARACTER(len=100) :: istate_ncfile     ! File holding initial state estimate
  CHARACTER(len=100) :: year              ! Year of initial state estimate
  REAL, ALLOCATABLE :: var2d(:,:)         ! Array for reading state variables


  member2d:DO member = 1 , dim_ens
     ! ******************************************
     ! *** Open file containing initial state ***
     ! ******************************************

     yr = year_init + member - 1
     WRITE(year, '(i4)') yr
     istate_ncfile= trim(istate_dir)//'iced.'//trim(year)//'-01-01-00000.nc'
     s = 1
     stat(s) = NF90_OPEN(istate_ncfile , NF90_NOWRITE, ncid_in)

     WRITE (*,'(/9x, a, 3x, a)') '2D initial state estimate file:', istate_ncfile

     DO i = 1, s
        IF (stat(i) .NE. NF90_NOERR) THEN
           WRITE(*,'(/9x, a, 3x, a)') &
                'NetCDF error in opening initial state file:', istate_ncfile
           CALL abort_parallel()
        END IF
     END DO

     ! ******************************************
     ! *** Read file containing initial state ***
     ! ******************************************

     ALLOCATE (var2d(nx_global, ny_global))

     ! Calculate offsets in case not already calculated
     CALL calc_2d_offset()

     ! Loop over all state variables using state variable .NC id array
     id2d: DO idx = 1, size(id2d_list)
        s=1
        stat(s) = NF90_INQ_VARID(ncid_in, trim(id2d_list(idx)), id_2dvar)

        ! Read state variable data from file
        pos = (/ 1, 1 /)
        cnt = (/ nx_global , ny_global /)
        s = s + 1
        stat(s) = NF90_GET_VAR(ncid_in, id_2dvar, var2d, start=pos, count=cnt)

        DO i = 1, s
           IF (stat(i) .NE. NF90_NOERR) THEN
              WRITE(*,'(/9x, a, 3x, a)') &
                   'NetCDF error in reading initial state file, var=', &
                   trim(id2d_list(idx))
              CALL abort_parallel()
           END IF
        END DO

        ! Write fields into state vector
        DO j = 1,ny_global
           DO i = 1,nx_global
              ens_p(i+(j-1)*nx_global + var2d_offset(idx), member) = var2d(i,j)
           END DO
        END DO
     END DO id2d

     DEALLOCATE(var2d)

     ! *******************************************
     ! *** Close file containing initial state ***
     ! *******************************************

     s = 1
     stat(s) = NF90_CLOSE(ncid_in)

     DO i = 1, s
        IF (stat(i) .NE. NF90_NOERR) THEN
           WRITE(*,'(/9x, a, 3x, a)') &
                'NetCDF error in closing initial state file:', istate_ncfile
           CALL abort_parallel()
        END IF
     END DO
  END DO member2d

END SUBROUTINE fill2d_ensarray

SUBROUTINE fill3d_ensarray(dim_p, dim_ens, ens_p)

! !DESCRIPTION:
! Fill ensemble array with 3d state variables from initial state file.

! !USES:
  USE netcdf
  USE ice_calendar, &
       ONLY: year_init
  USE mod_assimilation, &
       ONLY: istate_dir
  USE ice_domain_size, &
       ONLY: nx_global, ny_global, ncat
  USE mod_parallel_pdaf, &
       ONLY: abort_parallel

  IMPLICIT NONE

! !ARGUMENTS
  INTEGER, INTENT(in) :: dim_p                     ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                   ! Size of ensemble
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)     ! PE-local state ensemble

! *** local variables ***
  INTEGER :: s, i, j, k, idx, member,yr ! Counters
  INTEGER :: stat(20000)                ! Status flag for NetCDF commands
  INTEGER :: pos(3),cnt(3)              ! Vectors for 3D reading fields
  INTEGER :: ncid_in                    ! ID for NetCDF file
  INTEGER :: id_3dvar                   ! IDs for fields
  CHARACTER(len=100) :: istate_ncfile     ! File holding initial state estimate
  CHARACTER(len=100) :: year              ! Year of initial state estimate
  REAL, ALLOCATABLE :: var3d(:,:,:)       ! Array for reading state variables


  member3d:DO member = 1 , dim_ens
     ! ******************************************
     ! *** Open file containing initial state ***
     ! ******************************************

     yr = year_init + member - 1
     WRITE(year, '(i4)') yr
     istate_ncfile= trim(istate_dir)//'iced.'//trim(year)//'-01-01-00000.nc'
     s = 1
     stat(s) = NF90_OPEN(istate_ncfile , NF90_NOWRITE, ncid_in)

     WRITE (*,'(/9x, a, 3x, a)') '3D initial state estimate file:', istate_ncfile

     DO i = 1, s
        IF (stat(i) .NE. NF90_NOERR) THEN
           WRITE(*,'(/9x, a, 3x, a)') &
                'NetCDF error in opening initial state file:', istate_ncfile
           CALL abort_parallel()
        END IF
     END DO

     ! ******************************************
     ! *** Read file containing initial state ***
     ! ******************************************

     ALLOCATE (var3d(nx_global, ny_global, ncat))

     ! Calculate offsets in case not already calculated
     CALL calc_3d_offset()

     ! Loop over all state variables using state variable .NC id array
     id3d: DO idx = 1, size(id3d_list)
        s=1
        stat(s) = NF90_INQ_VARID(ncid_in, trim(id3d_list(idx)), id_3dvar)

        ! Read state variable data from file
        pos = (/ 1, 1, 1 /)
        cnt = (/ nx_global , ny_global , ncat /)
        s = s + 1
        stat(s) = NF90_GET_VAR(ncid_in, id_3dvar, var3d, start=pos, count=cnt)

        DO i = 1, s
           IF (stat(i) .NE. NF90_NOERR) THEN
              WRITE(*,'(/9x, a, 3x, a)') &
                   'NetCDF error in reading initial state file, var=', &
                   trim(id3d_list(idx))
              CALL abort_parallel()
           END IF
        END DO

        ! Write fields into state vector:
        DO k = 1,ncat
           DO j = 1,ny_global
              DO i = 1,nx_global
                 ens_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                      var3d_offset(idx), member) = var3d(i,j,k)
              END DO
           END DO
        END DO
     END DO id3d

     DEALLOCATE(var3d)

     ! *******************************************
     ! *** Close file containing initial state ***
     ! *******************************************

     s = 1
     stat(s) = NF90_CLOSE(ncid_in)

     DO i = 1, s
        IF (stat(i) .NE. NF90_NOERR) THEN
           WRITE(*,'(/9x, a, 3x, a)') &
                'NetCDF error in closing initial state file:', istate_ncfile
           CALL abort_parallel()
        END IF
     END DO
  END DO member3d

END SUBROUTINE fill3d_ensarray

SUBROUTINE fill2d_statevector(dim_p, state_p)

! !DESCRIPTION:
! Fill state vector with 2d state variables.

! !USES:
  USE netcdf
  USE ice_domain_size, &
       ONLY: nx_global, ny_global
  USE ice_state, &
       ONLY: uvel, vvel
  USE ice_flux, &
       ONLY: stressp_1, stressp_2, stressp_3, stressp_4, stressm_1, &
       stressm_2, stressm_3, stressm_4, stress12_1, stress12_2, &
       stress12_3, stress12_4, sst
#ifdef USE_STRESS
  USE ice_dyn_eap, &
       ONLY: a11_1, a11_2, a11_3, a11_4, a12_1, a12_2, a12_3, a12_4
#endif

  IMPLICIT NONE

! !ARGUMENTS
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state

! *** local variables ***
  INTEGER :: i, j        ! Counters


  ! Calculate offsets in case not already calculated
  CALL calc_2d_offset()

   ! Fill state vector with 2d variables
  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + uvel_offset) = uvel(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + vvel_offset) = vvel(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + stressp_1_offset) = stressp_1(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + stressp_2_offset) = stressp_2(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + stressp_3_offset) = stressp_3(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + stressp_4_offset) = stressp_4(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + stressm_1_offset) = stressm_1(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + stressm_2_offset) = stressm_2(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + stressm_3_offset) = stressm_3(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + stressm_4_offset) = stressm_4(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + stress12_1_offset) = stress12_1(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + stress12_2_offset) = stress12_2(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + stress12_3_offset) = stress12_3(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + stress12_4_offset) = stress12_4(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + sst_offset) = sst(i+1,j+1,1)
     END DO
  END DO

#ifdef USE_STRESS
  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + a11_1_offset) = a11_1(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + a11_2_offset) = a11_2(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + a11_3_offset) = a11_3(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + a11_4_offset) = a11_4(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + a12_1_offset) = a12_1(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + a12_2_offset) = a12_2(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + a12_3_offset) = a12_3(i+1,j+1,1)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        state_p(i+(j-1)*nx_global + a12_4_offset) = a12_4(i+1,j+1,1)
     END DO
  END DO
#endif

END SUBROUTINE fill2d_statevector

SUBROUTINE fill3d_statevector(dim_p, state_p)

! !DESCRIPTION:
! Fill state vector with 3d state variables.

! !USES:
  USE netcdf
  USE ice_domain_size, &
       ONLY: nx_global, ny_global, ncat
  USE ice_state, &
       ONLY: aicen, vicen, vsnon, trcrn, nt_Tsfc, nt_iage, &
       nt_FY, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, &
       nt_sice, nt_qice, nt_qsno


  IMPLICIT NONE

! !ARGUMENTS
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state

! *** local variables ***
  INTEGER :: i, j, k        ! Counters


  ! Calculate offsets in case not already calculated
  CALL calc_3d_offset()

  ! Fill state vector with 3d variables
  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                aicen_offset) = aicen(i+1,j+1,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                vicen_offset) = vicen(i+1,j+1,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                vsnon_offset) = vsnon(i+1,j+1,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                Tsfcn_offset) = trcrn(i+1,j+1,nt_Tsfc,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                iage_offset) = trcrn(i+1,j+1,nt_iage,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                FY_offset) = trcrn(i+1,j+1,nt_FY,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                alvl_offset) = trcrn(i+1,j+1,nt_alvl,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                vlvl_offset) = trcrn(i+1,j+1,nt_vlvl,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                apnd_offset) = trcrn(i+1,j+1,nt_apnd,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                hpnd_offset) = trcrn(i+1,j+1,nt_hpnd,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                ipnd_offset) = trcrn(i+1,j+1,nt_ipnd,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice001_offset) = trcrn(i+1,j+1,nt_sice,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice001_offset) = trcrn(i+1,j+1,nt_qice,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice002_offset) = trcrn(i+1,j+1,nt_sice+1,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice002_offset) = trcrn(i+1,j+1,nt_qice+1,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice003_offset) = trcrn(i+1,j+1,nt_sice+2,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice003_offset) = trcrn(i+1,j+1,nt_qice+2,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice004_offset) = trcrn(i+1,j+1,nt_sice+3,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice004_offset) = trcrn(i+1,j+1,nt_qice+3,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice005_offset) = trcrn(i+1,j+1,nt_sice+4,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice005_offset) = trcrn(i+1,j+1,nt_qice+4,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice006_offset) = trcrn(i+1,j+1,nt_sice+5,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice006_offset) = trcrn(i+1,j+1,nt_qice+5,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice007_offset) = trcrn(i+1,j+1,nt_sice+6,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice007_offset) = trcrn(i+1,j+1,nt_qice+6,k,1)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qsno001_offset) = trcrn(i+1,j+1,nt_qsno,k,1)
        END DO
     END DO
  END DO

END SUBROUTINE fill3d_statevector

SUBROUTINE distrib2d_statevector(dim_p, state_p)

! !DESCRIPTION:
! Distribute state vector 2d state variables.

! !USES:
  USE netcdf
  USE ice_domain_size, &
       ONLY: nx_global, ny_global
  USE ice_state, &
       ONLY: uvel, vvel
  USE ice_flux, &
       ONLY: stressp_1, stressp_2, stressp_3, stressp_4, stressm_1, &
       stressm_2, stressm_3, stressm_4, stress12_1, stress12_2, &
       stress12_3, stress12_4, sst
#ifdef USE_STRESS
  USE ice_dyn_eap, &
       ONLY: a11_1, a11_2, a11_3, a11_4, a12_1, a12_2, a12_3, a12_4
#endif

  IMPLICIT NONE

! !ARGUMENTS
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state

! *** local variables ***
  INTEGER :: i, j        ! Counters


  ! Calculate offsets in case not already calculated
  CALL calc_2d_offset()

  ! Distribute state vector 2d variables
  DO j = 1,ny_global
     DO i = 1,nx_global
        uvel(i+1,j+1,1) = state_p(i+(j-1)*nx_global + uvel_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        vvel(i+1,j+1,1) = state_p(i+(j-1)*nx_global + vvel_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressp_1(i+1,j+1,1) = state_p(i+(j-1)*nx_global + stressp_1_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressp_2(i+1,j+1,1) = state_p(i+(j-1)*nx_global + stressp_2_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressp_3(i+1,j+1,1) = state_p(i+(j-1)*nx_global + stressp_3_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressp_4(i+1,j+1,1) = state_p(i+(j-1)*nx_global + stressp_4_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressm_1(i+1,j+1,1) = state_p(i+(j-1)*nx_global + stressm_1_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressm_2(i+1,j+1,1) = state_p(i+(j-1)*nx_global + stressm_2_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressm_3(i+1,j+1,1) = state_p(i+(j-1)*nx_global + stressm_3_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressm_4(i+1,j+1,1) = state_p(i+(j-1)*nx_global + stressm_4_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stress12_1(i+1,j+1,1) = state_p(i+(j-1)*nx_global + stress12_1_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stress12_2(i+1,j+1,1) = state_p(i+(j-1)*nx_global + stress12_2_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stress12_3(i+1,j+1,1) = state_p(i+(j-1)*nx_global + stress12_3_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stress12_4(i+1,j+1,1) = state_p(i+(j-1)*nx_global + stress12_4_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        sst(i+1,j+1,1) = state_p(i+(j-1)*nx_global + sst_offset)
     END DO
  END DO

#ifdef USE_STRESS
  DO j = 1,ny_global
     DO i = 1,nx_global
        a11_1(i+1,j+1,1) = state_p(i+(j-1)*nx_global + a11_1_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a11_2(i+1,j+1,1) = state_p(i+(j-1)*nx_global + a11_2_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a11_3(i+1,j+1,1) = state_p(i+(j-1)*nx_global + a11_3_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a11_4(i+1,j+1,1) = state_p(i+(j-1)*nx_global + a11_4_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a12_1(i+1,j+1,1) = state_p(i+(j-1)*nx_global + a12_1_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a12_2(i+1,j+1,1) = state_p(i+(j-1)*nx_global + a12_2_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a12_3(i+1,j+1,1) = state_p(i+(j-1)*nx_global + a12_3_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a12_4(i+1,j+1,1) = state_p(i+(j-1)*nx_global + a12_4_offset)
     END DO
  END DO
#endif

END SUBROUTINE distrib2d_statevector

SUBROUTINE distrib3d_statevector(dim_p, state_p)

  ! !DESCRIPTION:
  ! Distribute state vector 3d state variables.

  ! !USES:
  USE netcdf
  USE ice_domain_size, &
       ONLY: nx_global, ny_global, ncat
  USE ice_state, &
       ONLY: aicen, vicen, vsnon, trcrn, nt_Tsfc, nt_iage, &
       nt_FY, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, &
       nt_sice, nt_qice, nt_qsno, aice, aice0, vice
  USE ice_constants, &
       ONLY: c0, puny


  IMPLICIT NONE

  ! !ARGUMENTS
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state

  ! *** local variables ***
  INTEGER :: i, j, k        ! Counters
  REAL :: catsum            ! Sum of categories
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?

  ! Calculate offsets in case not already calculated
  CALL calc_3d_offset()

  ! ************************************
  ! Distribute state vector 3d variables
  ! ************************************

  ! If routine is called for first time then don't need to worry about PDAF
  ! creating/destroying ice.
  IF(firsttime) THEN
     DO k =1,ncat
        DO j = 1,ny_global
           DO i = 1,nx_global
              aicen(i+1,j+1,k,1) = &
                   state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                   aicen_offset)
           END DO
        END DO
     END DO

     DO k =1,ncat
        DO j = 1,ny_global
           DO i = 1,nx_global
              vicen(i+1,j+1,k,1) = &
                   state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                   vicen_offset)
           END DO
        END DO
     END DO

     DO k =1,ncat
        DO j = 1,ny_global
           DO i = 1,nx_global
              vsnon(i+1,j+1,k,1) = &
                   state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                   vsnon_offset)
           END DO
        END DO
     END DO

     ! Update logical switch
     firsttime = .FALSE.

  ! Prevent PDAF from creating/destroying ice.
  ELSE
     DO k =1,ncat
        DO j = 1,ny_global
           DO i = 1,nx_global
              ! Don't let PDAF create ice
              IF (aicen(i+1,j+1,k,1) < puny) THEN
                 state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                      aicen_offset) = aicen(i+1,j+1,k,1)
                 ! Don't let PDAF destroy ice
              ELSE IF (state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global+aicen_offset) < &
                   puny) THEN
                 state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                      aicen_offset) = aicen(i+1,j+1,k,1)
                 ! Otherwise update CICE with PDAF
              ELSE
                 aicen(i+1,j+1,k,1) = &
                      state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                      aicen_offset)
              END IF
           END DO
        END DO
     END DO

     DO k =1,ncat
        DO j = 1,ny_global
           DO i = 1,nx_global
              ! If aicen category is 0, set vicen category to 0
              IF (aicen(i+1,j+1,k,1) < puny) THEN
                 vicen(i+1,j+1,k,1) = c0
                 ! Don't let PDAF create ice
              ELSE IF (vicen(i+1,j+1,k,1) < puny) THEN
                 state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                      vicen_offset) = vicen(i+1,j+1,k,1)
                 ! Don't let PDAF destroy ice
              ELSE IF (state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global+vicen_offset) < &
                   puny) THEN
                 state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                      vicen_offset) = vicen(i+1,j+1,k,1)
                 ! Otherwise update CICE with PDAF
              ELSE
                 vicen(i+1,j+1,k,1) = &
                      state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                      vicen_offset)
              END IF
           END DO
        END DO
     END DO

     DO k =1,ncat
        DO j = 1,ny_global
           DO i = 1,nx_global
              ! If aicen category is 0, set vsnon category to 0
              IF (aicen(i+1,j+1,k,1) < puny) THEN
                 vsnon(i+1,j+1,k,1) = c0
                 ! Don't let PDAF create ice
              ELSE IF (vsnon(i+1,j+1,k,1) < puny) THEN
                 state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                      vsnon_offset) = vsnon(i+1,j+1,k,1)
                 ! Don't let PDAF destroy ice
              ELSE IF (state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global+vsnon_offset) < &
                   puny) THEN
                 state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                      vsnon_offset) = vsnon(i+1,j+1,k,1)
                 ! Otherwise update CICE with PDAF
              ELSE
                 vsnon(i+1,j+1,k,1) = &
                      state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                      vsnon_offset)
              END IF
           END DO
        END DO
     END DO
  END IF

! *******************************************************
! Update remaining state variables in same way regardless
! of whether routine is called for first time
! *******************************************************

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_Tsfc,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                Tsfcn_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_iage,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                iage_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_FY,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                FY_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_alvl,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                alvl_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_vlvl,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                vlvl_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_apnd,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                apnd_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_hpnd,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                hpnd_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_ipnd,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                ipnd_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice001_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_qice,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice001_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice+1,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice002_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_qice+1,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice002_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice+2,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice003_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_qice+2,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice003_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice+3,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice004_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_qice+3,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice004_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice+4,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice005_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_qice+4,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice005_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice+5,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice006_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_qice+5,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice006_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice+6,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice007_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_qice+6,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qice007_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_qsno,k,1) = &
                state_p(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                qsno001_offset)
        END DO
     END DO
  END DO

END SUBROUTINE distrib3d_statevector

SUBROUTINE statevar_brutemod()

! !DESCRIPTION:
  ! Temporary location for brute force modifications to state variables.
  ! Also contains updates of some state variable derived quantities.

! !USES:
  USE ice_domain_size, &
       ONLY: nx_global, ny_global, ncat
  USE ice_state, &
       ONLY: aicen
  USE ice_constants, &
         ONLY: c1

  IMPLICIT NONE

! !ARGUMENTS

! *** local variables ***
  INTEGER :: i, j           ! Counters
  REAL :: catsum            ! Sum of categories

  ! Confirm that probability distribution does not sum to greater than one.
  DO j = 1,ny_global
     DO i = 1,nx_global
        catsum = SUM(aicen(i+1,j+1,1:ncat,1))
        IF (catsum > c1) THEN
           aicen(i+1,j+1,1:ncat,1) = aicen(i+1,j+1,1:ncat,1)/catsum
        END IF
     END DO
  END DO

END SUBROUTINE statevar_brutemod

END MODULE mod_statevector
