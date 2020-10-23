!$Id: prepoststep_ens_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_ens_pdaf --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
! 
! The routine is called for global filters (e.g. SEIK)
! before the analysis and after the ensemble transformation.
! For local filters (e.g. LSEIK) the routine is called
! before and after the loop over all local analysis
! domains.
! The routine provides full access to the state 
! estimate and the state ensemble to the user.
! Thus, user-controlled pre- and poststep 
! operations can be performed here. For example 
! the forecast and the analysis states and ensemble
! covariance matrix can be analyzed, e.g. by 
! computing the estimated variances. 
! For the offline mode, this routine is the place
! in which the writing of the analysis ensemble
! can be performed.
!
! If a user considers to perform adjustments to the 
! estimates (e.g. for balances), this routine is 
! the right place for it.
!
! Implementation for the 2D offline example
! with parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: screen, filtertype, subtype, forget, local_range, &
       locweight, srange, rms_obs, delt_obs
  USE mod_parallel_pdaf, &
       ONLY: mype_filter
    USE output_netcdf_asml, &
       ONLY: init_netcdf_asml, write_netcdf_asml, close_netcdf_asml
  USE ice_calendar, &
       ONLY: dt, npt, idate0, time
  USE ice_domain_size, &
       ONLY: nx_global, ny_global, ncat

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step (negative for call after forecast)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is initialised and can be used freely here (not for SEEK!)
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_X_update       (as U_prepoststep)
!EOP

! *** local variables ***
  INTEGER :: i, j, member, domain     ! Counters
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  REAL :: invdim_ens                   ! Inverse ensemble size
  CHARACTER(len=2) :: stepstr         ! String for time step
  CHARACTER(len=3) :: anastr          ! String for call type (initial, forecast, analysis)


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype_filter == 0) THEN
     IF (firsttime) THEN
        WRITE (*, '(8x, a)') 'Analyze initial state ensemble'
        anastr = 'ini'
        CALL init_netcdf_asml(idate0, dt, nx_global, ny_global, ncat, &
             filtertype, subtype, dim_ens, forget, local_range, locweight, &
             srange, rms_obs, delt_obs, npt)
        CALL close_netcdf_asml()
     ELSE
        IF (step<0) THEN
           WRITE (*, '(8x, a)') 'Analyze and write forecasted state ensemble'
           anastr = 'for'
        ELSE
           WRITE (*, '(8x, a)') 'Analyze and write assimilated state ensemble'
           anastr = 'ana'
        END IF
     END IF
  END IF

! *****************************************
! *** Calculations before analysis step ***
! *****************************************

!!$  IF (step < 0) THEN
!!$    
!!$  END IF

! ****************************************
! *** Calculations after analysis step ***
! ****************************************

!!$  IF (step > 0) THEN
!!$
!!$  END IF

! **************************************
! *** Begin statistical calculations ***
! **************************************

  ! Initialize numbers
  invdim_ens    = 1.0 / REAL(dim_ens)

  ! *** Compute mean state
  IF (mype_filter == 0) WRITE (*, '(8x, a)') '--- compute ensemble mean'

  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i, member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)

! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  ! *****
  ! Empty
  ! *****

! *****************
! *** Screen IO ***
! *****************

  ! *****
  ! Empty
  ! *****

! *******************
! *** File output ***
! *******************

  ! *** Now write analysis ensemble ***
  WRITE (*, '(8x, a)') '--- write ensemble and state estimate'

  ! Set string for time step
  IF (step>=0) THEN
     WRITE (stepstr, '(i2.2)') step
  ELSE
     WRITE (stepstr, '(i2.2)') -step
  END IF

  CALL write_netcdf_asml(anastr, step, time, dim_p, nx_global, ny_global,&
       ncat, state_p, dim_ens, ens_p)

  CALL close_netcdf_asml()

! ********************
! *** Finishing up ***
! ********************

  ! Deallocate observation arrays
  CALL deallocate_obs_pdafomi(step)

  firsttime = .FALSE.

END SUBROUTINE prepoststep_ens_pdaf
