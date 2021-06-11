!$Id: assimilate_pdaf.F90 1411 2013-09-25 14:04:41Z lnerger $
!BOP
!
! !ROUTINE: assimilate_pdaf - Routine to control perform analysis step
!
! !INTERFACE:
SUBROUTINE assimilate_pdaf()

! !DESCRIPTION:
! This routine is called during the model integrations at each time 
! step. It check whether the forecast phase is completed. If so, 
! PDAF_put_state_X is called to perform the analysis step.
!
! !REVISION HISTORY:
! 2013-08 - Lars Nerger - Initial code for NEMO
! Later revisions - see svn log
!
  ! !USES:
  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAFomi_assimilate_local, PDAF_get_localfilter
  USE mod_parallel_pdaf, &     ! Parallelization variables
       ONLY: mype_world, abort_parallel
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: filtertype
    USE ice_calendar, &
	 ONLY: nyr,month

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: step
! CAlls: PDAF_assimilate_X
!EOP

! Local variables
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: localfilter          ! Flag for domain-localized filter (1=true)

  ! External subroutines
  ! Interface between model and PDAF, and prepoststep
  EXTERNAL :: collect_state_pdaf, &  ! Collect a state vector from model fields
       distribute_state_pdaf, &      ! Distribute a state vector to model fields
       next_observation_pdaf, &      ! Provide time step of next observation
       prepoststep_ens_pdaf          ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, & ! Provide number of local analysis domains
       init_dim_l_pdaf, &            ! Initialize state dimension for local analysis domain
       g2l_state_pdaf, &             ! Get state on local analysis domain from global state
       l2g_state_pdaf                ! Update global state from state on local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: &
       init_dim_obs_pdafomi, &     ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &           ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi        ! Get dimension of obs. vector for local analysis domain


! *********************************
! *** Call assimilation routine ***
! *********************************

  ! Check  whether the filter is domain-localized
  CALL PDAF_get_localfilter(localfilter)

  IF (localfilter==1) THEN
     ! here define dates for changing forgetting factor
     !IF (nyr /= 1) THEN
     !   IF (month == 1) THEN
     !      CALL PDAF_reset_forget(0.99)
     !   END IF
     !   IF (month == 2) THEN
     !      CALL PDAF_reset_forget(0.99)
     !   END IF
     !   IF (month == 3) THEN
     !      CALL PDAF_reset_forget(0.99)
     !   END IF
     !   IF (month == 4) THEN
     !      CALL PDAF_reset_forget(0.95)
     !   END IF
     !   IF (month == 5) THEN
     !      CALL PDAF_reset_forget(0.95)
     !   END IF
     !   IF (month == 6) THEN
     !      CALL PDAF_reset_forget(0.95)
     !   END IF
     !   IF (month == 7) THEN
     !      CALL PDAF_reset_forget(0.95)
     !   END IF
     !   IF (month == 8) THEN
     !      CALL PDAF_reset_forget(0.95)
     !   END IF
     !    IF (month == 9) THEN
     !      CALL PDAF_reset_forget(0.95)
     !   END IF
     !   IF (month == 10) THEN
     !      CALL PDAF_reset_forget(0.99)
     !   END IF
     !   IF (month == 11) THEN
     !      CALL PDAF_reset_forget(0.99)
     !   END IF
     !   IF (month == 12) THEN
     !      CALL PDAF_reset_forget(0.99)
     !   END IF
     !END IF
     CALL PDAFomi_assimilate_local(collect_state_pdaf,&
          distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi,&
          prepoststep_ens_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,&
          init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf,&
          next_observation_pdaf, status_pdaf)
  ELSE
     WRITE (*,'(a)') 'ERROR - global filter not implemented, stopping.'
     CALL  abort_parallel()
  END IF


! Check for errors during execution of PDAF
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF_put_state - stopping! (PE ', mype_world,')'
     CALL  abort_parallel()
  END IF

END SUBROUTINE assimilate_pdaf
