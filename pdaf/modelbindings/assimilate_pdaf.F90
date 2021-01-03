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
  USE mod_parallel_pdaf, &
       ONLY: mype_world, abort_parallel
  USE mod_assimilation, &
       ONLY: filtertype, dim_state_p
  USE mod_iau, &
       ONLY: apply_inc, check_iau_apply, state_inc, iau_switch, &
       iau_apply

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: step
! Calls: PDAF_get_localfilter
! Calls: PDAFomi_assimilate_local
! Calls: distribute_inc (if IAU enabled)
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

  ! Check whether to apply updates for IAU
  IF (iau_switch) THEN
     ! Determine whether IAU applied on this step or not
     ! NOTE: Distinct to check_iau_compute!
     CALL check_iau_apply(iau_apply)
     IF (iau_apply) THEN
        CALL apply_inc(dim_state_p, state_inc)
     END IF
  END IF

END SUBROUTINE assimilate_pdaf
