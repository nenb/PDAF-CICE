!$Id: callback_obs_pdafomi.F90 496 2020-06-09 15:26:17Z lnerger $
!> callback_obs_pdafomi
!!
!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specifc routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!!
!! **Adding an observation type:**
!! When adding an observation type, one has to add one module
!! obs_TYPE_pdafomi (based on the template obs_TYPE_pdafomi_TEMPLATE.F90).
!! In addition one has to add a call to the different routines include
!! in this file. It is recommended to keep the order of the calls
!! consistent over all files. 
!! 
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_TYPE.
!!
SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)

  ! Include functions for different observations
  USE obs_ice_concen_pdafomi, ONLY: assim_ice_concen, init_dim_obs_ice_concen
  USE obs_ice_hi_m_pdafomi, ONLY: assim_ice_hi_m, init_dim_obs_ice_hi_m
  USE obs_ice_hi_dist_pdafomi, ONLY: assim_ice_hi_dist, init_dim_obs_ice_hi_dist

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step      !< Current time step
  INTEGER, INTENT(out) :: dim_obs !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_ice_concen ! Observation dimensions
  INTEGER :: dim_obs_ice_hi_m ! Observation dimensions
  INTEGER :: dim_obs_ice_hi_dist ! Observation dimensions

! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_ice_concen = 0
  dim_obs_ice_hi_m = 0
  dim_obs_ice_hi_dist = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  IF (assim_ice_concen) CALL init_dim_obs_ice_concen(step, dim_obs_ice_concen)
  IF (assim_ice_hi_m) CALL init_dim_obs_ice_hi_m(step, dim_obs_ice_hi_m)
  IF (assim_ice_hi_dist) CALL init_dim_obs_ice_hi_dist(step, dim_obs_ice_hi_dist)

  dim_obs = dim_obs_ice_concen + dim_obs_ice_hi_m + dim_obs_ice_hi_dist

END SUBROUTINE init_dim_obs_pdafomi


!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_TYPE.
!!
SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

  ! Include functions for different observations
  USE obs_ice_concen_pdafomi, ONLY: obs_op_ice_concen
  USE obs_ice_hi_m_pdafomi, ONLY: obs_op_ice_hi_m
  USE obs_ice_hi_dist_pdafomi, ONLY: obs_op_ice_hi_dist

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs            !< Dimension of full observed state
  REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
  REAL, INTENT(inout) :: ostate(dim_obs)  !< PE-local full observed state


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

  ! The order of these calls is not relevant as the setup
  ! of the overall observation vector is defined by the
  ! order of the calls in init_dim_obs_pdafomi
  CALL obs_op_ice_concen(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_ice_hi_m(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_ice_hi_dist(dim_p, dim_obs, state_p, ostate)

END SUBROUTINE obs_op_pdafomi


!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

  ! Include functions for different observations
  USE obs_ice_concen_pdafomi, ONLY: init_dim_obs_l_ice_concen
  USE obs_ice_hi_m_pdafomi, ONLY: init_dim_obs_l_ice_hi_m
  USE obs_ice_hi_dist_pdafomi, ONLY: init_dim_obs_l_ice_hi_dist


  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs  !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Call init_dim_obs_l specific for each observation
  CALL init_dim_obs_l_ice_concen(domain_p, step, dim_obs, dim_obs_l)
  CALL init_dim_obs_l_ice_hi_m(domain_p, step, dim_obs, dim_obs_l)
  CALL init_dim_obs_l_ice_hi_dist(domain_p, step, dim_obs, dim_obs_l)

END SUBROUTINE init_dim_obs_l_pdafomi


!-------------------------------------------------------------------------------

!> Call-back routine for deallocate_obs
!!
!! This routine calls the routine PDAFomi_deallocate_obs
!! for each observation type
!!
SUBROUTINE deallocate_obs_pdafomi(step)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_deallocate_obs
  ! Include observation types (rename generic name)
  USE obs_ice_concen_pdafomi, ONLY: obs_ice_concen => thisobs
  USE obs_ice_hi_m_pdafomi, ONLY: obs_ice_hi_m => thisobs
  USE obs_ice_hi_dist_pdafomi, ONLY: obs_ice_hi_dist => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step   !< Current time step


! *************************************
! *** Deallocate observation arrays ***
! *************************************

  CALL PDAFomi_deallocate_obs(obs_ice_concen)
  CALL PDAFomi_deallocate_obs(obs_ice_hi_m)
  CALL PDAFomi_deallocate_obs(obs_ice_hi_dist)

END SUBROUTINE deallocate_obs_pdafomi
