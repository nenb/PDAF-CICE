!$Id: init_n_domains_pdaf.F90 332 2019-12-30 09:37:03Z lnerger $
!> /brief  Set number of local analysis domains
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF
!!
!! The routine is called in PDAF_X_update 
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains. 
!! It has to set the number of local analysis 
!! domains for the PE-local domain.
!!
!! Implementation for the 2D online example
!! without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)

  USE ice_domain_size, &
       ONLY: nx_global, ny_global

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step        !< Current time step
  INTEGER, INTENT(out) :: n_domains_p !< PE-local number of analysis domains


! ************************************
! *** Initialize number of domains ***
! ************************************

  ! Use horizontal localization
  n_domains_p = nx_global*ny_global

END SUBROUTINE init_n_domains_pdaf
