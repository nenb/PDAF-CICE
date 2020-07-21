!$Id: init_dim_l_pdaf.F90 336 2020-01-21 13:23:07Z lnerger $
!>  Set dimension of local model state
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during analysis step
!! in PDAF_X_update in the loop over all local
!! analysis domains. It has to set the dimension
!! of the local model  state on the current analysis
!! domain.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! 2013-02 - Lars Nerger - Initial code
!! Later revisions - see repository log
!!
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: coords_l
  USE ice_domain_size, &
       ONLY: nx_global, ny_global
  USE ice_grid, &
       ONLY: tlon, tlat
  USE mod_statevector, &
       ONLY: calc_local_dim
  USE ice_constants, &
       ONLY: pi

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    !< Local state dimension

  ! Local variables
  INTEGER :: i,j   ! Model grid indices


! ****************************************
! *** Initialize local state dimension ***
! ****************************************

  ! Local dimension is total number of state variables
  ! NOTE: Each category for a 3D state variable is
  ! defined as a *new* state variable.
  dim_l = calc_local_dim()


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! Coordinates are defined using T longitude/latitude grid values

  ! First, compute (i,j) grid coordinates of local domain
  j = INT(CEILING(REAL(domain_p)/REAL(nx_global)))
  i = INT(domain_p) - (j-1)*REAL(nx_global)

  ! Now, convert to T longitude/latitude grid values.
  ! NOTE: tlon and tlat are in radians, and there are ghost cells.
  coords_l(1)=tlon(i+1,j+1,1)
  coords_l(2)=tlat(i+1,j+1,1)

END SUBROUTINE init_dim_l_pdaf
