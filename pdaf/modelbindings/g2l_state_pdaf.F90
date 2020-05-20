!$Id: g2l_state_pdaf.F90 332 2019-12-30 09:37:03Z lnerger $
!>  Restrict a model state to a local analysis domain
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during the loop over all
!! local analysis domains in PDAF_X_update
!! before the analysis on a single local analysis 
!! domain.  It has to initialize elements of the 
!! state vector for the local analysis domains from
!! the PE-local full state vector.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

  USE mod_statevector, &
       ONLY: var2d_offset, var3d_offset
  USE ice_domain_size, &
       ONLY: nx_global, ny_global, ncat

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  REAL, INTENT(in)    :: state_p(dim_p) !< PE-local full state vector 
  REAL, INTENT(out)   :: state_l(dim_l) !< State vector on local analysis domain

  ! Local variables
  INTEGER :: i        ! Counter
  INTEGER :: a, b, c  ! Variables for 3D state variable index


! *************************************
! *** Initialize local state vector ***
! *************************************
  
  ! The local state vector elements are the 2D state variables as well as
  ! each different category for the 3D state variables
  DO i = 1, dim_l
     ! First, fill the local state vector with 2D state variables.
     IF( i .LE. size(var2d_offset) ) THEN
        state_l(i) = state_p( var2d_offset(i) + domain_p )
     ELSE
        ! Now, fill the local state vector with 3D state variables.
        ! Each 3D state variable has dimension nx*ny*ncat. We first break
        ! this into chunks of size nx*ny. Then we identify each chunk as a
        ! different element of the local state vector.
        a = i - size(var2d_offset)
        b = INT( CEILING( REAL(a)/REAL(ncat) ) )
        c = MOD(a-1,ncat)
        state_l(i) = state_p( var3d_offset(b) + (nx_global*ny_global*c) + domain_p )
     END IF
  END DO

END SUBROUTINE g2l_state_pdaf
