MODULE PDAFomi_obs_op_extra
!-------------------------------------------------------------------------------
! Extra observation operators required for CICE

  USE PDAFomi_obs_f, ONLY: obs_f

CONTAINS

  SUBROUTINE PDAFomi_obs_op_f_ice_thickness(thisobs, nrows, state_p, obs_f_all, offset_obs)

    USE ice_domain_size, &
         ONLY: ncat
    USE ice_constants, &
         ONLY: c0, puny

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be averaged
    REAL, INTENT(in)    :: state_p(:)      !< PE-local model state (dim_p)
    REAL, INTENT(inout) :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)
    INTEGER, INTENT(inout) :: offset_obs   !< Offset of current observation in overall observation vector

! *** Local variables ***
    INTEGER :: i, row                  ! Counter
    REAL, ALLOCATABLE :: ostate1_p(:)  ! temporary quantity aicen*vicen
    REAL, ALLOCATABLE :: ostate2_p(:)  ! temporary quantity aice


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    obs_exist:IF (thisobs%dim_obs_p>0) THEN
       ALLOCATE(ostate1_p(thisobs%dim_obs_p))
       ALLOCATE(ostate2_p(thisobs%dim_obs_p))
       ! Initialize observed part of state vector by summing categories
       ! and forming quotient -> hi = vicen / aicen
       DO i = 1, thisobs%dim_obs_p
          ostate1_p(i) = c0
          ostate2_p(i) = c0
          DO row = 1, nrows
             ostate1_p(i) = ostate1_p(i) + &
                  state_p(thisobs%id_obs_p(row,i)) * state_p(thisobs%id_obs_p(row+nrows,i))
             ostate2_p(i) = ostate2_p(i) + state_p(thisobs%id_obs_p(row,i))
          END DO
          IF (ostate1_p(i) > puny .AND. ostate2_p(i) > puny) THEN
             obs_f_all(offset_obs+i) = ostate1_p(i) / ostate2_p(i)
          ELSE
             obs_f_all(offset_obs+i) = c0
          END IF
       END DO

       ! *** Store offset (mandatory!)
       thisobs%off_obs_f = offset_obs

       ! Increment offset in observaton vector
       offset_obs = offset_obs + thisobs%dim_obs_f

       DEALLOCATE(ostate1_p,ostate2_p)
    END IF obs_exist

  END SUBROUTINE PDAFomi_obs_op_f_ice_thickness

END MODULE PDAFomi_obs_op_extra
