MODULE PDAFomi_obs_op_extra
!-------------------------------------------------------------------------------
! Extra observation operators required for CICE

CONTAINS

  SUBROUTINE PDAFomi_obs_op_f_ice_thickness(dim_p, nobs_f_all, nobs_p_one, nobs_f_one, &
       id_obs_p_one, state_p, obs_f_all, offset_obs)

    USE ice_domain_size, &
         ONLY: ncat
    USE ice_constants, &
         ONLY: c0, puny

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p               !< PE-local state dimension
    INTEGER, INTENT(in) :: nobs_f_all          !< Length of obs. vector for all observations
    INTEGER, INTENT(in) :: nobs_p_one          !< PE-local number observations of current observation type
    INTEGER, INTENT(in) :: nobs_f_one          !< Full number observations of current observation type
    INTEGER, INTENT(in) :: id_obs_p_one(:, :)  !< Index of current observations in PE-local state vector (nrows, nobs_p_one)
    REAL, INTENT(in)    :: state_p(:)          !< PE-local model state (dim_p)
    REAL, INTENT(inout) :: obs_f_all(:)        !< Full observed state for all observation types (nobs_f_all)
    INTEGER, INTENT(inout) :: offset_obs       !< Offset of current observation in overall observation vector

! *** Local variables ***
    INTEGER :: i, k                    ! Counter
    REAL, ALLOCATABLE :: ostate1_p(:)  ! local observed part of aicen state vector
    REAL, ALLOCATABLE :: ostate2_p(:)  ! local observed part of vicen state vector


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    obs_exist:IF (nobs_p_one>0) THEN
       ALLOCATE(ostate1_p(nobs_p_one))
       ALLOCATE(ostate2_p(nobs_p_one))

       ! Initialize observed part of state vector by summing categories
       ! and forming quotient -> hi = vicen / aicen
       DO i = 1, nobs_p_one
          ostate1_p(i) = c0
          ostate2_p(i) = c0
          DO k = 1, ncat
             ! aicen
             ostate1_p(i) = ostate1_p(i) + state_p(id_obs_p_one(k,i))
             ! vicen
             ostate2_p(i) = ostate2_p(i) + state_p(id_obs_p_one(k+ncat,i))
          END DO
          IF (ostate1_p(i) > puny .AND. ostate2_p(i) > puny) THEN
             obs_f_all(offset_obs+i) = ostate2_p(i) / ostate1_p(i)
          ELSE
             obs_f_all(offset_obs+i) = c0
          END IF
       ENDDO

       ! Increment offset in observaton vector
       offset_obs = offset_obs + nobs_f_one

       DEALLOCATE(ostate1_p, ostate2_p)
    END IF obs_exist

  END SUBROUTINE PDAFomi_obs_op_f_ice_thickness

END MODULE PDAFomi_obs_op_extra
