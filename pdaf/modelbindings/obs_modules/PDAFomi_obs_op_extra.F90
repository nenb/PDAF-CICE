MODULE PDAFomi_obs_op_extra
!-------------------------------------------------------------------------------
! Extra observation operators required for CICE

  USE PDAFomi_obs_f, ONLY: obs_f, PDAFomi_gather_obsstate_f

CONTAINS

  SUBROUTINE PDAFomi_obs_op_f_ice_thickness(thisobs, nrows, state_p, obs_f_all, offset_obs)

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
    REAL, ALLOCATABLE :: ostate_p(:)       ! local observed part of state vector
    REAL, ALLOCATABLE :: ostate1_p(:)  ! temporary quantity aicen*vicen
    REAL, ALLOCATABLE :: ostate2_p(:)  ! temporary quantity aice


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_f_ice_thickness - thisobs%id_obs_p is not allocated'
       END IF

       ! *** PE-local: Initialize observed part state vector by averaging

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
          ALLOCATE(ostate1_p(thisobs%dim_obs_p))
          ALLOCATE(ostate2_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
          ALLOCATE(ostate1_p(1))
          ALLOCATE(ostate2_p(1))
       END IF

       DO i = 1, thisobs%dim_obs_p
          ostate1_p(i) = c0
          ostate2_p(i) = c0
          DO row = 1, nrows
             ostate1_p(i) = ostate1_p(i) + &
                  state_p(thisobs%id_obs_p(row+nrows,i))
             ostate2_p(i) = ostate2_p(i) + state_p(thisobs%id_obs_p(row,i))
          END DO
          IF (ostate1_p(i) > puny .AND. ostate2_p(i) > puny) THEN
             ostate_p(i) = (ostate1_p(i) / ostate2_p(i))*ostate2_p(i)
          ELSE
             ostate_p(i) = c0
          END IF
       ENDDO

       ! *** Store offset (mandatory!)
       thisobs%off_obs_f = offset_obs

       ! *** Global: Gather full observed state vector
       CALL PDAFomi_gather_obsstate_f(thisobs, ostate_p, obs_f_all, offset_obs)

       ! *** Clean up
       DEALLOCATE(ostate_p, ostate1_p,ostate2_p)

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_f_ice_thickness

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Ice Concentration Observation Operator !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE PDAFomi_obs_op_f_ice_concen(thisobs, nrows, state_p, obs_f_all, offset_obs)

    USE ice_constants, &
         ONLY: c0, puny

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be summed
    REAL, INTENT(in)    :: state_p(:)      !< PE-local model state (dim_p)
    REAL, INTENT(inout) :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)
    INTEGER, INTENT(inout) :: offset_obs   !< Offset of current observation in overall observation vector

! *** Local variables ***
    INTEGER :: i, row                  ! Counter
    REAL, ALLOCATABLE :: ostate_p(:)   ! local observed part of state vector
    REAL, ALLOCATABLE :: ostate1_p(:)  ! temporary quantity aicen


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_f_ice_concen - thisobs%id_obs_p is not allocated'
       END IF

       ! *** PE-local: Initialize observed part state vector by averaging

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
          ALLOCATE(ostate1_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
          ALLOCATE(ostate1_p(1))
       END IF

       ! Initialize observed part of state vector by summing categories
       ! and forming quotient -> hi = vicen / aicen
       DO i = 1, thisobs%dim_obs_p
          ostate1_p(i) = c0
          DO row = 1, nrows
             ostate1_p(i) = ostate1_p(i) + &
                  state_p(thisobs%id_obs_p(row,i))
          END DO
          IF (ostate1_p(i) > puny) THEN
             ostate_p(i) = ostate1_p(i)
          ELSE
             ostate_p(i) = c0
          END IF
       ENDDO

       ! *** Store offset (mandatory!)
       thisobs%off_obs_f = offset_obs

       ! *** Global: Gather full observed state vector
       CALL PDAFomi_gather_obsstate_f(thisobs, ostate_p, obs_f_all, offset_obs)

       ! *** Clean up
       DEALLOCATE(ostate_p, ostate1_p)

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_f_ice_concen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Monthly Mean Total Ice Thickness Observation Operator !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE PDAFomi_obs_op_f_ice_hi_m(thisobs, nrows, state_p, obs_f_all, offset_obs)

    USE ice_constants, &
         ONLY: c0, puny

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be summed - not used here!!!
    REAL, INTENT(in)    :: state_p(:)      !< PE-local model state (dim_p)
    REAL, INTENT(inout) :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)
    INTEGER, INTENT(inout) :: offset_obs   !< Offset of current observation in overall observation vector
    
! *** Local variables ***
    INTEGER :: i                       ! Counter
    REAL, ALLOCATABLE :: ostate_p(:)   ! local observed part of state vector
    REAL, ALLOCATABLE :: ostate1_p(:)  ! temporary quantity aicen


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_f_ice_hi_m - thisobs%id_obs_p is not allocated'
       END IF

       ! *** PE-local: Initialize observed part state vector by averaging

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
          ALLOCATE(ostate1_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
          ALLOCATE(ostate1_p(1))
       END IF

       ! We already have hi_m in the state vector so obs_op is the identity
       DO i = 1, thisobs%dim_obs_p
          ostate_p(i) = state_p(thisobs%id_obs_p(1,i))
       ENDDO

       ! *** Store offset (mandatory!)
       thisobs%off_obs_f = offset_obs

       ! *** Global: Gather full observed state vector
       CALL PDAFomi_gather_obsstate_f(thisobs, ostate_p, obs_f_all, offset_obs)

       ! *** Clean up
       DEALLOCATE(ostate_p, ostate1_p)

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_f_ice_hi_m

END MODULE PDAFomi_obs_op_extra
