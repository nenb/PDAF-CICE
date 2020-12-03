MODULE PDAFomi_obs_op_extra
!-------------------------------------------------------------------------------
! Extra observation operators required for CICE

  USE PDAFomi_obs_f, ONLY: obs_f, PDAFomi_gather_obsstate

CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Ice Concentration Observation Operator !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE PDAFomi_obs_op_ice_concen(thisobs, nrows, state_p, obs_f_all)

    USE ice_constants, &
         ONLY: c0, puny

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be summed
    REAL, INTENT(in)    :: state_p(:)      !< PE-local model state (dim_p)
    REAL, INTENT(inout) :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)

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
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_ice_concen - thisobs%id_obs_p is not allocated'
       END IF

       ! *** PE-local: Initialize observed part state vector by averaging

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
          ALLOCATE(ostate1_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
          ALLOCATE(ostate1_p(1))
       END IF

       ! Construct observed part of state vector by summing categories
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

       ! *** Global: Gather full observed state vector
       CALL PDAFomi_gather_obsstate(thisobs, ostate_p, obs_f_all)

       ! *** Clean up
       DEALLOCATE(ostate_p, ostate1_p)

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_ice_concen

END MODULE PDAFomi_obs_op_extra
