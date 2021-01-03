MODULE mod_cleanup

! !DESCRIPTION:
! Cleanup routines called from finalize_pdaf
! !USES:
  IMPLICIT NONE

CONTAINS
!------------------------------------------------------------------------------
  SUBROUTINE PDAF_CICE_deallocate()

  ! !DESCRIPTION:
  ! Deallocate arrays allocated for PDAF-CICE

  ! !USES:
    USE mod_iau, &
         ONLY: state_inc, iau_switch

    IMPLICIT NONE


    IF (iau_switch) THEN
       DEALLOCATE(state_inc)
    END IF

  END SUBROUTINE PDAF_CICE_deallocate

END MODULE mod_cleanup
