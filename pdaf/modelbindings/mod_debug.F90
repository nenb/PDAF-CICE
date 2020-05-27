MODULE mod_debug

! !DESCRIPTION:
! This module provides variables & routines for
! debugging PDAF-CICE.

! !USES:
  IMPLICIT NONE
  SAVE

! Coordiantes of gridpoint for debugging
  INTEGER :: i_gp_debug = 26
  INTEGER :: j_gp_debug = 12

CONTAINS

  SUBROUTINE update_output (i_gp, j_gp, dim_p, state_p)

! !DESCRIPTION:
! This routine outputs original CICE state variable values along
! with PDAF updated state variable values at a single gridpoint.

!   !USES:
    USE mod_statevector, &
         ONLY: calc_2d_offset, calc_3d_offset, aicen_offset, vicen_offset
    USE ice_domain_size, &
         ONLY: nx_global, ny_global, ncat
    USE ice_state, &
         ONLY: uvel, vvel, aicen, vicen, vsnon, trcrn, nt_Tsfc, nt_iage, &
       nt_FY, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, &
       nt_sice, nt_qice, nt_qsno
    USE ice_flux, &
         ONLY: stressp_1, stressp_2, stressp_3, stressp_4, stressm_1, &
         stressm_2, stressm_3, stressm_4, stress12_1, stress12_2, &
         stress12_3, stress12_4, sst
#ifdef USE_STRESS
    USE ice_dyn_eap, &
         ONLY: a11_1, a11_2, a11_3, a11_4, a12_1, a12_2, a12_3, a12_4
#endif

    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(in) :: i_gp   ! i-coordinate of gridpoint in statevector
    INTEGER, INTENT(in) :: j_gp   ! j-coordinate of gridpoint in statevector
    INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
    REAL, INTENT(in)    :: state_p(dim_p)          ! PE-local model state

    ! Local variables
    INTEGER :: k  ! Counter
    LOGICAL :: exist


    ! Calculate offsets in case not already calculated
    CALL calc_2d_offset()
    CALL calc_3d_offset()

    ! INTRODUCE SOME SCREEN TEXT EXPLAINING WHAT IS OUTPUT. ALSO OUTPUT
    ! GRIDPOINT VALUE FOR DEBUG OUTPUT AND WHAT ENSEMBLE MEMBER IS PRESENT.

    INQUIRE(file="debug.txt", exist=exist)
    IF (exist) THEN
       OPEN(11, file="debug.txt", status="old", position="append", action="write")
    ELSE
       OPEN(11, file="debug.txt", status="new", action="write")
    END IF

    DO k = 1, ncat
       WRITE(11,'(/9x, a, i1.1, a, 3x, f10.7)') 'aicen original (k=', k,'):', &
            aicen(i_gp+1,j_gp+1,k,1)
       WRITE(11,'(/9x, a, i1.1, a, 3x, f10.7)') 'aicen update(k=', k,'):', &
            state_p(i_gp+(j_gp-1)*nx_global+(k-1)*nx_global*ny_global + aicen_offset)
       WRITE(11,'(/9x, a, i1.1, a, 3x, f10.7)') 'vicen original (k=', k,'):', &
            vicen(i_gp+1,j_gp+1,k,1)
       WRITE(11,'(/9x, a, i1.1, a, 3x, f10.7)') 'vicen update(k=', k,'):', &
            state_p(i_gp+(j_gp-1)*nx_global+(k-1)*nx_global*ny_global + vicen_offset)
    END DO

    CLOSE(11)

  END SUBROUTINE update_output

  SUBROUTINE modif_update_output (i_gp, j_gp, dim_p, state_p)

! !DESCRIPTION:
  ! This routine outputs modified PDAF update values for the different
  ! state variables at a single gridpoint (e.g. update values with negative
  ! sea-ice concentration are modified) at a single gridpoint.

!   !USES:
    USE mod_statevector, &
         ONLY: calc_2d_offset, calc_3d_offset, aicen_offset, vicen_offset
    USE ice_domain_size, &
         ONLY: nx_global, ny_global, ncat
    USE ice_state, &
         ONLY: uvel, vvel, aicen, vicen, vsnon, trcrn, nt_Tsfc, nt_iage, &
       nt_FY, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, &
       nt_sice, nt_qice, nt_qsno
    USE ice_flux, &
         ONLY: stressp_1, stressp_2, stressp_3, stressp_4, stressm_1, &
         stressm_2, stressm_3, stressm_4, stress12_1, stress12_2, &
         stress12_3, stress12_4, sst
#ifdef USE_STRESS
    USE ice_dyn_eap, &
         ONLY: a11_1, a11_2, a11_3, a11_4, a12_1, a12_2, a12_3, a12_4
#endif

    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(in) :: i_gp   ! i-coordinate of gridpoint in statevector
    INTEGER, INTENT(in) :: j_gp   ! j-coordinate of gridpoint in statevector
    INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
    REAL, INTENT(in)    :: state_p(dim_p)          ! PE-local model state

    ! Local variables
    INTEGER :: k  ! Counter

    ! Calculate offsets in case not already calculated
    CALL calc_2d_offset()
    CALL calc_3d_offset()

    DO k = 1, ncat
       WRITE(*,'(/9x, a, i1.1, a, 3x, f10.7)') 'aicen modified update(k=', k,'):', &
            state_p(i_gp+(j_gp-1)*nx_global+(k-1)*nx_global*ny_global + aicen_offset)
       WRITE(*,'(/9x, a, i1.1, a, 3x, f10.7)') 'vicen modified update(k=', k,'):', &
            state_p(i_gp+(j_gp-1)*nx_global+(k-1)*nx_global*ny_global + vicen_offset)
    END DO

  END SUBROUTINE modif_update_output

  ! TO DO
  ! - i_gp and j_gp as namelist values and remove hardcoded values in this module; also link to PDAF-OMI debug gridpoint
  ! - write modified update output routine (see above paragraph)
  ! - introduce namelist debug flag for switching these routines on/off
  ! - determine whether screen output is helpful for individ

END MODULE mod_debug
