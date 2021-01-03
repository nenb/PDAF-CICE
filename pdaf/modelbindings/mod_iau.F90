MODULE mod_iau

! !DESCRIPTION:
! This module provides variables & routines for
! incremental analysis updates (IAU).

! !USES:
  USE ice_domain_size, ONLY: nx_global, ny_global, ncat
  USE mod_statevector, &
       ONLY: uvel_offset, vvel_offset, stressp_1_offset, &
       stressp_2_offset, stressp_3_offset, stressp_4_offset, &
       stressm_1_offset, stressm_2_offset, stressm_3_offset, &
       stressm_4_offset, stress12_1_offset, stress12_2_offset, &
       stress12_3_offset, stress12_4_offset, sst_offset, &
       hi_m_offset, aicen_offset, vicen_offset, vsnon_offset, &
       Tsfcn_offset, iage_offset, FY_offset, alvl_offset, &
       vlvl_offset, apnd_offset, hpnd_offset, ipnd_offset, &
       sice001_offset, qice001_offset, sice002_offset, &
       qice002_offset, sice003_offset, qice003_offset, &
       sice004_offset, qice004_offset, sice005_offset, &
       qice005_offset, sice006_offset, qice006_offset, &
       sice007_offset, qice007_offset, qsno001_offset, &
       hi_dist_offset, calc_2d_offset, calc_3d_offset
#ifdef USE_STRESS
  USE mod_statevector, &
       ONLY: a11_1_dim_state, a11_2_dim_state, a11_3_dim_state, &
       a11_4_dim_state, a12_1_dim_state, a12_2_dim_state, &
       a12_3_dim_state, a12_4_dim_state
#endif

  IMPLICIT NONE

  SAVE

  ! Define variables to store analysis increments
  REAL,  ALLOCATABLE    :: state_inc(:)    ! Vector holding increments

  ! Define flags for control of IAU
  LOGICAL :: iau_switch    ! Control if IAU performed or not
  LOGICAL :: iau_apply    ! Control whether IAU applied on given timestep
  LOGICAL :: iau_compute   ! Control whether IAU computed on given timestep

  ! Define variables for IAU config
  INTEGER :: iau_timesteps ! Number of timesteps IAU is spread over


CONTAINS

SUBROUTINE init_statevector_inc(dim_p, state_inc, state_p)

! !DESCRIPTION:
! Initialise the statevector increment with the forecast.

  IMPLICIT NONE

  ! !ARGUMENTS
  INTEGER, INTENT(in) :: dim_p             ! PE-local state dimension
  REAL, INTENT(inout) :: state_inc(dim_p)  ! PE-local model state increment
  REAL, INTENT(in) :: state_p(dim_p)       ! PE-local model state


  state_inc = state_p

END SUBROUTINE init_statevector_inc

SUBROUTINE compute_statevector_inc(dim_p, state_inc, state_p)

  ! !DESCRIPTION:
  ! Compute the statevector increment between forecast and analysis.

  ! !USES:
  IMPLICIT NONE

  ! !ARGUMENTS
  INTEGER, INTENT(in) :: dim_p             ! PE-local state dimension
  REAL, INTENT(inout) :: state_inc(dim_p)  ! PE-local model state increment
  REAL, INTENT(in) :: state_p(dim_p)       ! PE-local model state


  ! Compute the increment
  ! NOTE: state_inc initialised with forecast in init_statevector_inc
  state_inc = (state_p - state_inc)/iau_timesteps

END SUBROUTINE compute_statevector_inc

SUBROUTINE apply_inc(dim_p, state_inc)

! !DESCRIPTION:
! Apply the statevector increment to model
!
  ! !USES:
  USE mod_statevector, &
       ONLY: physics_check
  USE ice_blocks, &
       ONLY: nx_block, ny_block
  USE ice_itd, &       ! Update CICE aggregate quantities
       ONLY: aggregate
  USE ice_grid, &
       ONLY: tmask
  USE ice_domain, &
       ONLY: nblocks
  USE ice_domain_size, &
       ONLY: nx_global, ny_global, ncat, max_ntrcr
  USE ice_state        ! Variables required for aggregate subroutine

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p            ! PE-local state dimension
  REAL, INTENT(inout) :: state_inc(dim_p) ! PE-local state vector increment

! *** local variables ***
  INTEGER :: iblk                     ! Counter


! *****************************************
! *** Apply increments to model fields  ***
!******************************************

  CALL apply2d_inc(dim_p, state_inc)
  CALL apply3d_inc(dim_p, state_inc)

! *****************************************
! *** Adjustments after state increment ***
! *****************************************

  ! Check that PDAF updates satisfy physical laws.
  ! Modify updates that do not satisfy physical laws.
  CALL physics_check()

  ! Update aggregate quantities from CICE. Should be
  ! called AFTER physics_check.
  DO iblk = 1, nblocks
     !-------------------------------------------------------------
     ! aggregate tracers
     !-------------------------------------------------------------
     CALL aggregate (nx_block, ny_block, &
          aicen(:,:,:,iblk),  &
          trcrn(:,:,:,:,iblk),&
          vicen(:,:,:,iblk),  &
          vsnon(:,:,:,iblk),  &
          aice (:,:,  iblk),  &
          trcr (:,:,:,iblk),  &
          vice (:,:,  iblk),  &
          vsno (:,:,  iblk),  &
          aice0(:,:,  iblk),  &
          tmask(:,:,  iblk),  &
          max_ntrcr,          &
          trcr_depend)
  END DO

END SUBROUTINE apply_inc

SUBROUTINE apply2d_inc(dim_p, state_inc)

! !DESCRIPTION:
! Apply statevector increment to 2d state variables.

  ! !USES:
  USE mod_statevector, &
       ONLY: hi_m
  USE netcdf
  USE ice_state, &
       ONLY: uvel, vvel
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
  INTEGER, INTENT(in) :: dim_p              ! PE-local state dimension
  REAL, INTENT(inout) :: state_inc(dim_p)   ! PE-local model state increment

! *** local variables ***
  INTEGER :: i, j        ! Counters


  ! Calculate offsets in case not already calculated
  CALL calc_2d_offset()

  ! Distribute state vector 2d variables
  DO j = 1,ny_global
     DO i = 1,nx_global
        uvel(i+1,j+1,1) = uvel(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + uvel_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        vvel(i+1,j+1,1) = vvel(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + vvel_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressp_1(i+1,j+1,1) = stressp_1(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + stressp_1_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressp_2(i+1,j+1,1) = stressp_2(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + stressp_2_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressp_3(i+1,j+1,1) = stressp_3(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + stressp_3_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressp_4(i+1,j+1,1) = stressp_4(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + stressp_4_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressm_1(i+1,j+1,1) = stressm_1(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + stressm_1_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressm_2(i+1,j+1,1) = stressm_2(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + stressm_2_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressm_3(i+1,j+1,1) = stressm_3(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + stressm_3_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stressm_4(i+1,j+1,1) = stressm_4(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + stressm_4_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stress12_1(i+1,j+1,1) = stress12_1(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + stress12_1_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stress12_2(i+1,j+1,1) = stress12_2(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + stress12_2_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stress12_3(i+1,j+1,1) = stress12_3(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + stress12_3_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        stress12_4(i+1,j+1,1) = stress12_4(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + stress12_4_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        sst(i+1,j+1,1) = sst(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + sst_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        hi_m(i+1,j+1) = hi_m(i+1,j+1) + &
             state_inc(i+(j-1)*nx_global + hi_m_offset)
     END DO
  END DO

#ifdef USE_STRESS
  DO j = 1,ny_global
     DO i = 1,nx_global
        a11_1(i+1,j+1,1) = a11_1(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + a11_1_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a11_2(i+1,j+1,1) = a11_2(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + a11_2_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a11_3(i+1,j+1,1) = a11_3(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + a11_3_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a11_4(i+1,j+1,1) = a11_4(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + a11_4_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a12_1(i+1,j+1,1) = a12_1(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + a12_1_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a12_2(i+1,j+1,1) = a12_2(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + a12_2_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a12_3(i+1,j+1,1) = a12_3(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + a12_3_offset)
     END DO
  END DO

  DO j = 1,ny_global
     DO i = 1,nx_global
        a12_4(i+1,j+1,1) = a12_4(i+1,j+1,1) + &
             state_inc(i+(j-1)*nx_global + a12_4_offset)
     END DO
  END DO
#endif

END SUBROUTINE apply2d_inc

SUBROUTINE apply3d_inc(dim_p, state_inc)

  ! !DESCRIPTION:
  ! Apply statevector increment to 3d state variables.

  ! !USES:
  USE mod_statevector, &
       ONLY: hi_dist
  USE ice_state, &
       ONLY: aicen, vicen, vsnon, trcrn, nt_Tsfc, nt_iage, &
       nt_FY, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, &
       nt_sice, aice, aice0, vice
  USE ice_constants, &
       ONLY: c0, puny


  IMPLICIT NONE

  ! !ARGUMENTS
  INTEGER, INTENT(in) :: dim_p             ! PE-local state dimension
  REAL, INTENT(inout) :: state_inc(dim_p)  ! PE-local model state increment

  ! *** local variables ***
  INTEGER :: i, j, k        ! Counters
 
  ! Calculate offsets in case not already calculated
  CALL calc_3d_offset()

  ! *****************************************************
  ! Apply increment to enthalpy state vector 3d variables
  ! *****************************************************

  ! ***********   WARNING   *****************************
  !
  ! Routine to apply increment to enthalpy variables must
  ! be called BEFORE other variables.
  !
  ! *****************************************************

  CALL apply_enthalpies_inc(dim_p, state_inc)

  ! ***********************************************
  ! Apply increment to non-enthalpy state variables
  ! ***********************************************

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           aicen(i+1,j+1,k,1) = aicen(i+1,j+1,k,1) + &
                state_inc(i+(j-1)*nx_global + (k-1)*nx_global*ny_global + &
                aicen_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           vicen(i+1,j+1,k,1) = vicen(i+1,j+1,k,1) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                vicen_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           vsnon(i+1,j+1,k,1) = vsnon(i+1,j+1,k,1) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                vsnon_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice,k,1) = trcrn(i+1,j+1,nt_sice,k,1) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice001_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice+1,k,1) = trcrn(i+1,j+1,nt_sice+1,k,1) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice002_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice+2,k,1) = trcrn(i+1,j+1,nt_sice+2,k,1) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice003_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice+3,k,1) = trcrn(i+1,j+1,nt_sice+3,k,1) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice004_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice+4,k,1) = trcrn(i+1,j+1,nt_sice+4,k,1) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice005_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice+5,k,1) = trcrn(i+1,j+1,nt_sice+5,k,1) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice006_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_sice+6,k,1) = trcrn(i+1,j+1,nt_sice+6,k,1) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                sice007_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_Tsfc,k,1) = trcrn(i+1,j+1,nt_Tsfc,k,1) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                Tsfcn_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_iage,k,1) = trcrn(i+1,j+1,nt_iage,k,1) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                iage_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           trcrn(i+1,j+1,nt_FY,k,1) = trcrn(i+1,j+1,nt_FY,k,1) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                FY_offset)
        END DO
     END DO
  END DO

  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           hi_dist(i+1,j+1,k) = hi_dist(i+1,j+1,k) + &
                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                hi_dist_offset)
        END DO
     END DO
  END DO

  ! ****************************************************************
  ! Following state variables not incremented to avoid CICE crashes!
  ! ****************************************************************

!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_alvl,k,1) = trcrn(i+1,j+1,nt_alvl,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                alvl_offset)
!        END DO
!     END DO
!  END DO
!
!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_vlvl,k,1) = trcrn(i+1,j+1,nt_vlvl,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                vlvl_offset)
!        END DO
!     END DO
!  END DO

!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_apnd,k,1) = trcrn(i+1,j+1,nt_apnd,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                apnd_offset)
!        END DO
!     END DO
!  END DO

!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_hpnd,k,1) = trcrn(i+1,j+1,nt_hpnd,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                hpnd_offset)
!        END DO
!     END DO
!  END DO

!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_ipnd,k,1) = trcrn(i+1,j+1,nt_ipnd,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                ipnd_offset)
!        END DO
!     END DO
!  END DO


!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_qice,k,1) = trcrn(i+1,j+1,nt_qice,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                qice001_offset)
!        END DO
!     END DO
!  END DO


!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_qice+1,k,1) = trcrn(i+1,j+1,nt_qice+1,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                qice002_offset)
!        END DO
!     END DO
!  END DO


!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_qice+2,k,1) = trcrn(i+1,j+1,nt_qice+2,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                qice003_offset)
!        END DO
!     END DO
!  END DO

!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_qice+3,k,1) = trcrn(i+1,j+1,nt_qice+3,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                qice004_offset)
!        END DO
!     END DO
!  END DO


!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_qice+4,k,1) = trcrn(i+1,j+1,nt_qice+4,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                qice005_offset)
!        END DO
!     END DO
!  END DO


!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_qice+5,k,1) = trcrn(i+1,j+1,nt_qice+5,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                qice006_offset)
!        END DO
!     END DO
!  END DO


!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_qice+6,k,1) = trcrn(i+1,j+1,nt_qice+6,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                qice007_offset)
!        END DO
!     END DO
!  END DO

!  DO k =1,ncat
!     DO j = 1,ny_global
!        DO i = 1,nx_global
!           trcrn(i+1,j+1,nt_qsno,k,1) = trcrn(i+1,j+1,nt_qsno,k,1) + &
!                state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
!                qsno001_offset)
!        END DO
!     END DO
!  END DO

END SUBROUTINE apply3d_inc

SUBROUTINE apply_enthalpies_inc(dim_p, state_inc)

! !DESCRIPTION:
! Apply increment to enthalpy variables when PDAF creates ice.

! !USES:
  USE ice_domain_size, only: nilyr
  USE ice_state, &
       ONLY: aicen, trcrn, nt_qice, nt_qsno, nt_Tsfc, nt_sice
  USE ice_constants
  USE ice_zbgc_shared, only: min_salin

  IMPLICIT NONE

! !ARGUMENTS
  INTEGER, INTENT(in) :: dim_p             ! PE-local state dimension
  REAL, INTENT(inout) :: state_inc(dim_p)  ! PE-local model state increment

  ! *** local variables ***
  INTEGER :: i, j, k, m     ! Counters
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  REAL :: Ti                ! Enthalpy temp
  REAL :: Tmltz             ! Melting temperature profile of ice
  REAL :: nsal = 0.407      ! salinity constant
  REAL :: msal = 0.573      ! salinity constant
  REAL :: saltmax = 3.2     ! max salinity at ice base for BL99 (ppt)
  REAL :: zn                ! normalized ice thickness
  REAL :: sal               ! salinity profile
  REAL :: slope             ! temperature profile slope
  REAL :: Tf                ! freezing temperature at bottom of ice (top ocean)


  DO k =1,ncat
     DO j = 1,ny_global
        DO i = 1,nx_global
           IF (aicen(i+1,j+1,k,1) == c0 .AND. state_inc(i+(j-1)*nx_global + &
                (k-1)*nx_global*ny_global + aicen_offset) > puny) THEN
              Tf = -1.8
              Ti = min(-puny, trcrn(i+1,j+1,nt_Tsfc,k,1))
              trcrn(i+1,j+1,nt_qsno,k,1) = -rhos*Lfresh
              DO m=1, nilyr
                 ! FIND LINEAR TEMPERATURE PROFILE
                 slope = Tf - Ti
                 Ti = Ti + slope*(m-p5)/nilyr
                 ! FIND T_M (MELTING TEMPERATURE PROFILE)
                 zn = (m-p5)/nilyr
                 sal=(saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
                 sal= max(sal, min_salin)
                 Tmltz = -sal*depressT
                 !Tmltz = -1.8
                 ! CALC ENTHALPY AND ASSIGN
                 trcrn(i+1,j+1,nt_qice+m,k,1) = &
                      -(rhoi * (cp_ice*(Tmltz-Ti) &
                      + Lfresh*(c1-Tmltz/Ti) - cp_ocn*Tmltz))
              END DO
           END IF
           IF (aicen(i+1,j+1,k,1) == c0 .AND. state_inc(i+(j-1)*nx_global + &
                (k-1)*nx_global*ny_global + aicen_offset) > 0.1) THEN
              WRITE(*,*) 'WARNING: creating significant ice at grid cell'
              WRITE(*,*) 'i, j, ncat, new aicen:', i, j, k, &
                   state_inc(i+(j-1)*nx_global+(k-1)*nx_global*ny_global + &
                   aicen_offset)
           END IF
        END DO
     END DO
  END DO

END SUBROUTINE apply_enthalpies_inc

SUBROUTINE check_iau_compute(iau_compute)

! !DESCRIPTION:
! Routine to determine whether IAU increments should be computed.
!
! !USES:
  USE ice_calendar, &
       ONLY: mday, istep, dt
  USE ice_constants, &
       ONLY: secday

  IMPLICIT NONE

  ! !ARGUMENTS:
  LOGICAL, INTENT(inout) :: iau_compute ! Flag for computing IAU

  ! local variables
  INTEGER :: elapsed_days      ! Number of days since beginning of run


  iau_compute=.FALSE.

  elapsed_days = INT((istep * dt) / secday)

  ! Compute IAU on first day of every month
  IF (mday == 1) THEN
     ! Don't compute IAU on very first day of simulation
     IF (elapsed_days > 0) THEN
        iau_compute = .TRUE.
     END IF
  END IF

END SUBROUTINE check_iau_compute

SUBROUTINE check_iau_apply(iau_apply)

! !DESCRIPTION:
! Routine to determine whether IAU increments should be applied.
!
! !USES:
  USE ice_calendar, &
       ONLY: mday, istep, dt
  USE ice_constants, &
       ONLY: secday

  IMPLICIT NONE

  ! !ARGUMENTS:
  LOGICAL, INTENT(inout) :: iau_apply ! Flag for applying IAU

  ! local variables
  INTEGER, SAVE :: counter = 0 ! Count number of timesteps IAU applied
  INTEGER :: elapsed_days      ! Number of days since beginning of run
  LOGICAL, SAVE :: inc_computed = .FALSE. ! Flag whether IAU computed


  iau_apply=.FALSE.

  ! Count number of days since run started
  elapsed_days = INT((istep * dt) / secday)

  ! Check to see if IAU has been computed (computed on first day
  ! of month except first month)
  IF (mday == 1) THEN
     IF (elapsed_days > 0) THEN
        inc_computed = .TRUE.
     END IF
  END IF

  ! If IAU computed, check whether to apply
  IF (inc_computed) THEN
     ! Only apply IAU for iau_timesteps
     IF (counter < iau_timesteps) THEN
        iau_apply = .TRUE.
        counter = counter + 1
     ELSE
        ! Reset flag and counter
        inc_computed = .FALSE.
        counter = 0
     END IF
  END IF

END SUBROUTINE check_iau_apply

END MODULE mod_iau
