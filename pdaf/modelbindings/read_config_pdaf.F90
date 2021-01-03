!$Id: read_config_pdaf.F90 3 2013-09-05 10:28:51Z lnerger $
!BOP
!
! !ROUTINE: read_config_pdaf - Read configuration for PDAF
!
! !INTERFACE:
SUBROUTINE read_config_pdaf()

! !DESCRIPTION:
! This routine read the namelist file with
! parameters controlling data assimilation with
! PDAF.

! !USES:
  USE mod_parallel_pdaf, &
       ONLY: mype_world
  USE mod_assimilation, &
       ONLY: filtertype, subtype, dim_ens, delt_obs, screen, &
       type_forget, forget, local_range, locweight, srange, &
       istate_dir
  USE mod_iau, &
       ONLY: iau_switch, iau_timesteps
  USE obs_ice_concen_pdafomi, &
       ONLY: rms_ice_concen
  USE output_netcdf_asml, &
       ONLY: file_asml

  IMPLICIT NONE
!EOP

! Local variables
  CHARACTER(len=100) :: nmlfile   ! name of namelist file

  NAMELIST /pdaf_nml/ istate_dir, filtertype, subtype, dim_ens,&
       delt_obs, screen, type_forget, forget,&
       local_range, locweight, srange, rms_ice_concen,&
       iau_switch, iau_timesteps, file_asml

! ****************************************************
! ***   Initialize PDAF parameters from namelist   ***
! ****************************************************

! *** Read namelist file ***
  nmlfile ='namelist.pdaf'

  OPEN (30,file=nmlfile)
  READ (30,NML=pdaf_nml)
  CLOSE (30)

! *** Print configuration variables ***
  showconf: IF (mype_world .EQ. 0) THEN

    WRITE (*,'(/1x,a)') '-- Overview of PDAF configuration --'
    WRITE (*,'(3x,a)') 'PDAF [pdaf_nml]:'
    WRITE (*,'(5x,a,i10)')    'filtertype   ', filtertype
    WRITE (*,'(5x,a,i10)')    'subtype      ', subtype
    WRITE (*,'(5x,a,i10)')    'dim_ens      ', dim_ens
    WRITE (*,'(5x,a,i10)')    'delt_obs     ', delt_obs
    WRITE (*,'(5x,a,i10)')    'screen       ', screen
    WRITE (*,'(5x,a,i10)')    'type_forget  ', type_forget
    WRITE (*,'(5x,a,f10.2)')  'forget       ', forget
    WRITE (*,'(5x,a,es10.2)') 'local_range  ', local_range
    WRITE (*,'(5x,a,i10)')    'locweight    ', locweight
    WRITE (*,'(5x,a,es10.2)') 'srange       ', srange
    WRITE (*,'(5x,a,es10.2)') 'rms_ice_concen    ', rms_ice_concen
    WRITE (*,'(5x,a,a)')      'istate_dir   ', istate_dir
    WRITE (*,'(5x,a,i10)')    'iau_timesteps     ', iau_timesteps
    WRITE (*,'(5x,a,l)')      'iau_switch   ', iau_switch
    WRITE (*,'(5x,a,a)')      'file_asml    ', file_asml
    WRITE (*,'(1x,a)') '-- End of PDAF configuration overview --'

  END IF showconf

END SUBROUTINE read_config_pdaf
