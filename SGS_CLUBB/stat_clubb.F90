#ifdef CLUBB
! $Id: stat_clubb.F90 1555 2014-08-23 15:07:29Z dschanen@uwm.edu $
module stat_clubb

  implicit none

  public :: stats_clubb_update

#ifdef SILHS
  public stats_clubb_silhs_update
#endif

  private

  contains 
!---------------------------------------------------------------------------------------------------
  subroutine stats_clubb_update( upwp, vpwp, up2, vp2, wprtp, wpthlp, &
    wp2, wp3, rtp2, thlp2, rtpthlp, cloud_frac, rcm, um, vm )

! Description:
!   Update statistics for CLUBB variables
!
! References:
!   None
!---------------------------------------------------------------------------------------------------
  use grid, only: nx, ny, nzm, nz, dimx1_s, dimx2_s, dimy1_s, dimy2_s

  use hbuffer, only: hbuf_put, hbuf_avg_put

  ! Modules from CLUBB
  use clubb_api_module, only: &
    core_rknd, & ! Constant
    lin_interpolate_two_points_api, & ! Procedure(s)
    gr

  implicit none

  real(kind=core_rknd), dimension(nx, ny, nz), intent(in) :: &
    upwp,        &! u'w'                          [m^2/s^2]
    vpwp,        &! u'w'                          [m^2/s^2]
    up2,         &! u'^2                          [m^2/s^2]
    vp2,         &! v'^2                          [m^2/s^2]
    wprtp,       &! w' r_t'                       [(m kg)/(s kg)]
    wpthlp,      &! w' th_l'                      [(m K)/s]
    wp2,         &! w'^2                          [m^2/s^2]
    rtp2,        &! r_t'^2                        [(kg/kg)^2]
    thlp2,       &! th_l'^2                       [K^2]
    rtpthlp,     &! r_t' th_l'                    [(kg K)/kg]
    cloud_frac,  &! Cloud Fraction                [-]
    rcm           ! Cloud water                   [kg/kg]

  ! w'^3 is requires additional ghost points on the x and y dimension
  real(kind=core_rknd), dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nz), intent(in) :: &
    wp3,&    ! w'^3                       [m^3/s^3]
    um, &    ! x-wind                     [m/s]
    vm       ! y-wind                     [m/s]

  ! Local variables
  real, dimension(nzm) :: &
    upwp_avg,   &
    vpwp_avg,   &
    up2_avg,    &
    vp2_avg,    &
    wprtp_avg,  &
    wpthlp_avg, &
    wp2_avg,    &
    thlp2_avg,  &
    rtp2_avg,   &
    rtpthlp_avg !,&
!   sigma_sqd_w_avg, &
!   Kh_zt_avg,  &
!   tau_zm_avg

  real :: factor_xy

  integer :: i, j, k

  !---------------------------------------------------------
  ! CLUBB variables
  ! Notes: The variables located on the vertical velocity levels 
  ! must be interpolated for the stats grid, which is on the pressure levels.
  ! -dschanen 21 Jul 2008
  factor_xy = 1. / real( nx*ny )

  upwp_avg   = 0.0
  vpwp_avg   = 0.0
  vp2_avg    = 0.0
  up2_avg    = 0.0
  wprtp_avg  = 0.0
  wpthlp_avg = 0.0
  wp2_avg    = 0.0

  thlp2_avg   = 0.0
  rtp2_avg    = 0.0
  rtpthlp_avg = 0.0

  ! Here we omit the ghost point, since the SAM stats don't have one
  do i = 1, nx
    do j = 1, ny
      do k = 1, nzm
        upwp_avg(k) = upwp_avg(k) + real( lin_interpolate_two_points_api( &
          gr%zt(k+1), gr%zm(k+1), gr%zm(k), upwp(i,j,k+1), upwp(i,j,k) ) )
        vpwp_avg(k) = vpwp_avg(k) + real( lin_interpolate_two_points_api( &
          gr%zt(k+1), gr%zm(k+1), gr%zm(k), vpwp(i,j,k+1), vpwp(i,j,k) ) )
        vp2_avg(k) = vp2_avg(k) + real( lin_interpolate_two_points_api( &
          gr%zt(k+1), gr%zm(k+1), gr%zm(k), vp2(i,j,k+1), vp2(i,j,k) ) )
        up2_avg(k) = up2_avg(k) + real( lin_interpolate_two_points_api( &
          gr%zt(k+1), gr%zm(k+1), gr%zm(k), up2(i,j,k+1), up2(i,j,k) ) )
        wprtp_avg(k) = wprtp_avg(k) + real( lin_interpolate_two_points_api( &
          gr%zt(k+1), gr%zm(k+1), gr%zm(k), wprtp(i,j,k+1), wprtp(i,j,k) ) )
        wpthlp_avg(k) = wpthlp_avg(k) + real( lin_interpolate_two_points_api( &
          gr%zt(k+1), gr%zm(k+1), gr%zm(k), wpthlp(i,j,k+1), wpthlp(i,j,k) ) )
        wp2_avg(k) = wp2_avg(k) + real( lin_interpolate_two_points_api( &
          gr%zt(k+1), gr%zm(k+1), gr%zm(k), wp2(i,j,k+1), wp2(i,j,k) ) )
        rtp2_avg(k) = rtp2_avg(k) + real( lin_interpolate_two_points_api( &
          gr%zt(k+1), gr%zm(k+1), gr%zm(k), rtp2(i,j,k+1), rtp2(i,j,k) ) )
        thlp2_avg(k) = thlp2_avg(k) + real( lin_interpolate_two_points_api( &
          gr%zt(k+1), gr%zm(k+1), gr%zm(k), thlp2(i,j,k+1), thlp2(i,j,k) ) )
        rtpthlp_avg(k) = rtpthlp_avg(k) + real( lin_interpolate_two_points_api( &
          gr%zt(k+1), gr%zm(k+1), gr%zm(k), rtpthlp(i,j,k+1), rtpthlp(i,j,k) ) )
      end do ! k = 1..nzm
    end do ! j = 1..ny
  end do ! i = 1..nx

  ! Velocity grid variables
  call hbuf_put('UPWP', upwp_avg, factor_xy)
  call hbuf_put('VPWP', vpwp_avg, factor_xy)
  call hbuf_put('VP2', vp2_avg, factor_xy)
  call hbuf_put('UP2', up2_avg, factor_xy)
  call hbuf_put('WPRTP', wprtp_avg, factor_xy)
  call hbuf_put('WPTHLP', wpthlp_avg, factor_xy)
  call hbuf_put('WP2', wp2_avg, factor_xy)
  call hbuf_put('RTP2', rtp2_avg, factor_xy)
  call hbuf_put('THLP2', thlp2_avg, factor_xy)
  call hbuf_put('RTPTHLP', rtpthlp_avg, factor_xy)

  ! CLUBB thermodynamic grid varibles (SAM pressure levels + ghost point)
  call hbuf_avg_put('CLD_FRAC', real( cloud_frac(1:nx,1:ny,2:nz) ), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('RCM', real( rcm(1:nx,1:ny,2:nz) ), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('UM', real( um(1:nx,1:ny,2:nz) ), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('VM', real( vm(1:nx,1:ny,2:nz) ), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('WP3', real( wp3(1:nx,1:ny,2:nz) ), 1,nx, 1,ny, nzm, 1.)

  return
  end subroutine stats_clubb_update

#ifdef SILHS
!---------------------------------------------------------------------------------------------------
  subroutine stats_clubb_silhs_update( )

! Description:
!   Update statistics for CLUBB SILHS variables
!
! References:
!   None
!---------------------------------------------------------------------------------------------------
    use grid, only: nx, ny, nzm, nz

    use hbuffer, only: hbuf_put, hbuf_avg_put

    use microphysics, only: &
      nmicro_fields, mkname, index_water_vapor

    ! Modules from CLUBB
    use clubb_api_module, only: &
      core_rknd, & ! Constant
      lin_interpolate_two_points_api, & ! Procedure(s)
      gr, &
      d_variables, &
      iiPDF_chi, iiPDF_w, &
      iiPDF_rr, iiPDF_rs, iiPDF_ri, &
      iiPDF_Nr, iiPDF_Ns, iiPDF_Ni, iiPDF_Ncn, &
      iirrm, iiNrm, iirsm, iirim, & ! Variables
      iiNsm, iiNim

    use clubb_silhs_vars, only: &
      LH_microphys_calls


    use clubb_silhs_vars, only: &
      LH_rt, LH_t, X_nl_all_levs, LH_sample_point_weights, LH_t_avg_tndcy, &
      LH_micro_field_avg_tndcy

    implicit none

    ! Local Variables
    real, dimension(nx,ny,nzm) :: &
      LH_rt_weighted, &
      LH_t_weighted

    real, dimension(nx,ny,nzm,d_variables) :: &
      X_nl_all_levs_weighted

    character(len=8) :: stat_name
    integer :: indx, ivar, k

    ! ---- Begin Code ----

    ! Determine cloud weighted sample averages
    LH_rt_weighted   = 0.
    LH_t_weighted    = 0.
    X_nl_all_levs_weighted = 0.

    do indx = 1, LH_microphys_calls
      do k = 1, nzm
        LH_rt_weighted(:,:,k) = LH_rt_weighted(:,:,k) &
          + LH_rt(:,:,k,indx) * LH_sample_point_weights(:,:,indx)
        LH_t_weighted(:,:,k)  = LH_t_weighted(:,:,k) &
          + LH_t(:,:,k,indx) * LH_sample_point_weights(:,:,indx)

        do ivar = 1, d_variables
          X_nl_all_levs_weighted(:,:,k,ivar) = X_nl_all_levs_weighted(:,:,k,ivar) &
            + X_nl_all_levs(:,:,k,indx,ivar) * LH_sample_point_weights(:,:,indx)
        end do

      end do ! k = 1..nzm
    end do ! indx = 1..LH_microphys_calls

    LH_rt_weighted = LH_rt_weighted / real( LH_microphys_calls )
    LH_t_weighted = LH_t_weighted / real( LH_microphys_calls )
    X_nl_all_levs_weighted = X_nl_all_levs_weighted / real( LH_microphys_calls )

    call hbuf_avg_put( 'LH_RT', LH_rt_weighted, 1,nx, 1,ny, nzm, 1. )
    call hbuf_avg_put( 'LH_TL',  LH_t_weighted, 1,nx, 1,ny, nzm, 1. )

    do ivar = 1, d_variables
      if ( ivar  == iiPDF_chi ) then
        stat_name = "LH_S_MEL"
      else if ( ivar == iiPDF_w ) then
        stat_name = "LH_W"
      else if ( ivar == iiPDF_rr ) then
        stat_name = "LH_RR"
      else if ( ivar == iiPDF_rs ) then
        stat_name = "LH_RS"
      else if ( ivar == iiPDF_ri ) then
        stat_name = "LH_RI"
      else if ( ivar == iiPDF_Nr ) then
        stat_name = "LH_NR"
      else if ( ivar == iiPDF_Ns ) then
        stat_name = "LH_NS"
      else if ( ivar == iiPDF_Ni ) then
        stat_name = "LH_NI"
      else if ( ivar == iiPDF_Ncn ) then
        stat_name = "LH_NC"
      end if ! ivar

      call hbuf_avg_put( stat_name, X_nl_all_levs_weighted(:,:,:,ivar), 1,nx, 1,ny, nzm, 1. )
    end do

    ! Tendency averages

    call hbuf_avg_put( 'LH_TL_MC', real( LH_t_avg_tndcy ), &
                       1,nx, 1,ny, nzm, 1. )

    do ivar = 1, nmicro_fields
      if ( ivar == index_water_vapor ) then
        stat_name = 'LH_RT_MC'
      else if ( ivar == iirrm ) then
        stat_name = 'LH_RR_MC'
      else if ( ivar == iirsm ) then
        stat_name = 'LH_RS_MC'
      else if ( ivar == iirim ) then
        stat_name = 'LH_RI_MC'
      else if ( ivar == iiNim ) then
        stat_name = 'LH_NI_MC'
      else if ( ivar == iiNrm ) then
        stat_name = 'LH_NR_MC'
      else if ( ivar == iiNsm ) then
        stat_name = 'LH_NS_MC'
      else
        stat_name = ''
      end if
      if ( stat_name /= '' ) then
        call hbuf_avg_put( stat_name, &
                           real( LH_micro_field_avg_tndcy(:,:,:,ivar) ), &
                           1,nx, 1,ny, nzm, 1. )
      end if
    end do

    return
  end subroutine stats_clubb_silhs_update
#endif /* SILHS */

end module stat_clubb
#endif /* CLUBB */
