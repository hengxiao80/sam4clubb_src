#ifdef CLUBB
!-------------------------------------------------------------------------------
! $Id: clubb_sgs.F90 1629 2015-03-30 19:11:32Z raut@uwm.edu $
module clubb_sgs

! Description:
!   Contains function and subroutines for interfacing with the UW Milwaukee
!   Single-Column Model and also the CLUBB-SILHS subcolumn generator.

! References:
!   See DOC/CLUBB/clubb_doc/CLUBBeqns.pdf in this directory.
!-------------------------------------------------------------------------------

  use clubb_api_module, only: &
    setup_clubb_core_api, & ! Procedure(s0
    advance_clubb_core_api, &
    cleanup_clubb_core_api, &

    time_precision, & ! Constant(s)
    core_rknd

  use domain, only: &
    nsubdomains_x, & ! Constant(s)
    nsubdomains_y

  implicit none

  private

  public :: clubb_sgs_setup, advance_clubb_sgs, clubb_sgs_cleanup, &
            apply_clubb_sgs_tndcy_scalars, apply_clubb_sgs_tndcy_mom

  logical, private :: lstats_clubb

  integer, dimension(nsubdomains_x*nsubdomains_y), private :: &
    sample_nodes

#ifdef SILHS
  integer, private, save :: lh_iter = 0
#endif /* SILHS */
  contains
!-------------------------------------------------------------------------------
  subroutine clubb_sgs_setup( dt_clubb, z, rho, zi, rhow, tv0, tke )

! Description:
!   Initialize UWM CLUBB.

! References:
!   None
!-------------------------------------------------------------------------------

    ! From the CLUBB directory
    use clubb_api_module, only: &
      clubb_no_error, &
      set_clubb_debug_level_api, & ! Subroutine

      nparams, & ! Constants

      em_min, w_tol_sqd, rt_tol, thl_tol, &
      zero_threshold, fstderr, fstdout, pi, &

      zt2zm_api, &
      gr, & ! Derived type

      read_parameters_api, & ! Subroutine

      stats_init_api, & ! Subroutine

      l_use_boussinesq, & ! Variables
      l_tke_aniso, &

      iirrm, iiNrm, iirsm, iirim, &
      iirgm, iiNsm, iiNim, iiNgm, &
      l_mix_rat_hm, l_frozen_hm

    ! From the SAM directory
    use grid, only: rank, nx, ny, nz, nzm, dx, dy, time, case, caseid, & ! Variable(s)
      nrestart, dimx1_s, dimx2_s, dimy1_s, dimy2_s, ntracers

    use params, only: latitude0, longitude0 ! Variable(s)

    use domain, only: nx_gl, ny_gl, YES3D ! Constant(s)

    use params, only: lcond, cp ! Constants

    use sgs_params, only: doclubb_sfc_fluxes  ! Variable(s)

    use microphysics, only: &
      mkname, nmicro_fields ! Variable(s)

#ifdef SILHS
    use silhs_api_module, only: &
      l_lh_cloud_weighted_sampling

    use clubb_api_module, only: &
      d_variables, & ! Variable
      setup_corr_varnce_array_api, & ! Procedure
      l_fix_chi_eta_correlations,   &
      hydromet_list, & ! Variable(s)
      hydromet_tol, &
      hmp2_ip_on_hmm2_ip, &
      Ncnp2_on_Ncnm2, &
!     cm3_per_m3, &
      rr_tol,     &
      ri_tol,     &
      rs_tol,     &
      rg_tol,     &
      Nr_tol,     &
      Ni_tol,     &
      Ns_tol,     &
      Ng_tol,     &
      Nc_tol

#endif /* SILHS */

    use clubbvars, only: &
      upwp,   &! u'w'.                 [m^2/s^2]
      vpwp,   &! u'w'.                 [m^2/s^2]
      up2,    &! u'^2                  [m^2/s^2]
      vp2,    &! v'^2                  [m^2/s^2]
      wprtp,  &! w' r_t'.              [(m kg)/(s kg)]
      wpthlp, &! w' th_l'.             [(m K)/s]
      wprcp,  &! w' r_c'               [(kg/kg) m/s]
      wp2,    &! w'^2.                 [m^2/s^2]
      rtp2,   &! r_t'^2.               [(kg/kg)^2]
      thlp2,  &! th_l'^2.              [K^2]
      rtpthlp,&! r_t' th_l'.           [(kg K)/kg]
      wp3      ! w'^3.                 [m^3/s^3]

    use clubbvars, only: &
      tracer_tndcy, & ! Time tendency of the SAM set of tracers
      t_tndcy,  & ! CLUBB contribution to moist static energy  [K/s]
      qc_tndcy, & ! CLUBB contribution to liquid water         [kg/kg/s]
      qv_tndcy, & ! CLUBB contribution to vapor water          [kg/kg/s]
      u_tndcy,  & ! CLUBB contribution to x-wind               [m/s^2]
      v_tndcy     ! CLUBB contribution to y-wind               [m/s^2]

    use clubbvars, only: &
      sclrp2,      & ! Passive scalar variance.       [{units vary}^2]
      sclrpthlp,   & ! Passive scalar covariance.     [{units vary}^2]
      sclrprtp,    & ! Passive scalar covariance.     [{units vary}^2]
      wpsclrp        ! w'sclr'                        [units vary m/s]

    use clubbvars, only: &
      rho_ds_zm,       & ! Dry, static density on momentum levels      [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermodynamic levels [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density on momentum levels [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density on thermo. levels  [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levels  [K]
      thv_ds_zt          ! Dry, base-state theta_v on thermo. levels   [K]

    use clubbvars, only: &
      sclr_tol,   & ! Tolerance on high-order scalars
      edsclr_dim, & ! Number of eddy-diffusivity scalars
      sclr_dim      ! Numer of high-order scalars

    use clubbvars, only: &
      tndcy_precision ! Precision of CLUBB's contribution to the tendencies of mean variables

#ifdef SILHS
  use clubb_silhs_vars, only: &
    lh_sequence_length, &
    lh_microphys_type, &
    lh_microphys_disabled, &
    lh_microphys_non_interactive, &
    lh_microphys_calls
#endif

#ifdef CLUBB
#ifdef SILHS
    use clubb_silhs_vars, only: &
      lh_rt, &
      lh_t, &
      X_nl_all_levs, &
      lh_sample_point_weights, &
      X_mixt_comp_all_levs, &
      micro_field_prior, &
      lh_micro_field_sum_tndcy, &
      lh_micro_field_avg_tndcy
#endif
    use clubb_api_module, only: &
      setup_pdf_indices_api, &
      genrand_init_api, &
      genrand_intg, &
      hmp2_ip_on_hmm2_ip_ratios_type
#endif /* CLUBB */

    implicit none

    ! Constant parameters
    logical, parameter :: &
      l_uv_nudge       = .false.,  & ! Use u/v nudging (not used)
      l_implemented    = .true.      ! Implemented in a host model (always true)

    integer, parameter :: &
      grid_type    = 2, &  ! The 2 option specifies stretched thermodynamic levels
      iunit = 50           ! Fortran I/O unit

    character(len=6), parameter :: &
      saturation_equation = "flatau" ! Flatau polynomial approximation for SVP

#ifdef SILHS
    character(len=*), parameter :: &
      input_file_cloud = "/silhs_corr_matrix_cloud.in", &
      input_file_below = "/silhs_corr_matrix_below.in"
#endif
    real(kind=core_rknd), parameter :: &
      theta0   = 300._core_rknd, &! Reference temperature                     [K]
      ts_nudge = 86400._time_precision ! Time scale for u/v nudging (not used)     [s]

    ! Input Variables
    real(kind=time_precision), intent(in) :: &
      dt_clubb ! SAM-CLUBB subcycled model timestep   [s]

    real, dimension(nzm), intent(in) :: &
      z,   & ! Thermodynamic/Scalar grid in SAM         [m]
      rho    ! Thermodynamic/Scalar density in SAM      [kg/m^3]

    real, dimension(nz), intent(in) :: &
      zi, & ! Momentum/Vertical Velocity grid in SAM    [m]
      rhow  ! Momentum/Vertical Velocity density in SAM [kg/m^3]

    real, dimension(nzm), intent(in) :: &
      tv0 ! Virtual potential temperature from SAM      [K]

    real, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), intent(in) :: &
      tke   ! SGS TKE                   [m^2/s]

    ! Local Variables
    real(kind=core_rknd), dimension(nparams) :: &
      clubb_params ! These adjustable CLUBB parameters (C1, C2 ...)

    ! 1D variables with ghost points at the lowest level
    real(kind=core_rknd), dimension(nz) :: &
      zt, & ! Thermodynamic grid        [m]
      zm, & ! Momentum grid             [m]
      em    ! Turbulent kinetic energy  [-]

    logical :: l_stats  ! Stats enabled (T/F)

    real(kind=time_precision) :: &
      stats_tsamp, & ! Sampling interval for a single column of CLUBB data  [s]
      stats_tout     ! Output interval for a single column of CLUBB data    [s]

    character(len=10)  :: stats_fmt     ! Format of stats output (netCDF/GrADS)
    character(len=250) :: fname_prefix  ! Prefix for stats filename


    real(kind=core_rknd), dimension(nx) :: &
      rlon ! Longitude for stats [degrees E]

    real(kind=core_rknd), dimension(ny) :: &
      rlat ! Latitude for stats [degrees N]

    integer :: &
      err_code,  &    ! Code for when CLUBB fails
      i, j, ig, jg, & ! Loop indices
      ilen            ! Length of a string

    integer :: hydromet_dim 
    logical :: l_host_applies_sfc_fluxes ! Whether the host model applies the surface fluxes

    integer :: indx
    logical :: l_silhs_out

#ifdef SILHS
    integer(kind=genrand_intg) :: lh_seed
    type(hmp2_ip_on_hmm2_ip_ratios_type) :: hmp2_ip_on_hmm2_ip_ratios
#endif

    namelist /stats_setting/ l_stats, stats_fmt, stats_tsamp, stats_tout, &
      sample_nodes

#ifdef SILHS
    namelist /clubb_silhs/ lh_microphys_type, lh_microphys_calls, &
      lh_sequence_length, lh_seed, l_fix_chi_eta_correlations, &
      l_lh_cloud_weighted_sampling, hmp2_ip_on_hmm2_ip_ratios, Ncnp2_on_Ncnm2
#endif

!-------------------------------------------------------------------------------
!  SAM uses an Arakawa C type grid for the 3D quantities.  The UWM SCM has an
!  additional `ghost' point on the lowest pressure/thermodynamic level.
!  i.e.
!
!  SAM vert. vel. grid          UWM SCM moment. grid
!
!      Dimension    Elevation         Dimension       Elevation 
!  . . . (nz  ) . . zi(nz  )    . . . (gr%nz  ) . . gr%zm(gr%nz  ) . . .
!  . . . (nz-1) . . zi(nz-1)    . . . (gr%nz-1) . . gr%zm(gr%nz-1) . . .
!           |          |                   |             |
!  . . . (1   ) . . zi(1   )    . . . (1        ) . . gr%zm(1      )   . . .
!
!  In SAM the lowest grid point on the vertical velocity levels (or `interface' 
!  levels) is always 0 meters.  The UWM SCM supports an arbitrary starting
!  point for the momentum grid, but this code assumes 0 meters.
!
!  SAM pressure grid            UWM SCM thermo. grid
!
!      Dimension    Elevation         Dimension       Elevation 
!  . . . (nz-1) . . z(nz-1)     . . . (gr%nz  ) . . gr%zt(gr%nz  ) . . .
!  . . . (nz-2) . . z(nz-2)     . . . (gr%nz-1) . . gr%zt(gr%nz-1) . . .
!           |          |                   |                |
!  . . . (1   ) . . z(1   )     . . . (2        ) . . gr%zt(2        ) . . .
!  / / /  N/A   / /   N/A       / / / (1        ) / / gr%zt(1        ) / / /
!
!  Note that the lowest SCM point is below ground.
!-------------------------------------------------------------------------------

    !----- Begin Code -----

    ! Set the ghost point to be the distance between the first interface level,
    ! which is always zero, and the first pressure level.
    zt(1)    = real( -z(1), kind=core_rknd ) ! [m]
    ! All other pressure levels are determined by the host model
    zt(2:nz) = real( z(1:nzm), kind=core_rknd ) ! [m]

    zm = real( zi, kind=core_rknd )

    ! Set the SCM parameters (C2, etc. ) based on default values
    !call read_parameters_api( -99, "", clubb_params )

    ! Set the SCM parameters (C2, etc. ) based on a namelist
    call read_parameters_api( iunit, "CLUBB_PARAMETERS/tunable_parameters.in", clubb_params )

    ! Set the debug level.  Level 2 has additional computational expense since
    ! it checks the array variables in CLUBB for invalid values.
    call set_clubb_debug_level_api( 2 )


    ! Sanity check
    if ( sclr_dim > 0 .and. edsclr_dim > 0 ) then
      write(fstderr,*) "Only one scalar scheme can be enabled at one time"
      call task_abort()
    end if

    ! This is the tolerance on total water in the CLUBB SCM
    ! Other tracers will need this value set according to their order of 
    ! magnitude and the units they are in. Keep in mind that the variable
    ! sclrp2 will be clipped to a minimum value of sclr_tol^2
    sclr_tol(1:sclr_dim) = 1.e-8_core_rknd ! total water is in kg/kg

    ! Determine whether clubb is applying the surface flux or the host model 
    ! from the namelist variable doclubb_sfc_fluxes
    l_host_applies_sfc_fluxes = .not. doclubb_sfc_fluxes

    ! Determine total number of sample variates other than t, rt, and w.
    indx = 0
    do i = 1, nmicro_fields
      select case ( trim( mkname(i) ) )
      case ( 'QR', 'QP' )
        indx = indx + 1
        iirrm = indx

      case ( 'QI' )
        indx = indx + 1
        iirim = indx

      case ( 'QS' )
        indx = indx + 1
        iirsm = indx

      case ( 'QG' )
        ! This is not currently sampled, but we need the index to copy the
        ! mean from saved microphysics field
        indx = indx + 1
        iirgm = indx

      case ( 'CONP', 'NR' )
        indx = indx + 1
        iiNrm = indx

      case ( 'NI' )
        indx = indx + 1
        iiNim = indx

      case ( 'NS' )
        indx = indx + 1
        iiNsm = indx

      case ( 'NG' )
        indx = indx + 1
        ! See note above for QG.
        iiNgm = indx

      end select
    end do ! 1..nmicro_fields

    hydromet_dim = indx

    call setup_clubb_core_api     &
         ( nz, theta0, ts_nudge, & ! In
           hydromet_dim,  sclr_dim, &  !  In
           sclr_tol, edsclr_dim, clubb_params, & ! In
           l_host_applies_sfc_fluxes, & ! In
           l_uv_nudge, saturation_equation,  & ! In
           l_implemented, grid_type, zm(2), zm(1), zm(nz), & ! In
           zm(1:nz), zt(1:nz), & ! In
           zm(1), & ! In
           err_code )

    if ( err_code /= CLUBB_no_error ) then
      write(fstderr,*) "Initialization of CLUBB failed"
      call task_abort()
    end if

    ! Initialize stats_setting
    l_stats     = .false.
    stats_fmt   = "grads"
    stats_tsamp = 60._time_precision
    stats_tout  = 60._time_precision
    sample_nodes(:) = -1 ! Which nodes are outputting CLUBB stats columns

    ! Figure out which node and points we're sampling
    open(unit=iunit, file="clubb_stats")
    read(unit=iunit, nml=stats_setting)
    close(unit=iunit)

    if ( is_a_sample_node( rank ) .and. l_stats ) then

      ! Figure out the position on the global grid for the purposes of computing lon/lat
      call task_rank_to_index( rank, ig, jg )

      ! These formulas come from from setgrid.f90; note that latitude0 and longitude0
      ! must be set to the correct values in the prm file to get the correct position
      ! on the globe.
      do j=1,ny
        rlat(j) = real( latitude0, kind=core_rknd )+real( dy, kind=core_rknd )* &
                  real( j+jg-(ny_gl+YES3D-1)/2-1, kind=core_rknd )* & 
                    2.5e-8_core_rknd*360._core_rknd
      end do
      do i=1,nx
        rlon(i) = real( longitude0, kind=core_rknd )+real( dx, kind=core_rknd )/ &
                  cos( real( latitude0, kind=core_rknd ) *pi/180._core_rknd )* &
                    real(i+ig-nx_gl/2-1, core_rknd)*2.5e-8_core_rknd*360._core_rknd
      end do

      ! Use the case and caseid variables to set the filename
      fname_prefix = trim( case )//"_"//trim( caseid )
      ilen = len( trim( fname_prefix ) )
      fname_prefix = trim( fname_prefix )//"_node_0000"
      write(unit=fname_prefix(ilen+7:ilen+10),fmt='(i4.4)') rank

      l_silhs_out = .false.
      ! Use a bogus date, since SAM does not track the year, and it would require
      ! some work to convert the `day' variable to MMDD format
      call stats_init_api( iunit, fname_prefix, "./OUT_STAT/", l_stats, &
                       stats_fmt, stats_tsamp, stats_tout, "clubb_stats", &
                       nz, nx, ny, zt, zm, nz, zt, nz, zm, 1, 4, 1900, &
                       rlon(:), rlat(:), time, dt_clubb, l_silhs_out )

      ! If CLUBB stats are on for this node, toggle a flag in this module
      write(fstdout,*) "CLUBB stats enabled"
      lstats_clubb = .true.
    else
      lstats_clubb = .false.

    end if ! is_a_sample_node and l_stats

#ifdef SILHS
    ! Default values for namelist parameters
    lh_microphys_type = lh_microphys_non_interactive
    lh_microphys_calls = 2
    lh_sequence_length = 1
    lh_seed = 5489_genrand_intg
    l_fix_chi_eta_correlations = .true.
    l_lh_cloud_weighted_sampling = .true.

    ! Read the namelist from the prm file
    open(unit=iunit, file=trim( case )//"/prm")
    read(unit=iunit, nml=clubb_silhs, iostat=i)
    close(unit=iunit)

    if ( lh_microphys_type /= lh_microphys_disabled ) then

      ! These are needed for the PDF code -dschanen 17 June 2014
      allocate( hydromet_list(hydromet_dim) )
      allocate( hydromet_tol(hydromet_dim) )
      allocate( l_mix_rat_hm(hydromet_dim) ) 
      allocate( l_frozen_hm(hydromet_dim) ) 
      allocate( hmp2_ip_on_hmm2_ip(hydromet_dim) )
       
      if ( iirrm > 0 ) then
        hydromet_list(iirrm) = "rrm"
        hydromet_tol(iirrm)  = rr_tol
        l_mix_rat_hm(iirrm)  = .true.
        l_frozen_hm(iirrm)   = .false.
        hmp2_ip_on_hmm2_ip(iirrm) = hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip
      end if

      if ( iirim > 0 ) then
        hydromet_list(iirim) = "rim"
        hydromet_tol(iirim)  = ri_tol
        l_mix_rat_hm(iirim)  = .true.
        l_frozen_hm(iirim)   = .true.
        hmp2_ip_on_hmm2_ip(iirim) = hmp2_ip_on_hmm2_ip_ratios%rip2_ip_on_rim2_ip
      end if

      if ( iirsm > 0 ) then
        hydromet_list(iirsm) = "rsm"
        hydromet_tol(iirsm)  = rs_tol
        l_mix_rat_hm(iirsm)  = .true.
        l_frozen_hm(iirsm)   = .true.
        hmp2_ip_on_hmm2_ip(iirsm) = hmp2_ip_on_hmm2_ip_ratios%rsp2_ip_on_rsm2_ip
      end if

      if ( iirgm > 0 ) then
        hydromet_list(iirgm) = "rgm"
        hydromet_tol(iirgm)  = rg_tol
        l_mix_rat_hm(iirgm)  = .true.
        l_frozen_hm(iirgm)   = .true.
        hmp2_ip_on_hmm2_ip(iirgm) = hmp2_ip_on_hmm2_ip_ratios%rgp2_ip_on_rgm2_ip
      end if

      if ( iiNrm > 0 ) then
        hydromet_list(iiNrm) = "Nrm"
        hydromet_tol(iiNrm)  = Nr_tol
        l_mix_rat_hm(iiNrm)  = .false.
        l_frozen_hm(iiNrm)   = .false.
        hmp2_ip_on_hmm2_ip(iiNrm) = hmp2_ip_on_hmm2_ip_ratios%Nrp2_ip_on_Nrm2_ip
      end if

      if ( iiNim > 0 ) then
        hydromet_list(iiNim) = "Nim"
        hydromet_tol(iiNim)  = Ni_tol
        l_mix_rat_hm(iiNim)  = .false.
        l_frozen_hm(iiNim)   = .true.
        hmp2_ip_on_hmm2_ip(iiNim) = hmp2_ip_on_hmm2_ip_ratios%Nip2_ip_on_Nim2_ip
      end if

      if ( iiNsm > 0 ) then
        hydromet_list(iiNsm) = "Nsm"
        hydromet_tol(iiNsm)  = Ns_tol
        l_mix_rat_hm(iiNsm)  = .false.
        l_frozen_hm(iiNsm)   = .true.
        hmp2_ip_on_hmm2_ip(iiNsm) = hmp2_ip_on_hmm2_ip_ratios%Nsp2_ip_on_Nsm2_ip
      end if
      
      if ( iiNgm > 0 ) then
        hydromet_list(iiNgm) = "Ngm"
        hydromet_tol(iiNgm)  = Ng_tol
        l_mix_rat_hm(iiNgm)  = .false.
        l_frozen_hm(iiNgm)   = .true.
        hmp2_ip_on_hmm2_ip(iiNgm) = hmp2_ip_on_hmm2_ip_ratios%Ngp2_ip_on_Ngm2_ip
      end if
      
      call setup_pdf_indices_api( hydromet_dim, iirrm, iiNrm, & ! In
                                  iirim, iiNim, iirsm, iiNsm, & ! In
                                  iirgm, iiNgm ) ! In

      ! Determine d_variables and other LH indices by reading in the correlation
      ! files and from indexes determined above
      call setup_corr_varnce_array_api &
           ( trim( case )//input_file_cloud, trim( case )//input_file_below, & ! In
             iunit ) ! In

      ! Allocate based on lh_microphys_calls and d_variables
      allocate( lh_rt(nx,ny,nz,lh_microphys_calls), lh_t(nx,ny,nzm,lh_microphys_calls), &
                X_nl_all_levs(nx,ny,nz,lh_microphys_calls,d_variables), &
                X_mixt_comp_all_levs(nx,ny,nz,lh_microphys_calls), &
                lh_sample_point_weights(nx,ny,lh_microphys_calls), &
                micro_field_prior(nx,ny,nzm,nmicro_fields), &
                lh_micro_field_sum_tndcy(nx,ny,nzm,nmicro_fields), &
                lh_micro_field_avg_tndcy(nx,ny,nzm,nmicro_fields) )

    end if ! lh_microphys_type /= disabled

#endif /* SILHS */
    ! If this is restart run, just return at this point and do not re-initialize
    !  any variables as we would a run starting from the beginning.

    if ( nrestart /= 0 ) return

#ifdef SILHS
    ! Here we make the seed a function of processor number, so we don't create
    ! an artificial correlation between the samples within the subdomains.
    ! Note: this does not allow reproducibility when changing the number of processors
    !       after a restart.
    ! -dschanen 10 Jul 2013
    lh_seed = lh_seed + int( rank, kind = genrand_intg )
    call genrand_init_api( put=lh_seed )
#endif

    if ( sclr_dim > 0 ) then
      sclrp2    = 0._core_rknd
      sclrprtp  = 0._core_rknd
      sclrpthlp = 0._core_rknd
      wpsclrp   = 0._core_rknd
    end if

    ! Initialize CLUBB's tendencies to 0
    t_tndcy  = 0._tndcy_precision
    qc_tndcy = 0._tndcy_precision
    qv_tndcy = 0._tndcy_precision
    u_tndcy  = 0._tndcy_precision
    v_tndcy  = 0._tndcy_precision

    if ( ntracers > 0 ) then
      tracer_tndcy = 0._tndcy_precision
    end if

    ! SAM's dynamical core is anelastic, so l_use_boussineq should probably be
    ! set to false generally, as it is by default in the CLUBB SCM.
    if ( l_use_boussinesq ) then
      rho_ds_zm(:) = 1._core_rknd
      rho_ds_zt(:) = 1._core_rknd
      ! Set the value of dry, base-state theta_v.
      thv_ds_zm(:) = theta0
      thv_ds_zt(:) = theta0
    else 
      ! Set variables for the use of the anelastic equation set in CLUBB.
      ! Set the value of dry, static, base-state density.
      rho_ds_zm(:) = real( rhow(:), kind=core_rknd )
      rho_ds_zt(2:nz) = real( rho(1:nzm), kind=core_rknd )
      rho_ds_zt(1) = LIN_EXT( rho_ds_zt(3), rho_ds_zt(2), gr%zt(3), gr%zt(2), gr%zt(1) )
      ! Set the value of dry, base-state theta_v.
      thv_ds_zt(2:nz) = real( tv0(1:nzm), kind=core_rknd )
      thv_ds_zt(1) = real( tv0(1), kind=core_rknd )
      thv_ds_zm(:) = zt2zm_api( thv_ds_zt )
    end if
    ! Set the value of inverse dry, static, base-state density based on the
    ! value of dry, static, base-state density.
    invrs_rho_ds_zm(:) = 1.0_core_rknd / rho_ds_zm(:)
    invrs_rho_ds_zt(:) = 1.0_core_rknd / rho_ds_zt(:)

    ! Determine the initial value of some variables as in WRF-CLUBB

    wprtp(:,:,:)       = 0._core_rknd ! w'rt'
    wpthlp(:,:,:)      = 0._core_rknd ! w'thl'
    wprcp(:,:,:)       = 0._core_rknd ! w'rc'
    wp3(:,:,:)         = 0._core_rknd ! w'^3
    wp2(:,:,:)         = w_tol_sqd    ! w'^2
    up2(:,:,:)         = w_tol_sqd    ! u'^2
    vp2(:,:,:)         = w_tol_sqd    ! v'^2
    rtp2(:,:,:)        = rt_tol**2    ! rt'^2
    thlp2(:,:,:)       = thl_tol**2   ! thl'^2
    rtpthlp(:,:,:)     = 0._core_rknd ! rt'thl'
    upwp(:,:,:)        = 0._core_rknd ! u'w'
    vpwp(:,:,:)        = 0._core_rknd ! v'w'

    do i=1, nx, 1
      do j=1, ny, 1

        ! Extrapolate intial SGS TKE and use it to compute wp2
        ! This value is going to depend on initial noise and whether 
        ! Smagorinksy diffusion is enabled
        em(2:nz) = real( tke(i,j,1:nzm), kind=core_rknd )
        em(1)    = LIN_EXT( em(3), em(2), gr%zt(3), gr%zt(2), gr%zt(1) )
        em(1:nz) = max( zt2zm_api( em(1:nz) ), em_min )

!       em(:) = 1.0 ! Use this value for comparing DYCOMS II RF02 to the CLUBB SCM.

        !!!! Initialize w'^2 based on initial SGS TKE !!!!

        if ( l_tke_aniso ) then

          ! SGS TKE:  em = (1/2) * ( w'^2 + u'^2 + v'^2 )
          ! Evenly divide SGS TKE into its component
          ! contributions (w'^2, u'^2, and v'^2).

          wp2(i,j,1:nz) = (2._core_rknd/3._core_rknd) * em(1:nz)
          up2(i,j,1:nz) = (2._core_rknd/3._core_rknd) * em(1:nz)
          vp2(i,j,1:nz) = (2._core_rknd/3._core_rknd) * em(1:nz)

        else

          ! Assume isotropy for initialization of wp2
          ! SGS TKE:  em = (3/2) * w'^2

          wp2(i,j,1:nz) = (2._core_rknd/3._core_rknd) * em(1:nz)

        end if

      end do ! j=1..ny
    end do ! i=1..nx

    return
  end subroutine clubb_sgs_setup

!-------------------------------------------------------------------------------
  subroutine advance_clubb_sgs( dt_clubb, itime, stats_nsamp, stats_nout, &
                                rho, rhow, wsub, u, v, w, qpl, qci, qpi, &
                                t, qv, qcl, &
                                thlp2_forcing ) ! In/Out

! Description:
!   Advance Cloud Layers Unified By Binormals one timestep.

! References:
!   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!     Method and Model Description'' Golaz, et al. (2002)
!   JAS, Vol. 59, pp. 3540--3551.
!-------------------------------------------------------------------------------

    ! From SAM
    use grid, only: &
      nx, ny, nxp1, nyp1, nz, nzm,&! Local grid dimensions
      nx_gl, ny_gl,   &! Global grid dimensions
      dimx1_s, dimx2_s, dimy1_s, dimy2_s,& ! Scalars dimensions
      dimx1_u, dimx2_u, dimy1_u, dimy2_u,& ! U wind dimensions
      dimx1_v, dimx2_v, dimy1_v, dimy2_v,& ! V wind dimensions
      dimx1_w, dimx2_w, dimy1_w, dimy2_w,& ! W wind dimensions
      YES3D, rank, pres, dompi, dx, dy,  & 
      ntracers 

    use params, only: cp, lfus, lsub, &
      ug, vg ! ug and vg are scalars, not arrays

    use sgs_params, only: doclubb ! Variable(s)

    use vars, only: &
      fcory, fluxbt, fluxbq, fluxbu, fluxbv, gamaz, prespot ! Variables

    use microphysics, only: nmicro_fields

    use clubbvars, only: &
      upwp,        &! u'w'.                 [m^2/s^2]
      vpwp,        &! u'w'.                 [m^2/s^2]
      up2,         &! u'^2                  [m^2/s^2]
      vp2,         &! v'^2                  [m^2/s^2]
      wprtp,       &! w' r_t'.              [(m kg)/(s kg)]
      wpthlp,      &! w' th_l'.             [(m K)/s]
      wprcp,       &! w' r_c'.              [(kg/kg) m/s]
      wp2,         &! w'^2.                 [m^2/s^2]
      rtp2,        &! r_t'^2.               [(kg/kg)^2]
      thlp2,       &! th_l'^2.              [K^2]
      rtpthlp,     &! r_t' th_l'.           [(kg K)/kg]
      rcm,         &! Cloud water           [kg/kg]
      cloud_frac,  &! Cloud Fraction.       [-]
      rcm_in_layer,&! rcm in cloud layer    [kg/kg]
      cloud_cover, &! Cloud Cover           [-]
      wp3,         &! w'^3.                 [m^3/s^3]
      um,          &! x-wind                [m/s]
      vm            ! y-wind                [m/s]

    use clubbvars, only: &
      sclrp2,      & ! Passive scalar variance.       [{units vary}^2]
      sclrpthlp,   & ! Passive scalar covariance.     [{units vary}^2]
      sclrprtp,    & ! Passive scalar covariance.     [{units vary}^2]
      wpsclrp        ! w'sclr'                        [units vary m/s]

    use clubbvars, only: &
      u_tndcy,&  ! CLUBB contribution to the x wind
      v_tndcy,&  ! CLUBB contribution to the y wind
      qv_tndcy,& ! CLUBB contribution to vapor water mixing ratio
      qc_tndcy,& ! CLUBB contribution to liquid water mixing ratio
      t_tndcy    ! CLUBB contribution to moist static energy

    use clubbvars, only: &
      tracer_tndcy ! CLUBB contribution to a set of tracers

    use clubbvars, only: &
      rho_ds_zm,       & ! Dry, static density on momentum levels      [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermodynamic levels [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density on momentum levels [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density on thermo. levels  [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levels  [K]
      thv_ds_zt          ! Dry, base-state theta_v on thermo. levels   [K]

    use clubbvars, only: &
      sclr_dim, & ! Constant(s)
      edsclr_dim

    use clubbvars, only: &
      tndcy_precision ! Constant(s)

    use tracers, only: &
      fluxbtr, & ! Variable(s)
      tracer

    use calc_vars_util, only: &
        t2thetal, & ! Procedure(s)
        thetal2t

    ! From CLUBB
    use clubb_api_module, only: &
      clubb_no_error, &    ! Constant
      clubb_at_least_debug_level_api, & ! Function

      zm2zt_api, zt2zm_api, & ! Functions
      gr, & ! Derived type

      l_stats, l_stats_samp, & ! Logicals

      stats_begin_timestep_api, stats_end_timestep_api, & ! Subroutines

      pdf_parameter, & ! Derived type

      fstderr, & ! Constant

      clubb_i, clubb_j, & ! Variable(s)

      hydromet_dim, & ! Variable(s)

      iirrm, iiNrm, iirsm, iirim, & ! Variable(s)
      iiNsm, iiNim, iiNgm, iirgm

#ifdef SILHS
    use microphysics, only: &
      conc, micro_field, nmicro_fields, mkname ! Variable(s)
#else
    use microphysics, only: &
      micro_field, nmicro_fields, mkname ! Variable(s)
#endif /* SILHS */

#ifdef SILHS

    use clubb_silhs_vars, only: &
      lh_microphys_calls, &
      lh_microphys_type, &
      lh_microphys_disabled, &
      lh_sequence_length

    use clubb_api_module, only: &

      ! Note that as a 1D variable this will vary from column to column and will
      ! be overwritten for each i,y corresponding to SAM's x and y coordinates
      Lscale, & ! Variable(s)

      corr_array_n_cloud, & ! Variable(s)
      corr_array_n_below, &
      d_variables, &

      iiPDF_chi, iiPDF_w, &
      iiPDF_rr, iiPDF_rs, iiPDF_ri, &
      iiPDF_Nr, iiPDF_Ns, iiPDF_Ni, iiPDF_Ncn, &

      stats_accumulate_hydromet_api, &
      
      hydromet_pdf_parameter,   & ! Type(s)

      setup_pdf_parameters_api, &    ! Procedure(s)

      dp ! Constant

    use silhs_api_module, only: &
      lh_subcolumn_generator_api, & ! Procedure
      clip_transform_silhs_output_api, & ! Procedure
      lh_clipped_variables_type ! Type

    use clubb_silhs_vars, only: &
      lh_rt, & ! Variable(s)
      lh_t, &
      X_nl_all_levs, &
      lh_sample_point_weights, &
      X_mixt_comp_all_levs

#endif /* SILHS */

    implicit none

    ! Parameters
    logical, parameter :: & 
      l_implemented = .true., & ! CLUBB is implemented in a host model, so this is true
      l_advect      = .false. ! Whether to advect around the high-order moments

    real(kind=core_rknd), parameter, dimension(nz) :: &
      zero = 0.0_core_rknd ! Field of zeros

    ! Input
    real(kind=time_precision), intent(in) :: &
      dt_clubb ! Timestep size for CLUBB [s]

    integer, intent(in) :: &
      itime, &       ! Current timestep
      stats_nsamp, & ! Sampling interval for a single column of CLUBB data [timestep]   
      stats_nout     ! Output interval for a single column of CLUBB data [timestep]

    real, intent(in), dimension(nzm) :: &
      rho  ! Air density        [kg/m^3]

    real, intent(in), dimension(nz) :: &
      wsub,&! Imposed vertical velocity       [m/s]
      rhow  ! Density on vert velocity grid   [kg/m^3]

    real, intent(in), dimension(dimx1_u:dimx2_u,dimy1_u:dimy2_u,nzm) :: &
      u  ! u wind    [m/s]

    real, intent(in), dimension(dimx1_v:dimx2_v,dimy1_v:dimy2_v,nzm) :: &
      v  ! v wind    [m/s]

    real, intent(in), dimension(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) :: &
      w ! Vertical wind  [m/s]

    real, intent(in), dimension(nx,ny,nzm) :: &
      qpl,& ! Liquid water mixing ratio (precipitation) [kg/kg]
      qci,& ! Cloud ice water mixing ratio              [kg/kg]
      qpi   ! Snow + graupel mixing ratio (precip)      [kg/kg]

    real, intent(in), dimension(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm) :: &
      t     ! Moist static energy           [K]

    real, intent(in), dimension(nx,ny,nzm) :: &
      qv, & ! Water vapor mixing ratio                  [kg/kg]
      qcl   ! Liquid water mixing ratio (condensate)    [kg/kg]

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(nzm) :: &
      thlp2_forcing   ! <th_l'^2> forcing (momentum levels)   [K^2/s]

    ! Local Variables
    real(kind=core_rknd) :: &
      wpthlp_sfc, &! w' theta_l' at surface    [(m K)/s]
      wprtp_sfc,  &! w' r_t' at surface        [(kg m)/( kg s)]
      upwp_sfc,   &! u'w' at surface           [m^2/s^2]
      vpwp_sfc     ! v'w' at surface           [m^2/s^2]

    real(kind=core_rknd), dimension(nz) :: &
      thlm,    &! Liquid water potential temperature (theta_l)  [K]
      rtm,     &! Total water mixing ratio                      [kg/kg] 
      p_in_Pa, &! Pressure                                      [Pa] 
      rho_zt,  &! Density on pressure levels                    [kg/m^3]
      rho_zm,  &! Density on momentum levels                    [kg/m^3]
      exner,   &! Exner function                                [-]
      wm_zm,   &! Imposed subs. + perturbation w on vertical vel. levels    [m/s]
      wm_zt,   &! Imposed subs. + perturbation w on pressure levels         [m/s]
      rfrzm     ! Total ice-phase water mixing ratios           [kg/kg]

    real( kind = core_rknd ), dimension(nz, hydromet_dim) :: &
      wp2hmp,  &
      rtphmp,  &
      thlphmp

    real, dimension(nz) :: &
      dum     ! Dummy array for advection

    real(kind=core_rknd), allocatable, dimension(:,:) :: &
      sclrm,          & ! Array for high order passive scalars
      sclrm_forcing,  & ! Large-scale forcing array for passive scalars
      edsclrm,        & ! Array for eddy passive scalars
      edsclrm_forcing   ! Large-scale forcing array for eddy passive scalars

    real(kind=core_rknd), allocatable, dimension(:) :: &
      wpedsclrp_sfc, & ! Array for passive scalar surface flux
      wpsclrp_sfc      ! Array for high order scalar surface flux

    ! Thermo grid versions of variables on the momentum grid
    real, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nz) :: &
      wp2_zt, rtp2_zt, thlp2_zt, rtpthlp_zt, &
      wprtp_zt, wpthlp_zt, up2_zt, vp2_zt, &
      um_r4, vm_r4, um_old, vm_old ! wind arrays

    real(kind=tndcy_precision), dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nz) :: &
      um_change, vm_change ! Change in u/v      [m/s^2]

    type(pdf_parameter), allocatable, dimension(:) :: &
      pdf_params ! PDF parameters       [units vary]

    real(kind=core_rknd), dimension(nz,hydromet_dim) :: &
      hydromet,    & ! Mean value of hydrometeor fields            [units vary]
      wphydrometp    ! Covariances between w and the hydrometeors  [units vary]

#ifdef SILHS
    real(kind=core_rknd), dimension(nz) :: &
      Ncnm    ! Mean simplified cloud nuclei conc.; Nc=Ncn*H(chi) [num/kg]

    real( kind = core_rknd ), dimension(nz, hydromet_dim) :: &
      hydrometp2  ! Overall variance of the hydrometeors   [units vary]

    real(kind=core_rknd), dimension(nz) :: &
      Nc_in_cloud        ! Mean (in-cloud) cloud droplet conc.  [num/kg]

    real(kind=core_rknd), dimension(nz,lh_microphys_calls) :: &
      lh_thl ! Liquid potential temperature subcolumns  [K]

    type(hydromet_pdf_parameter), dimension(nz) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters      [units vary]

    real( kind = core_rknd ), dimension(d_variables,d_variables,nz) :: &
      corr_array_1, &        ! Correlation matrix for the first pdf component    [-]
      corr_array_2, &        ! Correlation matrix for the second pdf component   [-]
      corr_cholesky_mtx_1, & ! Transposed correlation cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed correlation cholesky matrix, 2nd comp. [-]

    real( kind = core_rknd ), dimension(d_variables,nz) :: &
      mu_x_1,    & ! Mean array for the 1st PDF component                 [units vary]
      mu_x_2,    & ! Mean array for the 2nd PDF component                 [units vary]
      sigma_x_1, & ! Standard deviation array for the 1st PDF component   [units vary]
      sigma_x_2    ! Standard deviation array for the 2nd PDF component   [units vary]

    type(lh_clipped_variables_type), dimension(nz,lh_microphys_calls) :: &
      lh_clipped_vars

    ! Whether to call Ncn_to_Nc or not. This call has caused trouble with the MG 1.0 microphysics,
    ! since the Nc-tendency is not updated. See cpt:comment:177:ticket:198 and cpt:ticket:201.
    logical :: l_use_Ncn_to_Nc = .true. 
#endif /* SILHS */

    real(kind=core_rknd), dimension(nz) :: &
      ice_supersat_frac, &
      radf

    ! Horizontal grid spacings (i.e., dx and dy), used for computing Lscale_max
    real(kind=core_rknd) :: host_dx, host_dy ! [m]

    integer :: err_code

    ! Array indices
    integer :: i, j, k, ig, jg, ip1, jp1, jm1, indx

!-------------------------------------------------------------------------------

    !----- Begin Code -----

    call t_startf( 'advance_clubb' ) ! For timing

    ! Initialize err_code to CLUBB_no_error.  In the event of the singular 
    ! matrix, etc. the variable will be set to the appropriate error code 
    ! within advance_clubb_core_api
    err_code = CLUBB_no_error

    host_dx = real( dx, kind=core_rknd )
    host_dy = real( dy, kind=core_rknd )

    ! Feed nothing into radf (set it to zero)
    radf(1:nz) = 0.0_core_rknd

    ! Density is in correct units
    rho_zt(2:nz) = real( rho(1:nzm), kind=core_rknd )
    rho_zt(1) = LIN_EXT( rho_zt(3), rho_zt(2), gr%zt(3), gr%zt(2), gr%zt(1) )

    rho_zm(1:nz) = real( rhow(1:nz), kind=core_rknd )

    ! Compute and extrapolate Exner function
    exner(2:nz) = 1.0_core_rknd / real( prespot(1:nzm), kind=core_rknd )
    exner(1)    = 1.0_core_rknd / LIN_EXT( exner(3), exner(2), gr%zt(3), gr%zt(2), gr%zt(1) )

    ! Allocate passive scalar arrays
    allocate( wpsclrp_sfc(sclr_dim), sclrm(nz,sclr_dim), &
              sclrm_forcing(nz,sclr_dim) )
    allocate( wpedsclrp_sfc(edsclr_dim), edsclrm(nz,edsclr_dim), &
              edsclrm_forcing(nz,edsclr_dim) )

    ! Allocate variables for the PDF closure scheme
    allocate( pdf_params(1:nz) )
 
    do i = 1, nx, 1
      do j = 1, ny, 1

        ip1 = min( nxp1, i+1 )  ! This is redundant, but we include it for safety
        jp1 = min( nyp1, j+1 )  ! This prevents an array out of bounds error
                                !   for dvdt in a 2D simulation

        ! Average u-wind (east-west wind) to scalar points.
        um_r4(i,j,2:nz) = 0.5 * ( u(i,j,1:nzm) + u(ip1,j,1:nzm) ) + ug
        um_r4(i,j,1) = um_r4(i,j,2)

        ! Average v-wind (north-south wind) to scalar points.
        vm_r4(i,j,2:nz) = 0.5 * ( v(i,j,1:nzm) + v(i,jp1,1:nzm) ) + vg
        vm_r4(i,j,1) = vm_r4(i,j,2)
      end do
    end do

    ! Adjust the ghost points to allow for interpolation back on to 
    !   the u & v grid points
    if ( dompi ) then
      call task_exchange( um_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                          nzm, 3,3,3,3, ntracers+nmicro_fields+19)
      call task_exchange( vm_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                          nzm, 3,3,3,3, ntracers+nmicro_fields+20)
    else
      call bound_exchange( um_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                           nzm, 3,3,3,3, ntracers+nmicro_fields+19)
      call bound_exchange( vm_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                           nzm, 3,3,3,3, ntracers+nmicro_fields+20)
    end if
    ! Lower Boundary condition on u/v
    um_r4(:,:,1) = um_r4(:,:,2)
    vm_r4(:,:,1) = vm_r4(:,:,2)

    ! Preserve value of u and v to calculate total change from CLUBB
    um_old = um_r4
    vm_old = vm_r4

    ! Copy the SAM precision values into CLUBB precision arrays
    um = real( um_r4, kind=core_rknd )
    vm = real( vm_r4, kind=core_rknd )

    do i=1, nx, 1

      do j=1, ny, 1
        
        ! Sample from a single column
        if ( is_a_sample_node( rank ) .and. lstats_clubb ) then
          ! Set the x and y within advance_clubb_core_api for stats purposes
          clubb_i = i
          clubb_j = j
          call stats_begin_timestep_api(itime,stats_nsamp,stats_nout )
        else
          l_stats_samp = .false.
        end if

        ! The 2-D flux arrays are already in the correct units

        ! The 2-D flux arrays are already in the correct units
        wprtp_sfc  = real( fluxbq(i,j), kind=core_rknd )  ! [m kg/kg s] 
        wpthlp_sfc = real( fluxbt(i,j), kind=core_rknd )  ! [m K/s]
! Vince Larson set sfc momentum flux constant, as a temporary band-aid.  
! 25 Feb 2008.
        ! These are set for the purposes of computing sfc_var, but this value is
        ! not applied to the value of u and v in SAM.
        upwp_sfc = real( fluxbu(i,j), kind=core_rknd )
        vpwp_sfc = real( fluxbv(i,j), kind=core_rknd )
! End of Vince Larson's change

        ! Set the surface flux of the two scalar types to the tracer flux at the
        ! bottom of the domain, and set edsclrm to the tracer
        do indx = 1, edsclr_dim, 1
          wpedsclrp_sfc(indx) = real( fluxbtr(i,j,indx), kind=core_rknd )
          edsclrm(2:nz,indx)  = real( tracer(i,j,1:nzm,indx), kind=core_rknd )
          edsclrm(1,indx) = real( LIN_EXT( edsclrm(3,indx), edsclrm(2,indx), &
                                     gr%zt(3), gr%zt(2), gr%zt(1) ), kind=core_rknd )

          edsclrm_forcing(1:nz,indx) = 0.0_core_rknd
        end do

        do indx = 1, sclr_dim, 1
          wpsclrp_sfc(indx) = real( fluxbtr(i,j,indx), kind=core_rknd )
          sclrm(2:nz,indx)  = real( tracer(i,j,1:nzm,indx), kind=core_rknd )
          sclrm(1,indx) = LIN_EXT( sclrm(3,indx), sclrm(2,indx), &
                                 gr%zt(3), gr%zt(2), gr%zt(1) )
          sclrm_forcing(1:nz,indx) = 0.0_core_rknd
        end do

  
        ! Check for negative values of water vapor being fed from SAM into CLUBB
        if ( clubb_at_least_debug_level_api( 2 ) ) then
          do k=1,nzm
            if ( qv(i,j,k) < 0. ) then
              write(fstderr,*) 'SAM has fed into CLUBB negative rv at grid point i,j,k =', &
                i, j, k
            end if
          end do         

          ! Check for negative values of cloud water being fed from SAM into CLUBB
          do k=1,nzm
            if ( qcl(i,j,k)  < 0. ) then
              write(fstderr,*) 'SAM has fed into CLUBB negative qcl at grid point i,j.k =', &
                i, j, k
            end if
          end do
        end if ! clubb_at_least_debug_level_api( 2 )

        ! Total water. Since the SCM does not account for ice, we sum only the
        ! non-precipitating liquid and vapor

        ! Total water is the sum of non-precipitating liquid + vapor
        rtm(2:nz) = real( qv(i,j,1:nzm) + qcl(i,j,1:nzm), kind=core_rknd )
        rtm(1)    = rtm(2)

        ! Cloud water is total non-precipitating liquid
        rcm(i,j,2:nz) = real( qcl(i,j,1:nzm), kind=core_rknd )
        rcm(i,j,1)    = 0.0_core_rknd ! No below ground cloud water

        ! Note: t is moist static energy, which is not quite the same as liquid
        ! potential temperature.
        thlm(2:nz) = real( t2thetal( t(i,j,1:nzm), gamaz(1:nzm), &
                                     qpl(i,j,1:nzm), qci(i,j,1:nzm), &
                                     qpi(i,j,1:nzm), prespot(1:nzm) ), &
                           kind = core_rknd )
        thlm(1)    = thlm(2)

        ! The w variable requires no extrapolation

    ! Vince Larson added option for l_advect = .true. .  13 Mar 2008.
    !   SAM's subroutine 'subsidence' imposes wsub on t, q, u, and v.
    !   SAM advects all means using u, v, w.
    !   When implemented in a host model, CLUBB imposes wm_zm/wm_zt on higher-order
    ! moments but not means. 
    !   (l_advect=.true.) advects all higher-order moments using u, v, w.
         if ( l_advect )  then
            wm_zt(1) = 0._core_rknd
            wm_zt(2:nz) = real( wsub(1:nzm), kind=core_rknd ) ! Use this if l_advect = .true.
            wm_zm = zt2zm_api( wm_zt )
         else ! l_advect = .false.
            ! Higher-order moments are advected vertically but not horizontally.
            ! In principle, this could lead to undesirable accumulation.
            wm_zt(1) = 0._core_rknd ! Set ghost point to 0.
            wm_zt(2:nz) = real( wsub(1:nzm), kind=core_rknd ) ! wsub is on the t-levels
            wm_zm(1:nz) = zt2zm_api( wm_zt ) ! Interpolate imposed subsidence to m-levels

            ! Resolved vertical velocity is on the momentum levels
            wm_zm(1:nz) = wm_zm(1:nz) + real( w(i,j,1:nz), kind=core_rknd ) 
            ! Interpolate resolved w to t-levels
            wm_zt(1:nz) = wm_zt + zm2zt_api( real( w(i,j,1:nz), kind=core_rknd ) )
         end if
    ! End Vince Larson's commenting

        ! Add in pressure perturbation, extrapolate, & convert from mb to Pa.
        ! Vince Larson of UWM removed perturbation pressure to avoid 
        !      negative pressure at domain top in ARM9707.  22 Dec 2007.
    !   pr(2:nz) = 100. *  ( pres(1:nzm) + p(i,j,1:nzm) )
    !   pr(1)    = 100. * LIN_EXT( pres(2)+p(i,j,2), pres(1)+p(i,j,1), &
    !                              gr%zt(3), gr%zt(2), gr%zt(1) )
        P_in_Pa(2:nz) = 100._core_rknd *  real( pres(1:nzm), kind=core_rknd )
        P_in_Pa(1)    = LIN_EXT( P_in_Pa(3), P_in_Pa(2), &
                                 gr%zt(3), gr%zt(2), gr%zt(1) )

        !  End Vince Larson's change.

        ! Sum all forms of ice
        rfrzm(2:nz) = real( qpi(i,j,1:nzm) + qci(i,j,1:nzm), kind=core_rknd )
        rfrzm(1) = 0._core_rknd

        wp2hmp  = 0.0_core_rknd
        rtphmp  = 0.0_core_rknd
        thlphmp = 0.0_core_rknd


        hydromet = 0._core_rknd

#ifdef SILHS
        ! Ncnm is not mean cloud droplet concentration.  This needs to be
        ! corrected at some point.
        Ncnm(2:nz) = real( conc(i,j,1:nzm), kind=core_rknd )
#endif /* SILHS */

      do indx = 1, nmicro_fields
        select case ( trim( mkname(indx) ) )
        case ( 'QR', 'QP' )
          if ( iirrm > 0 ) then
            hydromet(2:nz,iirrm) = real( micro_field(i,j,:,indx), core_rknd )
          end if

        case ( 'QI' )
          if ( iirim > 0 ) then
            hydromet(2:nz,iirim)  = real( micro_field(i,j,:,indx), core_rknd )
          end if

        case ( 'QS' )
          if ( iirsm > 0 ) then
            hydromet(2:nz,iirsm) = real( micro_field(i,j,:,indx), core_rknd )
          end if

        case ( 'QG' )
          if ( iirgm > 0 ) then
            hydromet(2:nz,iirgm) = real( micro_field(i,j,:,indx), core_rknd )
          end if

        case ( 'CONP', 'NR' )
          if ( iiNrm > 0 ) then
            hydromet(2:nz,iiNrm) = real( micro_field(i,j,:,indx), core_rknd )
          end if

        case ( 'NI' )
          if ( iiNim > 0 ) then
            hydromet(2:nz,iiNim) = real( micro_field(i,j,:,indx), core_rknd )
          end if

        case ( 'NS' )
          if ( iiNsm > 0 ) then
            hydromet(2:nz,iiNsm) = real( micro_field(i,j,:,indx), core_rknd )
          end if

        case ( 'NG' )
          ! Note: graupel is not a part of X_nl_all_levs. These lines are
          ! strictly for the purpose of outputting graupel from a single column
          if ( iirgm > 0 ) then
            hydromet(2:nz,iirgm) = real( micro_field(i,j,:,indx), core_rknd )
          end if

         end select
       end do ! 1..nmicro_fields

        ! Note:  we could set wphydrometp in WRF-CLUBB accoring to:
        !        <w'hm'> = - K_r * d<hm>/dz.
        !        We will have to base the calculation off of the values of
        !        hydromet from a single timestep, unlike CLUBB standalone that
        !        has the benefit of using a Crank-Nicholson time-stepping scheme.
        wphydrometp = 0.0_core_rknd

        ! Call the single column model, CLUBB
        call advance_clubb_core_api &
             ( l_implemented, dt_clubb, real( fcory(j), kind=core_rknd ), gr%zm(1), &
               hydromet_dim, & ! In
               zero(:), zero(:), zero(:), zero(:), & ! In
               sclrm_forcing, edsclrm_forcing, zero(:), & ! In
               zero(:), zero(:), thlp2_forcing(:), & ! In
               zero(:), wm_zm(:), wm_zt(:), & ! In
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, & ! In
               wpsclrp_sfc, wpedsclrp_sfc, &  ! In
               P_in_Pa(:), rho_zm(:), rho_zt(:), exner(:), & ! In
               rho_ds_zm(:), rho_ds_zt(:), invrs_rho_ds_zm(:), & ! In
               invrs_rho_ds_zt(:), thv_ds_zm(:), thv_ds_zt(:), hydromet, & ! In
               rfrzm(:), radf, wphydrometp, wp2hmp, rtphmp, thlphmp,  & ! In
               host_dx, host_dy, & ! In
               um(i,j,:), vm(i,j,:), upwp(i,j,:), vpwp(i,j,:), up2(i,j,:), vp2(i,j,:), & ! In/out
               thlm(:), rtm(:), wprtp(i,j,:), wpthlp(i,j,:), & ! In/out
               wp2(i,j,:), wp3(i,j,:), rtp2(i,j,:), thlp2(i,j,:), rtpthlp(i,j,:), & ! In/out
               sclrm, sclrp2(i,j,:,:), sclrprtp(i,j,:,:), sclrpthlp(i,j,:,:), & ! In/out
               wpsclrp(i,j,:,:), edsclrm, err_code, & ! In/out
               rcm(i,j,:), wprcp(i,j,:), cloud_frac(i,j,:), ice_supersat_frac, & ! Out
               rcm_in_layer(i,j,:), cloud_cover(i,j,:), pdf_params ) ! Out
#ifdef SILHS
        do k = 1, gr%nz
          if ( pdf_params(k)%mixt_frac > 1._core_rknd .or. &
               pdf_params(k)%mixt_frac < 0._core_rknd ) then

             write(fstderr,*) "Error in gaus_mixt_points:  mixture " &
                                // "fraction, mixt_frac, does not lie in [0,1]."
            
             call task_abort( )

           end if ! mixt_frac test

           !print *, k, pdf_params(k)%mixt_frac 
        end do ! 1 .. gr%nz

        if ( lh_microphys_type /= lh_microphys_disabled ) then

          Nc_in_cloud = -999.0_core_rknd ! TODO Add this as an option

          ! hydrometp2 is intent (inout); set to 0 for now and handle it later.
          hydrometp2 = 0.0_core_rknd

          ! Setup the PDF parameters.
          call setup_pdf_parameters_api( nz, d_variables, dt_clubb, &                ! Intent(in)
                                         Nc_in_cloud, rcm, cloud_frac(i,j,:), &      ! Intent(in)
                                         ice_supersat_frac, hydromet, wphydrometp, & ! Intent(in)
                                         corr_array_n_cloud, corr_array_n_below, &   ! Intent(in)
                                         pdf_params, l_stats_samp, &                 ! Intent(in)
                                         hydrometp2, &                               ! Intent(inout)
                                         mu_x_1, mu_x_2, &                           ! Intent(out)
                                         sigma_x_1, sigma_x_2, &                     ! Intent(out)
                                         corr_array_1, corr_array_2, &               ! Intent(out)
                                         corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! Intent(out)
                                         hydromet_pdf_params )                       ! Intent(out)

          ! Call the modified subcolumn generator
          call lh_subcolumn_generator_api &
               ( lh_iter, d_variables, lh_microphys_calls, lh_sequence_length, gr%nz, & ! In
                 pdf_params, gr%dzt, rcm(i,j,:), Lscale, & ! In
                 rho_ds_zt, mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, & ! In
                 real( corr_cholesky_mtx_1, kind = dp ), & ! In
                 real( corr_cholesky_mtx_2, kind = dp ), & ! In
                 hydromet_pdf_params, & ! In
                 X_nl_all_levs(i,j,:,:,:), X_mixt_comp_all_levs(i,j,:,:), & ! Out
                 lh_sample_point_weights(i,j,:)  ) ! Out

          call clip_transform_silhs_output_api &
               ( gr%nz, lh_microphys_calls, d_variables, X_mixt_comp_all_levs(i,j,:,:), & ! In
                 X_nl_all_levs(i,j,:,:,:), pdf_params, l_use_Ncn_to_Nc, &
                 lh_clipped_vars )

          lh_rt(i,j,:,:) = lh_clipped_vars(:,:)%rt
         ! Convert the thetal sample points into moist static energy sample points
          lh_t(i,j,:,:) = convert_thl_to_t_LH( lh_clipped_vars%thl, gamaz, prespot, &
                                               X_nl_all_levs(i,j,:,:,:) )

          ! Increment the iteration count for the purpose of knowing whether to repeat
          lh_iter = lh_iter + 1 

          if ( is_a_sample_node( rank ) ) then 
            call stats_accumulate_hydromet_api( hydromet, rho_ds_zt ) ! In
          end if

        end if ! lh_microphys_type /= lh_microphys_disabled
#endif /* SILHS */

        if ( is_a_sample_node( rank ) ) then
          call stats_end_timestep_api( )
        end if

        ! Check if a critical error has occured within the CLUBB model
        if ( err_code /= clubb_no_error ) then
          call task_rank_to_index( rank, ig, jg )
          write(fstderr,*) "Task #:", rank
          write(fstderr,*) "Single-column model failed at: ", "nx=", i, ";", "ny=", j, ";"
          write(fstderr,*) "x global=", i+ig, ";", "y global=", j+jg, ";"
          call task_abort( )
        end if

        ! If we're not doing a doclubbnoninter run, then we feed the results back
        !   into the 3D SAM model arrays.  Here we compute the total tendency to
        !   allow for subcycling and save compute time.
        if ( doclubb ) then

          ! Check for negative values of water vapor
          if ( clubb_at_least_debug_level_api( 2 ) ) then
            do k=1,nz
              if ( ( rtm(k) - rcm(i,j,k) ) < 0._core_rknd ) then
                write(fstderr,*) 'CLUBB has produced negative rvm at grid level k=', k
              end if
            end do
          end if ! clubb_at_least_debug_level_api( 2 )

          ! Re-compute vapor for total water and liquid from CLUBB
          !qv(i,j,1:nzm) = rtm(2:nz) - rcm(i,j,2:nz)
          qv_tndcy(i,j,1:nzm) = &
            ( rtm(2:nz) - rcm(i,j,2:nz) - real( qv(i,j,1:nzm), kind=core_rknd ) ) / dt_clubb

          if ( clubb_at_least_debug_level_api( 2 ) ) then
            ! Check for negative values of cloud water
            do k=1,nz
              if ( rcm(i,j,k)  < 0._core_rknd ) then
                write(fstderr,*) 'CLUBB has produced negative rcm at grid level k=', k
              end if
            end do
          end if ! clubb_at_least_debug_level_api( 2 )

          ! Re-compute qcl based on new rcm 
          !qcl(i,j,1:nzm) = rcm(i,j,2:nz)
          ! Compute tendency of total water due to CLUBB
          qc_tndcy(i,j,1:nzm) = ( rcm(i,j,2:nz) - real( qcl(i,j,1:nzm), kind=core_rknd ) ) &
                              / dt_clubb

          ! Compute moist static energy based on new thetal
!         t(i,j,1:nzm) = thetal2t( real( thlm(2:nz) ), gamaz(1:nzm), &
!                                  qpl(i,j,1:nzm), qci(i,j,1:nzm), &
!                                  qpi(i,j,1:nzm), prespot(1:nzm) )

          ! Compute tendency of moist static energy due to CLUBB
          ! Note that this formula assumes qci/qpl/qpi won't change rapidly in
          ! the time between successive clubb calls in order to avoid calling 
          ! thetal2t on at every SAM timestep -dschanen 27 Oct 08
          t_tndcy(i,j,1:nzm) &
          = ( real( thetal2t( real( thlm(2:nz) ), gamaz(1:nzm), &
                              qpl(i,j,1:nzm), qci(i,j,1:nzm), qpi(i,j,1:nzm), &
                              prespot(1:nzm) ), kind = core_rknd ) &
              - real( t(i,j,1:nzm), kind=core_rknd ) )  / dt_clubb

          do indx = 1, edsclr_dim 
            tracer_tndcy(i,j,1:nzm,indx) = &
              ( edsclrm(2:nz,indx) - real( tracer(i,j,1:nzm,indx), kind=core_rknd ) ) &
               / dt_clubb
          end do

          do indx = 1, sclr_dim
            tracer_tndcy(i,j,1:nzm,indx) = &
              ( sclrm(2:nz,indx) - real( tracer(i,j,1:nzm,indx), kind=core_rknd ) ) / dt_clubb
          end do

        end if ! doclubb

      end do ! j

    end do ! i

    ! De-allocate temporary arrays.  This is just in case the compiler isn't
    ! 100% Fortran 95 compliant and doesn't de-allocate this memory when it
    ! leaves the scope of advance_clubb_sgs
    deallocate( wpsclrp_sfc, sclrm )
    deallocate( wpedsclrp_sfc, edsclrm )
    deallocate( pdf_params )

    ! Copy back the value from the CLUBB precision um and vm
    um_r4 = real( um )
    vm_r4 = real( vm )

    if ( doclubb ) then

      ! Adjust the ghost points to allow for interpolation back onto the u & v grid
      if ( dompi ) then
        call task_exchange( um_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                            nzm, 3,3,3,3, ntracers+nmicro_fields+19)
        call task_exchange( vm_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                            nzm, 3,3,3,3, ntracers+nmicro_fields+20)
      else
        call bound_exchange( um_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                             nzm, 3,3,3,3, ntracers+nmicro_fields+19)
        call bound_exchange( vm_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                             nzm, 3,3,3,3, ntracers+nmicro_fields+20)
      end if

      ! Compute the total change in u due to the CLUBB part of the code
      um_change = real( um_r4 - um_old, kind=tndcy_precision ) / dt_clubb
      vm_change = real( vm_r4 - vm_old, kind=tndcy_precision ) / dt_clubb

      ! Average the contributions of CLUBB to the wind back on to the u and v grid
      ! This has shown to make the model unstable at fine horizontal resolution.
      ! To interpolate across subdomain boundaries requires that we 
      !  transfer information using MPI (via task_exchange).
      do i=1, nx, 1
        do j=1, ny, 1
          jm1 = max( dimy1_s, j-1 ) ! For the 2D case vm wind

          ! The horiztontal grid in SAM is always evenly spaced, so we just use
          ! 0.5 *( x(n-1)+x(n) ) to interpolate back to the u,v point on the Arakawa C grid
          u_tndcy(i,j,1:nzm) = &
            0.4_tndcy_precision * & ! This is a made up coefficient to reduce numerical instability
            0.5_tndcy_precision * &
            real( um_change(i,j,2:nz) + um_change(i-1,j,2:nz), kind=tndcy_precision )
          v_tndcy(i,j,1:nzm) = &
            0.4_tndcy_precision * & ! This is a made up coefficient to reduce numerical instability
            0.5_tndcy_precision * &
            real( vm_change(i,j,2:nz) + vm_change(i,jm1,2:nz), kind=tndcy_precision )

        end do ! j
      
      end do ! i

    end if ! doclubb
  

! Vince Larson attempted to advect higher-order moments horizontally.  
!     26 Feb 2008.

! Horizontal advection of higher-order moments.

! The following method has the drawback of requiring two interpolations, 
!    which unnecesarily smooths the fields in the vertical.
! In preparation for advection, interpolate to thermodynamic (scalar) vertical gridpoints.
!   (wp3 is already on the thermodynamic gridpoints.)


!print*, 'Before advection, wp2(nx,ny,:) =', wp2(nx,ny,:)
! For now we default to not doing this, because the interpolation seems to cause
! and artificial rise in fields such as moisture at a coarse model resolution.
! -dschanen 29 Apr 2008
  if ( l_advect ) then

    do i=1, nx, 1
      do j=1, ny, 1

        wp2_zt(i,j,:)     = real( zm2zt_api( wp2(i,j,:) ) )
        up2_zt(i,j,:)     = real( zm2zt_api( up2(i,j,:) ) )
        vp2_zt(i,j,:)     = real( zm2zt_api( vp2(i,j,:) ) )
        rtp2_zt(i,j,:)    = real( zm2zt_api( rtp2(i,j,:) ) )
        thlp2_zt(i,j,:)   = real( zm2zt_api( thlp2(i,j,:) ) )
        rtpthlp_zt(i,j,:) = real( zm2zt_api( rtpthlp(i,j,:) ) )
        wprtp_zt(i,j,:)   = real( zm2zt_api( wprtp(i,j,:) ) )
        wpthlp_zt(i,j,:)  = real( zm2zt_api( wpthlp(i,j,:) ) )

      end do ! j
    end do   ! i

    if ( dompi ) then
  
      call task_exchange( wp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+10 )
      call task_exchange( rtp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+11 )
      call task_exchange( thlp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+12 )
      call task_exchange( rtpthlp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+13 )
      call task_exchange( wprtp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+14 )
      call task_exchange( wpthlp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+15 )
      call task_exchange( wp3(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+16 )
      call task_exchange( up2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+17 )
      call task_exchange( vp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+18 )
    else

      call bound_exchange( wp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+10 )
      call bound_exchange( rtp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+11 )
      call bound_exchange( thlp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+12 )
      call bound_exchange( rtpthlp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+13 )
      call bound_exchange( wprtp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+14 )
      call bound_exchange( wpthlp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+15 )
      call bound_exchange( wp3(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+16 )
      call bound_exchange( up2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+17 )
      call bound_exchange( vp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+18 )

    end if

    ! Now call the standard SAM advection subroutine for scalars
    call advect_scalar( wp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( wp3(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),  &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( rtp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( thlp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( rtpthlp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( wprtp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( wpthlp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( up2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( vp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

!print*, 'After advect, wp2_zt(dimx2_s,dimy2_s,:) =', wp2_zt(dimx2_s,dimy2_s,:)
!
!do i=dimx1_s, dimx2_s, 1
!  do j=dimy1_s, dimy2_s, 1
!    if ( any ( rtp2_zt(i,j,:) < 0.0 ) ) then
!    print*, 'After advect, rtp2_zt at ', i, j, " = ",  rtp2_zt(i,j,:)
!    end if
!  end do ! i
!end do ! j
! Now interpolate back to momentum gridpoints.
!   (wp3 is already on the thermodynamic gridpoints.)
! do i=dimx1_s, dimx2_s, 1
!   do j=dimy1_s, dimy2_s, 1
    do i=1, nx, 1
      do j=1, ny, 1
     
        wp2(i,j,:)     = zt2zm_api( real( wp2_zt(i,j,:), kind=core_rknd ) )
        up2(i,j,:)     = zt2zm_api( real( up2_zt(i,j,:), kind=core_rknd ) )
        vp2(i,j,:)     = zt2zm_api( real( vp2_zt(i,j,:), kind=core_rknd ) )
        rtp2(i,j,:)    = zt2zm_api( real( rtp2_zt(i,j,:), kind=core_rknd ) )
        thlp2(i,j,:)   = zt2zm_api( real( thlp2_zt(i,j,:), kind=core_rknd ) )
        rtpthlp(i,j,:) = zt2zm_api( real( rtpthlp_zt(i,j,:), kind=core_rknd ) )
        wprtp(i,j,:)   = zt2zm_api( real( wprtp_zt(i,j,:), kind=core_rknd ) )
        wpthlp(i,j,:)  = zt2zm_api( real( wpthlp_zt(i,j,:), kind=core_rknd ) )

      end do ! j
    end do   ! i

    ! Clip variances where the top point is negative
    where ( wp2(:,:,nz) < 0._core_rknd ) wp2(:,:,nz) = 0._core_rknd
    where ( up2(:,:,nz) < 0._core_rknd ) up2(:,:,nz) = 0._core_rknd
    where ( vp2(:,:,nz) < 0._core_rknd ) vp2(:,:,nz) = 0._core_rknd
    where ( rtp2(:,:,nz) < 0._core_rknd ) rtp2(:,:,nz) = 0._core_rknd
    where ( thlp2(:,:,nz) < 0._core_rknd ) thlp2(:,:,nz) = 0._core_rknd

    ! Clip variances where the bottom point is negative
    where ( wp2(:,:,1) < 0._core_rknd ) wp2(:,:,1) = 0._core_rknd
    where ( up2(:,:,1) < 0._core_rknd ) up2(:,:,1) = 0._core_rknd
    where ( vp2(:,:,1) < 0._core_rknd ) vp2(:,:,1) = 0._core_rknd
    where ( rtp2(:,:,1) < 0._core_rknd ) rtp2(:,:,1) = 0._core_rknd
    where ( thlp2(:,:,1) < 0._core_rknd ) thlp2(:,:,1) = 0._core_rknd


!do i=1, nx, 1
!  do j=1, ny, 1
!    if ( any ( rtp2(i,j,:) < 0.0 ) ) then
!    print*, 'After interp, rtp2 at ', i, j, " = ",  rtp2(i,j,:)
!    end if
!  end do ! i
!end do ! j
!
!print*, 'After interp back, wp2(nx,ny,:) =', wp2(nx,ny,:)
!! End of Vince Larson's changes.
    end if ! ladvect

    call t_stopf('advance_clubb') ! For timing

    return
  end subroutine advance_clubb_sgs

!-------------------------------------------------------------------------------
  subroutine apply_clubb_sgs_tndcy_mom( dudt, dvdt )

! Description:
!   Applies the tendency of CLUBB's contribution to the horizontal wind fields.
!
! References:
!   None
!-------------------------------------------------------------------------------


    use grid, only: &
      nx, nxp1, ny, nyp1, dimx1_s, dimx2_s, dimy1_s, dimy2_s, nz, nzm, na, &
      rank

    use clubbvars, only: &
      u_tndcy, & ! CLUBB contribution to the x wind
      v_tndcy  ! CLUBB contribution to the y wind

    implicit none

    intrinsic :: any

    ! In variables
    real, intent(inout), dimension(nxp1,ny,nzm,3) :: &
      dudt ! u wind tendency [m/s^2]

    real, intent(inout), dimension(nx,nyp1,nzm,3) :: &
      dvdt ! v wind tendency [m/s^2]

    ! --- Begin Code --- 

    call t_startf('apply_clubb_sgs_tndcy_mom') ! For timing

    ! Since dudt/dvdt are already time tendencies, we just add the contribution 
    ! to the existing SAM contribution
    dudt(1:nx,1:ny,1:nzm,na) = dudt(1:nx,1:ny,1:nzm,na) + real( u_tndcy(1:nx,1:ny,1:nzm) )
    dvdt(1:nx,1:ny,1:nzm,na) = dvdt(1:nx,1:ny,1:nzm,na) + real( v_tndcy(1:nx,1:ny,1:nzm) )

    call t_stopf('apply_clubb_sgs_tndcy_mom') ! For timing

    return
  end subroutine apply_clubb_sgs_tndcy_mom

!-------------------------------------------------------------------------------
  subroutine apply_clubb_sgs_tndcy_scalars( dt, t, qv, qcl )

! Description:
!   Applies the tendency of CLUBB's contribution over nclubb subcycles of SAM's
!   dynamical core timestep (i.e. if dt == 10. and nclubb = 6, then only 1/6 of
!   CLUBB's contribution will be applied).  This subroutine adjust variables not
!   related to SAM's wind.
!
! References:
!   None
!-------------------------------------------------------------------------------
    use grid, only: &
      nx, nxp1, ny, nyp1, dimx1_s, dimx2_s, dimy1_s, dimy2_s, nz, nzm, na, &
      rank

    use vars, only: rho

    use domain, only: &
      ntracers

    use tracers, only: &
      tracer

    use clubbvars, only: &
      t_tndcy, & ! CLUBB contribution to moist static energy
      qc_tndcy,& ! CLUBB contribution to liquid water mixing ratio
      qv_tndcy   ! CLUBB contribution to vapor water mixing ratio

    use clubbvars, only: &
      tracer_tndcy

    use clubbvars, only: &
      sclr_dim, & ! Constant(s)
      edsclr_dim

    use clubbvars, only: &
     rho_ds_zt, & ! Variable(s)
     rho_ds_zm

    use clubb_api_module, only: &
      clubb_at_least_debug_level_api, & ! Procedure(s)
      fill_holes_vertical_api, &
      fstderr ! Constant

    implicit none

    intrinsic :: any

    ! In variables
    real(kind=time_precision), intent(in) :: &
      dt ! Timestep [s]

    ! In/Out variables
    real, intent(inout), dimension(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm) :: &
      t     ! Moist static energy           [K]

    real, intent(inout), dimension(nx,ny,nzm) :: &
      qv, & ! Water vapor mixing ratio                  [kg/kg]
      qcl   ! Liquid water mixing ratio (condensate)    [kg/kg]

    ! Local Variables
    real(kind=core_rknd), dimension(nz) :: tmpqv, tmpqcl

    real(kind=core_rknd) :: threshold ! Threshold on clipping [units vary]

    integer :: i, j, ig, jg

    ! --- Begin Code --- 

    call t_startf('apply_clubb_sgs_tndcy_scalar') ! For timing

    tmpqv  = 0.0_core_rknd
    tmpqcl = 0.0_core_rknd

    ! Add clubb tendency to qv, qc, t, and tracers
    do i = 1, nx, 1
      do j = 1, ny, 1

        t(i,j,1:nzm) = t(i,j,1:nzm) + real( dt*t_tndcy(i,j,1:nzm) )

        tmpqv(2:nz)  = real( qv(i,j,1:nzm), kind=core_rknd )  + dt*qv_tndcy(i,j,1:nzm)
        tmpqcl(2:nz) = real( qcl(i,j,1:nzm), kind=core_rknd ) + dt*qc_tndcy(i,j,1:nzm)

        if ( edsclr_dim > 0 .or. sclr_dim > 0 ) then
          tracer(i,j,1:nzm,1:ntracers) = tracer(i,j,1:nzm,1:ntracers) &
            + real( dt*tracer_tndcy(i,j,1:nzm,1:ntracers) )
        end if

        ! Apply hole-filling scheme to qv as needed
        threshold = 0._core_rknd
        if ( any( tmpqv(2:nz) < threshold ) ) then

          ! CLUBB's tendency in this column will produce a negative vapor water,
          ! so we apply hole-filling
          if ( clubb_at_least_debug_level_api( 1 ) ) then
            call task_rank_to_index( rank, ig, jg )
            write(fstderr,*) "Task #:", rank
            write(fstderr,*) "Applying hole-filling scheme to vapor water mixing ratio at:", &
              "nx=", i, ";", "ny=", j, ";"
            write(fstderr,*) "x global=", i+ig, ";", "y global=", j+jg, ";"
          end if

          call fill_holes_vertical_api( 2, threshold, "zt", rho_ds_zt, rho_ds_zm, tmpqv )

        end if

        ! Update qv
        qv(i,j,1:nzm) = real( tmpqv(2:nz) )

        threshold = 0._core_rknd
        ! Apply hole-filling scheme to qcl as needed
        if ( any( tmpqcl(2:nz) < threshold ) ) then

          ! CLUBB's tendency in this column will produce a negative cloud water,
          ! so we apply hole-filling
          if ( clubb_at_least_debug_level_api( 1 ) ) then
            call task_rank_to_index( rank, ig, jg )
            write(fstderr,*) "Task #:", rank
            write(fstderr,*) "Applying hole-filling scheme to cloud water mixing ratio at:", &
              "nx=", i, ";", "ny=", j, ";"
            write(fstderr,*) "x global=", i+ig, ";", "y global=", j+jg, ";"
          end if

          call fill_holes_vertical_api( 2, threshold, "zt", rho_ds_zt, rho_ds_zm, tmpqcl )

        end if

        ! Update qcl
        qcl(i,j,1:nzm) = real( tmpqcl(2:nz) )

      end do ! j = 1, ny
    end do ! i = 1, nx

    call t_stopf('apply_clubb_sgs_tndcy_scalar') ! For timing

    return
  end subroutine apply_clubb_sgs_tndcy_scalars

!-------------------------------------------------------------------------------
  subroutine clubb_sgs_cleanup( )
!   Description:
!     De-allocate memory and exit.
!
!   References:
!     None
!-------------------------------------------------------------------------------
    use grid, only: rank

    use clubb_api_module, only: &
      stats_finalize_api

    implicit none

    !----- Begin Code -----

    call cleanup_clubb_core_api( .true. )

    if ( is_a_sample_node( rank ) ) then
      call stats_finalize_api( )
    end if

    return
  end subroutine clubb_sgs_cleanup

!-------------------------------------------------------------------------------
  FUNCTION LIN_EXT( var_high, var_low, height_high, height_low, height_ext )

! Author: Brian M. Griffin,  UW Milwaukee

! References: None

! Description:
! This function computes a linear extension of the value of variable.
! Given two known values of a variable at two height values, the value
! of that variable at a height outside of those two height levels 
! (rather than a height between those two height levels) is computed.
!
! Here is a diagram:
!
!  -------------------------------- Height to be extended to; linear extension
!
!  ################################ Height high, know variable value
!
!
!
!  ################################ Height low, know variable value
!
!
!
!  -------------------------------- Height to be extended to; linear extension
!
!
! FORMULA:
!
! variable(@ Height extension) =
!
! [ (variable(@ Height high) - variable(@ Height low)) / (Height high - Height low) ]
! * (Height extension - Height high)  +  variable(@ Height high)
!-------------------------------------------------------------------------------

    IMPLICIT NONE

    ! Input Variables
    REAL(kind=core_rknd), INTENT(IN):: var_high
    REAL(kind=core_rknd), INTENT(IN):: var_low
    REAL(kind=core_rknd), INTENT(IN):: height_high
    REAL(kind=core_rknd), INTENT(IN):: height_low
    REAL(kind=core_rknd), INTENT(IN):: height_ext

    ! Output Variable
    REAL(kind=core_rknd):: lin_ext

    !----- Begin Code -----

    lin_ext = ( var_high - var_low ) / ( height_high - height_low ) &
         * ( height_ext - height_high ) + var_high

    RETURN
  END FUNCTION LIN_EXT

  !-----------------------------------------------------------------------------
  logical function is_a_sample_node( rank )

  ! Description:
  !   Determine if we're output single-columns stats from this node.
  !
  ! References:
  !   None
  !-----------------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: any, spread, size

    ! Input Variable
    integer, intent(in) :: rank

    integer :: iter

    ! ---- Begin Code ----

    ! Initialize
    is_a_sample_node = .false.

    ! Determine if we're sampling a column of stats from this node
    do iter = 1, size( sample_nodes )
      if ( sample_nodes(iter) == rank ) then
        is_a_sample_node = .true.
        exit
      end if
    end do

    return
  end function is_a_sample_node
 !-----------------------------------------------------------------------------

#ifdef SILHS
!------------------------------------------------------------------------------
  function convert_thl_to_t_LH( lh_thl, gamaz, prespot, X_nl_all_levs ) &
    result( lh_t )

! Description:
!   This function converts from 'θ_l' (CLUBB's conserved temperature variable) 
!   to the conserved temperature quantity denoted in SAM as 't' (moist static 
!   energy).  This is necessary because SILHS generates the sample points in
!   terms of former and SAM and it's microphysics schemes need the latter.
!
! References:
!   None
!------------------------------------------------------------------------------
    use grid, only: nzm, nz ! Constant(s)

    use calc_vars_util, only: &
      thetal2t  ! Procedure(s)

    use clubb_api_module, only: &
      dp, &   ! Constant(s)
      core_rknd, &
      iiPDF_chi, &
      iiPDF_rr, &
      iiPDF_rs, &
      iiPDF_ri

    use clubb_api_module, only: &
      d_variables ! Variable(s)

    use clubb_silhs_vars, only: &
      lh_microphys_calls ! Variable(s)

    implicit none

    ! Input Variables
    real(kind=core_rknd), dimension(nz,lh_microphys_calls), intent(in) :: &
      lh_thl  ! Sample of thetal  [K]

    real, dimension(nzm), intent(in) :: &
      gamaz,  & ! grav/Cp*z         [m]
      prespot   ! 1/exner           [-]

    real(kind=dp), dimension(nz,lh_microphys_calls,d_variables), intent(in) :: &
      X_nl_all_levs   ! All lognormal variates [units vary]

    ! Output Variables
    real(kind=core_rknd), dimension(nzm,lh_microphys_calls) :: &
      lh_t ! Latin hypercube samples of moist static energy  [K]

    ! Local variables
    real, dimension(nzm,lh_microphys_calls) :: &
      qcl ! Liquid water  [kg/kg]

    real, dimension(nzm,lh_microphys_calls) :: &
      qpl, qci, qpi ! Rain, ice, and snow mixing ratio [kg/kg]

    integer :: indx ! Loop index

    ! ---- Begin Code ----

    qcl = 0.0
    qpl = 0.0
    qci = 0.0
    qpi = 0.0

    if ( iiPDF_chi > 0 ) qcl = real( max( X_nl_all_levs(2:nz,:,iiPDF_chi), &
                                          0._dp ) )
    if ( iiPDF_rr > 0 )  qpl = real( X_nl_all_levs(2:nz,:,iiPDF_rr) )
    if ( iiPDF_ri > 0 )  qci = real( X_nl_all_levs(2:nz,:,iiPDF_ri) )

    ! Note: this assumes no graupel samples
    if ( iiPDF_rs > 0 )  qpi = real( X_nl_all_levs(2:nz,:,iiPDF_rs) )

    forall ( indx=1:lh_microphys_calls )
      lh_t(:,indx) &
      = real( thetal2t( real( lh_thl(2:nz,indx) ), gamaz, qpl(:,indx), &
                        qci(:,indx), qpi(:,indx), prespot ), &
              kind = core_rknd )
    end forall

    ! Error checking
    if ( any(  lh_t < 0._dp ) ) then

      print *, "lh_t = "
      write(6,'(40f10.3)') lh_t
      print *, "thetal = "
      write(6,'(40f10.3)') lh_thl
      print *, "qcl = "
      write(6,'(40f10.3)') qcl
      print *, "qpl = "
      write(6,'(40f10.3)') qpl
      print *, "qci = "
      write(6,'(40f10.3)') qci
      print *, "qpi = "
      write(6,'(40f10.3)') qpi
      write(0,*) "Error converting from thetal to moist static energy in SILHS"

      call task_abort( )

    end if

    return
  end function convert_thl_to_t_LH
#endif /* SILHS */

end module clubb_sgs

#endif /*CLUBB*/
