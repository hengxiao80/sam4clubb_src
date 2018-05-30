module rad_simple_uwm
#ifdef CLUBB
!$Id $
  contains

  subroutine rad_simple_driver
! Adapted from rad_simple.f90 for DYCOMS II RF01/RF02 and the Altostatocumulous
! cases by Dave Schanen UWM.

    use grid
    use vars
    use params
    use rad, only: qrad

    use clubb_api_module, only: &
      lin_interpolate_two_points_api, lin_interpolate_on_grid_api, & ! From CLUBB
      time_precision, core_rknd ! Constant(s)

    implicit none

    ! External
    intrinsic :: trim

    ! Constant Parameters
    logical, parameter :: &
      l_center = .true. ! Use centered difference approx. in sunray_sw

    ! Local Variables
    integer i,j,k

    real :: Frad_LW(nz),FTHRL(nzm) ! Radiative flux, heating rate
    real :: Frad_SW(nz), Frad_SW_flip(nz), FTHRS(nzm) ! Short-wave radiative flux, heating rate

    real, dimension(nz) :: LWP, Heaviside
    real(kind=core_rknd), dimension(nz) :: zw_flip, zt_flip
    real, dimension(nzm) :: dzt_flip
    real(kind=core_rknd) :: z_i, ls_div, kap, F0, F1, omega, gc, alvdr, eff_drop_radius, amu0, Fs0
    real(kind=core_rknd), dimension(12) :: amu0_values, Fs_values
    integer :: year, nparam
    logical :: l_heaviside, l_fix_cos_solar_zen

    real(kind=time_precision) :: time_in_sec

    real(kind=core_rknd) :: rtm_at_k, rtm_at_km1, z_at_k, z_at_km1

    if(.not.dolongwave) return

    select case ( trim( case ) )
    case ( "DYCOMS_RF01" )
      ls_div = 3.75e-6
      kap = 85.
      F0 = 70.
      F1 = 22.
      l_heaviside = .false.

      ! The following we only need for DYCOMS II RF01 with doshortwave = .true.
      ! (not for the GCSS specified DYCOMS II RF01 setup)
      eff_drop_radius = 1.0e-5
      alvdr = 0.1
      gc = 0.85
      omega = 0.992
      Fs_values = 1212.75

      nparam = 1
      Fs_values = 1212.75
      year = 2001
      l_fix_cos_solar_zen = .false.

    case ( "DYCOMS_RF02" )
      ls_div = 3.75e-6
      kap = 85.
      F0 = 70.
      F1 = 22.
      l_heaviside = .true.

    case ( "CLEX9_OCT14", "CLEX9_NOV02" )
      ls_div = 3.75e-6
      kap = 94.2
      F0 = 104.
      F1 = 62.
      l_heaviside = .false.
      omega = 0.9965
      gc = 0.86
      alvdr = 0.1 ! Visible direct surface albedo
      eff_drop_radius = 1.e-5

      year = 2001
      l_fix_cos_solar_zen = .false.

      nparam = 10
      ! Cosine of the solar zenith angle        [-]
      amu0_values = (/0.0, 0.01, 0.1, 0.2, 0.3, 0.4, &
        0.5, 0.6, 0.7, 0.8, 0.9, 1.0 /)

      ! The incident of incoming SW insolation at cloud top the
      ! direction of the incoming beam (not the vertical)   [W/m^2]
      Fs_values = (/0.0, 715.86, 1073.577, 1165.0905, 1204.7033, 1227.6898, &
        1243.1772, 1254.5893, 1263.5491, 1270.8668, 1277.0474, 1282.3994 /)

    case ( "NOV11_ALTOCU" )

      ls_div = 3.75e-6
      l_fix_cos_solar_zen = .true.
      amu0 = 0.4329
      kap = 94.2
      F0 = 104.
      F1 = 62.
      l_heaviside = .false.
      omega = 0.9965
      gc = 0.86
      alvdr = 0.1 ! Visible direct surface albedo
      eff_drop_radius = 1.e-5

      year = 1999 ! This shouldn't be needed for fix cosine solar zenith angle
      l_fix_cos_solar_zen = .false.

      nparam = 1
      ! Cosine of the solar zenith angle        [-]
      amu0_values = amu0

      ! The incident of incoming SW insolation at cloud top the
      ! direction of the incoming beam (not the vertical)   [W/m^2]
      Fs_values = 1212.75

    case ( "JUNE25_ALTOCU" )
      ls_div = 3.75e-6
      kap = 100.
      F0 = 107.
      F1 = 61.
      l_heaviside = .false.
      omega = 0.992
      gc = 0.86
      alvdr = 0.1 ! Visible direct surface albedo
      eff_drop_radius = 1.e-5

      year = 1996
      l_fix_cos_solar_zen = .false.

      nparam = 10
      ! Cosine of the solar zenith angle        [-]
      amu0_values = (/0.0, 0.01, 0.1, 0.2, 0.3, 0.4, &
        0.5, 0.6, 0.7, 0.8, 0.9, 1.0 /)

      ! The incident of incoming SW insolation at cloud top the
      ! direction of the incoming beam (not the vertical)   [W/m^2]
      Fs_values = (/0.0, 715.86, 1073.577, 1165.0905, 1204.7033, 1227.6898, &
        1243.1772, 1254.5893, 1263.5491, 1270.8668, 1277.0474, 1282.3994 /)


    case default
      write(0,*) "Error, don't know what to do for this case"
      call task_abort()

    end select

    ! Zero out statistics
    do k=1,nzm
      radlwdn(k) =0.
      radqrlw(k) =0.
      radswdn(k) =0.
      radqrsw(k) =0.
    enddo

    do i=1,nx
      do j=1,ny
        ! Radiation

        ! Compute liquid water path from top of the model
        ! We define liquid water path on momentum levels

        LWP(nz) = 0.0
        do k = nz-1, 1, -1
          LWP(k) = LWP(k+1) + rho(k)*qcl(i,j,k)*(dz*adzw(k))
        end do  ! k = nz..1

        if ( l_heaviside ) then
          ! Find the height of the isotherm rtm = 8.0 g/kg.

          k = 1 ! Start at the ground
          do while ( k <= nzm .and. qcl(i,j,k)+qv(i,j,k) > 8.0e-3 )
            k = k + 1
          end do
          if ( k == nzm+1 .or. k == 1 ) then
            write(0,*) "Identification of 8.0 g/kg level failed"
            write(0,*) "Subroutine: rad_simple_driver. " & 
              // "File:rad_simple_uwm.F90"
            write(0,*) "i j k = ", i, j, k

            ! Do a best guess as to the inversion (RF01)
            ! This is necessary because in SAM standalone a column can sometimes
            ! have a total water profile that doesn't decrease linearly with height,
            ! but when we average in the horizontal the column average does.
            ! -dschanen 16 Feb 2010
            do k = 1, nzm
              if ( z(k) >= 800. ) exit
            end do
!       write(0,'(A3,3A12)') " k ", "height", "qcl", "qv"
!       do k = 1, nzm
!         write(0,'(I3,3G12.3)') k, z(k), qcl(i,j,k), qv(i,j,k)
!       end do
!       call task_abort()
          end if

          rtm_at_k = qv(i,j,k)+qcl(i,j,k)
          z_at_k = z(k)
          rtm_at_km1 = qv(i,j,k-1)+qcl(i,j,k-1)
          z_at_km1 = z(k-1)
          z_i = lin_interpolate_two_points_api( 8.0e-3_core_rknd, rtm_at_k, rtm_at_km1, z_at_k, z_at_km1 )

          ! Compute the Heaviside step function for z - z_i.
          do k = 1, nz, 1
            if ( zi(k) - z_i  <  0.0 ) then
              Heaviside(k) = 0.0
            else if ( zi(k) - z_i  ==  0.0 ) then
              Heaviside(k) = 0.5
            else if ( zi(k) - z_i  >  0.0 ) then
              Heaviside(k) = 1.0
            end if
          end do
        else
          Heaviside(:) = 0.0
        end if

        ! Compute radiative flux profile (Frad).
        ! Radiative flux is defined on momentum levels.

        do k = 1, nz, 1

          Frad_LW(k) = F0 * EXP( -kap * LWP(k) ) & 
                     + F1 * EXP( -kap * (LWP(1) - LWP(k)) )

          if ( Heaviside(k) > 0.0 .and. l_heaviside ) then
            Frad_LW(k) = Frad_LW(k) & 
                    + rhow(k) * Cp * ls_div * Heaviside(k) & 
                      * ( 0.25 * ((zi(k)-z_i)**(4.0/3.0)) & 
                    + z_i * ((zi(k)-z_i)**(1.0/3.0)) )
          else if ( z_i > 0. .and. zi(k) > z_i ) then
            Frad_LW(k) = Frad_LW(k) & 
                       + rhow(k) * Cp * ls_div & 
                         * ( 0.25 * ((zi(k)-z_i)**(4.0/3.0)) & 
                       + z_i * ((zi(k)-z_i)**(1.0/3.0)) )
          end if

        end do ! 1..nz


        ! Compute the LW radiative heating rate.
        ! The radiative heating rate is defined on thermodynamic levels.

        do k = 1, nzm, 1
          FTHRL(k) = ( prespot(k) ) * ( -1.0/(Cp*rho(k)) ) & 
                   * ( Frad_LW(k+1) - Frad_LW(k) ) / ( dz*adz(k) )
        end do
        ! Compute cosine of the solar zenith angle for the SW calculation
        if ( .not. l_fix_cos_solar_zen ) then
          ! Determine the day based on the fractional day in SAM.
          time_in_sec = day * 86400._time_precision
          amu0 = cos_solar_zen( 1, 1, year, time_in_sec, &
                                latitude(i,j), longitude(i,j)  )
          amu0 = max( amu0, 0._core_rknd ) ! The sunray_sw code can't handle negative cosine.
        end if

        if ( doshortwave .and. amu0 > 0. ) then

          ! Determine the value for Fs0
          if ( nparam > 1 ) then
            call lin_interpolate_on_grid_api( nparam, amu0_values(1:nparam), &
                                       Fs_values(1:nparam), amu0, Fs0 )
          else
            Fs0 = Fs_values(1)
          end if

          ! Determine grid information for sunray_sw
          zw_flip = zi(nz:1:-1) ! Vertical velocity levels in SAM
          zt_flip(nz) = -z(1) ! This is the ghost point found in CLUBB  (below ground)
          do k = nzm, 1, -1
            zt_flip(k) = z(nz-k+1) ! Pressure levels from SAM
          end do
          dzt_flip(1:nzm) = adz(1:nzm) * dz ! Difference in level

          call sunray_sw( qcl(i,j,nzm:1:-1), rho(nzm:1:-1), amu0, dzt_flip, nzm, &
                          zw_flip, zt_flip, &
                          eff_drop_radius, alvdr, gc, Fs0, omega, l_center, &
                          Frad_SW_flip )
          do k = nz, 1, -1
            Frad_SW(k) = Frad_SW_flip(nz-k+1)
          end do
          do k = 1, nzm, 1
            FTHRS(k) = ( prespot(k) ) * ( -1.0/(Cp*rho(k)) ) & 
                     * ( Frad_SW(k+1) - Frad_SW(k) ) / ( dz*adz(k) )
          end do

        else
          FTHRS(:) = 0.
          Frad_SW(:) = 0.
        end if ! doshortwave
        do k=1,nzm
          t(i,j,k) = t(i,j,k) + ( FTHRL(k)+FTHRS(k) ) * dtn ! add radiative heating to sli
          radlwdn(k) = radlwdn(k) + Frad_LW(k)   ! net lw flux for stats
          radswdn(k) = radswdn(k) + Frad_SW(k)   ! net lw flux for stats
          radqrlw(k) = radqrlw(k) + FTHRL(k)  ! net lw heating for stats
          radqrsw(k) = radqrsw(k) + FTHRS(k) ! net sw heating for stats
          qrad(i,j,k) = FTHRL(k) + FTHRS(k) ! store radiative heating for 3d output/stepout
        enddo

      end do
    end do


  end subroutine rad_simple_driver

!-----------------------------------------------------------------------
  function cos_solar_zen & 
           ( day, month, year, current_time, lat_in_degrees, & 
             lon_in_degrees )

! Description:
!   A function based on coefficients from Liou and the Clayson and
!   Curry formula.  Approximates the cosine of the solar zenith
!   angle anywhere in the world based on current Greenwich mean
!   time and the latitude and longitude.

! References:
!   Clayson and Curry formula from C. A. Clayson and J. A. Curry ,
!   J. Geophys.
!   Res. Vol. 101, No. C12, Pages 28515-28528, 15 Dec. 1996.
!   Liou ``An Introduction to Atmospheric Radiation''
!     Table 2.2 and Eqn. 2.2.10
!-----------------------------------------------------------------------

    use clubb_api_module, only: &
      pi_dp, fstderr, sec_per_day, sec_per_hr,& ! Variable(s)
      radians_per_deg_dp, &

      time_precision, dp, & ! Variable(s)

      gregorian2julian_day_api, compute_current_date_api, leap_year_api ! Procedure(s)

    implicit none

    ! External
    intrinsic :: sin, cos, mod, abs, int

    ! Constant Parameters

    ! Liou's coefficients
    real( kind = dp ), parameter :: & 
      c0 =  0.006918_dp,   & ! [-]
      c1 = -0.399912_dp,   & ! [-]
      c2 = -0.006758_dp,   & ! [-]
      c3 = -0.002697_dp,   & ! [-]
      d1 =  0.070257_dp,   & ! [-]
      d2 =  0.000907_dp,   & ! [-]
      d3 =  0.000148_dp      ! [-]

    ! Input Variables
    integer, intent(in) ::  & 
      day,    & ! Day of month at model start
      month,  & ! Month of year at model start
      year      ! Year at model start

    real(kind=time_precision), intent(in) :: & 
      current_time   ! Current time since start date [s]

    real, intent(in) :: & 
      lat_in_degrees, & ! Latitude       [degrees_N]
      lon_in_degrees    ! Longitude      [degrees_E]

    ! Return Variable
    real :: cos_solar_zen

    ! Local Variables
    real( kind = dp ) :: & 
      t,  & 
      delta, & 
      zln, & 
      longang, & 
      latang, & 
      hour,  & 
      present_time

    integer :: & 
      jul_day, days_in_year, & 
      present_year, present_month, present_day

    call compute_current_date_api( day, month, &
                               year, & 
                               current_time, & 
                               present_day, present_month, & 
                               present_year, & 
                               present_time )

    jul_day = gregorian2julian_day_api( present_day, present_month, &
                               present_year )

    if ( leap_year_api( present_year ) ) then
      days_in_year = 366
    else
      days_in_year = 365
    end if

    !delta_in_degrees = delta*(180/pi_dp)

    ! Compute hour angle (old code)
    ! h = 2*pi_dp*t_since_noon/86400

    ! Determine the number of hours
    hour = present_time / sec_per_hr

    t = 2._dp*pi_dp*dble( jul_day-1 )/dble( days_in_year )

    delta = c0  & 
          + c1*cos( t ) + d1*sin( t ) & 
          + c2*cos( 2._dp*t ) + d2*sin( 2._dp*t ) & 
          + c3*cos( 3._dp*t ) + d3*sin( 3._dp*t )

! The angle  longang  is equivalent to the
! hour angle in the formula for cosZ .
! References: Source file zenith.f
!   from http://magic.gfdi.fsu.edu/seaflux/DIURNAL/README.txt
!   Clayson and Curry formula from C. A. Clayson and J. A. Curry ,
!   J. Geophys.
!   Res. Vol. 101, No. C12, Pages 28515-28528, 15 Dec. 1996 .

!   June 6, 2006

    select case ( int( hour ) )
    case ( 0:11 )
      zln = 180.00_dp - hour*15.00_dp ! Known magic number
    case ( 12:23 )
      zln = 540.00_dp - hour*15.00_dp ! Known magic number
    case default
      zln = 0.0_dp
      write(unit=fstderr,fmt=*) "Hour=", hour
      stop " > 24 hours in cosine solar zenith code"
    end select

    longang = abs( dble( lon_in_degrees ) - zln ) * radians_per_deg_dp
    latang  = dble( lat_in_degrees ) * radians_per_deg_dp


    ! Cosine of the solar zenith angle (sometimes denoted amu0).
    cos_solar_zen = real( sin( latang ) * sin( delta ) & 
                  + cos( latang ) * cos( delta ) * cos( longang ) )

    !write(*,'(a,f15.6)') "cosine solar zenith angle=", cos_solar_zen !%% debug

    return

  end function cos_solar_zen
!-----------------------------------------------------------------------
  subroutine sunray_sw( qc3, rbm, xi_abs, dsigm, kk, & 
                        coamps_zm, coamps_zt, & 
                        radius, A, gc, Fs0, omega, l_center, & 
                        Frad_SW )

! Description:
! for CLEX altocumulus case
! Written by Geert Lenderink for implementation of Shettle and
! Weinman's formulation for radiative flux into the EUROCS
! stratocumulus case.
!
! Adapted by Vince Larson, Chris Golaz, Adam Smith, Michael Falk, and
! others for COAMPS and CLUBB.
!
! Subroutine to compute shortwave radiative flux.
!
!
! The code for sunray_sw was slightly reconstructed in order to more
! closely follow Geert Lenderink's code for the Duynkerke et al.
! EUROCS case.  The original formulation used in that paper comes from
! Shettle and Weinman case.
!
! Tau is now computed inside this routine, as it's not needed in
! nov11_rad.  Tau is also computed for each layer instead of as the
! total optical depth from the top of the domain to the bottom of the
! layer being computed.  This makes tau zero everywhere outside of
! cloud.  F_diff and F_dir are also zero outside of cloud.
!
! Comments by Michael Falk, 26 January 2005.
!
! ADDITION TO COMMENT BY ADAM SMITH, 28 July 2006.
! We are attempting to compare CLUBB results with the COAMPS-LES model.
! To keep the two codes as similar as possible, this subroutine is a
! near-duplicate of the rad_lwsw subroutine  used in COAMPS_LES.  I
! have  modified variable declarations as needed to match CLUBB
! standards, and sections that do not apply to CLUBB have been removed.
!
! For a full explanation of how this subroutine is implemented in CLUBB,
! please see the nov11_altocu_tndcy or jun25_altocu_tndcy subroutines
! within gcss.F.

! REFERENCES:
! see subroutine rad_lwsw.
!-----------------------------------------------------------------------

    use clubb_api_module, only: &
      Cp, rho_lw, pi, & ! Variable(s)
      lin_interpolate_two_points_api, & ! Procedure(s)
      core_rknd ! Constant(s)

    implicit none

! Input variables

    integer, intent(in) :: kk

    real, dimension(kk), intent(in) ::  & 
      qc3, & 
      rbm, & 
      dsigm

    real(kind=core_rknd), dimension(kk+1), intent(in) :: & 
      coamps_zm,  & ! Altitude of momentum levels w/ COAMPS grid indices      [m]
      coamps_zt  ! Altitude of thermodynamic levels w/ COAMPS grid indices [m]

    real(kind=core_rknd), intent(in) ::  & 
      xi_abs, & 
      radius,  & 
      A,  & 
      gc,  & 
      Fs0,  & 
      omega

    logical, intent(in) ::  & 
      l_center

! Output variables
    real, dimension(kk+1), intent(out) ::  & 
      Frad_SW


! Local Variables
    real, dimension(kk+1) :: & 
      tau,     & ! Optical depth of an incremental layer.        [-]
      F_diff,  & ! Diffuse component of SW radiation             [W/m^2]
      F_dir      ! Diffuse component of LW radiation             [W/m^2]

    real(kind=core_rknd), dimension(kk+1) :: & 
      taude      ! Delta-Eddington transformation of tau.        [-]

    real :: taupath, tauc, t1, t2, t3, c1, c2, omegade, & 
         x1, x2, x3, rk, rk2, xi_abs2, rp, alpha, beta, rtt, & 
         exmu0, expk, exmk, xp23p, xm23p, ap23b, taucde

    integer :: k

    real :: ff, gcde

!-----------------------------------------------------------------------
!  CONSTANTS/PARAMETERS
!
!  values added by Michael Falk and Adam Smith
!
! ff     : gc^2, denoted "g^2" in Duynkerke.            Unit: none
! gcde   : Delta-Eddington transformation of gc.
!          Notated g-prime in Duynkerke eqn.18.         Unit: None
!
!-----------------------------------------------------------------------

    ff = gc*gc
    gcde = gc/(1.0+gc)

!-----------------------------------------------------------------------
!
!  Computation of tau and omega variables
!
!
! tauc   : column total optical depth                   Unit: none
! taupath: column total Delta-Eddington optical depth.  Unit: none
! omega  : single-scattering albedo                     Unit: none
! omegade: D-E omega-- from Duynkerke eqn.18.           Unit: none
! taucde : D-E tauc -- from Duynkerke eqn.18.           Unit: none
! taude  : D-E tau  -- from Duynkerke eqn.18.           Unit: none
!
! Comments by Michael Falk, 26 January 2005
!
!-----------------------------------------------------------------------


    tauc = 0.0
    do k=1, kk

      tau(k) = 1.5 * qc3(k) * rbm(k) * dsigm(k)  & !/ aoz(i,j)
             / radius / rho_lw
      tauc = tauc + tau(k)
    enddo
    tau(kk+1) = tau(kk)

    omegade = (1.0-ff)*omega/(1.0-omega*ff)
    taucde = (1.0-omega*ff)*tauc

    do k=1, kk+1
      taude(k) = (1.0-omega*ff)*tau(k)
    enddo


!-----------------------------------------------------------------------
!
!  Computation of constants for radiative transfer equations
!
!  These variables come from Duynkerke eqn.20, which, with slight
!  modifications, matches Shettle and Weinman between eqns.12b and 13.
!  Duynkerke uses slightly different formulations than Shettle and
!  Weinman:
!
! F0(S&W)    = F0(Duynkerke)*pi.
! alpha(S&W) = alpha(Duynkerke)*F0(S&W).
! beta(S&W)  = beta(Duynkerke)*F0(S&W).
! c1(S&W)    = c1(Duynkerke)*F0(S&W).
! c2(S&W)    = c2(Duynkerke)*F0(S&W).
! x3(S&W)    = x3(Duynkerke)*F0(S&W).
! y3(S&W)    = y3(Duynkerke)*F0(S&W).
!
! F0 is divided out of each term in several equations, and then
! reintroduced when the actual radiative flux is computed.  The
! computations here follow Lenderink's original sunray_sw code which
! uses the Duynkerke formulations for these variables.
!
! x1     : term 1 in k equation                         Unit: none
! x2     : term 2 in k equation                         Unit: none
! rk     : k equation                                   Unit: none
! rk2    : k equation squared                           Unit: none
! x3     : term in denominator of alpha and beta        Unit: none
! rp     : p equation                                   Unit: none
! alpha  : alpha equation                               Unit: none
! beta   : beta equation                                Unit: none
!
! The following variables are used by Lenderink to solve for
! Duynkerke's parameters C1 and C2.  They are all dimensionless.
!
! rtt    : 2/3.
! exmu0  : exponential term used in S&W eqn.14-- originally from
!          eqn 1 (also Salby 9.35) in the source function for SW
!          radiation
! expk   : one of the coefficients of C2 on the left side of Shettle
!          and Weinman eqn.14, originally from eqn.12a and eqn.12b
! exmk   : reciprocal of expk, one of the coefficients of C1 on the
!          left side of Shettle and Weinman eqn.14, originally from
!          eqn.12a and eqn.12b
! xp23p  : coefficient of C1 - left side of Shettle and Weinman eqn.13
! xm23p  : coefficient of C2 - left side of Shettle and Weinman eqn.13
! ap23b  : right side of Shettle and Weinman eqn.13
! t1     : the other coefficient of C1 on left side of Shettle and
!          Weinman eqn.14
! t2     : the other coefficient of C2 on left side of Shettle and
!          Weinman eqn.14
! t3     : the coefficient of exmu0 on the right side of Shettle and
!          Weinman eqn.14
!
! Comments by Michael Falk, 26 January 2005
!
!-----------------------------------------------------------------------

    x1 = 1.0-omegade*gcde
    x2 = 1.0-omegade
    rk = sqrt( 3.0*x2*x1 )
    xi_abs2 = xi_abs*xi_abs
    rk2 = rk*rk
    x3 = 4.0*(1.0-rk2*xi_abs2)
    rp = sqrt( 3.0*x2/x1 )
    alpha = 3.0*omegade*xi_abs2*(1.0+gcde*x2)/x3
    beta = 3.0*omegade*xi_abs*(1.0+3.0*gcde*xi_abs2*x2)/x3

    rtt = 2.0/3.0
    exmu0 = exp( -taucde/xi_abs )
    expk = exp( rk*taucde )
    exmk = 1.0/expk
    xp23p = 1.0+rtt*rp
    xm23p = 1.0-rtt*rp
    ap23b = alpha+rtt*beta

    t1 = 1.-A-rtt*(1.0+A)*rp
    t2 = 1.-A+rtt*(1.0+A)*rp
    t3 = (1.-A)*alpha-rtt*(1.+A)*beta+A*xi_abs


!-----------------------------------------------------------------------
!
! Shettle and Weinman 13 and 14, adapted into Duynkerke, give two
! equations and two unknowns which can be solved by linear combination
! (C2) and then substitution (C1).
!
! Comments by Michael Falk, 26 January 2005
!
!-----------------------------------------------------------------------

    c2 = (xp23p*t3*exmu0-t1*ap23b*exmk) & 
       / (xp23p*t2*expk-xm23p*t1*exmk)
    c1 = (ap23b-c2*xm23p)/xp23p


!-----------------------------------------------------------------------
!
!  Computation of diffuse and direct shortwave flux
!
! F_diff is the first term in Duynkerke eqn.19.  The F0 and pi
! constants which are divided out in the Duynkerke formulation are
! reintroduced here.
!
! Duynkerke eqn.19's F_diff term comes from Shettle and Weinman eqn.8,
! where F_diff = F(upward)-F(downward).  Then F_diff = (-4/3)*pi*I1,
! where I1 is solved in Shettle and Weinman eqn.12b.  Capital P in
! Shettle and Weinman eqn.12b should actually be a lowercase p.
!
! F_dir is the second term in Duynkerke eqn.19.
!
! The negative sign for F_diff and F_dir is related to the definition
! of which direction is a positive flux.
!
! Comments by Michael Falk, 26 January 2005
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! Computation of shortwave fluxes on staggered grid
!
! For a full explanation see the "Computation of radiative fluxes on
! staggered grid" section above.  For Frad_SW to be correctly computed
! on w levels, the non-constant component of F_diff and F_dir, we
! compute taupath on w levels as a centered difference between tau
! values on mass levels.
!
! --------taupath-->--Frad_SW----------------    k = 1  (w level)
!        /                   \        |-dwm
! -taude----------------------radht----------    k = 1  (mass level)
!        \                   /        |-dmw
! --------taupath-->--Frad_SW----------------    k = 2  (w level)
!        /                   \
! -taude----------------------radht----------    k = 2  (mass level)
!
! Vince Larson changed the F variables to w levels  03 Feb 2005
! Michael Falk changed the loop to start at k=2 and then solved
! separately for k=1 so the array didn't go out of bounds.
!
! This code makes the same assumption as above that dwm=dmw.
!
! Comments by Michael Falk, 16 February 2005.                          c
!                                                                      c
! ADDITIONAL NOTE: The CLUBB parameterization is now set up to be
!                  compatible with the use of a stretched
!                  (or unevenly-spaced) grid, as well as with the use
!                  of an evenly-spaced grid.  Therefore, dwm is not
!                  necessarily equal to dmw.  Interpolation functions
!                  are used to compute any weighted averages, rather
!                  than using general numbers such as (1/2), which is
!                  compatible only with an evenly-spaced grid.
!                  Brian Griffin; May 10, 2008.
!
!-----------------------------------------------------------------------

    if ( l_center ) then
      taupath = 0.5*taude(1)
    else
      taupath = 0.
    endif

    F_diff(1) = (-4.0/3.0) * Fs0 & 
              * (  & 
                 rp * & 
                     (  & 
                        c1*exp( -rk*taupath ) & 
                      - c2*exp( rk*taupath ) & 
                     ) & 
                 - beta*exp( -taupath/xi_abs )  & 
                )
    F_dir(1) = -Fs0*xi_abs*exp( -taupath/xi_abs )
    Frad_SW(1) = F_diff(1) + F_dir(1)

    do k = 2, kk+1

      if ( l_center ) then
        taupath = taupath  & 
                + lin_interpolate_two_points_api( coamps_zm(k), coamps_zt(k-1),  &
                           coamps_zt(k), taude(k-1), taude(k) )
      else
        taupath = taupath + taude(k)
      endif


      F_diff(k) = (-4.0/3.0) * Fs0 & 
                * (  & 
                   rp * & 
                       (  & 
                          c1*exp( -rk*taupath ) & 
                        - c2*exp( rk*taupath ) & 
                       ) & 
                   - beta*exp( -taupath/xi_abs )  & 
                  )
      F_dir(k) = -Fs0*xi_abs*exp( -taupath/xi_abs )
      Frad_SW(k) = F_diff(k) + F_dir(k)

    enddo ! k=2..kk+1

    return
  end subroutine sunray_sw
#endif /* CLUBB */
end module rad_simple_uwm
