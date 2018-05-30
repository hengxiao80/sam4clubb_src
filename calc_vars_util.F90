!-------------------------------------------------------------------------------
! $Id: calc_vars_util.F90 1458 2014-07-23 22:00:19Z bmg2@uwm.edu $
#if defined( CLUBB ) || defined( UWM_STATS )
!-------------------------------------------------------------------------------
module calc_vars_util

public :: t2thetal, &
          thetal2t

contains

  !-----------------------------------------------------------------------------
  elemental function t2thetal( t, gamaz, qpl, qci, qpi, prespot ) result( thl )

    ! Description:
    ! Convert moist static energy to liquid water potential temperature.
    !
    ! The equation for moist (liquid/ice water) static energy, h_L, is:
    !
    ! h_L = c_p * T + g * z - L_v * ( q_c + q_r ) - L_s * ( q_i + q_s + q_g );
    !
    ! where T is absolute temperature, c_p is specific heat of air at a constant
    ! pressure, g is gravity, z is altitude, L_v is the latent heat of
    ! vaporization, L_s is the latent heat of sublimation, q_c, q_r, q_i, q_s,
    ! and q_g are the mixing ratios of cloud water, rain water, ice, snow, and
    ! graupel, respectively.  This can be rewritten as:
    !
    ! h_L = c_p * T + g * z - L_v * ( qcl + qpl ) - L_s * ( qci + qpi );
    !
    ! where qcl is non-precipitating liquid water (cloud water) mixing ratio,
    ! qpl is precipitating liquid water (rain water) mixing ratio, qci is
    ! non-precipitating frozen water (ice) mixing ratio, and qpi is
    ! precipitating frozen water (snow + graupel) mixing ratio.
    ! 
    ! The SAM model has a variable "t", which is used for moist static energy,
    ! which is defined as:
    !
    ! t = h_L / c_p.
    !
    ! The equation for h_L is rewritten in terms of "t" as:
    !
    ! t = T + ( g / c_p ) * z
    !     - ( L_v / c_p ) * ( qcl + qpl ) - ( L_s / c_p ) * ( qci + qpi ).
    !
    ! The equation for liquid water potential temperature, theta_l, is:
    !
    ! theta_l = theta - ( L_v / ( c_p * exner ) ) * q_c;
    !
    ! where theta is potential temperaure and exner is:
    !
    ! exner = ( p / p0 )^(R_d/c_p);
    !
    ! where p is pressure, p0 is a reference pressure of 1.0 x 10^5 Pa., and R_d
    ! is the gas constant for dry air.  Potential temperature is related to
    ! absolute temperature by:
    !
    ! T = theta * exner;
    !
    ! so the equation for theta_l can be rewritten in terms of T:
    !
    ! theta_l = ( T - ( L_v / c_p ) * q_c ) / exner.
    !
    ! The equation can be written for theta_l in terms of t as:
    !
    ! theta_l = ( t - ( g / c_p ) * z
    !             + ( L_v / c_p ) * qpl + ( L_s / c_p ) * ( qci + qpi ) )
    !           / exner.
    
    ! References:
    ! None
    !---------------------------------------------------------------------------

    use params, only: &
        fac_cond, & ! Variables
        fac_sub

    implicit none

    ! Input Variables
    real, intent(in) :: &
      t,       & ! Moist (liquid/ice) static energy (divided by Cp)   [K]
      gamaz,   & ! grav/Cp*z                                          [K]
      qpl,     & ! Rain water mixing ratio                            [kg/kg]
      qci,     & ! Cloud ice mixing ratio                             [kg/kg]
      qpi,     & ! Snow+Graupel mixing ration                         [kg/kg]
      prespot    ! Exner^-1                                           [-]

    ! Result
    real :: thl  ! Liquid water potential temperature                 [K]


    ! Compute thetal (don't include ice because CLUBB doesn't) 
    thl = ( t - gamaz + fac_cond * qpl + fac_sub * ( qci + qpi ) ) * prespot


    return

  end function t2thetal

  !-----------------------------------------------------------------------------
  elemental function thetal2t( thl, gamaz, qpl, qci, qpi, prespot ) result( t )

    ! Description:
    ! Convert moist static energy to liquid water potential temperature.
    !
    ! The equation for moist (liquid/ice water) static energy, h_L, is:
    !
    ! h_L = c_p * T + g * z - L_v * ( q_c + q_r ) - L_s * ( q_i + q_s + q_g );
    !
    ! where T is absolute temperature, c_p is specific heat of air at a constant
    ! pressure, g is gravity, z is altitude, L_v is the latent heat of
    ! vaporization, L_s is the latent heat of sublimation, q_c, q_r, q_i, q_s,
    ! and q_g are the mixing ratios of cloud water, rain water, ice, snow, and
    ! graupel, respectively.  This can be rewritten as:
    !
    ! h_L = c_p * T + g * z - L_v * ( qcl + qpl ) - L_s * ( qci + qpi );
    !
    ! where qcl is non-precipitating liquid water (cloud water) mixing ratio,
    ! qpl is precipitating liquid water (rain water) mixing ratio, qci is
    ! non-precipitating frozen water (ice) mixing ratio, and qpi is
    ! precipitating frozen water (snow + graupel) mixing ratio.
    ! 
    ! The SAM model has a variable "t", which is used for moist static energy,
    ! which is defined as:
    !
    ! t = h_L / c_p.
    !
    ! The equation for h_L is rewritten in terms of "t" as:
    !
    ! t = T + ( g / c_p ) * z
    !     - ( L_v / c_p ) * ( qcl + qpl ) - ( L_s / c_p ) * ( qci + qpi ).
    !
    ! The equation for liquid water potential temperature, theta_l, is:
    !
    ! theta_l = theta - ( L_v / ( c_p * exner ) ) * q_c;
    !
    ! where theta is potential temperaure and exner is:
    !
    ! exner = ( p / p0 )^(R_d/c_p);
    !
    ! where p is pressure, p0 is a reference pressure of 1.0 x 10^5 Pa., and R_d
    ! is the gas constant for dry air.  Potential temperature is related to
    ! absolute temperature by:
    !
    ! T = theta * exner;
    !
    ! so the equation for theta_l can be rewritten in terms of T:
    !
    ! theta_l = ( T - ( L_v / c_p ) * q_c ) / exner.
    !
    ! The equation can be written for t in terms of theta_l as:
    !
    ! t = theta_l * exner + ( g / c_p ) * z
    !     - ( L_v / c_p ) * qpl - ( L_s / c_p ) * ( qci + qpi ).
    
    ! References:
    ! None
    !---------------------------------------------------------------------------

    use params, only: &
        fac_cond, & ! Variables
        fac_sub

    implicit none

    ! Input Variables
    real, intent(in) :: &
      thl,     & ! Liquid water potential temperature    [K]
      gamaz,   & ! grav/Cp*z                             [K]
      qpl,     & ! Rain water mixing ratio               [kg/kg]
      qci,     & ! Cloud ice mixing ratio                [kg/kg]
      qpi,     & ! Snow+Graupel mixing ration            [kg/kg]
      prespot    ! Exner^-1                              [-]

    ! Result
    real :: t    ! Moist (liquid/ice) static energy (divided by Cp)   [K]


    ! Compute moist static energy 
    t = thl / prespot + gamaz - fac_cond * qpl - fac_sub * ( qci + qpi )


    return

  end function thetal2t

!-------------------------------------------------------------------------------

end module calc_vars_util

#endif /*CLUBB .or. UWM_STATS*/
