#ifdef UWM_STATS
module compute_chi_module
!===============================================================================

  ! Description:
  ! This module contains functions used to compute chi (s_mellor) and
  ! eta (t_mellor). In addition, the calculation of the beta function is
  ! included.

  ! References:
  !    None
  !-------------------------------------------------------------------------
  
  implicit none

  public :: compute_chi_eta, &
            compute_beta 

contains
  
  !=============================================================================
  subroutine compute_chi_eta( thetal, rt, pres_mb, prespot, & ! Intent(in)
                              chi, eta )                      ! Intent(out)

    ! Description:
    ! Function to compute and return chi (s_mellor) and eta (t_mellor).

    ! References:
    ! Larson, V.E., et al., 2005:  Supplying local microphysics
    !   parameterizations with information about subgrid variability:  Latin
    !   hypercube sampling.  J. Atmos. Sci., 62, 4010-4026.
    !-----------------------------------------------------------------------

    use grid, only: &
        nx,  & ! Variable(s)
        ny,  &
        nzm

    use params, only: &
        cp,    & ! Specific heat of air                     [J/(kg dry air * K)]
        lcond    ! Latent heat of condensation/vaporization [J/(kg water)]

    use compute_correlation_module, only: &
        mean_xy_domain  ! Procedure(s)

    implicit none

    ! Input Variable(s)
    real, dimension(nx,ny,nzm), intent(in) :: &
      thetal, & ! Liquid water potential temperature     [K]
      rt        ! Total water mixing ratio (cloud+vapor) [kg/kg]

    real, dimension(nzm), intent(in) :: &
      pres_mb, & ! Pressure                              [mb]
      prespot    ! 1/exner                               [-]

    ! Output Variable(s)
    real, dimension(nx,ny,nzm), intent(out) :: &
      chi, & ! Extended liquid water mixing ratio        [kg/kg]
      eta    ! Coordinate orthogonal to chi              [kg/kg]

    ! Local Variable(s)
    real, dimension(nzm) :: &
      mean_rt,  & ! Mean total water mixing ratio           [kg/kg]
      mean_thl    ! Mean liquid water potential temperature [K]

    real :: &
      mean_Tl,        & ! Mean liquid water temperature           [K]
      rsat_mean_Tl_p, & ! Saturation mixing ratio at the mean Tl  [kg/kg]
      beta_mean_Tl,   & ! Beta function at the mean Tl            [1/(kg/kg)]
      mean_chi,       & ! Mean extended liquid water mixing ratio [kg/kg]
      crt,            & ! Coefficient of rt in chi/eta equations  [-]
      cthl              ! Coefficient of thl in chi/eta equations [(kg/kg)/K]

    real, external :: qsatw

    integer :: i, j, k  ! Loop indices


    ! Find the horizontal mean values of rt and thetal.
    mean_rt  = mean_xy_domain( nzm, rt )
    mean_thl = mean_xy_domain( nzm, thetal )

    do k = 1, nzm, 1

       ! Calculate the horizontal mean value of liquid water temperature, Tl.
       mean_Tl = mean_thl(k) / prespot(k)

       ! Calculate saturation mixing ratio with respect to liquid water at the
       ! mean Tl.
       rsat_mean_Tl_p = qsatw( mean_Tl, pres_mb(k) )

       ! Calculate the value of the beta function at the mean Tl.
       beta_mean_Tl = compute_beta( mean_Tl )

       ! Calculate the mean value of extended liquid water mixing ratio, chi.
       mean_chi = ( mean_rt(k) - rsat_mean_Tl_p ) &
                  / ( 1.0 + beta_mean_Tl * rsat_mean_Tl_p )

       ! Calculate the coefficient of rt in the equations for chi and eta.
       crt = 1.0 / ( 1.0 + beta_mean_Tl * rsat_mean_Tl_p )

       ! Calculate the coefficient of thl in the equations for chi and eta.
       cthl &
       = ( ( 1.0 + beta_mean_Tl * mean_rt(k) ) * beta_mean_Tl * rsat_mean_Tl_p &
           / ( 1.0 + beta_mean_Tl * rsat_mean_Tl_p )**2 ) &
         * ( cp / lcond ) / prespot(k)

       do i = 1, nx, 1
          do j = 1, ny, 1

             ! Calculate the value of extended liquid water mixing ratio, chi.
             chi(i,j,k) = mean_chi + crt * ( rt(i,j,k) - mean_rt(k) ) &
                                   - cthl * ( thetal(i,j,k) - mean_thl(k) )

             ! Calculate the value of eta.
             eta(i,j,k) = crt * ( rt(i,j,k) - mean_rt(k) ) &
                          + cthl * ( thetal(i,j,k) - mean_thl(k) )

          enddo ! j = 1, ny, 1
       enddo ! i = 1, nx, 1

    enddo ! k = 1, nzm, 1


    return

  end subroutine compute_chi_eta

  !=============================================================================
  pure function compute_beta( Tl ) &
  result( beta )

    ! Description:
    ! Function to calculate the value of beta based on liquid water temperature.

    ! References:
    ! Lewellen, W. S. and S. Yoh, 1993:  Binormal model of ensemble partial
    !    cloudiness.  J. Atmos. Sci., 50, 1228-1237.
    !-----------------------------------------------------------------------

    use params, only: &
        cp,    & ! Specific heat of air                     [J/(kg dry air * K)]
        lcond, & ! Latent heat of condensation/vaporization [J/(kg water)]
        Rd_Rv, & ! Rd over Rv                               [-]
        rgas     ! Gas constant for dry air                 [J/(kg dry air * K)]

    implicit none

    ! Input Variable
    real, intent(in) :: &
      Tl    ! Liquid water temperature         [K]
    
    ! Return Variable
    real &
      beta    ! Value of the beta function     [(kg dry air)/(kg water)]


    ! Calculate the value of beta.
    ! Note:  the units of beta are kg (dry air) / kg (water), which are inverse
    !        units of a mixing ratio.
    beta = Rd_Rv * ( lcond / ( rgas * Tl ) ) * ( lcond / ( cp * Tl ) )

    
    return

  end function compute_beta

!===============================================================================
  
end module compute_chi_module
#endif /*UWM_STATS*/
