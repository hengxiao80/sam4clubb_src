#ifdef UWM_STATS
module compute_correlation_module
!===============================================================================

  ! Description:
  ! This module contains a function that is used to compute the correlation
  ! coefficient of two variables, and a number of functions that are used to
  ! compute mean values, variances, and covariances, either over the entire
  ! horizontal domain or within precipitation.

  implicit none

  public :: compute_correlation,     &
            mean_xy_domain,          &
            variance_xy_domain,      &
            covariance_xy_domain,    &
            mean_ip_xy_domain,       &
            variance_ip_xy_domain,   &
            covariance_ip_xy_domain, &
            sum_value_xy_domain,     &
            sum_value_xy

  ! Value to use for an undefined correlation.
  ! This value needs to be (much) greater than 1 or (much) less than -1.
  real, parameter, public :: undef_corr = -999.0

  ! Number of valid correlation values (for each correlation type and grid
  ! level) sampled during an statistical output period.
  integer, allocatable, public :: corravg_count(:,:)

  ! The value of average_type used for correlations.
  integer, parameter, public :: corr_avg = -1

  private

contains

  !=============================================================================
  function compute_correlation( num_z, variance_x, variance_y, covar_xy ) &
  result( corr_xy )

    ! Description:
    ! Function to compute and return a vertical profile of the correlation
    ! coefficient between variables x and y.

    ! References:
    !    None
    !-----------------------------------------------------------------------

    implicit none
    
    ! Input Variables
    integer, intent(in) :: &
      num_z    ! Num. of z-dir. grid pts. in vertical range

    real, intent(in), dimension(num_z) :: &
      variance_x, & ! Variance of x
      variance_y, & ! Variance of y
      covar_xy      ! Covariance between x and y
    
    ! Return Variable
    real, dimension(num_z) :: &
      corr_xy    ! Correlation between x and y

    ! Local Variables
    real :: &
      stdev_x, & ! Standard deviation of x
      stdev_y    ! Standard deviation of y

    integer :: k ! Loop index

    
    do k = 1, num_z, 1

       if ( variance_x(k) > 0.0 .and. variance_y(k) > 0.0 ) then

          ! Both x and y vary at this level.  Calculate the correlation.
      
          ! Compute the standard deviation of x and y.
          stdev_x = sqrt( variance_x(k) )
          stdev_y = sqrt( variance_y(k) )
      
          ! Compute the correlation of x and y for the current model level.
          corr_xy(k) = covar_xy(k) / ( stdev_x * stdev_y )
      
          ! Assertion check
          if ( corr_xy(k) > 1.0 ) then
             if ( corr_xy(k) < 1.0005 ) then
                ! Correlation should be 1, but is off slightly due to numerical
                ! round-off error.
                corr_xy(k) = 1.0
             else ! corr_xy(k) > 1.0005
                print *, "In compute_correlation, the correlation of x and y" &
                         //" is greater than 1.0005."
                print *, "lev., corr_xy, covar_xy, variance_x, variance_y = ", &
                         k, corr_xy(k), covar_xy(k), &
                         variance_x(k), variance_y(k)
             endif ! corr_xy(k) < 1.0005
          elseif ( corr_xy(k) < -1.0 ) then
             if ( corr_xy(k) > -1.0005 ) then
                ! Correlation should be -1, but is off slightly due to numerical
                ! round-off error.
                corr_xy(k) = -1.0
             else ! corr_xy(k) < -1.0005
                print *, "In compute_correlation, the correlation of x and y" &
                         //" is less than -1.0005."
                print *, "lev., corr_xy, covar_xy, variance_x, variance_y = ", &
                         k, corr_xy(k), covar_xy(k), &
                         variance_x(k), variance_y(k)
             endif ! corr_xy(k) < -1.0005
          endif ! corr_xy(k) > 1.0 or corr_xy(k) < 1.0

       else ! variance_x(k) = 0 or variance_y(k) = 0

          ! Either x or y (or both) is constant at this level.
          ! The correlation is undefined.
          corr_xy(k) = undef_corr

       endif ! variance_x(k) > 0 and variance_y(k) > 0
    
    enddo ! k = 1, num_z, 1


    return  

  end function compute_correlation

  !=============================================================================
  function mean_xy_domain( num_z, variable )

    ! Description:
    ! Calculates the mean value of a variable over the entire horizontal model
    ! domain for each vertical level in the range.

    !-----------------------------------------------------------------------

    use grid, only: &
        nx,          & ! Num. x-dir. grid pts. (in subdom. when nsubdomains > 1)
        ny,          & ! Num. y-dir. grid pts. (in subdom. when nsubdomains > 1)
        nsubdomains    ! Number of subdomains

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      num_z    ! Num. of z-dir. grid pts. in vertical range

    real, dimension(nx,ny,num_z), intent(in) :: &
      variable  ! Variable

    ! Return Variable
    real, dimension(num_z) :: &
      mean_xy_domain  ! Mean value over entire domain

    ! Local Variables
    integer :: &
      num_grd_pts_dom    ! Total number of gridpoints in entire domain

    real, dimension(num_z) :: &
      domain_sum_variable    ! Sum total of the variable over entire domain

    integer :: k    ! Loop index

   
    ! Number of grid points in the entire horizontal domain.
    num_grd_pts_dom = nx * ny * nsubdomains

    ! Sum value of the variable over the entire horizontal domain.
    domain_sum_variable = sum_value_xy_domain( num_z, variable )

    ! Mean value of the variable for the entire horizontal domain.
    do k = 1, num_z, 1
       mean_xy_domain(k) = domain_sum_variable(k) / real( num_grd_pts_dom )
    enddo ! k = 1, num_z, 1


    return
    
  end function mean_xy_domain

  !=============================================================================
  function variance_xy_domain( num_z, variable, mean )

    ! Description:
    ! Calculates the variance of a variable over the entire horizontal model
    ! domain for each vertical level in the range.

    !-----------------------------------------------------------------------

    use grid, only: &
        nx,          & ! Num. x-dir. grid pts. (in subdom. when nsubdomains > 1)
        ny,          & ! Num. y-dir. grid pts. (in subdom. when nsubdomains > 1)
        nsubdomains    ! Number of subdomains

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_z    ! Num. of z-dir. grid pts. in vertical range

    real, dimension(nx,ny,num_z), intent(in) :: &
      variable  ! Variable

    real, dimension(num_z), intent(in) :: &
      mean    ! Mean value of variable

    ! Return Variable
    real, dimension(num_z) :: &
      variance_xy_domain  ! Variance over entire domain

    ! Local Variables
    integer :: &
      num_grd_pts_dom    ! Total number of gridpoints in entire domain

    real, dimension(nx,ny,num_z) :: &
      variable_p_sqd    ! Value of ( variable - mean )^2.

    real, dimension(num_z) :: &
      domain_sum_calc    ! Sum total of calculation over entire domain

    integer :: i, j, k    ! Loop indices

   
    ! Number of grid points in the entire horizontal domain
    num_grd_pts_dom = nx * ny * nsubdomains

    ! Square of the deviation of the variable from horizontal-domain mean value.
    do k = 1, num_z, 1
       do i = 1, nx, 1
          do j = 1, ny, 1
             variable_p_sqd(i,j,k) = ( variable(i,j,k) - mean(k) )**2
          enddo ! j = 1, ny, 1
       enddo ! i = 1, nx, 1
    enddo ! k = 1, num_z, 1

    ! Sum of variable'^2 over the entire horizontal domain.
    domain_sum_calc = sum_value_xy_domain( num_z, variable_p_sqd )

    ! Variance of the variable for the entire horizontal domain.
    do k = 1, num_z, 1
       variance_xy_domain(k) = domain_sum_calc(k) / real( num_grd_pts_dom )
    enddo ! k = 1, num_z, 1


    return
    
  end function variance_xy_domain

  !=============================================================================
  function covariance_xy_domain( num_z, variable1, variable2, &
                                 mean1, mean2 )

    ! Description:
    ! Calculates the covariance between two variables over the entire horizontal
    ! model domain for each vertical level in the range.

    !-----------------------------------------------------------------------

    use grid, only: &
        nx,          & ! Num. x-dir. grid pts. (in subdom. when nsubdomains > 1)
        ny,          & ! Num. y-dir. grid pts. (in subdom. when nsubdomains > 1)
        nsubdomains    ! Number of subdomains

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_z    ! Num. of z-dir. grid pts. in vertical range

    real, dimension(nx,ny,num_z), intent(in) :: &
      variable1, & ! Variable 1
      variable2    ! Variable 2

    real, dimension(num_z), intent(in) :: &
      mean1, & ! Mean value of Variable 1
      mean2    ! Mean value of Variable 2

    ! Return Variable
    real, dimension(num_z) :: &
      covariance_xy_domain  ! Covariance over entire domain

    ! Local Variables
    integer :: &
      num_grd_pts_dom    ! Total number of gridpoints in entire domain

    real, dimension(nx,ny,num_z) :: &
      var1p_var2p    ! Value of ( variable1 - mean1 ) * ( variable2 - mean2 ).

    real, dimension(num_z) :: &
      domain_sum_calc    ! Sum total of calculation over entire domain

    integer :: i, j, k    ! Loop indices

   
    ! Number of grid points in the entire horizontal domain
    num_grd_pts_dom = nx * ny * nsubdomains

    ! Product of the deviation of the variable from its horizontal-domain mean
    ! value for both variables.
    do k = 1, num_z, 1
       do i = 1, nx, 1
          do j = 1, ny, 1
             var1p_var2p(i,j,k) &
             = ( variable1(i,j,k) - mean1(k) ) * ( variable2(i,j,k) - mean2(k) )
          enddo ! j = 1, ny, 1
       enddo ! i = 1, nx, 1
    enddo ! k = 1, num_z, 1

    ! Sum of var1'var2' over the entire horizontal domain.
    domain_sum_calc = sum_value_xy_domain( num_z, var1p_var2p )

    ! Covariance between the variables for the entire horizontal domain.
    do k = 1, num_z, 1
       covariance_xy_domain(k) = domain_sum_calc(k) / real( num_grd_pts_dom )
    enddo ! k = 1, num_z, 1


    return
    
  end function covariance_xy_domain

  !=============================================================================
  function mean_ip_xy_domain( num_z, variable, precipitating, num_precip )

    ! Description:
    ! Calculates the mean of a variable in precipitation over the entire
    ! horizontal model domain for each vertical level in the range.

    !-----------------------------------------------------------------------

    use grid, only: &
        nx, & ! Number of x-dir. grid pts. (in subdomain when nsubdomains > 1)
        ny    ! Number of y-dir. grid pts. (in subdomain when nsubdomains > 1)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_z    ! Num. of z-dir. grid pts. in vertical range

    real, dimension(nx,ny,num_z), intent(in) :: &
      variable,      & ! Variable
      precipitating    ! Precipitation is found when precipitating = 1.0

    integer, dimension(num_z), intent(in) :: &
      num_precip    ! Number of grid points in the entire domain with
                    ! precipitation.

    ! Return Variable
    real, dimension(num_z) :: &
      mean_ip_xy_domain  ! Mean value in precip over entire domain

    ! Local Variables
    real, dimension(nx,ny,num_z) :: &
      variable_ip    ! Value of variable in precipitation

    real, dimension(num_z) :: &
      domain_sum_variable   ! Sum total of the variable over entire domain

    integer :: i, j, k    ! Loop indices


    ! Value of variable in precipitation.
    do k = 1, num_z, 1
       do i = 1, nx, 1
          do j = 1, ny, 1
             variable_ip(i,j,k) = variable(i,j,k) * precipitating(i,j,k)
          enddo ! j = 1, ny, 1
       enddo ! i = 1, nx, 1
    enddo ! k = 1, num_z, 1

    ! Sum value of the variable in precipitation over the entire horizontal
    ! domain.
    domain_sum_variable = sum_value_xy_domain( num_z, variable_ip )

    ! Mean value of the variable for the entire horizontal domain.
    do k = 1, num_z, 1
       if ( num_precip(k) > 0 ) then
          mean_ip_xy_domain(k) = domain_sum_variable(k) / real( num_precip(k) )
       else
          mean_ip_xy_domain(k) = 0.0
       endif
    enddo ! k = 1, num_z, 1


    return
    
  end function mean_ip_xy_domain

  !=============================================================================
  function variance_ip_xy_domain( num_z, variable, precipitating, &
                                  mean_ip, num_precip )

    ! Description:
    ! Calculates the variance of a variable in precipitation over the entire
    ! horizontal model domain for each vertical level in the range.

    !-----------------------------------------------------------------------

    use grid, only: &
        nx, & ! Number of x-dir. grid pts. (in subdomain when nsubdomains > 1)
        ny    ! Number of y-dir. grid pts. (in subdomain when nsubdomains > 1)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_z    ! Num. of z-dir. grid pts. in vertical range

    real, dimension(nx,ny,num_z), intent(in) :: &
      variable,      & ! Variable
      precipitating    ! Precipitation is found when precipitating = 1.0

    real, dimension(num_z), intent(in) :: &
      mean_ip    ! Mean value of variable in precipitation.

    integer, dimension(num_z), intent(in) :: &
      num_precip    ! Number of grid points in the entire domain with
                    ! precipitation.

    ! Return Variable
    real, dimension(num_z) :: &
      variance_ip_xy_domain  ! Variance in precip over entire domain

    ! Local Variables
    real, dimension(nx,ny,num_z) :: &
      variable_ip_p_sqd    ! Value of ( variable - mean )^2 in precipitation.

    real, dimension(num_z) :: &
      domain_sum_calc    ! Sum total of calculation over entire domain

    integer :: i, j, k    ! Loop indices

   
    ! Square of the deviation of the variable from horizontal-domain mean value
    ! in precipitation.
    do k = 1, num_z, 1
       do i = 1, nx, 1
          do j = 1, ny, 1
             variable_ip_p_sqd(i,j,k) &
             = precipitating(i,j,k) * ( variable(i,j,k) - mean_ip(k) )**2
          enddo ! j = 1, ny, 1
       enddo ! i = 1, nx, 1
    enddo ! k = 1, num_z, 1

    ! Sum of variable'^2 in precipitation over the entire horizontal domain.
    domain_sum_calc = sum_value_xy_domain( num_z, variable_ip_p_sqd )

    ! Variance of the variable in precipitation for the entire horizontal
    ! domain.
    do k = 1, num_z, 1
       if ( num_precip(k) > 0 ) then
          variance_ip_xy_domain(k) = domain_sum_calc(k) / real( num_precip(k) )
       else
          variance_ip_xy_domain(k) = 0.0
       endif
    enddo ! k = 1, num_z, 1


    return
    
  end function variance_ip_xy_domain

  !=============================================================================
  function covariance_ip_xy_domain( num_z, variable1, variable2, precipitating,&
                                    mean1_ip, mean2_ip, num_precip )

    ! Description:
    ! Calculates the covariance between two variables in precipitation over the
    ! entire horizontal model domain for each vertical level in the range.

    !-----------------------------------------------------------------------

    use grid, only: &
        nx, & ! Number of x-dir. grid pts. (in subdomain when nsubdomains > 1)
        ny    ! Number of y-dir. grid pts. (in subdomain when nsubdomains > 1)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_z    ! Num. of z-dir. grid pts. in vertical range

    real, dimension(nx,ny,num_z), intent(in) :: &
      variable1,     & ! Variable 1
      variable2,     & ! Variable 2
      precipitating    ! Precipitation is found when precipitating = 1.0

    real, dimension(num_z), intent(in) :: &
      mean1_ip, & ! Mean value of Variable 1 in precipitation
      mean2_ip    ! Mean value of Variable 2 in precipitation

    integer, dimension(num_z), intent(in) :: &
      num_precip    ! Number of grid points in the entire domain with
                    ! precipitation.

    ! Return Variable
    real, dimension(num_z) :: &
      covariance_ip_xy_domain  ! Covariance over entire domain

    ! Local Variables
    real, dimension(nx,ny,num_z) :: &
      var1p_var2p_ip    ! Value of ( variable1 - mean1 ) * ( variable2 - mean2 )
                        ! in precipitation.

    real, dimension(num_z) :: &
      domain_sum_calc    ! Sum total of calculation over entire domain

    integer :: i, j, k    ! Loop indices

   
    ! Product of the deviation of the variable from its horizontal-domain mean
    ! value for both variables in precipitation.
    do k = 1, num_z, 1
       do i = 1, nx, 1
          do j = 1, ny, 1
             var1p_var2p_ip(i,j,k) &
             = precipitating(i,j,k) &
               * ( variable1(i,j,k) - mean1_ip(k) ) &
               * ( variable2(i,j,k) - mean2_ip(k) )
          enddo ! j = 1, ny, 1
       enddo ! i = 1, nx, 1
    enddo ! k = 1, num_z, 1

    ! Sum of var1'var2' in precipitation over the entire horizontal domain.
    domain_sum_calc = sum_value_xy_domain( num_z, var1p_var2p_ip )

    ! Covariance between the variables in precipitation for the entire
    ! horizontal domain.
    do k = 1, num_z, 1
       if ( num_precip(k) > 0 ) then
          covariance_ip_xy_domain(k) &
          = domain_sum_calc(k) / real( num_precip(k) )
       else
          covariance_ip_xy_domain(k) = 0.0
       endif
    enddo ! k = 1, num_z, 1


    return
    
  end function covariance_ip_xy_domain

  !=============================================================================
  function sum_value_xy_domain( num_z, value )

    ! Description:
    ! Calculates the sum of a value (variable or quantity) over the entire
    ! horizontal model domain for each vertical level in the range.

    !-----------------------------------------------------------------------

    use grid, only: &
        nx,    & ! Num. of x-dir. grid pts. (in subdomain when nsubdomains > 1)
        ny,    & ! Num. of y-dir. grid pts. (in subdomain when nsubdomains > 1)
        dompi    ! Flag for multitasking (true when nsubdomains > 1)

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      num_z    ! Num. of z-dir. grid pts. in vertical range

    real, dimension(nx,ny,num_z), intent(in) :: &
      value  ! Value (variable or quantity)

    ! Return Variable
    real, dimension(num_z) :: &
      sum_value_xy_domain  ! Sum of variable or quantity over entire domain

    ! Local Variables
    real, dimension(num_z) :: &
      sum_value  ! Running sum of the variable or quantity

   
    ! Sum of the variable or quantity in x and y over the domain (over the
    ! subdomain when nsubdomains > 1) for each vertical level.
    sum_value = sum_value_xy( nx, ny, num_z, value )

    ! When the number of subdomains (nsubdomains) > 1, find the overall sum of
    ! the variable or quantity in x and y over the entire domain by summing all
    ! the partial sums from each subdomain.
    if ( dompi ) then
       call task_sum_real( sum_value, sum_value_xy_domain, num_z )
    else
       sum_value_xy_domain = sum_value
    endif ! dompi


    return

  end function sum_value_xy_domain

  !=============================================================================
  function sum_value_xy( num_x, num_y, num_z, value )

    ! Description:
    ! Calculates the sum of a value (variable or quantity) over the specified
    ! horizontal model domain for each vertical level in the range.

    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_x, & ! Num. of x-dir. grid pts. in specified horizontal model domain
      num_y, & ! Num. of y-dir. grid pts. in specified horizontal model domain
      num_z    ! Num. of z-dir. grid pts. in vertical range

    real, dimension(num_x,num_y,num_z), intent(in) :: &
      value  ! Value (variable or quantity)

    ! Return Variable
    real, dimension(num_z) :: &
      sum_value_xy  ! Sum of variable or quantity over specified domain

    ! Local Variables
    integer :: i, j, k  ! Loop indices


    ! Sum of the variable or quantity in x and y over the specified horizontal
    ! model domain for each vertical level.
    do k = 1, num_z, 1

       sum_value_xy(k) = 0.0

       do i = 1, num_x, 1
          do j = 1, num_y, 1

             sum_value_xy(k) = sum_value_xy(k) + value(i,j,k)

          enddo ! j = 1, num_y, 1
       enddo ! i = 1, num_x, 1

    enddo ! k = 1, num_z, 1


    return

  end function sum_value_xy

  !=============================================================================

end module compute_correlation_module
#endif /*UWM_STATS*/
