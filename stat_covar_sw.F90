#ifdef UWM_STATS
! w and scalar covariance budget:
subroutine stat_covar_sw( s_curr, s_prev, w_curr, w_prev, swle )

  use grid

  implicit none

  ! Input variables
  real, dimension(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm), intent(in) :: &
    s_curr, & ! Current value of scalar variable                  [scalar units]
    s_prev    ! Previous value of scalar variable (before update) [scalar units]

  real, dimension(nx,ny,nz), intent(in) :: &
    w_curr, & ! Current vertical velocity                         [m/s]
    w_prev    ! Previous vertical velocity                        [m/s]

  ! Output Variables
  real, dimension(nzm), intent(out) :: &
    swle    ! Contribution of change in variable to s'w'  [(scalar units) m/s^2]

  ! Local Variables
  real, dimension(nx,ny,nzm) :: &
    w_curr_t, & ! Current vertical velocity (interp.)   [m/s]
    w_prev_t    ! Previous vertical velocity (interp.)  [m/s]

  real, dimension(nzm) :: &
    s_curr0, & ! Current average of scalar variable   [scalar units]
    s_prev0    ! Previous average of scalar variable  [scalar units]

  real, dimension(nzm) :: &
    w_curr0, & ! Current average of vertical velocity (interp.)   [m/s]
    w_prev0    ! Previous average of vertical velocity (interp.)  [m/s]

  real :: coef ! Inverse time step duration                       [1/s]

  integer :: i, j, k  ! Loop indices


  ! Current mean of the scalar variable.
  call averageXY_MPI( s_curr, dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, s_curr0 )

  ! Previous mean of the scalar variable.
  call averageXY_MPI( s_prev, dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, s_prev0 )

  ! Interpolate current and previous values of vertical velocity to the interior
  ! (thermodynamic) grid.
  do i = 1, nx, 1
     do j = 1, ny, 1
        do k = 1, nzm-1, 1
           w_curr_t(i,j,k) = 0.5 * ( w_curr(i,j,k) + w_curr(i,j,k+1) )
           w_prev_t(i,j,k) = 0.5 * ( w_prev(i,j,k) + w_prev(i,j,k+1) )
        enddo ! k = 1, nzm-1, 1
        ! At level k = nzm, w at the top level (nz) is set to be 0, following
        ! other statistics subroutines.
        w_curr_t(i,j,nzm) = 0.5 * w_curr(i,j,nzm)
        w_prev_t(i,j,nzm) = 0.5 * w_prev(i,j,nzm)
     enddo
  enddo

  ! Current mean of vertical velocity (interpolated to thermodynamic levels).
  call averageXY_MPI( w_curr_t, 1, nx, 1, ny, nzm, w_curr0 )

  ! Previous mean of vertical velocity (interpolated to thermodynamic levels).
  call averageXY_MPI( w_prev_t, 1, nx, 1, ny, nzm, w_prev0 )

  coef = 1.0/dtn

  do k = 1, nzm, 1

     swle(k) = 0.0

     do j = 1, ny, 1
        do i = 1, nx, 1

           swle(k) = swle(k) &
                     + ( s_curr(i,j,k) - s_curr0(k) ) &
                       * ( w_curr_t(i,j,k) - w_curr0(k) ) &
                     - ( s_prev(i,j,k) - s_prev0(k) ) &
                       * ( w_prev_t(i,j,k) - w_prev0(k) )

        enddo
     enddo

     swle(k) = swle(k) * coef

  enddo


  return

end subroutine stat_covar_sw

#endif /*UWM_STATS*/
