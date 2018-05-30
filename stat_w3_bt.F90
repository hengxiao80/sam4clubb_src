#ifdef PNNL_STATS
!Heng Xiao w'w'w' budget
subroutine stat_w3_bt(w_curr, w_prev, w3le )

  use grid

  implicit none

  ! Input variables
  real, dimension(nx,ny,nz), intent(in) :: &
    w_curr, & ! Current vertical velocity                         [m/s]
    w_prev    ! Previous vertical velocity                        [m/s]

  ! Output Variables
  real, dimension(nzm), intent(out) :: &
    w3le    ! Contribution of change in variable to w'w'w'  [m^3/s^4]

  ! Local Variables
  real, dimension(nx,ny,nzm) :: &
    w_curr_t, & ! Current vertical velocity (interp.)   [m/s]
    w_prev_t    ! Previous vertical velocity (interp.)  [m/s]

  real, dimension(nzm) :: &
    w_curr0, & ! Current average of vertical velocity (interp.)   [m/s]
    w_prev0    ! Previous average of vertical velocity (interp.)  [m/s]

  real :: coef ! Inverse time step duration                       [1/s]

  integer :: i, j, k  ! Loop indices

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

  coef = 1.0/(dtn*nx*ny)

  do k = 1, nzm, 1

     w3le(k) = 0.0

     do j = 1, ny, 1
        do i = 1, nx, 1

           w3le(k) = w3le(k) &
                       + ( w_curr_t(i,j,k) - w_curr0(k) )**3 &
                       - ( w_prev_t(i,j,k) - w_prev0(k) )**3

        enddo
     enddo

     w3le(k) = w3le(k) * coef

  enddo

  return

end subroutine stat_w3_bt

#endif /*PNNL_STATS*/
