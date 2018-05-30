#ifdef SILHS
! $Id: clubb_silhs_vars.F90 1554 2014-08-22 22:10:27Z dschanen@uwm.edu $
module clubb_silhs_vars

  use grid, only: &
    nx, &
    ny,&
    nz,&
    nzm,&
    dimx1_s,&
    dimx2_s,&
    dimy1_s,&
    dimy2_s

  use clubb_api_module, only: &
    core_rknd, & ! CLUBB core real kind
    dp

  implicit none

  private ! Default scope

  ! These are various parameters that control how SILHS is called.
  integer, parameter, public :: &
    lh_microphys_interactive      = 1, &
    lh_microphys_non_interactive  = 2, &
    lh_microphys_disabled         = 3

  integer, public :: lh_microphys_calls, lh_sequence_length, lh_microphys_type

  ! Allocatable variables that can change in dimension at runtime
  real(kind=core_rknd), public, allocatable, dimension(:,:,:,:) :: &
    LH_rt, & ! Latin hypercube samples of total water         [kg/kg]
    LH_t     ! Latin hypercube samples of moist static energy [K]

  real(kind=dp), public, allocatable, dimension(:,:,:,:,:) :: &
    X_nl_all_levs ! Lognormally distributed hydrometeors [units vary]

  integer, public, allocatable, dimension(:,:,:,:) :: &
    X_mixt_comp_all_levs ! Which mixture component the sample is in

  real(kind=core_rknd), public, allocatable, dimension(:,:,:) :: &
    LH_sample_point_weights ! Weights for cloud weighted sampling

  ! Static variables
  real(kind=core_rknd), public, dimension(nx,ny,nzm) :: &
    LH_t_sum_tndcy,     & ! Sum of all t LH tendencies      [K/s]
    LH_t_avg_tndcy,     & ! Average of all t LH tendencies  [K/s]
    LH_qn_sum_tndcy,    & ! Sum of all qn  LH tendencies    [kg/kg/s]
    LH_qn_avg_tndcy       ! Average of all qn LH tendencies [kg/kg/s]

  real, public, dimension(nx,ny,nzm) :: &
    t_prior,  & ! Saved value of t                 [K]
    qn_prior    ! Saved value of liquid water      [kg/kg]

  real, public, dimension(nx,ny,nz) :: &
    w_prior  ! Saved value of w  [m/s]

  real, public, allocatable, dimension(:,:,:,:) :: &
    micro_field_prior ! Saved values of the micro_fields     [units vary]

  real(kind=core_rknd), public, allocatable, dimension(:,:,:,:) :: &
    LH_micro_field_sum_tndcy, & ! Sum of all micro_field tendencies     [units vary/s]
    LH_micro_field_avg_tndcy    ! Average of all micro_field tendencies [units vary/s]

end module clubb_silhs_vars
#endif /* SILHS */
