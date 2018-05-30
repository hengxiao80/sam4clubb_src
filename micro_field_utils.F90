module micro_field_utils
#ifdef SILHS

 use grid, only: &
   nx, ny, nz, nzm ! Constant(s)

  use vars, only: &
    t, & ! Variable(s)
    w   

  use microphysics, only: &
    micro_field, & ! Variable(s)
    flag_wmass, &
    nmicro_fields
  
  implicit none

  private

  public :: micro_save, micro_restore, micro_copy_in, accumulate_LH_tendency, &
    clip_micro_field

  real, dimension(nx,ny,nzm), save :: &
    conc_save, & ! Saved value of cloud droplet number concentration
    t_cached, &
    qn_cached

  real, allocatable, dimension(:,:,:,:), save :: &
    micro_field_cached

  contains

!--------------------------------------------------------------------------------------------------
  subroutine micro_save( t_prior, w_prior, qn_prior, micro_field_prior )

! Description:
!   A simple subroutine that saves the values of the fields that are altered by
!   the Latin Hypercube sampling code.
!
! References:
!   None
!--------------------------------------------------------------------------------------------------
    use microphysics, only: conc, qn

    implicit none

    real, intent(out), dimension(nx,ny,nzm) :: &
      t_prior, & ! Moist static energy      [K]
      qn_prior   ! Cloud water mixing ratio [kg/kg]

    real, intent(out), dimension(nx,ny,nz) :: &
      w_prior  ! Vertical velocity    [m/s]

    real, intent(out), dimension(nx,ny,nzm, nmicro_fields) :: &
      micro_field_prior ! Microphysics fields [units vary]

    integer :: i, j, k

    ! ---- Begin Code ----

    forall ( i=1:nx, j=1:ny, k=1:nzm )
      t_prior(i,j,k)   = t(i,j,k)

      micro_field_prior(i,j,k,:) = micro_field(i,j,k,:)
    end forall

    forall ( i=1:nx, j=1:ny, k=1:nz )
      w_prior(i,j,k) = w(i,j,k)
    end forall

    ! Note: Ideally we could generalize this for all microphysics schemes in a
    ! cleaner way, but here we save cloud droplet number concentration and cloud
    ! water mixing ratio in local arrays. -dschanen 24 Aug 2012
    conc_save = conc
    qn_prior   = qn
    
    return
  end subroutine micro_save

  subroutine micro_restore( t_prior, w_prior, qn_prior, micro_field_prior )

    use microphysics, only: conc, qn

    implicit none

    real, intent(in), dimension(nx,ny,nzm) :: &
      t_prior, & ! Moist static energy  [K]
      qn_prior   ! Cloud water mixing ratio [kg/kg]

    real, intent(in), dimension(nx,ny,nz) :: &
      w_prior   ! Vertical velocity    [m/s]

    real, intent(in), dimension(nx,ny,nzm, nmicro_fields) :: &
      micro_field_prior ! Microphysics fields [units vary]

    integer :: i, j, k

    ! ---- Begin Code ----
    forall ( i=1:nx, j=1:ny, k=1:nzm )
      t(i,j,k)   = t_prior(i,j,k)

      micro_field(i,j,k,:) = micro_field_prior(i,j,k,:)
    end forall

    forall ( i=1:nx, j=1:ny, k=1:nz )
      w(i,j,k) = w_prior(i,j,k)
    end forall

    conc = conc_save
    qn   = qn_prior

    return
  end subroutine micro_restore

!--------------------------------------------------------------------------------------------------
  subroutine micro_copy_in( LH_iter, LH_microphys_calls, LH_rt, LH_t, micro_field_mean, &
                            X_nl_all_levs )
! Description:
!   Copy the subcolumns generated with SILHS into the variables that SAM uses in
!   the microphysics code.  A similar subroutine would be required for the
!   radiation or a chemistry simulator.
!
! References:
!   None
!--------------------------------------------------------------------------------------------------
    use clubb_api_module, only: &
      d_variables, &
      iiPDF_chi, & ! Variables
      iiPDF_w, &
      iiPDF_rr, &
      iiPDF_rs, &
      iiPDF_ri, &
      iiPDF_rg,&
      iiPDF_Nr, &
      iiPDF_Ns, &
      iiPDF_Ni, &
      iiPDF_Ng, &
      iiPDF_Ncn, &

      dp, &
      core_rknd

    use microphysics, only: &
      qn, & ! Variables
      conc, &
      index_water_vapor, &
      mkname, &
      nmicro_fields

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      LH_iter, &         ! Current subcolumn we're indexing
      LH_microphys_calls ! Total subcolumns

    real(kind=core_rknd), dimension(nx,ny,nz,LH_microphys_calls), intent(in) :: &
      LH_rt ! Values of total water         [kg/kg]
    real(kind=core_rknd), dimension(nx,ny,nzm,LH_microphys_calls), intent(in) :: &
      LH_t  ! Values of moist static energy [K]

    real, dimension(nx,ny,nzm,nmicro_fields), intent(in) :: &
      micro_field_mean ! Mean value for the micro fields  [units vary]

    real(kind=dp), dimension(nx,ny,nz,LH_microphys_calls,d_variables), intent(in) :: &
      X_nl_all_levs ! Values of various sample fields     [units vary]

    ! Local Variables
    real(kind=dp), dimension(nx,ny,nzm) :: &
      rc_sample ! Value of liquid water

    integer :: i, j, k, ivar, indx

    ! Cloud Liquid
    rc_sample = max( X_nl_all_levs(:,:,2:nz,LH_iter,iiPDF_chi), 0._dp )

    ! In the case of SAM1MOM this variable will be overwritten with a different
    ! value
    qn = real( rc_sample )

    ! Cloud droplet concentration
    conc = real( X_nl_all_levs(:,:,2:nz,LH_iter,iiPDF_Ncn) )

    ! It appears that index_water_vapor is actually the index of total water for
    ! all 3 microphysics schemes -dschanen 4 Sep 2012
    micro_field(1:nx,1:ny,1:nzm,index_water_vapor) = LH_rt(:,:,2:nz,LH_iter)

    ! Vertical velocity
    ! TODO interpolate to zi levels
    w(1:nx,1:ny,1:nz) = X_nl_all_levs(:,:,1:nz,LH_iter,iiPDF_w)

    ! Temperature
    t(1:nx,1:ny,1:nzm) = LH_t(:,:,1:nzm,LH_iter)

    do indx = 1, nmicro_fields, 1
      select case ( trim( mkname(indx) ) )

      ! Rain  
      case ( 'QR', 'QP' )
        micro_field(1:nx,1:ny,1:nzm,indx) = real( X_nl_all_levs(:,:,2:nz,LH_iter,iiPDF_rr) )
      case ( 'NR' )
        micro_field(1:nx,1:ny,1:nzm,indx) = real( X_nl_all_levs(:,:,2:nz,LH_iter,iiPDF_Nr) )

      ! Snow
      case ( 'QS' )
        micro_field(1:nx,1:ny,1:nzm,indx) = real( X_nl_all_levs(:,:,2:nz,LH_iter,iiPDF_rs) )
      case ( 'NS' )
        micro_field(1:nx,1:ny,1:nzm,indx) = real( X_nl_all_levs(:,:,2:nz,LH_iter,iiPDF_Ns) )

      ! Ice
      case ( 'QI' )
        micro_field(1:nx,1:ny,1:nzm,indx) = real( X_nl_all_levs(:,:,2:nz,LH_iter,iiPDF_ri) )
      case ( 'NI' )
        micro_field(1:nx,1:ny,1:nzm,indx) = real( X_nl_all_levs(:,:,2:nz,LH_iter,iiPDF_Ni) )

      ! Graupel (disabled; uses the mean value from micro_field)
      case ( 'QG' )
        micro_field(1:nx,1:ny,1:nzm,indx) = micro_field_mean(1:nx,1:ny,1:nzm,indx) 
      case ( 'NG' )
        micro_field(1:nx,1:ny,1:nzm,indx) = micro_field_mean(1:nx,1:ny,1:nzm,indx) 

      end select
    end do

    if ( .not. allocated( micro_field_cached ) ) then
      allocate( micro_field_cached(nx,ny,nzm,nmicro_fields) )
    end if

    forall( i=1:nx, j=1:ny, k=1:nzm, ivar=1:nmicro_fields )
      micro_field_cached(i,j,k,ivar) = micro_field(i,j,k,ivar)
    end forall

    forall( i=1:nx, j=1:ny, k=1:nzm )
      qn_cached(i,j,k) = qn(i,j,k)
      t_cached(i,j,k)  = t(i,j,k)
    end forall

    return
  end subroutine micro_copy_in

!--------------------------------------------------------------------------------------------------
  subroutine accumulate_LH_tendency( dtn, t, qn, micro_field, &
                                     LH_sample_point_weight, &
                                     LH_t_sum_tndcy, LH_qn_sum_tndcy, &
                                     LH_micro_field_sum_tndcy )
! Description:
!   This subroutine will increment the _tndcy arrays by a weighted tendency
!   that has been computed within the microphysics scheme.
!
! References:
!   None
!--------------------------------------------------------------------------------------------------
    use clubb_api_module, only: &
      core_rknd ! Constant

    use grid, only: &
      dimx1_s, dimx2_s, dimy1_s, dimy2_s ! Constants

    implicit none

    ! Input Variables
    real, intent(in) :: &
      dtn ! SAM dynamical timestep      [s]

    real, intent(in), dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) :: &
      t  ! Moist static energy  [K]

    real, intent(in), dimension(nx,ny,nzm) :: &
      qn ! Cloud water mixing ratio     [kg/kg]

    real, intent(in), dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields) :: &
      micro_field ! Microphysics fields [units vary]

    real(kind=core_rknd), dimension(nx,ny), intent(in) :: &
      LH_sample_point_weight ! Weight for the subcolumn

    real(kind=core_rknd), dimension(nx,ny,nzm), intent(inout) :: &
      LH_t_sum_tndcy, & ! Sum of all t tendencies          [K/s]
      LH_qn_sum_tndcy   ! Sum of all qn tendencies         [kg/kg/s]

    real(kind=core_rknd), dimension(nx,ny,nzm,nmicro_fields), intent(inout) :: &
      LH_micro_field_sum_tndcy  ! Sum of all micro_field tendencies     [units vary/s]

    integer :: k, ivar

    ! ---- Begin Code ----
    do k = 1, nzm
      ! Moist static energy
      LH_t_sum_tndcy(:,:,k) = LH_t_sum_tndcy(:,:,k) &
         + LH_sample_point_weight * ( t(1:nx,1:ny,k) - t_cached(1:nx,1:ny,k) ) / dtn

      ! Liquid water
      LH_qn_sum_tndcy(:,:,k) = LH_qn_sum_tndcy(:,:,k) &
         + LH_sample_point_weight * ( qn(:,:,k) - qn_cached(:,:,k) ) / dtn

      ! All other fields
      do ivar = 1, nmicro_fields
        LH_micro_field_sum_tndcy(:,:,k,ivar)  = LH_micro_field_sum_tndcy(:,:,k,ivar) &
          + LH_sample_point_weight * ( micro_field(1:nx,1:ny,k,ivar) &
                                       - micro_field_cached(:,:,k,ivar) ) / dtn
      end do
    end do ! k = 1 .. nzm

    return
  end subroutine accumulate_LH_tendency

!-------------------------------------------------------------------------------
  subroutine clip_micro_field

! Description:
!   Clip all microphysics variables to prevent negative mixing ratios and number
!   concentrations.
!
! References:
!   None
!-------------------------------------------------------------------------------
    use microphysics, only: qn, index_water_vapor ! Variable(s)

    use clubb_api_module, only: &
      core_rknd, & ! Constant
      fill_holes_vertical_api, & ! Procedure
      zero_threshold, rt_tol  ! Constant(s)

    use clubbvars, only: &
      rho_ds_zt, rho_ds_zm  ! Variable(s)

    implicit none

    ! External
    intrinsic :: real

    ! Local Variables
    real(kind=core_rknd), dimension(nz) :: clip_field

    integer :: i, j, ivar

    ! ---- Begin Code ----

    do i = 1, nx
      do j = 1, ny

        do ivar = 1, nmicro_fields

          ! Set the field to be clipped
          clip_field(2:nz) = real( micro_field(i,j,1:nzm,ivar), kind=core_rknd )
          clip_field(1) = 0._core_rknd ! Set ghost point to zero (this point is not referenced)

          ! Use a conservative clipping scheme for the mixing ratio variables
          if ( flag_wmass(ivar) == 1 ) then
            if ( ivar == index_water_vapor ) then
              call fill_holes_vertical_api( 2, rt_tol, "zt", rho_ds_zt, rho_ds_zm, clip_field )
            else
              call fill_holes_vertical_api( 2, zero_threshold, "zt", rho_ds_zt, rho_ds_zm, clip_field )
            end if

            ! Apply hard clipping if the above doesn't work
            where ( clip_field < zero_threshold ) clip_field = zero_threshold

          else ! Number concentrations are not conserved so hard clipping is applied

            where ( clip_field < zero_threshold )
              clip_field = zero_threshold
            end where

          end if ! flag_wmass(ivar) == 1 (a zero indicates the variable is not a mixing ratio)

          micro_field(i,j,1:nzm,ivar) = real( clip_field(2:nz) )

        end do ! 1 .. nmicro_fields

        ! Clip cloud water mixing ratio (this doesn't occur in micro_field)
        clip_field(2:nz) = real( qn(i,j,1:nzm), kind=core_rknd )
        call fill_holes_vertical_api( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, clip_field )
        ! Hard clipping if the above doesn't work.
        where ( clip_field < zero_threshold ) clip_field = zero_threshold
        qn(i,j,1:nzm) = real( clip_field(2:nz) )

      end do ! 1 .. ny
    end do ! 1 .. nx
    
    return
  end subroutine clip_micro_field
!-------------------------------------------------------------------------------

#endif /* SILHS */
end module micro_field_utils
