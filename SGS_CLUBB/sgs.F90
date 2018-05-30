module sgs

! module for original SAM subgrid-scale SGS closure (Smagorinsky or 1st-order TKE)
! Modified to use the UW Milwaukee CLUBB scheme.
! Marat Khairoutdinov, 2012

  use grid, only: nx,nxp1,ny,nyp1,YES3D,nzm,nz,dimx1_s,dimx2_s,dimy1_s,dimy2_s
  use params, only: dosgs
  use sgs_params, only: &
    dosmagor, &
    doclubb, &
    doclubbnoninter, &
    doclubb_sfc_fluxes, &
    nclubb
#ifdef CLUBB
  use clubb_api_module, only: core_rknd
#endif
  implicit none

!----------------------------------------------------------------------
! Required definitions:

!!! prognostic scalar (need to be advected arround the grid):

  integer, parameter :: nsgs_fields = 1   ! total number of prognostic sgs vars

  real sgs_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nsgs_fields)

!!! sgs diagnostic variables that need to exchange boundary information (via MPI):

  integer, parameter :: nsgs_fields_diag = 2   ! total number of diagnostic sgs vars

! diagnostic fields' boundaries:
  integer, parameter :: dimx1_d=0, dimx2_d=nxp1, dimy1_d=1-YES3D, dimy2_d=nyp1

  real sgs_field_diag(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm, nsgs_fields_diag)

  logical:: advect_sgs = .false. ! advect prognostics or not, default - not (Smagorinsky)
  logical, parameter:: do_sgsdiag_bound = .true.  ! exchange boundaries for diagnostics fields

! SGS fields that output by default (if =1).
  integer, parameter :: flag_sgs3Dout(nsgs_fields) = (/0/)
  integer, parameter :: flag_sgsdiag3Dout(nsgs_fields_diag) = (/0,0/)

  real fluxbsgs (nx, ny, 1:nsgs_fields) ! surface fluxes
  real fluxtsgs (nx, ny, 1:nsgs_fields) ! top boundary fluxes

!!! these arrays may be needed for output statistics:

  real sgswle(nz,1:nsgs_fields)  ! resolved vertical flux
  real sgswsb(nz,1:nsgs_fields)  ! SGS vertical flux
  real sgsadv(nz,1:nsgs_fields)  ! tendency due to vertical advection
  real sgslsadv(nz,1:nsgs_fields)  ! tendency due to large-scale vertical advection
  real sgsdiff(nz,1:nsgs_fields)  ! tendency due to vertical diffusion

!------------------------------------------------------------------
! internal (optional) definitions:

! make aliases for prognostic variables:

  real tke(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! SGS TKE
  equivalence (tke(dimx1_s,dimy1_s,1),sgs_field(dimx1_s,dimy1_s,1,1))

! make aliases for diagnostic variables:

  real tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
  real tkh (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy conductivity
  equivalence (tk(dimx1_d,dimy1_d,1), sgs_field_diag(dimx1_d, dimy1_d,1,1))
  equivalence (tkh(dimx1_d,dimy1_d,1), sgs_field_diag(dimx1_d, dimy1_d,1,2))


  real grdf_x(nzm)! grid factor for eddy diffusion in x
  real grdf_y(nzm)! grid factor for eddy diffusion in y
  real grdf_z(nzm)! grid factor for eddy diffusion in z
! Local diagnostics:

  real tkesbbuoy(nz), tkesbshear(nz),tkesbdiss(nz), tkesbdiff(nz)

#ifdef CLUBB
! Added for the advance_clubb_core call
  real( kind = core_rknd) :: thlp2_forcing(nzm)

  real(kind=core_rknd) :: &
    rtm_integral_before(nx,ny), rtm_integral_after(nx,ny), rtm_flux_top, rtm_flux_sfc
  real(kind=core_rknd) :: &
    thlm_integral_before(nx,ny), thlm_integral_after(nx,ny), thlm_before(nzm), thlm_after(nzm), &
    thlm_flux_top, thlm_flux_sfc

  real(kind=core_rknd), dimension(nzm) :: &
    rtm_column ! Total water (vapor + liquid)     [kg/kg]

  integer :: &
    stats_nsamp, & ! Sampling interval for a single column of CLUBB data [timestep]
    stats_nout     ! Output interval for a single column of CLUBB data [timestep]
#endif

  CONTAINS

! required microphysics subroutines and function:
!----------------------------------------------------------------------
!!! Read microphysics options from prm (namelist) file

  subroutine sgs_setparm()

    use grid, only: case

    implicit none

    integer :: ios, ios_missing_namelist, place_holder

    NAMELIST /SGS_CLUBB/ &
      dosmagor, &           ! Diagnostic Smagorinsky closure
      doclubb, &            ! Enable the interactive CLUBB code
      doclubbnoninter, &    ! Enable non-interactive CLUBB code
      doclubb_sfc_fluxes, & ! Use CLUBB to apply the surface fluxes
      nclubb                ! Number clubb calls per dynamical timestep

    NAMELIST /BNCUIODSBJCB/ place_holder

    dosmagor = .true. ! default
    doclubb = .false.
    doclubbnoninter = .false.
    doclubb_sfc_fluxes = .false.
    nclubb = 1

    place_holder = -1 ! To avoid a compiler warning -dschanen 23 Aug 14
    !----------------------------------
    !  Read namelist for microphysics options from prm file:
    !------------
    open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')

    read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
    rewind(55) !note that one must rewind before searching for new namelists

    read (55,SGS_CLUBB,IOSTAT=ios)

    advect_sgs = .not.dosmagor

    if (ios.ne.0) then
      !namelist error checking
      if(ios.ne.ios_missing_namelist) then
        write(*,*) '****** ERROR: bad specification in SGS_CLUBB namelist'
        call task_abort()
      end if
    end if
    close(55)

#ifndef CLUBB
    !--------------------------
    ! do a CLUBB sanity check
    if ( doclubb .or. doclubbnoninter ) then
      write(0,*) "Cannot call CLUBB if -DCLUBB is not in FFLAGS"
      call task_abort()
    end if
#endif

    return
  end subroutine sgs_setparm

!----------------------------------------------------------------------
!!! Initialize sgs:


  subroutine sgs_init()

    use grid, only: nrestart, dx, dy, dz, adz, masterproc
    use params, only: LES
#ifdef CLUBB
    use grid, only: dt, z, zi
    use vars, only: rho, rhow, tv0
    use clubb_api_module, only: time_precision, stats_tsamp, stats_tout
    use clubb_sgs, only: clubb_sgs_setup
#endif

    implicit none
#ifdef CLUBB
    ! External
    intrinsic :: nint
#endif

    integer k

    if(nrestart.eq.0) then

      sgs_field = 0.
      sgs_field_diag = 0.

      fluxbsgs = 0.
      fluxtsgs = 0.

    end if

    if(masterproc) then
      if(dosmagor) then
        write(*,*) 'Smagorinsky SGS Closure'
      else
        write(*,*) 'Prognostic TKE 1.5-order SGS Closure'
      end if
      if ( doclubb ) then
        write(*,*) 'CLUBB Parameterization'
      end if
    end if

    if(LES) then
      do k=1,nzm
        grdf_x(k) = dx**2/(adz(k)*dz)**2
        grdf_y(k) = dy**2/(adz(k)*dz)**2
        grdf_z(k) = 1.
      end do
    else
      do k=1,nzm
        grdf_x(k) = min(16.,dx**2/(adz(k)*dz)**2)
        grdf_y(k) = min(16.,dy**2/(adz(k)*dz)**2)
        grdf_z(k) = 1.
      end do
    end if

    sgswle = 0.
    sgswsb = 0.
    sgsadv = 0.
    sgsdiff = 0.
    sgslsadv = 0.

#ifdef CLUBB
    !------------------------------------------------------------------
    ! Do initialization for UWM CLUBB
    !------------------------------------------------------------------

    if ( doclubb .or. doclubbnoninter ) then
      call clubb_sgs_setup( real( dt*real( nclubb ), kind=time_precision), &
                            z, rho, zi, rhow, tv0, tke )
    end if

    thlp2_forcing(:) = 0.0_core_rknd

    stats_nsamp = nint( stats_tsamp / real( dt, kind=time_precision ) )
    stats_nout  = nint( stats_tout / real( dt, kind=time_precision ) )

#endif
    return
  end subroutine sgs_init

!----------------------------------------------------------------------
!!! make some initial noise in sgs:
!
  subroutine setperturb_sgs(ptype)

    use vars, only: q0, z

    implicit none

    integer, intent(in) :: ptype
    integer i,j,k

    select case (ptype)

    case(0)

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(k.le.4.and..not.dosmagor) then
              tke(i,j,k)=0.04*(5.-real(k))
            endif
          end do
        end do
      end do

    case(1)

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(q0(k).gt.6.e-3.and..not.dosmagor) then
              tke(i,j,k)=1.
            endif
          end do
        end do
      end do

    case(2)

    case(3)   ! gcss wg1 smoke-cloud case

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(q0(k).gt.0.5e-3.and..not.dosmagor) then
              tke(i,j,k)=1.
            endif
          end do
        end do
      end do


    case(4)  ! gcss wg1 arm case

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(z(k).le.150..and..not.dosmagor) then
              tke(i,j,k)=0.15*(1.-z(k)/150.)
            endif
          end do
        end do
      end do


    case(5)  ! gcss wg1 BOMEX case

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(z(k).le.3000..and..not.dosmagor) then
              tke(i,j,k)=1.-z(k)/3000.
            endif
          end do
        end do
      end do

    case(6)  ! GCSS Lagragngian ASTEX


      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(q0(k).gt.6.e-3.and..not.dosmagor) then
              tke(i,j,k)=1.
            endif
          end do
        end do
      end do


    case default

    end select
    return
  end subroutine setperturb_sgs

!----------------------------------------------------------------------
!!! Estimate Courant number limit for SGS
!

  subroutine kurant_sgs(cfl)

    use grid, only: dt, dx, dy, dz, adz, adzw
    implicit none

    real, intent(out) :: cfl

    integer k
    real tkhmax(nz)

    do k = 1,nzm
      tkhmax(k) = maxval(tkh(1:nx,1:ny,k))
    end do

    cfl = 0.
    do k=1,nzm
      cfl = max(cfl,        &
         0.5*tkhmax(k)*grdf_z(k)*dt/(dz*adzw(k))**2, &
         0.5*tkhmax(k)*grdf_x(k)*dt/dx**2, &
         real(YES3D)*0.5*tkhmax(k)*grdf_y(k)*dt/dy**2)
    end do

    return
  end subroutine kurant_sgs


!----------------------------------------------------------------------
!!! compute sgs diffusion of momentum:
!
  subroutine sgs_mom()

#ifdef CLUBB
    use vars, only: dudt, dvdt ! Variable(s)
    use clubb_sgs, only: apply_clubb_sgs_tndcy_mom ! Procedure(s)
#endif /* CLUBB */

    implicit none

    call diffuse_mom()

#ifdef CLUBB
    if ( doclubb ) then
      call apply_clubb_sgs_tndcy_mom( dudt, dvdt ) ! in/out
    end if
#endif

  end subroutine sgs_mom

!----------------------------------------------------------------------
!!! compute sgs diffusion of scalars:
!
  subroutine sgs_scalars()

    use vars
    use microphysics
    use tracers
    use params, only: dotracers !, doclubb, doclubb_sfc_fluxes
#ifdef PNNL_STATS
    use calc_vars_util, only: t2thetal
#endif
#ifdef CLUBB
    use clubb_sgs, only: advance_clubb_sgs, clubb_sgs_setup, clubb_sgs_cleanup, &
      apply_clubb_sgs_tndcy_scalars, &
      apply_clubb_sgs_tndcy_mom  ! Subroutines
    use calc_vars_util, only: &
      t2thetal    ! Functions
    use clubbvars, only: edsclr_dim, sclr_dim, wprtp, wpthlp, rho_ds_zt, rho_ds_zm, &
      rtm_spurious_source, thlm_spurious_source
    use clubb_api_module, only: &
      time_precision, &
      vertical_integral_api, & ! Function
      calculate_spurious_source_api, &
      gr, &   ! Variable(s)
      stats_tsamp, &
      stats_tout
    use sgs_params, only: doclubb, doclubbnoninter, nclubb, doclubb_sfc_fluxes
#endif  /*CLUBB*/
    implicit none

    real dummy(nz)
    real f2lediff_xy(nz), f2lediss_xy(nz), fwlediff_xy(nz)
    real f2lediff_z(nz), f2lediss_z(nz), fwlediff_z(nz)
    real fluxbtmp(nx,ny), fluxttmp(nx,ny) !bloss
    integer i, j, k
#ifdef PNNL_STATS
    !MWSWong:convert HL to THL
    real thelprev(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real thelcurr(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real f0(nzm), df0(nzm)
    real r2dx,r2dy,r2dx0,r2dy0,r2dz
    integer kb,kc,jb,jc,n
    real factor_xy
    real q1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real qtogprev(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real qtogcurr(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real f0_qtog(nzm), df0_qtog(nzm)
    real f0_q0(nzm), df0_q0(nzm)
#endif /*PNNL_STATS*/


#ifdef PNNL_STATS
    if(dostatis) then
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            thelprev(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                        qci(i,j,k), qpi(i,j,k), prespot(k) )
            qtogprev(i,j,k) = 0.0
            do n = 1, nmicro_fields ! prevents accessing unreferenced memory
              qtogprev(i,j,k) = qtogprev(i,j,k) + flag_wmass(n)*micro_field(i,j,k,n)
            end do
            q1(i,j,k) = micro_field(i,j,k,index_water_vapor)
          end do
        end do
      end do
    end if
#endif /*PNNL_STATS*/

!      call diffuse_scalar(t,fluxbt,fluxtt,tdiff,twsb, &
!                           t2lediff,t2lediss,twlediff,.true.)
    f2lediff_xy = 0.0
    f2lediss_xy = 0.0
    fwlediff_xy = 0.0

    call diffuse_scalar_xy(t,fluxbt,fluxtt,tdiff,twsb, &
                         f2lediff_xy,f2lediss_xy,fwlediff_xy,.true.)
    f2lediff_z =0.0
    f2lediss_z =0.0
    fwlediff_z =0.0
#ifdef CLUBB
    ! Diffuse moist static energy in the vertical only if CLUBB is not being
    ! called
    if ( .not. doclubb ) then
      call diffuse_scalar_z(t,fluxbt,fluxtt,tdiff,twsb, &
                          f2lediff_z,f2lediss_z,fwlediff_z,.true.)
    else ! doclubb
      if ( doclubb_sfc_fluxes ) then
        ! The flux will be applied in advance_clubb_core, so the 2nd argument
        ! is zero.
        call fluxes_scalar_z(t,fzero,fluxtt,tdiff,twsb, & 
                             f2lediff_z,f2lediss_z,fwlediff_z,.true.)
      else
        call fluxes_scalar_z(t,fluxbt,fluxtt,tdiff,twsb, & 
                             f2lediff_z,f2lediss_z,fwlediff_z,.true.)
      end if
    end if
#else
    call diffuse_scalar_z(t,fluxbt,fluxtt,tdiff,twsb, &
                         f2lediff_z,f2lediss_z,fwlediff_z,.true.)
#endif

    t2lediff = f2lediff_xy + f2lediff_z
    t2lediss = f2lediss_xy + f2lediss_z
    twlediff = fwlediff_xy + fwlediff_z


    if(advect_sgs) then
!         call diffuse_scalar(tke,fzero,fzero,dummy,sgswsb, &
!                                    dummy,dummy,dummy,.false.)
      call diffuse_scalar_xy(tke,fzero,fzero,dummy,sgswsb, &
                                 dummy,dummy,dummy,.false.)
      call diffuse_scalar_z(tke,fzero,fzero,dummy,sgswsb, &
                                 dummy,dummy,dummy,.false.)
    end if


!
!    diffusion of microphysics prognostics:
!
    call micro_flux()

    total_water_evap = total_water_evap - total_water()

    do k = 1,nmicro_fields
      if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
       .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
       .or. doprecip.and.flag_precip(k).eq.1 ) then

        fluxbtmp(1:nx,1:ny) = fluxbmk(1:nx,1:ny,k)
        fluxttmp(1:nx,1:ny) = fluxtmk(1:nx,1:ny,k)
!           call diffuse_scalar(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
!                mkdiff(:,k),mkwsb(:,k), dummy,dummy,dummy,.false.)

        !--MWSWong: want to output q2 dissip/diftr terms
        f2lediff_xy = 0.0
        f2lediss_xy = 0.0
        fwlediff_xy = 0.0

        if(k.ne.index_water_vapor) then
          call diffuse_scalar_xy(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
             mkdiff(:,k),mkwsb(:,k), dummy,dummy,dummy,.false.)
        else !k==index_water_vapor
          call diffuse_scalar_xy(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
            mkdiff(:,k),mkwsb(:,k), f2lediff_xy,f2lediss_xy,fwlediff_xy,.true.)
        end if

        !-MWSWong: want to output q2 dissip/diftr terms
        f2lediff_z =0.0
        f2lediss_z =0.0
        fwlediff_z =0.0

        if(k.ne.index_water_vapor) then
          call diffuse_scalar_z(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
             mkdiff(:,k),mkwsb(:,k),dummy,dummy,dummy,.false.)
        else ! k==index_water_vapor
          if(.not. doclubb) then
            call diffuse_scalar_z(micro_field(:,:,:,k),fluxbtmp,fluxttmp, & 
                  mkdiff(:,k),mkwsb(:,k), f2lediff_z,f2lediss_z,fwlediff_z,.true.)
          else   ! doclubb
            call fluxes_scalar_z(micro_field(:,:,:,k),fluxbtmp,fluxttmp, & 
                  mkdiff(:,k),mkwsb(:,k), f2lediff_z,f2lediss_z,fwlediff_z,.true.)
          end if

          q2lediff = f2lediff_xy + f2lediff_z
          q2lediss = f2lediss_xy + f2lediss_z
          qwlediff = fwlediff_xy + fwlediff_z

        end if

      end if
    end do

    total_water_evap = total_water_evap + total_water()

#ifdef UWM_MISC
    ! Update microphysical variables qv, qcl, qpl, qci, and qpi after the
    ! micro_field array has been changed due to diffusion.
    call micro_diagnose()
#endif /*UWM_MISC*/

    ! diffusion of tracers:

    if(dotracers) then

      call tracers_flux()

      do k = 1,ntracers

        fluxbtmp = fluxbtr(:,:,k)
        fluxttmp = fluxttr(:,:,k)
!          call diffuse_scalar(tracer(:,:,:,k),fluxbtmp,fluxttmp, &
!               trdiff(:,k),trwsb(:,k), &
!               dummy,dummy,dummy,.false.)
        call diffuse_scalar_xy(tracer(:,:,:,k),fluxbtmp,fluxttmp, &
             trdiff(:,k),trwsb(:,k), &
             dummy,dummy,dummy,.false.)

#ifdef CLUBB
        ! Only diffuse the tracers if CLUBB is either disabled or using the
        ! eddy scalars code to diffuse them.
        if ( .not. doclubb .or. ( doclubb .and. edsclr_dim < 1 .and. sclr_dim < 1 ) ) then
          call diffuse_scalar_z(tracer(:,:,:,k),fluxbtmp,fluxttmp, &
               trdiff(:,k),trwsb(:,k), &
               dummy,dummy,dummy,.false.)
        end if
#else
        call diffuse_scalar_z(tracer(:,:,:,k),fluxbtmp,fluxttmp, &
             trdiff(:,k),trwsb(:,k), &
             dummy,dummy,dummy,.false.)
#endif
!!$          call diffuse_scalar(tracer(:,:,:,k),fluxbtr(:,:,k),fluxttr(:,:,k),trdiff(:,k),trwsb(:,k), &
!!$                           dummy,dummy,dummy,.false.)

      end do

    end if



#ifdef PNNL_STATS
!(following MICRO_M2005/microphysics.F90)
    ! MWSWong:THL and QTOG budget terms
    if(dostatis) then
      do k=1,nzm
        theldiff(k)=0.
        qtogdiff(k)=0.
        thelwlediff(k)=0.
        qtogwlediff(k)=0.
        do j=1,ny
          do i=1,nx
            thelcurr(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                        qci(i,j,k), qpi(i,j,k), prespot(k) )
            qtogcurr(i,j,k) = 0.0
            do n = 1, nmicro_fields ! prevents accessing unreferenced memory
              qtogcurr(i,j,k) = qtogcurr(i,j,k) + flag_wmass(n)*micro_field(i,j,k,n)
            end do
            theldiff(k) = theldiff(k)+( thelcurr(i,j,k)-thelprev(i,j,k))
            qtogdiff(k) = qtogdiff(k)+( qtogcurr(i,j,k)-qtogprev(i,j,k))
          end do
        end do
      end do


      call stat_varscalar(thelcurr,thelprev,f0,df0,thel2lediff)
      call setvalue(thelwlediff,nzm,0.)
      call stat_sw2(thelcurr,thelprev,thelwlediff)

      call stat_varscalar(qtogcurr,qtogprev,f0_qtog,df0_qtog,qtog2lediff)
      call setvalue(qtogwlediff,nzm,0.)
      call stat_sw2(qtogcurr,qtogprev,qtogwlediff)

      call averageXY_MPI(micro_field(:,:,:,index_water_vapor), &
                         dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,f0_q0)
      call averageXY_MPI(q1,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,df0_q0)

      do k=1,nzm
        qthellediff(k)=0.
        do j=1,ny
          do i=1,nx
            qthellediff(k) &
            = qthellediff(k) &
              + (thelcurr(i,j,k)-f0(k)) &
                * (micro_field(i,j,k,index_water_vapor)-f0_q0(k)) &
              - (thelprev(i,j,k)-df0(k)) * (q1(i,j,k)-df0_q0(k))
          end do
        end do
        qthellediff(k)=qthellediff(k)*(1./(dtn*nx*ny))
      end do


      factor_xy=1./float(nx*ny)
      r2dx0=1./(2.*dx)
      r2dy0=1./(2.*dy)
      do k=1,nzm
        thel2lediss(k)=0.
        qtog2lediss(k)=0.
        qthellediss(k)=0.
        kc=min(nzm,k+1)
        kb=max(1,k-1)
        r2dz=2./((kc-kb)*(adzw(k+1)+adzw(k))*dz)
        r2dx=r2dx0*sqrt((kc-kb)*dx*r2dz) ! grid anisotropy correction
        r2dy=r2dy0*sqrt((kc-kb)*dx*r2dz)
        do j=1,ny
          jc=j+YES3D
          jb=j-YES3D
          do i=1,nx
            thel2lediss(k)=thel2lediss(k)-tkh(i,j,k)*factor_xy*( &
                    ((thelcurr(i,j,kc)-f0(kc)-thelcurr(i,j,kb)+f0(kb))*r2dz)**2 ) !&
            !-tkh(i,j,k)*factor_xy*( &
            !    ((thelcurr(i+1,j,k)-thelcurr(i-1,j,k))*r2dx)**2+ &
            !    ((thelcurr(i,jc,k)-thelcurr(i,jb,k))*r2dy)**2 )
            qtog2lediss(k)=qtog2lediss(k)-tkh(i,j,k)*factor_xy*( &
                    ((qtogcurr(i,j,kc)-f0_qtog(kc)-qtogcurr(i,j,kb)+f0_qtog(kb))*r2dz)**2 ) !&
            !-tkh(i,j,k)*factor_xy*( &
            !   ((qtogcurr(i+1,j,k)-qtogcurr(i-1,j,k))*r2dx)**2+ &
            !   ((qtogcurr(i,jc,k)-qtogcurr(i,jb,k))*r2dy)**2 )
            qthellediss(k) &
            = qthellediss(k) &
              - tkh(i,j,k) &
                * ( factor_xy * ( thelcurr(i,j,kc)-f0(kc) - thelcurr(i,j,kb)+f0(kb)) &
                              * ( micro_field(i,j,kc,index_water_vapor)-f0_q0(kc) &
                                  - micro_field(i,j,kb,index_water_vapor)+f0_q0(kb)) *r2dz**2 ) !&
            !-tkh(i,j,k)*factor_xy*( &
            !   ((thelcurr(i+1,j,k)-thelcurr(i-1,j,k))*(micro_field(i+1,j,k,index_water_vapor)-micro_field(i-1,j,k,index_water_vapor)) )*r2dx**2+ &
            !   ((thelcurr(i,jc,k)-thelcurr(i,jb,k))*(micro_field(i,jc,k,index_water_vapor)-micro_field(i,jb,k,index_water_vapor)) )*r2dy**2 )
          end do
        end do
        !thel2lediss(k)=thel2lediss(k)*2.*factor_xy
        !qtog2lediss(k)=qtog2lediss(k)*2.*factor_xy
        !qthellediss(k)=qthellediss(k)*2.*factor_xy
        thel2lediss(k)=thel2lediss(k)*2.
        qtog2lediss(k)=qtog2lediss(k)*2.
        qthellediss(k)=qthellediss(k)*2.

        !write(6,*) 'SGS_CLUBB/sgs: thel2diss(k) t2lediss(k) ', thel2lediss(k), t2lediss(k), 'k', k, &
!					thelcurr(1,1,k), t(1,1,k)
!         write(6,*) 'SGS_CLUBB/sgs: qtog2diss(k) q2lediss(k) ', qtog2lediss(k), q2lediss(k), 'k', k, &
!					qtogcurr(1,1,k), micro_field(1,1,k,index_water_vapor)


      end do
    end if
#endif /*PNNL_STATS*/

#ifdef CLUBB
    if ( doclubb ) then

      ! Recalculate q, qv, qcl based on new micro_fields (updated by horizontal
      ! diffusion)
      call micro_update()

      ! Then Re-compute q/qv/qcl based on values computed in CLUBB
      call apply_clubb_sgs_tndcy_scalars &
           ( real( dtn, kind=time_precision), & ! in
             t, qv, qcl) ! in/out

      call micro_adjust( qv, qcl ) ! in

      ! Calculate the vertical integrals for RTM and THLM again so
      ! calculate whether CLUBB is a spurious source or sink of either.
      ! - nielsenb UWM 4 Jun 2010
      do i = 1,nx
        do j = 1,ny
          rtm_flux_top = rho_ds_zm(nz) * real( wprtp(i,j,nz), kind=core_rknd )
          rtm_flux_sfc = rho_ds_zm(1) * real( fluxbq(i,j), kind=core_rknd )
          rtm_column   = real( qv(i,j,1:nzm) + qcl(i,j,1:nzm), kind=core_rknd )
          rtm_integral_after(i,j) = vertical_integral_api( (nz - 2 + 1), rho_ds_zt(2:nz), &
                                        rtm_column, gr%invrs_dzt(2:nz) )

          rtm_spurious_source(i,j) = calculate_spurious_source_api( rtm_integral_after(i,j), &
                                                     rtm_integral_before(i,j), &
                                                     rtm_flux_top, rtm_flux_sfc, &
                                                     0.0_core_rknd, real( dtn, kind=core_rknd) )

          thlm_flux_top = rho_ds_zm(nz) * wpthlp(i,j,nz)
          thlm_flux_sfc = rho_ds_zm(1) * real( fluxbt(i,j), core_rknd )

          thlm_after = real( t2thetal( t(i,j,1:nzm), gamaz(1:nzm), &
                                       qpl(i,j,1:nzm), qci(i,j,1:nzm), &
                                       qpi(i,j,1:nzm), prespot(1:nzm) ), &
                             kind = core_rknd )

          thlm_integral_after(i,j) = vertical_integral_api( (nz - 2 + 1), rho_ds_zt(2:nz), &
                                                     thlm_after(1:nzm), gr%invrs_dzt(2:nz))

          thlm_spurious_source(i,j) = calculate_spurious_source_api( thlm_integral_after(i,j), &
                                                         thlm_integral_before(i,j), &
                                                         thlm_flux_top, thlm_flux_sfc, &
                                                         0.0_core_rknd, real( dtn, kind=core_rknd ))
        end do
      end do
      ! End spurious source calculation

    end if! doclubb
#endif


    return
  end subroutine sgs_scalars

!----------------------------------------------------------------------
!!! compute sgs processes (beyond advection):
!
  subroutine sgs_proc()

    use grid, only: nstep,dt,icycle
    use params, only: dosmoke

#ifdef CLUBB
    use clubb_sgs, only: advance_clubb_sgs, clubb_sgs_setup, clubb_sgs_cleanup, &
      apply_clubb_sgs_tndcy_scalars, &
      apply_clubb_sgs_tndcy_mom  ! Subroutines
    use microphysics, only: micro_update
    use calc_vars_util, only: &
      t2thetal    ! Functions
    use clubbvars, only: edsclr_dim, sclr_dim, wprtp, wpthlp, rho_ds_zt, rho_ds_zm, &
      rtm_spurious_source, thlm_spurious_source
    use clubb_api_module, only: &
      time_precision, &
      vertical_integral_api, & ! Function
      calculate_spurious_source_api, &
      gr, &   ! Variable(s)
      stats_tsamp, &
      stats_tout
    use sgs_params, only: doclubb, doclubbnoninter, nclubb, doclubb_sfc_fluxes
    use vars, only: qpl, qv, qcl, gamaz, prespot, t, qci, qpi, u, v, w, rho, rhow, wsub
    use grid, only: dtn
#endif  /*CLUBB*/

    implicit none
#ifdef CLUBB
    integer :: i, j
#endif
!    SGS TKE equation:

    if(dosgs) call tke_full()

#ifdef CLUBB
!----------------------------------------------------------
!     Do a timestep with CLUBB if enabled:

    if ( doclubb .or. doclubbnoninter ) then
      ! In case of ice fall, we recompute qci here for the
      ! single-moment scheme.  Also, subsidence, diffusion and advection have
      ! been applied to micro_field but not qv/qcl so they must be updated.
      call micro_update()
    end if ! doclubb .or. doclubbnoninter

    if ( doclubb ) then
      ! Calculate the vertical integrals for RTM and THLM so we can later
      ! calculate whether CLUBB is a spurious source or sink of either.
      ! - nielsenb UWM 4 Jun 2010
      do i = 1,nx
        do j = 1,ny
          rtm_column = real( qv(i,j,1:nzm) + qcl(i,j,1:nzm), kind=core_rknd )
          rtm_integral_before(i,j) = vertical_integral_api( (nz - 2 + 1), rho_ds_zt(2:nz), &
                                       rtm_column, gr%invrs_dzt(2:nz) )

          thlm_before = real( t2thetal( t(i,j,1:nzm), gamaz(1:nzm), &
                                        qpl(i,j,1:nzm), qci(i,j,1:nzm), &
                                        qpi(i,j,1:nzm), prespot(1:nzm) ), &
                              kind = core_rknd )

          thlm_integral_before(i,j) = vertical_integral_api( (nz - 2 + 1), rho_ds_zt(2:nz), &
                                                        thlm_before(1:nzm), gr%invrs_dzt(2:nz) )
        end do
      end do
      ! End vertical integral

    end if ! doclubb

    if ( doclubb .or. doclubbnoninter ) then

      ! We call CLUBB here because adjustments to the wind
      ! must occur prior to adams() -dschanen 26 Aug 2008
      ! Here we call clubb only if nstep divides the current timestep,
      ! or we're on the very first timestep
      if ( nstep == 1 .or. mod( nstep, nclubb ) == 0 ) then

        call advance_clubb_sgs &
             ( real( dtn*real( nclubb ), kind=time_precision), & ! in
               nstep, stats_nsamp, stats_nout, & ! in
               rho, rhow, wsub, u, v, w, qpl, qci, qpi, & ! in
               t, qv, qcl, &
               thlp2_forcing ) ! in/out
      end if ! nstep == 1 .or. mod( nstep, nclubb) == 0

    end if ! doclubb .or. doclubbnoninter

#endif

    return
  end subroutine sgs_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
  subroutine sgs_diagnose()
! None

  end subroutine sgs_diagnose


!----------------------------------------------------------------------
!!!! Collect sgs history statistics (vertical profiles)
!
  subroutine sgs_statistics()

    use hbuffer, only: hbuf_put, hbuf_avg_put
    use params, only : lcond
#ifdef CLUBB
    use sgs_params, only: doclubb, doclubbnoninter
    use clubbvars, only: upwp, vpwp, up2, vp2, wprtp, wpthlp, &
        wp2, wp3, rtp2, thlp2, rtpthlp, cloud_frac, rcm, um, vm
    use stat_clubb, only: stats_clubb_update ! Procedure(s)
#endif /* CLUBB */

    implicit none

    real factor_xy
    real tkz(nzm), tkhz(nzm)
    integer i,j,k !,n
! character(LEN=6) :: statname  !bloss: for conditional averages

    call t_startf ('sgs_statistics')

    factor_xy = 1./float(nx*ny)

    do k=1,nzm
      tkz(k) = 0.
      tkhz(k) = 0.
      do j=1,ny
        do i=1,nx
          tkz(k)=tkz(k)+tk(i,j,k)
          tkhz(k)=tkhz(k)+tkh(i,j,k)
        end do
      end do
    end do

    call hbuf_avg_put('TKES',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)

    call hbuf_put('TK',tkz,factor_xy)
    call hbuf_put('TKH',tkhz,factor_xy)

!---------------------------------------------------------
! SGS TKE Budget:

    call hbuf_put('ADVTRS',sgswle(:,1),factor_xy)
    call hbuf_put('BUOYAS',tkesbbuoy,factor_xy)
    call hbuf_put('SHEARS',tkesbshear,factor_xy)
    call hbuf_put('DISSIPS',tkesbdiss,factor_xy)
#ifdef CLUBB
        if ( doclubb .or. doclubbnoninter ) then
          call stats_clubb_update( upwp, vpwp, up2, vp2, wprtp, wpthlp, &
                wp2, wp3, rtp2, thlp2, rtpthlp, cloud_frac, rcm, um, vm )
        end if
#endif
    call t_stopf ('sgs_statistics')

    return
  end subroutine sgs_statistics

!----------------------------------------------------------------------
! called when stepout() called

  subroutine sgs_print()
#ifdef CLUBB
    use clubbvars, only: wp2, wp3, up2, vp2 ! Variable(s)
#endif
    implicit none

    call fminmax_print('tke:',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
    call fminmax_print('tk:',tk,0,nxp1,1-YES3D,nyp1,nzm)
    call fminmax_print('tkh:',tkh,0,nxp1,1-YES3D,nyp1,nzm)
#ifdef CLUBB /* Print some diagnostics for the state of CLUBB */
    call fminmax_print('wp2 sgs:',real( wp2 ),1,nx,1,ny,nz)
    call fminmax_print('wp3 sgs:',real( wp3 ),1,nx,1,ny,nz)
    call fminmax_print('up2 sgs:',real( up2 ),1,nx,1,ny,nz)
    call fminmax_print('vp2 sgs:',real( vp2 ),1,nx,1,ny,nz)
#endif /* CLUBB */

    return
  end subroutine sgs_print

!----------------------------------------------------------------------
!!! Initialize the list of sgs statistics
!
  subroutine sgs_hbuf_init(namelist,deflist,unitlist,status,average_type,count,sgscount)
#ifdef CLUBB
use sgs_params, only: doclubb, doclubbnoninter
#endif
    implicit none
    character(*) namelist(*), deflist(*), unitlist(*)
    integer status(*),average_type(*),count,sgscount
#ifdef CLUBB
    ! Local Variables
    integer :: &
      stat, & ! Whether the variable is enabled ( > 0 == enabled)
      stat_type ! Conditional stats?  All CLUBB stats are not.

    character(len=20) :: stat_name ! Name of the stat
    character(len=80) :: definition ! Description of the variable
    character(len=10) :: units    ! Units on the variable

    ! ---- Begin Code ----
    sgscount = 0 ! Note: sgs_count is named trcount in huffer.F90

    if ( doclubb .or. doclubbnoninter ) then

      ! Code template taken from hbuffer.F90 -dschanen 8 Sept 2013
      open(77,file='RUNDATA/clubb_lst', status='old',form='formatted')

333   continue
      read(77,err=333,end=444,fmt=*) stat,stat_type,stat_name,definition,units
      if ( stat.gt.0 ) then
        sgscount = sgscount + 1
        namelist(count+sgscount) = stat_name
        deflist(count+sgscount)  = definition
        unitlist(count+sgscount) = units
        status(count+sgscount)   = stat
        average_type(count+sgscount) = stat_type
      end if

      goto 333
444   continue

      count = count + sgscount ! Update total number of current stats

      close(unit=77)
    end if 

#endif /* CLUBB */

    return
  end subroutine sgs_hbuf_init


end module sgs
