module microphysics

! module for original SAM bulk microphysics
! Marat Khairoutdinov, 2006

use grid, only: nx,ny,nzm,nz, dimx1_s,dimx2_s,dimy1_s,dimy2_s ! subdomain grid information 
use params, only: doprecip, docloud
#ifdef CLUBB
use sgs_params, only: doclubb
#endif
use micro_params
implicit none

!----------------------------------------------------------------------
!!! required definitions:

integer, parameter :: nmicro_fields = 3   ! total number of prognostic water vars

!!! microphysics prognostic variables are storred in this array:

real micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields)

integer, parameter :: flag_wmass(nmicro_fields) = (/1,1,0/)
integer, parameter :: index_water_vapor = 1 ! index for variable that has water vapor
integer, parameter :: index_cloud_ice = -1   ! index for cloud ice (sedimentation)
#ifdef PNNL_STATS
integer :: &
  indx_qr = 2,  & ! Index for rain water mixing ratio (theta-l/qtog adv. budgets)
  indx_qi = -1, & ! Index for ice mixing ratio (theta-l/qtog adv. budgets)
  indx_qs = -1, & ! Index for snow mixing ratio (theta-l/qtog adv. budgets)
  indx_qg = -1    ! Index for graupel mixing ratio (theta-l/qtog adv. budgets)
#endif /* PNNL_STATS */
integer, parameter :: flag_precip(nmicro_fields) = (/0,1,1/)

! both variables correspond to mass, not number
integer, parameter :: flag_number(nmicro_fields) = (/0,0,1/)

! SAM1MOM 3D microphysical fields are output by default.
integer, parameter :: flag_micro3Dout(nmicro_fields) = (/0,0,1/)

real fluxbmk (nx, ny, 1:nmicro_fields) ! surface flux of tracers
real fluxtmk (nx, ny, 1:nmicro_fields) ! top boundary flux of tracers

!!! these arrays are needed for output statistics:

real mkwle(nz,1:nmicro_fields)  ! resolved vertical flux
real mkwsb(nz,1:nmicro_fields)  ! SGS vertical flux
real mkadv(nz,1:nmicro_fields)  ! tendency due to vertical advection
real mklsadv(nz,1:nmicro_fields)  ! tendency due to large-scale vertical advection
real mkdiff(nz,1:nmicro_fields)  ! tendency due to vertical diffusion

!======================================================================
! UW ADDITIONS

#ifdef UWM_MICRO_CHANGES
logical :: docloudsed  ! Flag for sedimentation of cloud water.
#endif
!bloss: arrays with names/units for microphysical outputs in statistics.
character*3, dimension(nmicro_fields) :: mkname
character*80, dimension(nmicro_fields) :: mklongname
character*10, dimension(nmicro_fields) :: mkunits
real, dimension(nmicro_fields) :: mkoutputscale
#ifdef UWM_STATS
real :: mstor(nzm, nmicro_fields)
#endif

! END UW ADDITIONS
!======================================================================

!------------------------------------------------------------------
! Optional (internal) definitions)

! make aliases for prognostic variables:
! note that the aliases should be local to microphysics

real q(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! total nonprecipitating water
real qp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! total precipitating water
real conp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! precipitating water number concentration
equivalence (q(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,1))
equivalence (qp(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,2))
equivalence (conp(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,3))
#ifdef SILHS
real conc(nx, ny, nzm)  ! Cloud droplet number concentration [#/kg]
#endif

real qn(nx,ny,nzm)  ! cloud condensate (liquid + ice)

real qpsrc(nz)  ! source of precipitation microphysical processes
real qpfall(nz) ! source of precipitating water due to fall out in a given level
real qpevp(nz)  ! sink of precipitating water due to evaporation

real vrain, vsnow, vgrau, crain, csnow, cgrau  ! precomputed coefs for precip terminal velocity

#ifdef UWM_STATS
! Statistics involving microphysics processes.
real qp_auto(nx,ny,nzm)  ! Autoconversion rate
real qp_accr(nx,ny,nzm)  ! Accretion rate
real qp_evap(nx,ny,nzm)  ! Evaporation rate
#endif /*UWM_STATS*/
CONTAINS

! required microphysics subroutines and function:
!----------------------------------------------------------------------
!!! Read microphysical options from prm file and allocate variables
!
subroutine micro_setparm()
  use vars
  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder
  
  !======================================================================
  ! UW ADDITION
  NAMELIST /MICRO_DRIZZLE/ &
#ifdef UWM_MICRO_CHANGES
       docloudsed, &    ! Flag for cloud water sedimentation
#endif
       Nc0 ! (cm-3) Prescribed cloud drop concentration 

  !bloss: Create dummy namelist, so that we can figure out error code
  !       for a mising namelist.  This lets us differentiate between
  !       missing namelists and those with an error within the namelist.
  NAMELIST /BNCUIODSBJCB/ place_holder

  Nc0 = 40. ! default cloud drop number concentration

#ifdef UWM_MICRO_CHANGES
  docloudsed = .true.  ! Cloud water sedimentation is turned on by default
#endif
  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 

  !bloss: get error code for missing namelist (by giving the name for
  !       a namelist that doesn't exist in the prm file).
  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  rewind(55) !note that one must rewind before searching for new namelists

  !bloss: read in MICRO_DRIZZLE namelist
  read (55,MICRO_DRIZZLE,IOSTAT=ios)

  if (ios.ne.0) then
     !namelist error checking
     if(ios.ne.ios_missing_namelist) then
        write(*,*) '****** ERROR: bad specification in MICRO_DRIZZLE namelist'
        call task_abort()
     elseif(masterproc) then
        write(*,*) '****************************************************'
        write(*,*) '***** No MICRO_DRIZZLE namelist in prm file ********'
        write(*,*) '****************************************************'
     end if
  end if
  close(55)

  if(masterproc) then
     write(*,*) 'Cloud droplet number concentration = ', Nc0, '/cm3'
#ifdef UWM_MICRO_CHANGES
     write(*,*) 'Cloud water sedimentation is allowed:', docloudsed
#endif
  end if
#ifdef UWM_STATS
  mstor = 0.0
#endif
  ! END UW ADDITION
  !======================================================================

end subroutine micro_setparm

!----------------------------------------------------------------------
!!! Initialize microphysics:

subroutine micro_init()

  use vars, only: q0
#ifdef SILHS
  use clubb_api_module, only: cm3_per_m3
  use vars, only: rho
#endif
  use grid, only: nrestart
#ifdef UWM_STATS
  integer k, n
#else
  integer k
  
#endif
  if(doprecip) call precip_init() 

  if(nrestart.eq.0) then

     micro_field = 0.
     do k=1,nzm
      q(:,:,k) = q0(k)
     end do
     qn = 0.
     fluxbmk = 0.
     fluxtmk = 0.

#ifdef CLUBB
     if(docloud.or.doclubb) then
#else
     if(docloud) then
#endif
       call cloud()
       call micro_diagnose()
#ifdef SILHS
       do k = 1, nzm
         conc(:,:,k) = Nc0 / rho(k) * cm3_per_m3
       end do
#endif
     end if

  end if

  mkwle = 0.
  mkwsb = 0.
  mkadv = 0.
  mkdiff = 0.

  qpsrc = 0.
  qpevp = 0.

  mkname(1) = 'QT'
  mklongname(1) = 'TOTAL WATER: VAPOR + CONDENSATE'
  mkunits(1) = 'g/kg'
  mkoutputscale(1) = 1.e3

  mkname(2) = 'QR'
  mklongname(2) = 'DRIZZLE MASS MIXING RATIO'
  mkunits(2) = 'g/kg'
  mkoutputscale(2) = 1.e3

  mkname(3) = 'CONP'
  mklongname(3) = 'DRIZZLE DROP CONCENTRATION'
  mkunits(3) = 'cm-3'
  mkoutputscale(3) = 1.e-6
#ifdef UWM_STATS
  do n=1, nmicro_fields
    do k=1, nzm
      mstor(k,n) = sum(micro_field(1:nx, 1:ny, k, n))
    end do
  end do
#endif

end subroutine micro_init

!----------------------------------------------------------------------
!!! fill-in surface and top boundary fluxes:
!
subroutine micro_flux()

  use vars, only: fluxbq, fluxtq

#ifdef CLUBB
  use params, only: doclubb_sfc_fluxes
  if ( doclubb .and. doclubb_sfc_fluxes ) then
    fluxbmk(:,:,index_water_vapor) = 0.0
  else
    fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
  end if
#else
  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
#endif
  fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)

end subroutine micro_flux

!----------------------------------------------------------------------
!!! compute local microphysics processes (bayond advection and SGS diffusion):
!
subroutine micro_proc()

   use grid, only: icycle
#ifdef UWM_STATS
   use grid, only: &
       nstep,        &
       nsave3D,      &
       nsave3Dstart, &
       nsave3Dend
#endif /* UWM_STATS */
   integer k

   ! Update bulk coefficient
   if(doprecip.and.icycle.eq.1) call precip_init() 

   if(docloud) then
     call cloud()
     if(doprecip) call precip_proc()
     call micro_diagnose()
   end if
#ifdef CLUBB
   if(doclubb) then
     if(doprecip) call precip_proc()
     call micro_diagnose()
   end if
#endif
#ifdef SILHS
   if ( .not. doclubb .and. .not. docloud ) then
     if(doprecip) call precip_proc()
   end if
#endif
#ifdef UWM_STATS
   if ( mod( nstep, nsave3D) .eq. 0 .and. nstep .ge. nsave3Dstart &
        .and. nstep .le. nsave3Dend ) then
      call write_3d_micro_fields()
   endif
#endif /* UWM_STATS */

end subroutine micro_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine micro_diagnose()
 
   use vars
   integer i,j,k

   do k=1,nzm
    do j=1,ny
     do i=1,nx
       qv(i,j,k) = q(i,j,k) - qn(i,j,k)
       qcl(i,j,k) = qn(i,j,k)
       qci(i,j,k) = 0.
       qpl(i,j,k) = qp(i,j,k)
       qpi(i,j,k) = 0.
     end do
    end do
   end do

end subroutine micro_diagnose

!----------------------------------------------------------------------
!!! function to compute terminal velocity for precipitating variables:

real function term_vel_qp(i,j,k,ind)
  
  integer, intent(in) ::  i,j,k ! current indexes
  real rvdr
  integer ind ! placeholder dummy variable

  term_vel_qp = 0.
  if(qp(i,j,k).gt.qp_threshold) then
      conp(i,j,k) = max(qp(i,j,k)*coefconpmin, conp(i,j,k))
      rvdr = (coefrv*qp(i,j,k)/conp(i,j,k))**0.333
      term_vel_qp= max(0.,1.2e4*rvdr-0.2)
  else
      qp(i,j,k)=0.
      conp(i,j,k) = 0.
  endif

end function term_vel_qp

real function term_vel_conp(i,j,k,ind)
  
  integer, intent(in) ::  i,j,k ! current indexes
  real rvdr
  integer ind ! placeholder dummy variable

  term_vel_conp = 0.
  if(qp(i,j,k).gt.qp_threshold) then
      conp(i,j,k) = max(qp(i,j,k)*coefconpmin, conp(i,j,k))
      rvdr = (coefrv*qp(i,j,k)/conp(i,j,k))**0.333
      term_vel_conp= max(0.,0.7e4*rvdr-0.1)
  else
      qp(i,j,k)=0.
      conp(i,j,k) = 0.
  endif

end function term_vel_conp

!----------------------------------------------------------------------
!!! compute sedimentation 
!
subroutine micro_precip_fall()
  
  use vars
  use params, only : pi

  real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real f0(nzm),df0(nzm)
  real omega(nx,ny,nzm)
  integer ind
  integer i,j,k
#ifdef PNNL_STATS
  real qtog1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real qtog_avg(nzm), qtog1_avg(nzm)
#endif /*PNNL_STATS*/

! Initialize arrays that accumulate surface precipitation flux

 if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
   do j=1,ny
    do i=1,nx
     precsfc(i,j)=0.
    end do
   end do
   do k=1,nzm
    precflux(k) = 0.
   end do
 end if

 do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
    qpfall(k)=0.
    tlat(k) = 0.
 end do

#ifdef UWM_STATS
 ! Note:  this was moved here from below in order to properly account for the
 !        sedimentation (fall) of rain water mixing ratio in the budgets for
 !        t'^2 and t'w'.
 if(dostatis) then
   do k=1,nzm
     do j=dimy1_s,dimy2_s
       do i=dimx1_s,dimx2_s
          df(i,j,k) = t(i,j,k)
       end do
     end do
   end do
#ifdef PNNL_STATS
   do k = 1, nzm
      do j = 1, ny
         do i = 1, nx
            ! The grand total (including precipitation) water mixing ratio for
            ! KK microphysics is qtog = qv + qcl + qpl (water vapor + cloud
            ! water + rain water), which is qtog = q + qp.
            qtog1(i,j,k) = q(i,j,k) + qp(i,j,k)
         enddo
       enddo
   enddo
#endif /*PNNL_STATS*/
 endif
#endif
 call precip_fall(qp, term_vel_qp, 0, omega, ind)
 call precip_fall(conp, term_vel_conp, 3, omega, ind)

! keep reasonable values:

do k=1,nzm
  do j=1,ny
      do i=1,nx
         if(qp(i,j,k).gt.qp_threshold) then
             conp(i,j,k) = max(qp(i,j,k)*coefconpmin, conp(i,j,k))
         else
             q(i,j,k) = q(i,j,k) +qp(i,j,k)
             qp(i,j,k)=0.
             conp(i,j,k) = 0.
         endif
      end do
  end do
end do

 if(dostatis) then
#if 0
! Note: The code below was commented out and moved up in order to properly account
!   for the sedimentation (fall) of rain water mixing ratio in the budgets
!   for the t'^2 and t'w'.
#endif /* false */
#ifndef UWM_STATS
   do k=1,nzm
     do j=dimy1_s,dimy2_s
       do i=dimx1_s,dimx2_s
          df(i,j,k) = t(i,j,k)
       end do
     end do
   end do
#endif
#ifdef PNNL_STATS
   do k = 1, nzm
      do j = 1, ny
         do i = 1, nx
            ! The grand total (including precipitation) water mixing ratio for
            ! KK microphysics is qtog = qv + qcl + qpl (water vapor + cloud
            ! water + rain water), which is qtog = q + qp.
            qtog(i,j,k) = q(i,j,k) + qp(i,j,k)
         enddo
       enddo
   enddo
#endif /*PNNL_STATS*/
   call stat_varscalar(t,df,f0,df0,t2leprec)
   call setvalue(twleprec,nzm,0.)
   call stat_sw2(t,df,twleprec)
#ifdef PNNL_STATS
   call stat_varscalar(qtog,qtog1,qtog_avg,qtog1_avg,qtog2leprec)
   call setvalue(qtogwleprec,nzm,0.)
   call stat_sw2(qtog,qtog1,qtogwleprec)
#endif /*PNNL_STATS*/
 endif

end subroutine micro_precip_fall

!----------------------------------------------------------------------
!!!! Collect microphysics history statistics (vertical profiles)
!
subroutine micro_statistics()
  
  use vars
  use hbuffer, only: hbuf_avg_put, hbuf_put
#ifdef UWM_STATS
  use params, only : lcond, LES
#else
  use params, only : lcond
#endif /*UWM_STATS*/
#ifdef UWM_STATS
  use calc_vars_util, only: &
      t2thetal  ! Procedure(s)

  use compute_correlation_module, only: &
      mean_xy_domain,          & ! Procedure(s)
      variance_xy_domain,      &
      covariance_xy_domain,    &
      mean_ip_xy_domain,       &
      variance_ip_xy_domain,   &
      covariance_ip_xy_domain, &
      sum_value_xy_domain
#endif /*UWM_STATS*/

  real tmp(2), factor_xy 
  real qcz(nzm), qiz(nzm), qrz(nzm), qsz(nzm), qgz(nzm), nrz(nzm), omg
#ifdef SILHS
  real ncz(nzm)
#endif
  integer i,j,k,n
  character(LEN=6) :: statname  !bloss: for conditional averages

#ifdef UWM_STATS
  !-------------------------------------------------------------
  ! Special Statistics -- Variances and Covariances
  !-------------------------------------------------------------
  real rt(nx,ny,nzm)
  real rc(nx,ny,nzm)
  real rr(nx,ny,nzm)
  real thl(nx,ny,nzm)
  real coef
  real Ncld(nx,ny,nzm)
  real Nr(nx,ny,nzm)
  real w_at_t(nx,ny,nzm)
  real precipitating(nx,ny,nzm)
  integer num_precip(nzm)
  ! Means
  real rtm(nzm)
  real rcm(nzm)
  real rrm(nzm)
  real thlm(nzm)
  real Ncm(nzm)
  real Nrm(nzm)
  real wm_at_t(nzm)
  real evapm(nzm)
  real autom(nzm)
  real accrm(nzm)
  real precip_frac(nzm)
  real rrm_ip(nzm)
  real Nrm_ip(nzm)
  real rtm_ip(nzm)
  real thlm_ip(nzm)
  real wm_ip_at_t(nzm)
  ! Variances
  real rtp2(nzm)
  real thlp2(nzm)
  real rcp2(nzm)
  real rrp2(nzm)
  real Nrp2(nzm)
  real Ncp2(nzm)
  real wp2(nzm)
  real rrp2_ip(nzm)
  real Nrp2_ip(nzm)
  real rtp2_ip(nzm)
  real thlp2_ip(nzm)
  real wp2_ip(nzm)
  ! Covariances
  real rtpthlp(nzm)
  real rtprrp(nzm)
  real thlprrp(nzm)
  real rcprrp(nzm)
  real rtpNrp(nzm)
  real thlpNrp(nzm)
  real rcpNrp(nzm)
  real rtpNcp(nzm)
  real thlpNcp(nzm)
  real rcpNcp(nzm)
  real rrpNrp(nzm)
  real wprtp(nzm)
  real wpthlp(nzm)
  real wprrp(nzm)
  real wpNrp(nzm)
  real wpNcp(nzm)
  real wpevapp(nzm)
  real rtpevapp(nzm)
  real thlpevapp(nzm)
  real wpautop(nzm)
  real rtpautop(nzm)
  real thlpautop(nzm)
  real wpaccrp(nzm)
  real rtpaccrp(nzm)
  real thlpaccrp(nzm)
  real rtpthlp_ip(nzm)
  real wprrp_ip(nzm)
  real wpNrp_ip(nzm)
  real rtprrp_ip(nzm)
  real thlprrp_ip(nzm)
  real rtpNrp_ip(nzm)
  real thlpNrp_ip(nzm)
  real rrpNrp_ip(nzm)

#endif /*UWM_STATS*/
  call t_startf ('statistics')

  factor_xy = 1./float(nx*ny)

  do k=1,nzm
      tmp(1) = dz/rhow(k)
      tmp(2) = tmp(1) / dtn
      mkwsb(k,1) = mkwsb(k,1) * tmp(1) * rhow(k) * lcond
      mkwle(k,1) = mkwle(k,1)*tmp(2)*rhow(k)*lcond + mkwsb(k,1)
#ifdef CLUBB
      if((docloud.or.doclubb).and.doprecip) then
#else
      if(docloud.and.doprecip) then
#endif /* CLUBB */
        mkwsb(k,2) = mkwsb(k,2) * tmp(1) * rhow(k) * lcond
        mkwle(k,2) = mkwle(k,2)*tmp(2)*rhow(k)*lcond + mkwsb(k,2)
      endif
#ifdef UWM_STATS     
      do n=1, nmicro_fields
        mstor(k,n) = SUM(micro_field(1:nx,1:ny,k,n))-mstor(k,n)
      end do
#endif
  end do

  call hbuf_put('QTFLUX',mkwle(:,1),factor_xy)
  call hbuf_put('QTFLUXS',mkwsb(:,1),factor_xy)
  call hbuf_put('QPFLUX',mkwle(:,2),factor_xy)
  call hbuf_put('QPFLUXS',mkwsb(:,2),factor_xy)

  do k=1,nzm
    qcz(k) = 0.
    qiz(k) = 0.
    qrz(k) = 0.
    qsz(k) = 0.
    qgz(k) = 0.
    do j=1,ny
    do i=1,nx
      qcz(k)=qcz(k)+qcl(i,j,k)
      qiz(k)=0.
      qrz(k)=qrz(k)+qpl(i,j,k)
      qsz(k)=0.
      qgz(k)=0.
    end do
    end do
  end do

  call hbuf_put('QC',qcz,1.e3*factor_xy)
  call hbuf_put('QR',qrz,1.e3*factor_xy)

  call hbuf_avg_put('CONP',conp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.e-6)
#ifdef SILHS
  call hbuf_avg_put('NC',conc,1,nx,1,ny,nzm,1.e-6)
#endif

  call hbuf_put('QTADV',mkadv(:,1)+qifall,factor_xy*86400000./dtn)
  call hbuf_put('QTDIFF',mkdiff(:,1),factor_xy*86400000./dtn)
  call hbuf_put('QTSINK',qpsrc,-factor_xy*86400000./dtn)
  call hbuf_put('QTSRC',qpevp,-factor_xy*86400000./dtn)
  call hbuf_put('QPADV',mkadv(:,2),factor_xy*86400000./dtn)
  call hbuf_put('QPDIFF',mkdiff(:,2),factor_xy*86400000./dtn)
  call hbuf_put('QPFALL',qpfall,factor_xy*86400000./dtn)
  call hbuf_put('QPSRC',qpsrc,factor_xy*86400000./dtn)
  call hbuf_put('QPEVP',qpevp,factor_xy*86400000./dtn)

  do n = 1,nmicro_fields
     if(flag_wmass(n).lt.1) then
        ! remove factor of rho from number concentrations
        mklsadv(1:nzm,n) = mklsadv(1:nzm,n)*rho(:)
#ifdef UWM_STATS
        mstor(1:nzm,n) = mstor(1:nzm,n) * rho(:)
#endif
     end if
     call hbuf_put(trim(mkname(n))//'LSADV', &
          mklsadv(:,n),mkoutputscale(n)*factor_xy*86400.)
  end do

  do n = 1,ncondavg

     do k=1,nzm
        qcz(k) = 0.
        qrz(k) = 0.
        nrz(k) = 0.
#ifdef SILHS
        ncz(k) = 0.
#endif
        do j=1,ny
           do i=1,nx
              qcz(k)=qcz(k)+qcl(i,j,k)*condavg_mask(i,j,k,n)
              qrz(k)=qrz(k)+qpl(i,j,k)*condavg_mask(i,j,k,n)
              ! drizzle number concentration
              nrz(k)=nrz(k)+rho(k)*micro_field(i,j,k,3)*condavg_mask(i,j,k,n)
#ifdef SILHS
              ncz(k)=ncz(k)+rho(k)*conc(i,j,k)*condavg_mask(i,j,k,n)
#endif
           end do
        end do
     end do

     call hbuf_put('QC' // TRIM(condavgname(n)),qcz,1.e3)
     call hbuf_put('QR' // TRIM(condavgname(n)),qrz,1.e3)
#ifdef UWM_MISC
     call hbuf_put('CONP' // TRIM(condavgname(n)),nrz,1.e-6)
#else /* This looks like a bug  -dschanen 22 Aug 2012 */
     call hbuf_put('CONP' // TRIM(condavgname(n)),qrz,1.e-6)
#endif
#ifdef SILHS
     call hbuf_put('NC' // TRIM(condavgname(n)),ncz,1.e-6)
#endif
#ifdef UWM_STATS
     call hbuf_put(trim(mkname(n))//'STO', &
          mstor(:,n), mkoutputscale(n)*factor_xy*86400./dtn)
#endif
  end do

#ifdef UWM_STATS
  !-------------------------------------------------------------
  ! Special Statistics -- Variances and Covariances
  !-------------------------------------------------------------

  ! Place variables in 3D arrays.
  do k = 1, nzm
     do i = 1, nx
        do j = 1, ny

           ! Moisture variables
           rt(i,j,k) = qv(i,j,k) + qcl(i,j,k)
           rc(i,j,k) = qcl(i,j,k)
           rr(i,j,k) = qpl(i,j,k)

           ! Thetal
           thl(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                  0.0, 0.0, prespot(k) )

           ! Drop and droplet concentrations
           Nr(i,j,k) = conp(i,j,k)

           if (LES) then
              coef=0.
           else
              coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
           endif
           if ( rc(i,j,k) .gt. coef ) then
              Ncld(i,j,k) = Nc0*1.0e+6
           else
              Ncld(i,j,k) = 0.0
           endif

           ! Vertical Velocity
           w_at_t(i,j,k) = 0.5 * ( w(i,j,k) + w(i,j,k+1) )

           ! Statistics for in-precip values.
           if ( qpl(i,j,k) > 0.0 ) then
              precipitating(i,j,k) = 1.0
           else
              precipitating(i,j,k) = 0.0
           endif

        enddo ! j = 1, ny
     enddo ! i = 1, nx
  enddo ! k = 1, nzm

  num_precip = int( sum_value_xy_domain( nzm, precipitating ) )

  ! Find the horizontal mean values of the variables.
  rtm     = mean_xy_domain( nzm, rt )
  rcm     = mean_xy_domain( nzm, rc )
  rrm     = mean_xy_domain( nzm, rr )
  thlm    = mean_xy_domain( nzm, thl )
  Nrm     = mean_xy_domain( nzm, Nr )
  Ncm     = mean_xy_domain( nzm, Ncld )
  wm_at_t = mean_xy_domain( nzm, w_at_t )

  evapm = mean_xy_domain( nzm, qp_evap )
  autom = mean_xy_domain( nzm, qp_auto )
  accrm = mean_xy_domain( nzm, qp_accr )

  precip_frac = mean_xy_domain( nzm, precipitating )

  ! Statistics for in-precip values.
  rrm_ip = mean_ip_xy_domain( nzm, rr, precipitating, num_precip )
  Nrm_ip = mean_ip_xy_domain( nzm, Nr, precipitating, num_precip )

  rtm_ip     = mean_ip_xy_domain( nzm, rt, precipitating, num_precip )
  thlm_ip    = mean_ip_xy_domain( nzm, thl, precipitating, num_precip )
  wm_ip_at_t = mean_ip_xy_domain( nzm, w_at_t, precipitating, num_precip )

  ! Find the horizontal variances of the variables.
  rtp2  = variance_xy_domain( nzm, rt, rtm )
  rcp2  = variance_xy_domain( nzm, rc, rcm )
  rrp2  = variance_xy_domain( nzm, rr, rrm )
  thlp2 = variance_xy_domain( nzm, thl, thlm )
  Nrp2  = variance_xy_domain( nzm, Nr, Nrm )
  Ncp2  = variance_xy_domain( nzm, Ncld, Ncm )
  wp2   = variance_xy_domain( nzm, w_at_t, wm_at_t )

  ! Statistics for in-precip values.
  rrp2_ip = variance_ip_xy_domain( nzm, rr, precipitating, rrm_ip, num_precip )
  Nrp2_ip = variance_ip_xy_domain( nzm, Nr, precipitating, Nrm_ip, num_precip )

  rtp2_ip  = variance_ip_xy_domain( nzm, rt, precipitating, &
                                    rtm_ip, num_precip )
  thlp2_ip = variance_ip_xy_domain( nzm, thl, precipitating, &
                                    thlm_ip, num_precip )
  wp2_ip   = variance_ip_xy_domain( nzm, w_at_t, precipitating, &
                                    wm_ip_at_t, num_precip )

  ! Find the horizontal covariances of the variables.
  rtpthlp = covariance_xy_domain( nzm, rt, thl, rtm, thlm )
  rtprrp  = covariance_xy_domain( nzm, rt, rr, rtm, rrm )
  thlprrp = covariance_xy_domain( nzm, thl, rr, thlm, rrm )
  rcprrp  = covariance_xy_domain( nzm, rc, rr, rcm, rrm )
  rtpNrp  = covariance_xy_domain( nzm, rt, Nr, rtm, Nrm )
  thlpNrp = covariance_xy_domain( nzm, thl, Nr, thlm, Nrm )
  rcpNrp  = covariance_xy_domain( nzm, rc, Nr, rcm, Nrm )
  rtpNcp  = covariance_xy_domain( nzm, rt, Ncld, rtm, Ncm )
  thlpNcp = covariance_xy_domain( nzm, thl, Ncld, thlm, Ncm )
  rcpNcp  = covariance_xy_domain( nzm, rc, Ncld, rcm, Ncm )
  rrpNrp  = covariance_xy_domain( nzm, rr, Nr, rrm, Nrm )
  wprtp   = covariance_xy_domain( nzm, w_at_t, rt, wm_at_t, rtm )
  wpthlp  = covariance_xy_domain( nzm, w_at_t, thl, wm_at_t, thlm )
  wprrp   = covariance_xy_domain( nzm, w_at_t, rr, wm_at_t, rrm )
  wpNrp   = covariance_xy_domain( nzm, w_at_t, Nr, wm_at_t, Nrm )
  wpNcp   = covariance_xy_domain( nzm, w_at_t, Ncld, wm_at_t, Ncm )

  wpevapp   = covariance_xy_domain( nzm, w_at_t, qp_evap, wm_at_t, evapm )
  rtpevapp  = covariance_xy_domain( nzm, rt, qp_evap, rtm, evapm )
  thlpevapp = covariance_xy_domain( nzm, thl, qp_evap, thlm, evapm )
  wpautop   = covariance_xy_domain( nzm, w_at_t, qp_auto, wm_at_t, autom )
  rtpautop  = covariance_xy_domain( nzm, rt, qp_auto, rtm, autom )
  thlpautop = covariance_xy_domain( nzm, thl, qp_auto, thlm, autom )
  wpaccrp   = covariance_xy_domain( nzm, w_at_t, qp_accr, wm_at_t, accrm )
  rtpaccrp  = covariance_xy_domain( nzm, rt, qp_accr, rtm, accrm )
  thlpaccrp = covariance_xy_domain( nzm, thl, qp_accr, thlm, accrm )

  ! Statistics for in-precip values.
  rtpthlp_ip = covariance_ip_xy_domain( nzm, rt, thl, precipitating, &
                                        rtm_ip, thlm_ip, num_precip )
  wprrp_ip   = covariance_ip_xy_domain( nzm, w_at_t, rr, precipitating, &
                                        wm_ip_at_t, rrm_ip, num_precip )
  wpNrp_ip   = covariance_ip_xy_domain( nzm, w_at_t, Nr, precipitating, &
                                        wm_ip_at_t, Nrm_ip, num_precip )
  rtprrp_ip  = covariance_ip_xy_domain( nzm, rt, rr, precipitating, &
                                        rtm_ip, rrm_ip, num_precip )
  thlprrp_ip = covariance_ip_xy_domain( nzm, thl, rr, precipitating, &
                                        thlm_ip, rrm_ip, num_precip )
  rtpNrp_ip  = covariance_ip_xy_domain( nzm, rt, Nr, precipitating, &
                                        rtm_ip, Nrm_ip, num_precip )
  thlpNrp_ip = covariance_ip_xy_domain( nzm, thl, Nr, precipitating, &
                                        thlm_ip, Nrm_ip, num_precip )
  rrpNrp_ip  = covariance_ip_xy_domain( nzm, rr, Nr, precipitating, &
                                        rrm_ip, Nrm_ip, num_precip )

  ! Output mean horizontal averages of the mean values (for cross-checking).
  call hbuf_put('RTM',rtm,1.0)
  call hbuf_put('RCM',rcm,1.0)
  call hbuf_put('RRM',rrm,1.0)
  call hbuf_put('THLM',thlm,1.0)
  call hbuf_put('WM',wm_at_t,1.0)
  call hbuf_put('WP2',wp2,1.0)
  call hbuf_put('WPRTP',wprtp,1.0)
  call hbuf_put('WPTHLP',wpthlp,1.0)
  ! Output mean horizontal averages of the statistical parameters.
  call hbuf_put('NCM',Ncm,1.0)
  call hbuf_put('NRM',Nrm,1.0)
  call hbuf_put('EVAPM',evapm,1.0)
  call hbuf_put('AUTOM',autom,1.0)
  call hbuf_put('ACCRM',accrm,1.0)
  call hbuf_put('RTP2',rtp2,1.0)
  call hbuf_put('RTPTHLP',rtpthlp,1.0)
  call hbuf_put('THLP2',thlp2,1.0)
  call hbuf_put('RCP2',rcp2,1.0)
  call hbuf_put('RTPRRP',rtprrp,1.0)
  call hbuf_put('THLPRRP',thlprrp,1.0)
  call hbuf_put('RRP2',rrp2,1.0)
  call hbuf_put('RCPRRP',rcprrp,1.0)
  call hbuf_put('RTPNRP',rtpNrp,1.0)
  call hbuf_put('THLPNRP',thlpNrp,1.0)
  call hbuf_put('NRP2',Nrp2,1.0)
  call hbuf_put('RCPNRP',rcpNrp,1.0)
  call hbuf_put('RTPNCP',rtpNcp,1.0)
  call hbuf_put('THLPNCP',thlpNcp,1.0)
  call hbuf_put('NCP2',Ncp2,1.0)
  call hbuf_put('RCPNCP',rcpNcp,1.0)
  call hbuf_put('RRPNRP',rrpNrp,1.0)
  call hbuf_put('WPRRP',wprrp,1.0)
  call hbuf_put('WPNRP',wpNrp,1.0)
  call hbuf_put('WPNCP',wpNcp,1.0)
  call hbuf_put('WPEVAPP',wpevapp,1.0)
  call hbuf_put('RTPEVAPP',rtpevapp,1.0)
  call hbuf_put('THLPEVAPP',thlpevapp,1.0)
  call hbuf_put('WPAUTOP',wpautop,1.0)
  call hbuf_put('RTPAUTOP',rtpautop,1.0)
  call hbuf_put('THLPAUTOP',thlpautop,1.0)
  call hbuf_put('WPACCRP',wpaccrp,1.0)
  call hbuf_put('RTPACCRP',rtpaccrp,1.0)
  call hbuf_put('THLPACCRP',thlpaccrp,1.0)
  ! Statistics for in-precip values.
  call hbuf_put('PREC_FRAC',precip_frac,1.0)
  call hbuf_put('RRM_IP',rrm_ip,1.0)
  call hbuf_put('NRM_IP',Nrm_ip,1.0)
  call hbuf_put('RTM_IP',rtm_ip,1.0)
  call hbuf_put('THLM_IP',thlm_ip,1.0)
  call hbuf_put('WM_IP',wm_ip_at_t,1.0)
  call hbuf_put('RRP2_IP',rrp2_ip,1.0)
  call hbuf_put('NRP2_IP',Nrp2_ip,1.0)
  call hbuf_put('RTP2_IP',rtp2_ip,1.0)
  call hbuf_put('THLP2_IP',thlp2_ip,1.0)
  call hbuf_put('WP2_IP',wp2_ip,1.0)
  call hbuf_put('RTPTHLP_IP',rtpthlp_ip,1.0)
  call hbuf_put('WPRRP_IP',wprrp_ip,1.0)
  call hbuf_put('WPNRP_IP',wpNrp_ip,1.0)
  call hbuf_put('RTPRRP_IP',rtprrp_ip,1.0)
  call hbuf_put('THLPRRP_IP',thlprrp_ip,1.0)
  call hbuf_put('RTPNRP_IP',rtpNrp_ip,1.0)
  call hbuf_put('THLPNRP_IP',thlpNrp_ip,1.0)
  call hbuf_put('RRPNRP_IP',rrpNrp_ip,1.0)

#endif /*UWM_STATS*/
  call t_stopf ('statistics')

end subroutine micro_statistics

!----------------------------------------------------------------------
! called when stepout() called

subroutine micro_print()
  call fminmax_print('conp:',conp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
#ifdef SILHS
  call fminmax_print('conc:',conc,1,nx,1,ny,nzm)
#endif
end subroutine micro_print

!----------------------------------------------------------------------
!!! Initialize the list of microphysics statistics 
!
subroutine micro_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)

  use vars


   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count,trcount
   integer ntr, n


   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QTFLUX'
   deflist(count) = 'Nonprecipitating water flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QTFLUXS'
   deflist(count) = 'Nonprecipitating-water flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QPFLUX'
   deflist(count) = 'Precipitating-water turbulent flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QPFLUXS'
   deflist(count) = 'Precipitating-water turbulent flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QC'
   deflist(count) = 'Liquid cloud water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'CONP'
   deflist(count) = 'Drizzle drop concentration'
   unitlist(count) = '1/cm3'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QR'
   deflist(count) = 'Rain water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

#ifdef SILHS
   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NC'
   deflist(count) = 'Cloud drop concentration'
   unitlist(count) = '1/cm3'
   status(count) = 1    
   average_type(count) = 0
#endif
   do n = 1,nmicro_fields
      count = count + 1
      trcount = trcount + 1
      namelist(count) = TRIM(mkname(n))//'LSADV'
      deflist(count) = 'Source of '//TRIM(mklongname(n))//' due to large-scale vertical advection'
      unitlist(count) = TRIM(mkunits(n))//'day'
      status(count) = 1    
      average_type(count) = 0
   end do

  !bloss: setup to add an arbitrary number of conditional statistics
   do n = 1,ncondavg

      count = count + 1
      trcount = trcount + 1
      namelist(count) = 'QC' // TRIM(condavgname(n))
      deflist(count) = 'Mean Liquid cloud water in ' // TRIM(condavglongname(n))
      unitlist(count) = 'g/kg'
      status(count) = 1    
      average_type(count) = n

      count = count + 1
      trcount = trcount + 1
      namelist(count) = 'QR' // TRIM(condavgname(n))
      deflist(count) = 'Mean Drizzle Mixing ratio in ' // TRIM(condavglongname(n))
      unitlist(count) = 'g/kg'
      status(count) = 1    
      average_type(count) = n

      count = count + 1
      trcount = trcount + 1
      namelist(count) = 'CONP' // TRIM(condavgname(n))
      deflist(count) = 'Mean drizzle drop concentration in ' // TRIM(condavglongname(n))
      unitlist(count) = '/cm3'
      status(count) = 1    
      average_type(count) = n

#ifdef SILHS
      count = count + 1
      trcount = trcount + 1
      namelist(count) = 'NC' // TRIM(condavgname(n))
      deflist(count) = 'Mean cloud drop concentration in ' // TRIM(condavglongname(n))
      unitlist(count) = '/cm3'
      status(count) = 1    
      average_type(count) = n
#endif
   end do

#ifdef UWM_STATS
   !-------------------------------------------------------------
   ! Special Statistics -- Variances and Covariances
   !-------------------------------------------------------------

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTM'
   deflist(count) = 'Mean total water mixing ratio'
   unitlist(count) = 'kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RCM'
   deflist(count) = 'Mean cloud water mixing ratio'
   unitlist(count) = 'kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RRM'
   deflist(count) = 'Mean rain water mixing ratio'
   unitlist(count) = 'kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'THLM'
   deflist(count) = 'Mean liquid water potential temperature'
   unitlist(count) = 'K'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WM'
   deflist(count) = 'Mean vertical velocity'
   unitlist(count) = 'm/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WP2'
   deflist(count) = 'Variance of vertical velocity'
   unitlist(count) = 'm2/s2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WPRTP'
   deflist(count) = 'Covariance between w and rt'
   unitlist(count) = 'kg/kg m/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WPTHLP'
   deflist(count) = 'Covariance between w and theta-l'
   unitlist(count) = 'K m/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NCM'
   deflist(count) = 'Mean cloud droplet concentration'
   unitlist(count) = 'm-3'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NRM'
   deflist(count) = 'Mean rain drop concentration'
   unitlist(count) = 'm-3'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'EVAPM'
   deflist(count) = 'Mean rain water evaporation rate'
   unitlist(count) = 'kg/kg/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'AUTOM'
   deflist(count) = 'Mean rain water autoconversion rate'
   unitlist(count) = 'kg/kg/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'ACCRM'
   deflist(count) = 'Mean rain water accretion rate'
   unitlist(count) = 'kg/kg/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTP2'
   deflist(count) = 'Variance of rt'
   unitlist(count) = 'kg2/kg2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'THLP2'
   deflist(count) = 'Variance of theta-l'
   unitlist(count) = 'K2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTPTHLP'
   deflist(count) = 'Covariance between rt and theta-l'
   unitlist(count) = 'K kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RCP2'
   deflist(count) = 'Variance of rc'
   unitlist(count) = 'kg2/kg2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTPRRP'
   deflist(count) = 'Covariance between rt and rr'
   unitlist(count) = 'kg2/kg2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'THLPRRP'
   deflist(count) = 'Covariance between theta-l and rr'
   unitlist(count) = 'K kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RRP2'
   deflist(count) = 'Variance of rr'
   unitlist(count) = 'kg2/kg2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RCPRRP'
   deflist(count) = 'Covariance between rc and rr'
   unitlist(count) = 'kg2/kg2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTPNRP'
   deflist(count) = 'Covariance between rt and Nr'
   unitlist(count) = 'm-3 kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'THLPNRP'
   deflist(count) = 'Covariance between theta-l and Nr'
   unitlist(count) = 'K m-3'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NRP2'
   deflist(count) = 'Variance of Nr'
   unitlist(count) = 'm-6'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RCPNRP'
   deflist(count) = 'Covariance between rc and Nr'
   unitlist(count) = 'm-3 kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTPNCP'
   deflist(count) = 'Covariance between rt and Nc'
   unitlist(count) = 'm-3 kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'THLPNCP'
   deflist(count) = 'Covariance between theta-l and Nc'
   unitlist(count) = 'K m-3'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NCP2'
   deflist(count) = 'Variance of Nc'
   unitlist(count) = 'm-6'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RCPNCP'
   deflist(count) = 'Covariance between rc and Nc'
   unitlist(count) = 'm-3 kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RRPNRP'
   deflist(count) = 'Covariance between rr and Nr'
   unitlist(count) = 'm-3 kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WPRRP'
   deflist(count) = 'Covariance between w and rr'
   unitlist(count) = 'kg/kg m/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WPNRP'
   deflist(count) = 'Covariance between w and Nr'
   unitlist(count) = 'm-2 s-1'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WPNCP'
   deflist(count) = 'Covariance between w and Nc'
   unitlist(count) = 'm-2 s-1'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WPEVAPP'
   deflist(count) = 'Covariance between w and rain water evaporation rate'
   unitlist(count) = 'm kg/kg/s2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTPEVAPP'
   deflist(count) = 'Covariance between rt and rain water evaporation rate'
   unitlist(count) = 'kg2/kg2/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'THLPEVAPP'
   deflist(count) = 'Covariance between thl and rain water evaporation rate'
   unitlist(count) = 'K kg/kg/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WPAUTOP'
   deflist(count) = 'Covariance between w and rain water autoconversion rate'
   unitlist(count) = 'm kg/kg/s2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTPAUTOP'
   deflist(count) = 'Covariance between rt and rain water autoconversion rate'
   unitlist(count) = 'kg2/kg2/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'THLPAUTOP'
   deflist(count) = 'Covariance between thl and rain water autoconversion rate'
   unitlist(count) = 'K kg/kg/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WPACCRP'
   deflist(count) = 'Covariance between w and rain water accretion rate'
   unitlist(count) = 'm kg/kg/s2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTPACCRP'
   deflist(count) = 'Covariance between rt and rain water accretion rate'
   unitlist(count) = 'kg2/kg2/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'THLPACCRP'
   deflist(count) = 'Covariance between thl and rain water accretion rate'
   unitlist(count) = 'K kg/kg/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'PREC_FRAC'
   deflist(count) = 'Precipitation Fraction'
   unitlist(count) = '-'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RRM_IP'
   deflist(count) = 'Mean rain water mixing ratio in precip'
   unitlist(count) = 'kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NRM_IP'
   deflist(count) = 'Mean rain drop concentration in precip'
   unitlist(count) = 'm-3'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTM_IP'
   deflist(count) = 'Mean total water mixing ratio in precip'
   unitlist(count) = 'kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'THLM_IP'
   deflist(count) = 'Mean liquid water potential temperature in precip'
   unitlist(count) = 'K'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WM_IP'
   deflist(count) = 'Mean vertical velocity in precip'
   unitlist(count) = 'm/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RRP2_IP'
   deflist(count) = 'Variance of rr in precip'
   unitlist(count) = 'kg2/kg2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NRP2_IP'
   deflist(count) = 'Variance of Nr in precip'
   unitlist(count) = 'm-6'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTP2_IP'
   deflist(count) = 'Variance of rt in precip'
   unitlist(count) = 'kg2/kg2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'THLP2_IP'
   deflist(count) = 'Variance of theta-l in precip'
   unitlist(count) = 'K2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WP2_IP'
   deflist(count) = 'Variance of vertical velocity in precip'
   unitlist(count) = 'm2/s2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTPTHLP_IP'
   deflist(count) = 'Covariance between rt and theta-l in precip'
   unitlist(count) = 'K kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WPRRP_IP'
   deflist(count) = 'Covariance between w and rr in precip'
   unitlist(count) = 'kg/kg m/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'WPNRP_IP'
   deflist(count) = 'Covariance between w and Nr in precip'
   unitlist(count) = 'm-2 s-1'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTPRRP_IP'
   deflist(count) = 'Covariance between rt and rr in precip'
   unitlist(count) = 'kg2/kg2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'THLPRRP_IP'
   deflist(count) = 'Covariance between theta-l and rr in precip'
   unitlist(count) = 'K kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RTPNRP_IP'
   deflist(count) = 'Covariance between rt and Nr in precip'
   unitlist(count) = 'm-3 kg/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'THLPNRP_IP'
   deflist(count) = 'Covariance between theta-l and Nr in precip'
   unitlist(count) = 'K m-3'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'RRPNRP_IP'
   deflist(count) = 'Covariance between rr and Nr in precip'
   unitlist(count) = 'kg/kg m-3'
   status(count) = 1    
   average_type(count) = 0
#endif /*UWM_STATS*/
end subroutine micro_hbuf_init


!-----------------------------------------------------------------------
! Supply function that computes total water in a domain:
!
real(8) function total_water()

  use vars, only : nstep,nprint,adz,dz,rho
  real(8) tmp
  integer i,j,k,m

  total_water = 0.
  do m=1,nmicro_fields
   if(flag_wmass(m).eq.1) then
    do k=1,nzm
      tmp = 0.
      do j=1,ny
        do i=1,nx
          tmp = tmp + micro_field(i,j,k,m)
        end do
      end do
      total_water = total_water + tmp*adz(k)*dz*rho(k)
    end do
   end if
  end do

end function total_water

! -------------------------------------------------------------------------------
! Dummy effective radius functions:

function Get_reffc() ! liquid water
  real, pointer, dimension(:,:,:) :: Get_reffc
end function Get_reffc

function Get_reffi() ! ice
  real, pointer, dimension(:,:,:) :: Get_reffi
end function Get_reffi

#ifdef CLUBB

!---------------------------------------------------------------------
subroutine micro_update()

! Description:
! This subroutine essentially does what micro_proc does but does not
! call any microphysics subroutines.  We need to do this for the 
! single-moment bulk microphysics (SAM1MOM) so that CLUBB gets a
! properly updated value of ice fed in. 
!
! -dschanen UWM
!---------------------------------------------------------------------

  ! Update the dynamical core variables (e.g. qv, qcl) with the value in
  ! micro_field.  Diffusion, advection, and other processes are applied to
  ! micro_field but not the variables in vars.f90
  call micro_diagnose()

  return
end subroutine micro_update

!---------------------------------------------------------------------
subroutine micro_adjust( new_qv, new_qc )
! Description:
!   Adjust total water in SAM based on values from CLUBB.
! References:
!   None
!---------------------------------------------------------------------

  implicit none

  real, dimension(nx,ny,nzm), intent(in) :: &
    new_qv, & ! Water vapor mixing ratio that has been adjusted by CLUBB [kg/kg]
    new_qc    ! Cloud water mixing ratio that has been adjusted by CLUBB [kg/kg]

  ! Total water mixing ratio
  q(1:nx,1:ny,1:nzm) = new_qv(1:nx,1:ny,1:nzm) &
                     + new_qc(1:nx,1:ny,1:nzm)

  ! Cloud water mixing ratio
  qn(1:nx,1:ny,1:nzm) = new_qc(1:nx,1:ny,1:nzm) 

  return
end subroutine micro_adjust

#endif /*CLUBB*/

#ifdef UWM_STATS
!---------------------------------------------------------------------------
subroutine write_3d_micro_fields()

! Description:
!   Produce 3D output files containing values of selected fields at 3D output
!   times.
!
!   This is largely 'write_3D_fields.F90', but reused here. 
!---------------------------------------------------------------------------

use vars, only: &
    lenstr,  & ! Variable(s)
    t,       &
    qv,      &
    qcl,     &
    w,       &
    qpl,     &
    gamaz,   &
    prespot

use domain, only: &
    nsubdomains_x, & ! Variable(s)
    nsubdomains_y

use grid, only:  &
    masterproc,  & ! Variable(s)
    output_sep,  &
    rank,        &
    nsubdomains, &
    nstep,       &
    RUN3D,       &
    save3Dbin,   &
    case,        &
    caseid,      &
    save3Dsep,   &
    nrestart,    &
    notopened3D, &
    nx,          &
    ny,          &
    nzm,         &
    z,           &
    pres,        &
    dx,          &
    dy,          &
    nstep,       &
    dt,          &
    day0,        &
    dompi,       &
    dogzip3D

use calc_vars_util, only: &
    t2thetal  ! Procedure(s)

use compute_chi_module, only: &
    compute_chi_eta  ! Procedure(s)

implicit none

character *120 filename
character *80 long_name
character *8 name
character *10 timechar
character *4 rankchar
character *5 sepchar
character *12 filetype
character *10 units
character *12 c_z(nzm), c_p(nzm), c_dx, c_dy, c_time
integer i, j, k, nfields, nfields1
real tmp(nx,ny,nzm)

real, dimension(nx,ny,nzm) :: &
  thl, & ! Liquid water potential temperature                 [K]
  rt,  & ! Total water mixing ratio                           [kg/kg]
  chi, & ! Extended liquid water mixing ratio                 [kg/kg]
  eta    ! Coordinate orthogonal to chi in PDF transformation [kg/kg]


call t_startf('3D_out')

nfields  = 11  ! number of 3D fields to save
nfields1 = 0   ! assertion check

if ( masterproc .or. output_sep ) then
  
   if ( output_sep ) then
      write(rankchar,'(i4)') rank
      sepchar="_"//rankchar(5-lenstr(rankchar):4)
   else
      sepchar=""
   endif
  
   write(rankchar,'(i4)') nsubdomains
   write(timechar,'(i10)') nstep
  
   do k = 1, 11-lenstr(timechar)-1
      timechar(k:k)='0'
   enddo

   if ( RUN3D ) then
      if ( save3Dbin ) then
         filetype = '_micro.bin3D'
      else
         filetype = '_micro.com3D'
      endif
    
      filename = './OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
                 rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)// &
                 filetype//sepchar
    
      open( 46, file=filename, status='unknown', form='unformatted' )

   else
    
      if ( save3Dbin ) then
         if ( save3Dsep ) then
            filetype = '_micro.bin3D'
         else
            filetype = '_micro.bin2D'
         endif
      else
         if ( save3Dsep ) then
            filetype = '_micro.com3D'
         else
            filetype = '_micro.com2D'
         endif
      endif
  
      if ( save3Dsep ) then
         filename = './OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
                    rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)// &
                    filetype//sepchar
         open( 46, file=filename, status='unknown', form='unformatted' )
      else
         filename = './OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
                    rankchar(5-lenstr(rankchar):4)//filetype//sepchar
         if ( nrestart .eq. 0 .and. notopened3D ) then
            open( 46, file=filename, status='unknown', form='unformatted' )
         else
            open( 46, file=filename, status='unknown', &
                  form='unformatted', position='append' )
         endif
         notopened3D=.false.
      endif ! save3Dsep

   endif ! RUN3D

   if ( masterproc ) then

      if ( save3Dbin ) then

         write(46) nx, ny, nzm, nsubdomains, &
                   nsubdomains_x, nsubdomains_y, nfields
         do k = 1, nzm
            write(46) z(k)
         enddo
         do k = 1, nzm
            write(46) pres(k)
         enddo
         write(46) dx
         write(46) dy
         write(46) nstep*dt/(3600.*24.)+day0

      else
      
         write(long_name,'(8i4)') nx, ny, nzm, nsubdomains, &
                                  nsubdomains_x, nsubdomains_y, nfields
         do k=1,nzm
            write(c_z(k),'(f12.3)') z(k)
         enddo
         do k=1,nzm
            write(c_p(k),'(f12.3)') pres(k)
         enddo
         write(c_dx,'(f12.0)') dx
         write(c_dy,'(f12.0)') dy
         write(c_time,'(f12.5)') nstep*dt/(3600.*24.)+day0
        
         write(46) long_name(1:32)
         write(46) c_time, c_dx, c_dy, (c_z(k),k=1,nzm), (c_p(k),k=1,nzm)

      endif ! save3Dbin

   endif ! masterproc

endif ! masterproc.or.output_sep

!--------------------------------------
! Output 3D fields
!--------------------------------------
do i = 1, nx, 1
   do j = 1, ny, 1
      do k = 1, nzm, 1

         ! Calculate rt
         rt(i,j,k) = qv(i,j,k) + qcl(i,j,k)

         ! Calculate thetal
         thl(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                0.0, 0.0, prespot(k) )

      enddo ! k = 1, nzm, 1
   enddo ! j = 1, ny, 1
enddo ! i = 1, nx, 1

! Calculate the values of chi and eta.
call compute_chi_eta( thl, rt, pres, prespot, &
                      chi, eta )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = 0.5 * ( w(i,j,k) + w(i,j,k+1) )
      enddo
   enddo
enddo
name = 'W'
long_name = 'Vertical Velocity'
units = 'm/s'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = rt(i,j,k)
      enddo
   enddo
enddo
name = 'RT'
long_name = 'Total water mixing ratio (vapor+cloud)'
units = 'kg/kg'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = thl(i,j,k)
      enddo
   enddo
enddo
name = 'THL'
long_name = 'Liquid water potential temperature'
units = 'K'
call compress3D( tmp, nx, ny, nzm, name, long_name,units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = chi(i,j,k)
      enddo
   enddo
enddo
name = 'CHI'
long_name = 'Extended liquid water mixing ratio, chi'
units = 'kg/kg'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = eta(i,j,k)
      enddo
   enddo
enddo
name = 'ETA'
long_name = 'Eta (orthogonal to chi in PDF trans.)'
units = 'kg/kg'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = qpl(i,j,k)
      enddo
   enddo
enddo
name = 'RR'
long_name = 'Rain water mixing ratio'
units = 'kg/kg'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = conp(i,j,k)
      enddo
   enddo
enddo
name = 'NR'
long_name = 'Rain drop concentration'
units = 'num/kg'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = qcl(i,j,k)
      enddo
   enddo
enddo
name = 'RC'
long_name = 'Cloud water mixing ratio'
units = 'kg/kg'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = qp_evap(i,j,k)
      enddo
   enddo
enddo
name = 'RR_EVAP'
long_name = 'Rain water evaporation rate'
units = 'kg/kg/s'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = qp_auto(i,j,k)
      enddo
   enddo
enddo
name = 'RR_AUTO'
long_name = 'Rain water autoconversion rate'
units = 'kg/kg/s'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = qp_accr(i,j,k)
      enddo
   enddo
enddo
name = 'RR_ACCR'
long_name = 'Rain water accretion rate'
units = 'kg/kg/s'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

call task_barrier()

if (nfields.ne.nfields1) then
   if (masterproc) print *, 'write_3D_micro_fields error:  nfields = ', &
                            nfields, '; nfields1 = ', nfields1
   call task_abort()
endif
if ( masterproc ) then
   close (46)
   if ( RUN3D .or. save3Dsep ) then
      if (dogzip3D) call systemf('gzip -f '//filename)
      print *, 'Writing 3D data. file:  '//filename
   else
      print *, 'Appending 3D data. file:  '//filename
   endif
endif

call t_stopf('3D_out')


return

end subroutine write_3d_micro_fields


#endif /* UWM_STATS */
end module microphysics



