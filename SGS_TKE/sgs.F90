module sgs

! module for original SAM subgrid-scale SGS closure (Smagorinsky or 1st-order TKE)
! Marat Khairoutdinov, 2012

use grid, only: nx,nxp1,ny,nyp1,YES3D,nzm,nz,dimx1_s,dimx2_s,dimy1_s,dimy2_s 
use params, only: dosgs
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

! --- Heng Xiao, 03/13/2024
! adding a switch for isotropic grdf
logical :: istropic_grdf = .false.
! --- Heng Xiao, 03/13/2024

logical:: dosmagor   ! if true, then use Smagorinsky closure

! Local diagnostics:

real tkesbbuoy(nz), tkesbshear(nz),tkesbdiss(nz), tkesbdiff(nz)

CONTAINS

! required microphysics subroutines and function:
!----------------------------------------------------------------------
!!! Read microphysics options from prm (namelist) file

subroutine sgs_setparm()

  use grid, only: case
  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder

  !======================================================================
  NAMELIST /SGS_TKE/ &
       dosmagor ! Diagnostic Smagorinsky closure

  NAMELIST /BNCUIODSBJCB/ place_holder

  dosmagor = .true. ! default 

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')

  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  rewind(55) !note that one must rewind before searching for new namelists

  read (55,SGS_TKE,IOSTAT=ios)

  advect_sgs = .not.dosmagor

  if (ios.ne.0) then
     !namelist error checking
     if(ios.ne.ios_missing_namelist) then
        write(*,*) '****** ERROR: bad specification in SGS_TKE namelist'
        call task_abort()
     end if
  end if
  close(55)

  ! END UW ADDITION
  !======================================================================

end subroutine sgs_setparm

!----------------------------------------------------------------------
!!! Initialize sgs:


subroutine sgs_init()

  use grid, only: nrestart, dx, dy, dz, adz, masterproc
  use params, only: LES
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
  end if

  if(LES) then
    do k=1,nzm
       ! Heng Xiao, 02/28/2024, 03/13/2024
       ! adding a switch for isotropic grdf
       if (istropic_grdf) then
         grdf_x(k) = 1.
         grdf_y(k) = 1.
       else
         grdf_x(k) = dx**2/(adz(k)*dz)**2
         grdf_y(k) = dy**2/(adz(k)*dz)**2
       end if
       ! Heng Xiao, 02/28/2024, 03/13/2024
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


end subroutine sgs_init

!----------------------------------------------------------------------
!!! make some initial noise in sgs:
!
subroutine setperturb_sgs(ptype)

use vars, only: q0, z
integer, intent(in) :: ptype
integer i,j,k

select case (ptype)

  case(0)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(k.le.4.and..not.dosmagor) then
            tke(i,j,k)=0.04*(5-k)
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
     YES3D*0.5*tkhmax(k)*grdf_y(k)*dt/dy**2)
end do

end subroutine kurant_sgs


!----------------------------------------------------------------------
!!! compute sgs diffusion of momentum:
!
subroutine sgs_mom()

   call diffuse_mom()

end subroutine sgs_mom

!----------------------------------------------------------------------
!!! compute sgs diffusion of scalars:
!
subroutine sgs_scalars()

  use vars
  use microphysics
  use tracers
  ! Heng Xiao: added for diagnosing SGS <thel'w'> and <qto'w'>
#ifdef PNNL_STATS
  use params, only: dotracers, rgas, cp
#else
  use params, only: dotracers
#endif
#ifdef PNNL_STATS
  use calc_vars_util, only: t2thetal
#endif
  implicit none

    real dummy(nz)
    real fluxbtmp(nx,ny), fluxttmp(nx,ny) !bloss
    integer k


#ifdef PNNL_STATS
    !MWSWong:convert HL to THL
    real thelprev(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real thelcurr(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real f0(nzm), df0(nzm)
    real r2dx,r2dy,r2dx0,r2dy0,r2dz
    integer i,j,kb,kc,jb,jc,n
    real factor_xy
    real q1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real qtogprev(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real qtogcurr(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real f0_qtog(nzm), df0_qtog(nzm)
    real f0_q0(nzm), df0_q0(nzm)


      if(dostatis) then
         do k = 1, nzm
            do j = 1, ny
               do i = 1, nx
                  thelprev(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                              qci(i,j,k), qpi(i,j,k), &
                                              prespot(k) )
                  qtogprev(i,j,k) = 0.0
                  do n = 1, nmicro_fields ! prevents accessing unreferenced memory 
                     qtogprev(i,j,k) = qtogprev(i,j,k) &
                                      + flag_wmass(n) * micro_field(i,j,k,n)
                  enddo
                  q1(i,j,k) = micro_field(i,j,k,index_water_vapor)
               enddo
            enddo
         enddo
      endif
#endif /*PNNL_STATS*/
      call diffuse_scalar(t,fluxbt,fluxtt,tdiff,twsb, &
                           t2lediff,t2lediss,twlediff,.true.)
    
      if(advect_sgs) then
         call diffuse_scalar(tke,fzero,fzero,dummy,sgswsb, &
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
#ifdef UWM_STATS
           if ( k == index_water_vapor ) then
              call diffuse_scalar(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
                   mkdiff(:,k),mkwsb(:,k), q2lediff,q2lediss,qwlediff,.true.)
           else
              call diffuse_scalar(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
                   mkdiff(:,k),mkwsb(:,k), dummy,dummy,dummy,.false.)
           endif
#else
           call diffuse_scalar(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
                mkdiff(:,k),mkwsb(:,k), dummy,dummy,dummy,.false.)
#endif /*UWM_STATS*/
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
          call diffuse_scalar(tracer(:,:,:,k),fluxbtmp,fluxttmp, &
               trdiff(:,k),trwsb(:,k), &
               dummy,dummy,dummy,.false.)
!!$          call diffuse_scalar(tracer(:,:,:,k),fluxbtr(:,:,k),fluxttr(:,:,k),trdiff(:,k),trwsb(:,k), &
!!$                           dummy,dummy,dummy,.false.)

        end do

      end if
#ifdef PNNL_STATS 
      ! MWSWong:THL and QTOG budget terms
      if(dostatis) then
         do k = 1, nzm
            theldiff(k)=0.
            qtogdiff(k)=0.
            thelwlediff(k)=0.
            qtogwlediff(k)=0.
            do j = 1, ny
               do i = 1, nx
                  thelcurr(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                              qci(i,j,k), qpi(i,j,k), &
                                              prespot(k) )
                  qtogcurr(i,j,k) = 0.0
                  do n = 1, nmicro_fields ! prevents accessing unreferenced memory 
                     qtogcurr(i,j,k) = qtogcurr(i,j,k) &
                                       + flag_wmass(n) * micro_field(i,j,k,n)
                  enddo
                  theldiff(k) = theldiff(k)+( thelcurr(i,j,k)-thelprev(i,j,k))
                  qtogdiff(k) = qtogdiff(k)+( qtogcurr(i,j,k)-qtogprev(i,j,k))
               enddo
            enddo
         enddo


         call stat_varscalar(thelcurr,thelprev,f0,df0,thel2lediff)
         call setvalue(thelwlediff,nzm,0.)
         call stat_sw2(thelcurr,thelprev,thelwlediff)

         call stat_varscalar(qtogcurr,qtogprev,f0_qtog,df0_qtog,qtog2lediff)
         call setvalue(qtogwlediff,nzm,0.)
         call stat_sw2(qtogcurr,qtogprev,qtogwlediff)

         call averageXY_MPI(micro_field(:,:,:,index_water_vapor), &
                            dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,f0_q0)
         call averageXY_MPI(q1,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,df0_q0)

         do k = 1, nzm
            qthellediff(k)=0.
            do j=1,ny
               do i=1,nx
                  qthellediff(k) &
                  = qthellediff(k) &
                    + ( thelcurr(i,j,k) - f0(k) ) &
                      * ( micro_field(i,j,k,index_water_vapor) - f0_q0(k) ) &
                    - ( thelprev(i,j,k) - df0(k) ) * ( q1(i,j,k) - df0_q0(k) )
               enddo
            enddo
            qthellediff(k)=qthellediff(k)*(1./(dtn*nx*ny))
         enddo


         factor_xy=1./float(nx*ny)
         r2dx0=1./(2.*dx)
         r2dy0=1./(2.*dy)
         do k = 1, nzm
            thel2lediss(k)=0.
            qtog2lediss(k)=0.
            qthellediss(k)=0.
            kc=min(nzm,k+1)
            kb=max(1,k-1)
            r2dz=2./((kc-kb)*(adzw(k+1)+adzw(k))*dz)
            r2dx=r2dx0*sqrt((kc-kb)*dx*r2dz) ! grid anisotropy correction
            r2dy=r2dy0*sqrt((kc-kb)*dx*r2dz)
            !Heng Xiao: for diagnosing SGS <thel'w'>
            thelws(k) = 0.
            qtows(k) = 0.
            do j = 1, ny
               jc=j+YES3D
               jb=j-YES3D
               do i = 1, nx
                  thel2lediss(k) &
                  = thel2lediss(k) &
                    - tkh(i,j,k) * factor_xy &
                      * ( ( ( thelcurr(i,j,kc) - f0(kc) &
                              - thelcurr(i,j,kb) + f0(kb) ) * r2dz )**2 ) !&
                  ! - tkh(i,j,k)*factor_xy*( &
                  !      ((thelcurr(i+1,j,k)-thelcurr(i-1,j,k))*r2dx)**2+ &
                  !      ((thelcurr(i,jc,k)-thelcurr(i,jb,k))*r2dy)**2 )
                  qtog2lediss(k) &
                  = qtog2lediss(k) &
                    - tkh(i,j,k) * factor_xy &
                      * ( ( ( qtogcurr(i,j,kc) - f0_qtog(kc) &
                              - qtogcurr(i,j,kb) + f0_qtog(kb) ) * r2dz )**2 ) !&
                  ! - tkh(i,j,k)*factor_xy*( &
                  !   ((qtogcurr(i+1,j,k)-qtogcurr(i-1,j,k))*r2dx)**2+ &
                  !   ((qtogcurr(i,jc,k)-qtogcurr(i,jb,k))*r2dy)**2 )
                  qthellediss(k) &
                  = qthellediss(k) &
                    - tkh(i,j,k) * ( factor_xy &
                      * ( thelcurr(i,j,kc) - f0(kc) &
                          - thelcurr(i,j,kb) + f0(kb) ) &
                      * ( micro_field(i,j,kc,index_water_vapor) - f0_q0(kc) &
                          - micro_field(i,j,kb,index_water_vapor) + f0_q0(kb) ) * r2dz**2 ) !&
                  ! - tkh(i,j,k)*factor_xy*( &
                  !   ((thelcurr(i+1,j,k)-thelcurr(i-1,j,k))*(micro_field(i+1,j,k,index_water_vapor)-micro_field(i-1,j,k,index_water_vapor))
                  !   )*r2dx**2+ &
                  !   ((thelcurr(i,jc,k)-thelcurr(i,jb,k))*(micro_field(i,jc,k,index_water_vapor)-micro_field(i,jb,k,index_water_vapor))
                  !   )*r2dy**2 )

                  !Heng Xiao: diagnosing the SGS <thel'w'>
                  !Note that the flux is on the interface levels not the mid-layer levels.
                  !The same is true for QTOFLX/QTOFLXS but it seems that the statistics
                  !output put them at the mid-layer levels (95% sure).
                  !The other fluxes (and other w related variables) calculated in statistics.F90
                  !are, however, interpolated to the mid-layer levels.
                  if (k .eq. 1) then
                    thelws(k) = thelws(k) + fluxbt(i,j)*(1000./pres0)**(rgas/cp)
                    qtows(k) = qtows(k) + 1.0e3*fluxbq(i,j)*(1000./pres0)**(rgas/cp) ! convert to g/kg
                  else
                    thelws(k) = thelws(k) - 0.5*(tkh(i,j,k)+tkh(i,j,k-1))* &
                                            (thelprev(i,j,k)-thelprev(i,j,k-1)) &
                                            /(adzw(k)*dz)
                    qtows(k) = qtows(k) - 0.5*(tkh(i,j,k)+tkh(i,j,k-1))* &
                                          (q1(i,j,k)-q1(i,j,k-1))*1.0e3 &
                                          /(adzw(k)*dz)
                  end if
               enddo
            enddo
            !thel2lediss(k)=thel2lediss(k)*2.*factor_xy
            !qtog2lediss(k)=qtog2lediss(k)*2.*factor_xy
            !qthellediss(k)=qthellediss(k)*2.*factor_xy
            thel2lediss(k)=thel2lediss(k)*2.
            qtog2lediss(k)=qtog2lediss(k)*2.
            qthellediss(k)=qthellediss(k)*2.

            !write(6,*) 'SGS_TKE/sgs: thel2diss(k) t2lediss(k) ', thel2lediss(k), t2lediss(k), 'k', k, &
!                                        thelcurr(1,1,k), t(1,1,k) 
!           write(6,*) 'SGS_TKE/sgs: qtog2diss(k) q2lediss(k) ', qtog2lediss(k), q2lediss(k), 'k', k, &
!                                       qtogcurr(1,1,k), micro_field(1,1,k,index_water_vapor)

         enddo
      endif
#endif /*PNNL_STATS*/



end subroutine sgs_scalars

!----------------------------------------------------------------------
!!! compute sgs processes (beyond advection):
!
subroutine sgs_proc()

   use grid, only: nstep,dt,icycle
   use params, only: dosmoke

!    SGS TKE equation:

     if(dosgs) call tke_full()

end subroutine sgs_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine sgs_diagnose()
! None 

end subroutine sgs_diagnose


!----------------------------------------------------------------------
!!!! Collect microphysics history statistics (vertical profiles)
!
subroutine sgs_statistics()
  
  use vars
  use hbuffer, only: hbuf_put, hbuf_avg_put
  use params, only : lcond

  real tmp(2), factor_xy 
  real tkz(nzm), tkhz(nzm)
  integer i,j,k,n
  character(LEN=6) :: statname  !bloss: for conditional averages

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

  call t_stopf ('sgs_statistics')

end subroutine sgs_statistics

!----------------------------------------------------------------------
! called when stepout() called

subroutine sgs_print()

 call fminmax_print('tke:',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tk:',tk,0,nxp1,1-YES3D,nyp1,nzm)
 call fminmax_print('tkh:',tkh,0,nxp1,1-YES3D,nyp1,nzm)

end subroutine sgs_print

!----------------------------------------------------------------------
!!! Initialize the list of sgs statistics 
!
subroutine sgs_hbuf_init(namelist,deflist,unitlist,status,average_type,count,sgscount)
character(*) namelist(*), deflist(*), unitlist(*)
integer status(*),average_type(*),count,sgscount

end subroutine sgs_hbuf_init


end module sgs



