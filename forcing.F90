
subroutine forcing
	
use vars
use params
use microphysics, only: micro_field, index_water_vapor, total_water
#ifdef PNNL_STATS
use microphysics, only: nmicro_fields, flag_wmass 
use calc_vars_util, only: t2thetal
#endif /*PNNL_STATS*/
use simple_ocean, only: sst_evolve

implicit none


integer i,j,k,n,nn,m,iz,iday0,iday
real coef, radtend, dayy
real tt(nzm,2),qq(nzm,2),uu(nzm,2),vv(nzm,2),ww(nzm,2)
real ratio1, ratio2, ratio_t1, ratio_t2
logical zgrid
#ifdef UWM_MISC
real qneg,qpoz, factor
#endif
#ifdef UWM_STATS
real t1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real q1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real t_avg(nzm), t1_avg(nzm)
real q_avg(nzm), q1_avg(nzm)
#ifdef PNNL_STATS
real thel1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real thel_avg(nzm), thel1_avg(nzm)
real qtog1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real qtog_avg(nzm), qtog1_avg(nzm)
#endif /*PNNL_STATS*/
#endif /*UWM_STATS*/


call t_startf ('forcing')


! if doseasons=.false. do perpetual forcing

 if(doseasons) then
   dayy = day
 else
   iday0 = day0
   iday = day
   dayy = day-iday
   dayy = iday0 + dayy
 end if


! ---------------------------------------------------------------
! Large-scale sounding:

  nn=1
  do i=1,nsnd-1
   if(day.gt.daysnd(i)) then
     nn=i
   endif
  end do

  do n=1,2
  
    m = nn+n-1
    if(zsnd(2,m).gt.zsnd(1,m)) then
      zgrid=.true.
    else if(psnd(2,m).lt.psnd(1,m)) then
      zgrid=.false.
    else
#ifdef UWM_MISC
     if(masterproc) then
        print*,'error in grid in snd'
        ! Added by dschanen UWM June 25 2008
        print *, "p= ", psnd
        print *, "z= ", zsnd
      end if
#else
      if(masterproc) print*,'error in grid in snd'
#endif      
      stop
    end if
    do iz = 1,nzm
      if(zgrid) then
        do i = 2,nzsnd
          if(z(iz).le.zsnd(i,m)) then
           coef = (z(iz)-zsnd(i-1,m))/(zsnd(i,m)-zsnd(i-1,m)) 
           tt(iz,n)=tsnd(i-1,m)+(tsnd(i,m)-tsnd(i-1,m))*coef
           tt(iz,n) = tt(iz,n) / prespot(iz)
           qq(iz,n)=qsnd(i-1,m)+(qsnd(i,m)-qsnd(i-1,m))*coef
           uu(iz,n)=usnd(i-1,m)+(usnd(i,m)-usnd(i-1,m))*coef
           vv(iz,n)=vsnd(i-1,m)+(vsnd(i,m)-vsnd(i-1,m))*coef
           goto 11
          endif
        end do
      else
        do i = 2,nzsnd
          if(pres(iz).ge.psnd(i,m)) then
           coef = (pres(iz)-psnd(i-1,m))/(psnd(i,m)-psnd(i-1,m))
           tt(iz,n)=tsnd(i-1,m)+(tsnd(i,m)-tsnd(i-1,m))*coef
           tt(iz,n) = tt(iz,n) /prespot(iz)
           qq(iz,n)=qsnd(i-1,m)+(qsnd(i,m)-qsnd(i-1,m))*coef
           uu(iz,n)=usnd(i-1,m)+(usnd(i,m)-usnd(i-1,m))*coef
           vv(iz,n)=vsnd(i-1,m)+(vsnd(i,m)-vsnd(i-1,m))*coef
           goto 11
          endif
        end do
      end if

      call atmosphere(z(iz-1)/1000.,ratio1,ratio2,ratio_t1)
      call atmosphere(z(iz)/1000.,ratio1,ratio2,ratio_t2)

      tt(iz,n)=ratio_t2/ratio_t1*tt(iz-1,n)
!      qq(iz,n)=max(0.,2.*qq(iz-1,n)-qq(iz-2,n))
      qq(iz,n) = qq(iz-1,n)*exp(-(z(iz)-z(iz-1))/3000.)
      uu(iz,n)=uu(iz-1,n)
      vv(iz,n)=vv(iz-1,n)

   11 continue

    end do ! iz 

  end do ! n

  coef=(day-daysnd(nn))/(daysnd(nn+1)-daysnd(nn))
  do k=1,nzm
    tg0(k)=tt(k,1)+(tt(k,2)-tt(k,1))*coef
    qg0(k)=qq(k,1)+(qq(k,2)-qq(k,1))*coef
    qg0(k)=qg0(k)*1.e-3
! Note that ug0 and vg0 maybe reset if dolargescale is true)
    ug0(k)=uu(k,1)+(uu(k,2)-uu(k,1))*coef - ug
    vg0(k)=vv(k,1)+(vv(k,2)-vv(k,1))*coef - vg
   end do

! ---------------------------------------------------------------
! Initialize tendencies:


do k=1,nzm
 ttend(k)=0.
#ifdef PNNL_STATS
 theltend(k)=0.
 qtogtend(k)=0.
 qtheltend(k)=0.
#endif /*PNNL_STATS*/
 qtend(k)=0.
end do


! Large-Scale Advection Forcing:


if(dolargescale.and.time.gt.timelargescale) then

  nn=1
  do i=1,nlsf-1
   if(day.gt.dayls(i)) then
     nn=i
   endif
  end do

  do n=1,2

    m = nn+n-1

    if(zls(2,m).gt.zls(1,m)) then
     zgrid=.true.
    else if(pls(2,m).lt.pls(1,m)) then
     zgrid=.false.
    else
     if(masterproc) print*,'error in grid in lsf'
     stop
    end if
    do iz = 1,nzm
     if(zgrid) then
      do i = 2,nzlsf
       if(z(iz).le.zls(i,m)) then
         coef = (z(iz)-zls(i-1,m))/(zls(i,m)-zls(i-1,m))
         tt(iz,n)=dtls(i-1,m)+(dtls(i,m)-dtls(i-1,m))*coef
         qq(iz,n)=dqls(i-1,m)+(dqls(i,m)-dqls(i-1,m))*coef
         uu(iz,n)=ugls(i-1,m)+(ugls(i,m)-ugls(i-1,m))*coef
         vv(iz,n)=vgls(i-1,m)+(vgls(i,m)-vgls(i-1,m))*coef
         ww(iz,n)=wgls(i-1,m)+(wgls(i,m)-wgls(i-1,m))*coef
         goto 12
       endif
      end do
     else
      do i = 2,nzlsf
       if(pres(iz).ge.pls(i,m)) then
         coef = (pres(iz)-pls(i-1,m))/(pls(i,m)-pls(i-1,m))
         tt(iz,n)=dtls(i-1,m)+(dtls(i,m)-dtls(i-1,m))*coef
         qq(iz,n)=dqls(i-1,m)+(dqls(i,m)-dqls(i-1,m))*coef
         uu(iz,n)=ugls(i-1,m)+(ugls(i,m)-ugls(i-1,m))*coef
         vv(iz,n)=vgls(i-1,m)+(vgls(i,m)-vgls(i-1,m))*coef
         ww(iz,n)=wgls(i-1,m)+(wgls(i,m)-wgls(i-1,m))*coef
         goto 12
       endif
      end do
     end if
     tt(iz,n)=0.
     qq(iz,n)=0.
     uu(iz,n)=uu(iz-1,n)
     vv(iz,n)=vv(iz-1,n)
     ww(iz,n)=0.
  12 continue

    end do

   end do ! n

#ifdef UWM_STATS
   if (dostatis) then
      do k=1,nzm
         do j=1,ny
            do i=1,nx
               t1(i,j,k) = t(i,j,k)
               q1(i,j,k) = micro_field(i,j,k,index_water_vapor)
#ifdef PNNL_STATS
               thel1(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                        qci(i,j,k), qpi(i,j,k), prespot(k) )
               qtog1(i,j,k) = 0.0
               do n = 1, nmicro_fields ! prevents accessing unreferenced memory 
                  qtog1(i,j,k) = qtog1(i,j,k) + flag_wmass(n)*micro_field(i,j,k,n)
               end do
#endif /*PNNL_STATS*/
               misc(i,j,k) = w(i,j,k)
            enddo
         enddo
      enddo ! k = 1, nzm
   endif
#endif /*UWM_STATS*/
   coef=(day-dayls(nn))/(dayls(nn+1)-dayls(nn))
   dosubsidence = .false.
   do k=1,nzm
    ttend(k)=tt(k,1)+(tt(k,2)-tt(k,1))*coef
    qtend(k)=qq(k,1)+(qq(k,2)-qq(k,1))*coef
    ug0(k)=uu(k,1)+(uu(k,2)-uu(k,1))*coef - ug
    vg0(k)=vv(k,1)+(vv(k,2)-vv(k,1))*coef - vg
    wsub(k)=ww(k,1)+(ww(k,2)-ww(k,1))*coef
    dosubsidence = dosubsidence .or. wsub(k).ne.0.
    do j=1,ny
     do i=1,nx
      t(i,j,k)=t(i,j,k)+ttend(k) * dtn
#ifdef UWM_MISC
      if ( do_spec_hum_to_mix_rat ) then
         ! The specified forcings were given in terms of specific humidity.
         ! They need to be put in terms of mixing ratio.
         micro_field(i,j,k,index_water_vapor) &
         = micro_field(i,j,k,index_water_vapor) &
           + ( 1.0 + micro_field(i,j,k,index_water_vapor) )**2 * qtend(k) * dtn
      else ! already in mixing ratio
         micro_field(i,j,k,index_water_vapor) = &
              micro_field(i,j,k,index_water_vapor) + qtend(k) * dtn
      endif ! do_spec_hum_to_mix_rat
#else
      micro_field(i,j,k,index_water_vapor) = &
                 max(0.,micro_field(i,j,k,index_water_vapor) + qtend(k) * dtn)
#endif /* UWM_MISC */
#ifdef PNNL_STATS
      thel(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                              qci(i,j,k), qpi(i,j,k), prespot(k) )
      qtog(i,j,k) = 0.0
      do n = 1, nmicro_fields ! prevents accessing unreferenced memory 
          qtog(i,j,k) = qtog(i,j,k) + flag_wmass(n)*micro_field(i,j,k,n)
      end do

      theltend(k) = theltend(k) + (thel(i,j,k)-thel1(i,j,k))/dtn/float(nx*ny)
      qtogtend(k) = qtogtend(k) + (qtog(i,j,k)-qtog1(i,j,k))/dtn/float(nx*ny)
#endif /*PNNL_STATS*/
     end do
    end do
   end do 

#ifdef PNNL_STATS
    ! Observed Large-Scale tendency to the covariance r'thl'
    call averageXY_MPI(micro_field(:,:,:,index_water_vapor),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,q_avg)
    call averageXY_MPI(q1,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,q1_avg)
    call averageXY_MPI(thel,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,thel_avg)
    call averageXY_MPI(thel1,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,thel1_avg)

    do k=1,nzm
       qtheltend(k)=0.
       do j=1,ny
         do i=1,nx
            qtheltend(k) = qtheltend(k) + (thel(i,j,k)-thel_avg(k))*(micro_field(i,j,k,index_water_vapor)-q_avg(k)) &
                          - (thel1(i,j,k)-thel1_avg(k))*(q1(i,j,k)-q1_avg(k))
         end do
       end do
       qtheltend(k)=qtheltend(k)/dtn/float(nx*ny)
    end do
#endif /*PNNL_STATS*/
#ifdef UWM_MISC

! fix water conservation here (Minghuai.Wang@pnnl.gov)
   do j=1, ny
    do i=1, nx
      qneg = 0.
      qpoz = 0.
      do k=1, nzm
        if(micro_field(i,j,k,index_water_vapor).lt.0.) then            
          qneg = qneg + micro_field(i,j,k,index_water_vapor)*adz(k)*dz*rho(k)
        else
          qpoz = qpoz + micro_field(i,j,k,index_water_vapor)*adz(k)*dz*rho(k)
        end if
      enddo
      if(qneg == 0 .and. qpoz == 0) then
        factor=1
      else
        factor = 1.+qneg/qpoz
      endif
      do k=1, nzm
        micro_field(i,j,k,index_water_vapor) = max(0.,micro_field(i,j,k,index_water_vapor)*factor)
      enddo 
    enddo
   enddo
#endif
   pres0 = pres0ls(nn)+(pres0ls(nn+1)-pres0ls(nn))*coef

   if(wgls_holds_omega) then
      ! convert omega (sitting in wsub) into large-scale vertical velocity.
      ! Note that omega was read in from SCAM IOP netcdf input file.
      do k = 1,nzm
         wsub(k) = -wsub(k)/rho(k)/ggr
      end do
   end if

   if(dosubsidence) call subsidence()

#ifdef UWM_STATS
   if (dostatis) then
      call stat_varscalar(t,t1,t_avg,t1_avg,t2leforc)
      call setvalue(twleforc,nzm,0.)
      call stat_sw2(t,t1,twleforc)
      call stat_varscalar(micro_field(:,:,:,index_water_vapor),q1,q_avg,q1_avg,q2leforc)
      call setvalue(qwleforc,nzm,0.)
      call stat_sw2(micro_field(:,:,:,index_water_vapor),q1,qwleforc)
#ifdef PNNL_STATS
     !MWSWong: THL and QTOG budgets
      call stat_varscalar(thel,thel1,thel_avg,thel1_avg,thel2leforc)
      call setvalue(thelwleforc,nzm,0.)
      call stat_sw2(thel,thel1,thelwleforc)
      call stat_varscalar(qtog,qtog1,qtog_avg,qtog1_avg,qtog2leforc)
      call setvalue(qtogwleforc,nzm,0.)
      call stat_sw2(qtog,qtog1,qtogwleforc)

     !MWSWong: budget for covariance r'thl'
      do k=1,nzm
        qthelleforc(k)=0.
        do j=1,ny
          do i=1,nx
             qthelleforc(k)=qthelleforc(k) + (thel(i,j,k)-thel_avg(k))*(micro_field(i,j,k,index_water_vapor)-q_avg(k)) &
                             - (thel1(i,j,k)-thel1_avg(k))*(q1(i,j,k)-q1_avg(k))
          end do
        end do
        qthelleforc(k)=qthelleforc(k)*(1./(dtn*nx*ny))
     end do

#endif /*PNNL_STATS*/
      ! Reset variable misc back to 0.
      misc(:,:,:) = 0.0
   endif

#endif /*UWM_STATS*/
end if 

!---------------------------------------------------------------------
! Prescribed Radiation Forcing:


if(doradforcing.and.time.gt.timelargescale) then

  nn=1
  do i=1,nrfc-1
   if(day.gt.dayrfc(i)) then
     nn=i
   endif
  end do

  do n=1,2

    m = nn+n-1 

    if(prfc(2,m).gt.prfc(1,m)) then
     zgrid=.true.
    else if(prfc(2,m).lt.prfc(1,m)) then
     zgrid=.false.
    else
     if(masterproc) print*,'error in grid in rad'
     stop
    end if
    do iz = 1,nzm
     if(zgrid) then
      do i = 2,nzrfc
       if(z(iz).le.prfc(i,m)) then
        tt(iz,n)=dtrfc(i-1,m)+(dtrfc(i,m)-dtrfc(i-1,m))/(prfc(i,m)-prfc(i-1,m)) &
                                                           *(z(iz)-prfc(i-1,m))
        goto 13
       endif
      end do
     else
      do i = 2,nzrfc
       if(pres(iz).ge.prfc(i,m)) then
        tt(iz,n)=dtrfc(i-1,m)+(dtrfc(i,m)-dtrfc(i-1,m))/(prfc(i,m)-prfc(i-1,m)) &
                                                           *(pres(iz)-prfc(i-1,m))
        goto 13
       endif
      end do
     end if
     tt(iz,n)=0.
  13 continue
    end do

  end do ! n

#ifdef UWM_STATS
  if (dostatis) then
     do k=1,nzm
        do j=1,ny
           do i=1,nx
              t1(i,j,k) = t(i,j,k)
#ifdef PNNL_STATS
!MWSWong:budget for THL 
              thel1(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                       qci(i,j,k), qpi(i,j,k), prespot(k) )
#endif /*PNNL_STATS*/
              misc(i,j,k) = w(i,j,k)
           enddo
        enddo
     enddo
  endif

#endif /*UWM_STATS*/
  coef=(day-dayrfc(nn))/(dayrfc(nn+1)-dayrfc(nn))
  do k=1,nzm
    radtend=tt(k,1)+(tt(k,2)-tt(k,1))*coef
    radqrlw(k)=radtend*float(nx*ny)
    radqrsw(k)=0.
    do j=1,ny
     do i=1,nx
       t(i,j,k)=t(i,j,k)+radtend*dtn 
     end do
    end do
  end do

#ifdef PNNL_STATS
!MWSWong:budget for THL 
  if (dostatis) then
     do k=1,nzm
        do j=1,ny
           do i=1,nx
              thel(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                      qci(i,j,k), qpi(i,j,k), prespot(k) )
           enddo
        enddo
     enddo
  endif

#endif /*PNNL_STATS*/
#ifdef UWM_STATS

   if (dostatis) then
      call stat_varscalar(t,t1,t_avg,t1_avg,t2lerad)
      call setvalue(twlerad,nzm,0.)
      call stat_sw2(t,t1,twlerad)
#ifdef PNNL_STATS
!MWSWong:add budget term(s) for THL
      call stat_varscalar(thel,thel1,thel_avg,thel1_avg,thel2lerad)
      call setvalue(thelwlerad,nzm,0.)
      call stat_sw2(thel,thel1,thelwlerad)
#endif /*PNNL_STATS*/
      ! Reset variable misc back to 0.
      misc(:,:,:) = 0.0
   endif
#endif /*UWM_STATS*/
endif


!----------------------------------------------------------------------------
! Surface flux forcing:

if(dosfcforcing.and.time.gt.timelargescale) then

   nn=1
   do i=1,nsfc-1
     if(day.gt.daysfc(i)) then
       nn=i
     endif
   end do

   coef=(day-daysfc(nn))/(daysfc(nn+1)-daysfc(nn))
   tabs_s=sstsfc(nn)+(sstsfc(nn+1)-sstsfc(nn))*coef
   fluxt0=(shsfc(nn)+(shsfc(nn+1)-shsfc(nn))*coef)/(rhow(1)*cp)
   fluxq0=(lhsfc(nn)+(lhsfc(nn+1)-lhsfc(nn))*coef)/(rhow(1)*lcond)
   tau0=tausfc(nn)+(tausfc(nn+1)-tausfc(nn))*coef

   do j=1,ny
    do i=1,nx
      sstxy(i,j) = tabs_s - t00
    end do
   end do

   if(dostatis) then
     sstobs = tabs_s  ! sst is not averaged over the sampling period
     lhobs = lhobs + fluxq0 * rhow(1)*lcond
     shobs = shobs + fluxt0 * rhow(1)*cp
   end if

endif

!-------------------------------------------------------------------------------

if(.not.dosfcforcing.and.dodynamicocean) then

   call sst_evolve()

endif

!-------------------------------------------------------------------------------
call t_stopf ('forcing')


end subroutine forcing
