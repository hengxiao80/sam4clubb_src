  
subroutine precip_proc

use vars
use microphysics
use micro_params
use params
#ifdef SILHS
use clubb_api_module, only: cm3_per_m3
#endif
#ifdef PNNL_STATS
use calc_vars_util, only: t2thetal
#endif /* PNNL_STATS */

implicit none

integer i,j,k
real auto, accr, dq, qsat
real qcc
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real f0(nzm),df0(nzm)
#ifdef PNNL_STATS
real thel1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real thel_avg(nzm),thel1_avg(nzm)
#endif /* PNNL_STATS */
      
call t_startf ('precip_proc')
     
if(dostatis) then
        
  do k=1,nzm
    do j=dimy1_s,dimy2_s
     do i=dimx1_s,dimx2_s
      df(i,j,k) = q(i,j,k)
     end do
    end do
  end do
         
#ifdef PNNL_STATS
  do k = 1, nzm
     do j = 1, ny
        do i = 1, nx
           thel1(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qp(i,j,k), &
                                    0.0, 0.0, prespot(k) )
        enddo
     enddo
  enddo

#endif /* PNNL_STATS */
endif


 
do k=1,nzm
 qpsrc(k)=0.
 qpevp(k)=0.
 do j=1,ny
  do i=1,nx	  
	  
#ifdef UWM_STATS
   ! Variables used to record statistics for microphysics process rates.
   qp_auto(i,j,k) = 0.0
   qp_accr(i,j,k) = 0.0
   qp_evap(i,j,k) = 0.0
#endif /*UWM_STATS*/
!-------     Autoconversion/accretion

   if(qn(i,j,k)+qp(i,j,k).gt.0.) then


         if(qn(i,j,k).gt.0.) then
 
           qcc = qn(i,j,k)

#ifdef SILHS
           auto = 1350.*qcc**1.47/( rho(k)*conc(i,j,k)/cm3_per_m3 )**1.79 
#else
           auto = 1350.*qcc**1.47/Nc0**1.79   ! Linearized drizzle autoconversion
                                                !(Khairoutdinov and Kogan 2000)
#endif
           accr = 67.*qcc**0.15*qp(i,j,k)**1.15 ! Linearized accretion

           qcc = qcc/(1.+dtn*(auto+accr))
           auto = auto*qcc
           accr = accr*qcc
           dq = min(qn(i,j,k),dtn*(auto+accr))
           qp(i,j,k) = qp(i,j,k) + dq
           conp(i,j,k) = conp(i,j,k) + max(0.,dq-dtn*accr)*coefconpmax ! all new drizzle drops have rd_min size
           q(i,j,k) = q(i,j,k) - dq
           qn(i,j,k) = qn(i,j,k) - dq
           qpsrc(k) = qpsrc(k) + dq
#ifdef UWM_STATS
           ! Statistics for microphysics process rates.
           qp_auto(i,j,k) = auto
           qp_accr(i,j,k) = accr
#endif /*UWM_STATS*/

         elseif(qp(i,j,k).gt.qp_threshold.and.qn(i,j,k).eq.0.) then

           qsat = qsatw(tabs(i,j,k),pres(k))
           dq = dtn * evapr1(k) * (qp(i,j,k)*conp(i,j,k)**2)**0.333 * (q(i,j,k) /qsat-1.) 
#ifdef UWM_STATS
           ! Statistics for microphysics process rates.
           qp_evap(i,j,k) = dq / dtn
#endif /*UWM_STATS*/
           dq = max(-0.5*qp(i,j,k),dq)
           conp(i,j,k) = conp(i,j,k) + dq/qp(i,j,k)*conp(i,j,k)
           qp(i,j,k) = qp(i,j,k) + dq
           q(i,j,k) = q(i,j,k) - dq
           qpevp(k) = qpevp(k) + dq

         else

           q(i,j,k) = q(i,j,k) + qp(i,j,k)
           qpevp(k) = qpevp(k) - qp(i,j,k)
           qp(i,j,k) = 0.
           conp(i,j,k) = 0.

         endif

    endif

    if(qp(i,j,k).lt.0.) then
      qp(i,j,k)=0.
      conp(i,j,k) = 0.
    end if

    conp(i,j,k) = max(qp(i,j,k)*coefconpmin, conp(i,j,k))  ! keep conp reasonable

  end do
 enddo
enddo
    


if(dostatis) then
                  
  call stat_varscalar(q,df,f0,df0,q2leprec)
  call setvalue(qwleprec,nzm,0.)
  call stat_sw2(q,df,qwleprec)

#ifdef PNNL_STATS
  do k = 1, nzm
     do j = 1, ny
        do i = 1, nx
           thel(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qp(i,j,k), &
                                   0.0, 0.0, prespot(k) )
        enddo
     enddo
  enddo

  call stat_varscalar(thel,thel1,thel_avg,thel1_avg,thel2leprec)
  call setvalue(thelwleprec,nzm,0.)
  call stat_sw2(thel,thel1,thelwleprec)

  do k = 1, nzm
     qthelleprec(k) = 0.0
     do j = 1, ny
        do i = 1, nx
           qthelleprec(k) &
           = qthelleprec(k) &
             + ( thel(i,j,k) - thel_avg(k) ) * ( q(i,j,k) - f0(k) ) &
             - ( thel1(i,j,k) - thel1_avg(k) ) * ( df(i,j,k) - df0(k) )
        enddo
     enddo
     qthelleprec(k) = qthelleprec(k) * ( 1.0 / ( dtn * nx * ny ) )
  enddo

#endif /* PNNL_STATS */
endif

call t_stopf ('precip_proc')

end subroutine precip_proc

