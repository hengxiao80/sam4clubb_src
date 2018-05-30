
subroutine radiation()

!	Radiation interface

#ifdef UWM_STATS
use vars
#endif /*UWM_STATS*/
#ifdef PNNL_STATS
!MWSWong: THL budget
use microphysics, only: micro_field, index_water_vapor
use calc_vars_util, only: t2thetal
#endif /*PNNL_STATS*/
use grid
use params, only: dosmoke, doradsimple
#ifdef ALTOCU_RAD
use rad_simple_uwm, only: rad_simple_driver
#endif
implicit none

#ifdef UWM_STATS
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real f0(nzm),df0(nzm)
integer i,j,k
#endif /*UWM_STATS*/
#ifdef PNNL_STATS
!MWSWong: THL budget
integer n
real df_thel(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real f0_thel(nzm),df0_thel(nzm),f0_q0(nzm)
#endif /*PNNL_STATS*/

call t_startf ('radiation')
	
#ifdef UWM_STATS
if (dostatisrad) then
   do k=1,nzm
      do j=1,ny
         do i=1,nx
            df(i,j,k) = t(i,j,k)
#ifdef PNNL_STATS
!MWSWong: THL budget
            df_thel(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                       qci(i,j,k), qpi(i,j,k), prespot(k) )
#endif /*PNNL_STATS*/
         enddo
      enddo
   enddo
endif
#endif /*UWM_STATS*/

if(doradsimple) then

!  A simple predefined radiation (longwave only)

    if(dosmoke) then
       call rad_simple_smoke()
    else
#ifdef ALTOCU_RAD /* Includes analytic shortwave */
       call rad_simple_driver
#else
       call rad_simple()
#endif

    end if
	 
else


! Call full radiation package:
 

    call rad_full()	
 
endif

#ifdef UWM_STATS
if (dostatisrad) then
#ifdef PNNL_STATS
!MWSWong: THL budget
   do k=1,nzm
      do j=1,ny
         do i=1,nx
            thel(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                    qci(i,j,k), qpi(i,j,k), prespot(k) )
         end do
      end do
   end do


   call stat_varscalar(thel,df_thel,f0_thel,df0_thel,thel2lerad)
   call setvalue(thelwlerad,nzm,0.)
   call stat_sw2(thel,df_thel,thelwlerad)

   ! QTHL budget
   call averageXY_MPI(micro_field(:,:,:,index_water_vapor),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,f0_q0)
   do k=1,nzm
      qthellerad(k)=0.
      do j=1,ny
        do i=1,nx
          qthellerad(k) &
          = qthellerad(k) &
            + (thel(i,j,k)-f0_thel(k)) &
              * (micro_field(i,j,k,index_water_vapor)-f0_q0(k)) &
            - (df_thel(i,j,k)-df0_thel(k)) &
              * (micro_field(i,j,k,index_water_vapor)-f0_q0(k))
        end do
      end do 
      qthellerad(k) = qthellerad(k)*(1./(dtn*nx*ny))
   end do
#endif /*PNNL_STATS*/
   call stat_varscalar(t,df,f0,df0,t2lerad)
   call setvalue(twlerad,nzm,0.)
   call stat_sw2(t,df,twlerad)
endif
#endif /*UWM_STATS*/

call t_stopf ('radiation')

end


