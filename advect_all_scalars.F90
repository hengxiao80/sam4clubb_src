subroutine advect_all_scalars()

  use vars
  use microphysics
  use sgs
  use tracers
  use params, only: dotracers
#ifdef PNNL_STATS
  use calc_vars_util, only: t2thetal
#endif /*PNNL_STATS*/
  implicit none
  real dummy(nz)
  integer k
#ifdef PNNL_STATS
  integer i,j
  real qpl_incl_bndry(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)
  real qci_incl_bndry(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)
  real qpi_incl_bndry(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)
  real thelcurr(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)
  real thelprev(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)
  real f0(nzm), df0(nzm), fff(nz)
  real coef, factor
  real qtogcurr(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)
  real qtogprev(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)
  real q1(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)
  real thel1(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)
  real qthel(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)
  real qthel1(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)
  real f0_qtog(nzm), df0_qtog(nzm), f0_q(nzm), df0_q(nzm)
  integer n

 !MWSWong:for THL budget 
     if(dostatis) then

        ! initialize thelprev and qtogprev this way to include all the
        ! information from boundary grid boxes, for thelprev and qtogprev will
        ! be advected horizontally as part of the stats code.

        if ( indx_qr > 0 ) then
           qpl_incl_bndry(:,:,:) = micro_field(:,:,:,indx_qr)
        else
           qpl_incl_bndry(:,:,:) = 0.0
        endif

        if ( indx_qi > 0 ) then
           qci_incl_bndry(:,:,:) = micro_field(:,:,:,indx_qi)
        else
           qci_incl_bndry(:,:,:) = 0.0
        endif

        if ( indx_qs > 0 .and. indx_qg > 0 ) then
           qpi_incl_bndry(:,:,:) &
           = micro_field(:,:,:,indx_qs) + micro_field(:,:,:,indx_qg)
        elseif ( indx_qs > 0 ) then
           qpi_incl_bndry(:,:,:) = micro_field(:,:,:,indx_qs)
        elseif ( indx_qg > 0 ) then
           qpi_incl_bndry(:,:,:) = micro_field(:,:,:,indx_qg)
        else
           qpi_incl_bndry(:,:,:) = 0.0
        endif

        do k = 1, nzm
           do j = dimy1_s, dimy2_s
              do i = dimx1_s, dimx2_s

                 thelprev(i,j,k) &
                 = t2thetal( t(i,j,k), gamaz(k), qpl_incl_bndry(i,j,k), &
                             qci_incl_bndry(i,j,k), qpi_incl_bndry(i,j,k), &
                             prespot(k) )

                 qtogprev(i,j,k) &
                 = micro_field(i,j,k,index_water_vapor) &
                   + qpl_incl_bndry(i,j,k) + qci_incl_bndry(i,j,k) &
                   + qpi_incl_bndry(i,j,k)

              enddo
           enddo
        enddo

        q1(:,:,:) = micro_field(:,:,:,index_water_vapor)  ! for q*thel budget
        thel1(:,:,:) = thelprev(:,:,:) !for q*thel budget

     endif
#endif /*PNNL_STATS*/


!---------------------------------------------------------
!      advection of scalars :

     call advect_scalar(t,tadv,twle,t2leadv,t2legrad,twleadv,.true.)
    
!
!    Advection of microphysics prognostics:
!

     do k = 1,nmicro_fields
        if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
         .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
         .or. doprecip.and.flag_precip(k).eq.1 ) &
#ifndef UWM_STATS /*Vanilla SAM advects all microphysics prognostics, but doesn't compute budgets for any of them */
           call advect_scalar(micro_field(:,:,:,k),mkadv(:,k),mkwle(:,k),dummy,dummy,dummy,.false.)
#else
        then ! in order to keep consistent with vanilla SAM, leave line continuation

          if (   k.eq.index_water_vapor) then
             !UWM wants to compute the total water vapor budget.
             call advect_scalar(micro_field(:,:,:,k),mkadv(:,k),mkwle(:,k),&
                              q2leadv,q2legrad,qwleadv,.true.)
          else
             !UWM does not want any other microphyscis prognostics to be included
             !in the budgets.
             call advect_scalar(micro_field(:,:,:,k),mkadv(:,k),mkwle(:,k),dummy,dummy,dummy,.false.)
          endif !index water vapor
        endif ! vanilla SAM line continuation
#endif
     end do

#ifdef UWM_MISC
     ! Update microphysical variables qv, qcl, qpl, qci, and qpi after the
     ! micro_field array has been changed due to advection.
     call micro_diagnose()
#endif /*UWM_MISC*/
#ifdef PNNL_STATS 
  ! MWSWong:convert HL to THL

     if(dostatis) then
         do k=1,nzm
           theladv(k)=0.
           qtogadv(k)=0.
           do j=1,ny
             do i=1,nx

              thel(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                          qci(i,j,k), qpi(i,j,k), prespot(k) )

              qtogcurr(i,j,k) = 0.0
              
              do n = 1, nmicro_fields ! prevents accessing unreferenced memory 
                 qtogcurr(i,j,k) = qtogcurr(i,j,k) + flag_wmass(n)*micro_field(i,j,k,n)
              end do

              theladv(k) = theladv(k) + (thel(i,j,k)-thelprev(i,j,k))
              qtogadv(k) = qtogadv(k) + (qtogcurr(i,j,k)-qtogprev(i,j,k))
              
             end do
           end do
         end do
         call stat_varscalar(thel,thelprev,f0,df0,thel2leadv)
         call stat_sw2(thel,thelprev,thelwleadv(1:nzm))
         call stat_varscalar(qtogcurr,qtogprev,f0_qtog,df0_qtog,qtog2leadv)
         call stat_sw2(qtogcurr,qtogprev,qtogwleadv(1:nzm))


         call averageXY_MPI(micro_field(:,:,:,index_water_vapor),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,f0_q)
         call averageXY_MPI(q1,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,df0_q)

        ! for qthl budget
         do k=1,nzm
            qthelleadv(k)=0.
            do j=1,ny
              do i=1,nx
                qthel(i,j,k)=(thel(i,j,k)-f0(k))*(micro_field(i,j,k,index_water_vapor)-f0_q(k))
                qthel1(i,j,k)=(thel1(i,j,k)-df0(k))*(q1(i,j,k)-df0_q(k))
                qthelleadv(k) = qthelleadv(k) + (qthel(i,j,k)-qthel1(i,j,k))/dtn/float(nx*ny)
              end do
            end do
         end do


!MWSWong: from ADV_***/scalar_advect.F90
!  Compute advection flux of variance for THL
     fff=0.
     do k=1,nzm
       do j=dimy1_s,dimy2_s
        do i=dimx1_s,dimx2_s
          thelprev(i,j,k) = (thelprev(i,j,k)-df0(k))**2
        end do
       end do
     end do
     coef = max(1.e-10,maxval(thelprev(dimx1_s:dimx2_s, dimy1_s:dimy2_s, 1:nzm))) !original
     thelprev(:,:,:) = thelprev(:,:,:) / coef
     if(RUN3D) then
       call advect_scalar3D(thelprev, u, v, w, rho, rhow, fff)
     else
       call advect_scalar2D(thelprev, u, w, rho, rhow, fff)   
     endif
#ifdef UWM_MISC  
     fff(:) = fff(:) * coef !weberjk, bugfix. fff needs to be unscaled.
#else
     thelprev(:,:,:) = thelprev(:,:,:) * coef
#endif  
     factor=dz/(nx*ny*dtn)
     do k = 1,nzm
      fff(k)=fff(k) * factor
     end do
     fff(nz)=0.
     do k = 1,nzm
      thel2legrad(k) = thel2leadv(k)
      thel2leadv(k)=-(fff(k+1)-fff(k))/(dz*adz(k)*rho(k))        
      thel2legrad(k)=thel2legrad(k)-thel2leadv(k)
     end do

!  Compute advection flux of variance for QTOG 
     fff=0.
     do k=1,nzm
       do j=dimy1_s,dimy2_s
        do i=dimx1_s,dimx2_s
          qtogprev(i,j,k) = (qtogprev(i,j,k)-df0_qtog(k))**2
        end do
       end do
     end do
     coef = max(1.e-10,maxval(qtogprev(dimx1_s:dimx2_s, dimy1_s:dimy2_s, 1:nzm)))
     qtogprev(:,:,:) = qtogprev(:,:,:) / coef

     if(RUN3D) then
       call advect_scalar3D(qtogprev, u, v, w, rho, rhow, fff)
     else
       call advect_scalar2D(qtogprev, u, w, rho, rhow, fff)      
     endif
#ifdef UWM_MISC  
     fff(:) = fff(:) * coef !weberjk, bugfix. fff needs to be unscaled.
#else
     qtogprev(:,:,:) = qtogprev(:,:,:) * coef
#endif  
     factor=dz/(nx*ny*dtn)
     do k = 1,nzm
      fff(k)=fff(k) * factor
     end do
     fff(nz)=0.
     do k = 1,nzm
      qtog2legrad(k) = qtog2leadv(k)
      qtog2leadv(k)=-(fff(k+1)-fff(k))/(dz*adz(k)*rho(k))       
      qtog2legrad(k)=qtog2legrad(k)-qtog2leadv(k)
     end do

     !  Compute advection flux of variance for Q*THL
     fff=0.
     do k=1,nzm
       do j=dimy1_s,dimy2_s
        do i=dimx1_s,dimx2_s
          qthel1(i,j,k) = (q1(i,j,k)-df0_q(k))*(thel1(i,j,k)-df0(k))
        end do
       end do
     end do
     coef = max(1.e-10,maxval(qthel1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, 1:nzm)))
     qthel1(:,:,:) = qthel1(:,:,:) / coef

     if(RUN3D) then
       call advect_scalar3D(qthel1, u, v, w, rho, rhow, fff)
     else
       call advect_scalar2D(qthel1, u, w, rho, rhow, fff)
     endif
#ifdef UWM_MISC  
     fff(:) = fff(:) * coef !weberjk, bugfix. fff needs to be unscaled.
#else
     qthel1(:,:,:) = qthel1(:,:,:) * coef
#endif  
     factor=dz/(nx*ny*dtn)
     do k = 1,nzm
      fff(k)=fff(k) * factor
     end do
     fff(nz)=0.
     do k = 1,nzm
      qthelgrad(k) = qthelleadv(k)
      qthelleadv(k)=-(fff(k+1)-fff(k))/(dz*adz(k)*rho(k))
      qthelgrad(k)=qthelgrad(k)-qthelleadv(k)
     end do

    end if   !dostat
#endif /*PNNL_STATS*/
!
!    Advection of sgs prognostics:
!

     if(dosgs.and.advect_sgs) then
       do k = 1,nsgs_fields
           call advect_scalar(sgs_field(:,:,:,k),sgsadv(:,k),sgswle(:,k),dummy,dummy,dummy,.false.)
       end do
     end if


!
!   Precipitation fallout:
!
    if(doprecip) then

       total_water_prec = total_water_prec + total_water()

       call micro_precip_fall()

       total_water_prec = total_water_prec - total_water()

#ifdef UWM_MISC
      ! Update microphysical variables qv, qcl, qpl, qci, and qpi after the
      ! micro_field array has been changed due to precipitation fallout.
      call micro_diagnose()
#endif /*UWM_MISC*/

    end if

 ! advection of tracers:

     if(dotracers) then

        do k = 1,ntracers
         call advect_scalar(tracer(:,:,:,k),tradv(:,k),trwle(:,k),dummy,dummy,dummy,.false.)
        end do

     end if

end subroutine advect_all_scalars
