#ifdef CLUBB
! $Id: main.F90 1559 2014-08-23 20:05:20Z dschanen@uwm.edu $
#endif
program crm

!       Main module.

use vars
use hbuffer
use microphysics
use sgs
#ifdef CLUBB
use sgs_params, only: doclubb, doclubbnoninter
use clubb_sgs, only: clubb_sgs_cleanup
#endif /* CLUBB */
use tracers
use movies, only: init_movies
#ifdef PNNL_STATS
use calc_vars_util, only: t2thetal
#endif /*PNNL_STATS*/
#ifdef SILHS
use clubb_api_module, only: &
  core_rknd ! Constants

use sgs_params, only: nclubb

use clubb_silhs_vars, only: &
  LH_rt, &  ! Variables
  LH_t, &
  X_nl_all_levs, &
  LH_sample_point_weights, &
  LH_t_sum_tndcy, & 
  LH_t_avg_tndcy, &
  LH_qn_sum_tndcy, & 
  LH_qn_avg_tndcy, &
  t_prior, &
  w_prior, &
  qn_prior, &
  micro_field_prior, &
  LH_micro_field_sum_tndcy, &
  LH_micro_field_avg_tndcy

use clubb_silhs_vars, only: &
  LH_microphys_type, & ! Variables
  LH_microphys_interactive, &
  LH_microphys_non_interactive, &
  LH_microphys_disabled, &
  LH_microphys_calls

use micro_field_utils, only: &
  micro_copy_in, & ! Procedure(s)
  micro_save, &
  micro_restore, &
  clip_micro_field, &
  accumulate_LH_tendency

!use microphysics, only: &
!  qn ! Variable(s)
#endif /* SILHS */
#ifdef UWM_STATS
use compute_correlation_module, only: &
    corravg_count
#endif /*UWM_STATS*/
implicit none

integer k, icyc, nn, nstatsteps

#ifdef SILHS
logical :: dostatis_save, doclubb_save
integer :: LH_iter
real :: dtn_save
#endif
#ifdef UWM_STATS

real t1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real q1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real w_prev(nx,ny,nz)
real w_curr(nx,ny,nz)
real t_avg(nzm), t1_avg(nzm)
real q_avg(nzm), q1_avg(nzm)
#endif /*UWM_STATS*/
#ifdef PNNL_STATS

!MWSWong:convert HL to THEL
real thel1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real thel_avg(nzm), thel1_avg(nzm)
real qtog1(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real qtog_avg(nzm), qtog1_avg(nzm)
#endif /*PNNL_STATS*/

#if defined(CLUBB) || defined(UWM_STATS)
!Array indicies for spurious RTM check (CLUBB)
!Array indicies (UWM_STATS)
integer :: i, j
#endif

#ifdef PNNL_STATS
integer :: n
#endif

!-------------------------------------------------------------------
! determine the rank of the current task and of the neighbour's ranks

call task_init() 
!------------------------------------------------------------------
! print time, version, etc

if(masterproc) call header()	
!------------------------------------------------------------------
! Initialize timing library.  2nd arg 0 means disable, 1 means enable

   call t_setoptionf (1, 0)
   call t_initializef ()

   call t_startf ('total')
   call t_startf ('initialize')
!------------------------------------------------------------------

call init()     ! initialize some statistics arrays
call setparm()	! set all parameters and constants

!------------------------------------------------------------------
! Initialize or restart from the save-dataset:

if(nrestart.eq.0) then
   day=day0 
   call setgrid() ! initialize vertical grid structure
   call setdata() ! initialize all variables
elseif(nrestart.eq.1) then
   call read_all()
   call setgrid() ! initialize vertical grid structure
   call diagnose()
   call sgs_init()
   call micro_init()  !initialize microphysics
elseif(nrestart.eq.2) then  ! branch run
   call read_all()
   call setgrid() ! initialize vertical grid structure
   call diagnose()
   call setparm() ! overwrite the parameters
   call sgs_init()
   call micro_init()  !initialize microphysics
   nstep = 0
   day0 = day
else
   print *,'Error: confused by value of NRESTART'
   call task_abort() 
endif

call init_movies()
call stat_2Dinit()
call tracers_init() ! initialize tracers
call setforcing()
if(masterproc) call printout()

!------------------------------------------------------------------
!  Initialize statistics buffer:

call hbuf_init()
	
!------------------------------------------------------------------
nstatis = nstat/nstatfrq
nstat = nstatis * nstatfrq
nstatsteps = 0
call t_stopf ('initialize')
!------------------------------------------------------------------
!   Main time loop    
!------------------------------------------------------------------

do while(nstep.lt.nstop.and.nelapse.gt.0) 
        
  nstep = nstep + 1
  time = time + dt
  day = day0 + nstep*dt/86400.
  nelapse = nelapse - 1
!------------------------------------------------------------------
!  Check if the dynamical time step should be decreased 
!  to handle the cases when the flow being locally linearly unstable
!------------------------------------------------------------------

  ncycle = 1

  call kurant()

  total_water_before = 0.
  total_water_after = 0.
  total_water_evap = 0.
  total_water_prec = 0.
  total_water_ls = 0.

  do icyc=1,ncycle

     icycle = icyc
     dtn = dt/ncycle
     dt3(na) = dtn
     dtfactor = dtn/dt

     if(mod(nstep,nstatis).eq.0.and.icycle.eq.ncycle) then
        nstatsteps = nstatsteps + 1
        dostatis = .true.
        if(masterproc) print *,'Collecting statistics...'
     else
        dostatis = .false.
     endif

     !bloss:make special statistics flag for radiation,since it's only updated at icycle==1.
     dostatisrad = .false.
     if(mod(nstep,nstatis).eq.0.and.icycle.eq.1) dostatisrad = .true.
#ifdef UWM_STATS

     if (dostatis) then
        do k=1,nzm
           do j=1,ny
              do i=1,nx
                 t1(i,j,k) = t(i,j,k)
                 q1(i,j,k) = micro_field(i,j,k,index_water_vapor)
#ifdef PNNL_STATS
                !MWSWong: THEL and QTOG budgets 
                 thel1(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                          qci(i,j,k), qpi(i,j,k), prespot(k) )
                 qtog1(i,j,k) = 0.0
                 do n = 1, nmicro_fields ! prevents accessing unreferenced memory 
                    qtog1(i,j,k) = qtog1(i,j,k) + flag_wmass(n)*micro_field(i,j,k,n)
                 end do

                !initalize
                 thel(i,j,k) = thel1(i,j,k) 
                 qtog(i,j,k) = qtog1(i,j,k) 
#endif /*PNNL_STATS*/
              enddo
           enddo
        enddo
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 w_prev(i,j,k) = w(i,j,k)
              enddo
           enddo
        enddo
     endif

#ifdef PNNL_STATS
     ! MWSWong: 21 August 2014
     ! Calculating thelstor using this time level's prespot (see thel
     ! calculation above)
     !
     ! Note: All other storage term calculations are based on their
     ! mean values,
     ! calculated in diagnose.F90 (e.g., t0, u0, thel0, etc.).
     ! However, thel0 does not use the updated prespot value, computed
     ! by pressz() at the end
     ! of diagnose(). Recalculating thelstor based on the *actual*
     ! time-level n thel here is
     ! therefore needed.
     do k=1,nzm
       thelstor(k)=0.
       do j=1,ny
         do i=1,nx
           thelstor(k)=thelstor(k)+thel(i,j,k)
         end do
       end do
       thelstor(k)=thelstor(k)/float(nx*ny)
     end do
        
     if (dompi) then
       call averageXY_MPI(thel,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,thelstor)
     end if
#endif /*PNNL_STATS*/

#endif /*UWM_STATS*/

!---------------------------------------------
!  	the Adams-Bashforth scheme in time

     call abcoefs()
 
!---------------------------------------------
!  	initialize stuff: 
	
     call zero()

!-----------------------------------------------------------

     total_water_before = total_water_before + total_water()

!-----------------------------------------------------------
!       Buoyancy term:

	     
     call buoyancy()

!------------------------------------------------------------

     total_water_ls =  total_water_ls - total_water()

!------------------------------------------------------------
!       Large-scale and surface forcing:

     call forcing()

!----------------------------------------------------------
!       Nadging to sounding:

     call nudging()

!----------------------------------------------------------
!   	spange-layer damping near the upper boundary:

     if(dodamping) call damping()

!----------------------------------------------------------

     total_water_ls =  total_water_ls + total_water()

!---------------------------------------------------------
!   Ice fall-out
   
      if(docloud) then
          call ice_fall()
      end if

!----------------------------------------------------------
!     Update scalar boundaries after large-scale processes:

     call boundaries(3)

!---------------------------------------------------------
!     Update boundaries for velocities:

      call boundaries(0)

!-----------------------------------------------
!     surface fluxes:

     if(dosurface) call surface()

!-----------------------------------------------------------
!  SGS physics:

     if (dosgs) call sgs_proc()

!----------------------------------------------------------
!     Fill boundaries for SGS diagnostic fields:

     call boundaries(4)
!-----------------------------------------------
!       advection of momentum:

     call advect_mom()

!----------------------------------------------------------
!	SGS effects on momentum:

     if(dosgs) call sgs_mom()

!-----------------------------------------------------------
!       Coriolis force:
	     
     if(docoriolis) call coriolis()
	 
!---------------------------------------------------------
!       compute rhs of the Poisson equation and solve it for pressure. 

     call pressure()

!---------------------------------------------------------
!       find velocity field at n+1/2 timestep needed for advection of scalars:
!  Note that at the end of the call, the velocities are in nondimensional form.
	 
     call adams()

!----------------------------------------------------------
!     Update boundaries for all prognostic scalar fields for advection:

     call boundaries(2)

!---------------------------------------------------------
!      advection of scalars :

     call advect_all_scalars()

!-----------------------------------------------------------
!    Convert velocity back from nondimensional form:

      call uvw()

!----------------------------------------------------------
!     Update boundaries for scalars to prepare for SGS effects:

     call boundaries(3)
   
!---------------------------------------------------------
!      SGS effects on scalars :

     if (dosgs) call sgs_scalars()

!-----------------------------------------------------------
!       Handle upper boundary for scalars

     if(doupperbound) call upperbound()

!-----------------------------------------------------------
!       Cloud condensation/evaporation and precipitation processes:

#if defined( CLUBB ) && !defined( SILHS )
      if(docloud.or.dosmoke.or.doclubb) call micro_proc()
#elif SILHS
      if ( doclubb .and. ( nstep == 1 .or. mod( nstep, nclubb ) == 0 ) ) then

        select case ( LH_microphys_type )

        case ( LH_microphys_interactive, LH_microphys_non_interactive )

          ! Backup micro_field
          call micro_save( t_prior, w_prior, qn_prior, micro_field_prior ) ! Out

          ! Zero out the tendency arrays 
          LH_t_sum_tndcy = 0._core_rknd
          LH_qn_sum_tndcy = 0._core_rknd
          LH_micro_field_sum_tndcy = 0._core_rknd

          ! Save the values of dostatis and doclubb, then overwrite them with
          ! .false.  We do this in order to avoid calling code that we don't
          ! want enabled when we're calling the microphysics with a sample point
          ! (e.g. micro_diagnose, the statistics, and the divide by cloud
          ! fraction code). We also subcycle the micro_proc call and so must
          ! save the value of dtn and reset at the end of this case statement.
          dostatis_save = dostatis
          doclubb_save = doclubb
          dtn_save = dtn ! Save the SAM timestep
          dostatis = .false.
          doclubb = .false.
          ! The new micro_proc timestep is the SAM timestep (dtn) times nclubb
          ! since we only call micro_proc when clubb and silhs are called.
          dtn = real( nclubb ) * dtn_save 

          ! Iterate over the total number of subcolumns and call the
          ! microphysics code on each
          do LH_iter = 1, LH_microphys_calls

            ! Copy the sample point into micro_field
            call micro_copy_in( LH_iter, LH_microphys_calls, LH_rt, LH_t, micro_field_prior, & ! In
                                X_nl_all_levs ) ! In


            ! Call the microphysics using the sample point
            call micro_proc()

            ! Sum the tendencies
            call accumulate_LH_tendency( dtn, t, qn, micro_field, & ! In
                                         LH_sample_point_weights(:,:,LH_iter), & ! In
                                         LH_t_sum_tndcy, LH_qn_sum_tndcy, & ! In/Out
                                         LH_micro_field_sum_tndcy ) ! In/Out

          end do ! 1..LH_microphys_calls

          ! Average tendencies
          LH_t_avg_tndcy = LH_t_sum_tndcy / real( LH_microphys_calls )
          LH_qn_avg_tndcy = LH_qn_sum_tndcy / real( LH_microphys_calls )
          LH_micro_field_avg_tndcy = LH_micro_field_sum_tndcy / real( LH_microphys_calls )

          ! Copy back the original value of t, w, and the micro fields
          call micro_restore( t_prior, w_prior, qn_prior, micro_field_prior ) ! In

          ! Revert doclubb and dostatis to their prior values
          dostatis = dostatis_save
          doclubb  = doclubb_save

          if ( LH_microphys_type == LH_microphys_interactive ) then

            ! Add tendencies to micro_field and moist static energy (note that we
            ! omit the ghost points on the horiztonal boundaries until later)

            micro_field(1:nx,1:ny,1:nzm,:) = micro_field(1:nx,1:ny,1:nzm,:) &
                                           + dtn * LH_micro_field_avg_tndcy(1:nx,1:ny,1:nzm,:)

            qn = qn + dtn * LH_qn_avg_tndcy
            ! Because of hole filling t may not be prefectly conserved with this
            ! code
            t(1:nx,1:ny,1:nzm) = t(1:nx,1:ny,1:nzm) &
                               + dtn * LH_t_avg_tndcy

            ! Clip all quantites within micro_field using the hole filling
            ! algorithm for mass variables and hard clipping for number
            ! concentrations.
            call clip_micro_field( )

            ! Update qcl, qv and other dynamical core variables
            call micro_update()

          else ! Non-iteractive
            ! Call the microphysics using the mean fields
            call micro_proc()
          end if

          ! Copy back the current value for dtn
          dtn = dtn_save

        case ( LH_microphys_disabled )

          ! Call the microphysics using the mean fields
          call micro_proc()

        case default

          if ( masterproc ) write(0,*) "Fatal error determining LH_microphys_type"
          call task_abort()

        end select ! LH_microphys_type

      else if ( docloud.or.dosmoke ) then
        call micro_proc()

      end if ! doclubb

#else
      if(docloud.or.dosmoke) call micro_proc()
#endif

!----------------------------------------------------------
!  Tracers' physics:

      call tracers_physics()

!-----------------------------------------------------------
!	Radiation

      if(dolongwave.or.doshortwave) then 
	call radiation()     
      end if

!-----------------------------------------------------------
!    Compute diagnostic fields:

      call diagnose()

!----------------------------------------------------------

      total_water_after = total_water_after + total_water()

#ifdef UWM_STATS
      if (dostatis) then
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  w_curr(i,j,k) = w(i,j,k)
               enddo
            enddo
         enddo
#ifdef PNNL_STATS
         do k=1,nzm
            do j=1,ny
              do i=1,nx
                 thel(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                         qci(i,j,k), qpi(i,j,k), prespot(k) )
                 qtog(i,j,k) = 0.0
                 do n = 1, nmicro_fields ! prevents accessing unreferenced memory 
                    qtog(i,j,k) = qtog(i,j,k) + flag_wmass(n)*micro_field(i,j,k,n)
                 end do
              end do
            end do
         end do
#endif /*PNNL_STATS*/
         call stat_varscalar( t, t1, t_avg, t1_avg, t2lebt )
         call setvalue( twlebt, nzm, 0. )
         call stat_covar_sw( t, t1, w_curr, w_prev, twlebt )
         call stat_varscalar( micro_field(:,:,:,index_water_vapor), q1, &
                              q_avg, q1_avg, q2lebt )
         call setvalue( qwlebt, nzm, 0. )
         call stat_covar_sw( micro_field(:,:,:,index_water_vapor), q1, &
                             w_curr, w_prev, qwlebt )
#ifdef PNNL_STATS
        !MWSWong:Convert HL to THEL 
         call stat_varscalar( thel, thel1, thel_avg, thel1_avg, thel2lebt )
         call setvalue( thelwlebt, nzm, 0. )
         call stat_covar_sw( thel, thel1, w_curr, w_prev, thelwlebt )
        !MWSWong:"grand total water variance budget" 
         call stat_varscalar( qtog, qtog1, qtog_avg, qtog1_avg, qtog2lebt )
         call setvalue(qtogwlebt, nzm, 0. )
         call stat_covar_sw( qtog, qtog1, w_curr, w_prev, qtogwlebt) 


        !local storage of the covariance r'thl'
        do k=1,nzm
           qthellebt(k)=0.
           do j=1,ny
              do i=1,nx
                 qthellebt(k) &
                 = qthellebt(k) &
                   + (thel(i,j,k)-thel_avg(k)) &
                     * (micro_field(i,j,k,index_water_vapor)-q_avg(k)) &
                   - (thel1(i,j,k)-thel1_avg(k)) * (q1(i,j,k)-q1_avg(k))
              end do
           end do
           qthellebt(k)=qthellebt(k)*(1./(dtn*nx*ny))
        end do 
        !Heng Xiao w'w'w' budget
        call setvalue( w3lebt, nzm, 0. )
        call stat_w3_bt(w_curr, w_prev, w3lebt)
#endif /*PNNL_STATS*/

      endif
#endif /*UWM_STATS*/
!----------------------------------------------------------
! Rotate the dynamic tendency arrays for Adams-bashforth scheme:

      nn=na
      na=nc
      nc=nb
      nb=nn

#ifdef PNNL_STATS
!----------------------------------------------------------
! MWSWong: set mstor to be the microphysical mixing ratios updated from the last icycle 
if (icyc.ne.ncycle) then
  qstor = q0
  tstor = t0
  ustor = u0
  vstor = v0
! ---MWSWong: adapted from pnnl_branch_r1061
  do k=1, nzm
    qi0(k) = sum(dble(qci(1:nx,1:ny,k)))/float(nx*ny)
    qrain0(k)=sum(qpl(1:nx,1:ny,k))/float(nx*ny)
    qice0(k)=sum(qci(1:nx,1:ny,k))/float(nx*ny) + sum(qpi(1:nx,1:ny,k))/float(nx*ny)
  end do
  qtostor = q0 - qi0
  thelstor = thel0 
!----

  do n=1, nmicro_fields
    do k=1, nzm
       mstor(k, n) = SUM(dble(micro_field(1:nx,1:ny,k,n)))  ! MWSWong: cast to double
    end do
  end do
end if ! icyc.ne.ncycle

#endif /*PNNL_STATS*/
   end do ! icycle	
          
!----------------------------------------------------------
!  collect statistics, write save-file, etc.

   call stepout(nstatsteps)
  
!----------------------------------------------------------
!----------------------------------------------------------

end do ! main loop

#ifdef CLUBB
! Deallocate CLUBB variables, etc.
! -UWM
if ( doclubb .or. doclubbnoninter ) call clubb_sgs_cleanup( )
#endif
#ifdef UWM_STATS
deallocate ( corravg_count )
#endif /*UWM_STATS*/
!----------------------------------------------------------
!----------------------------------------------------------

   call t_stopf('total')
   if(masterproc) call t_prf(rank)

#ifdef CLUBB
! Show that we have completed successfully. Used for scripting SAM runs.
! - nielsenb UWM 4/9/2008
print *,'Done!'

#endif
call task_stop()

end program crm
