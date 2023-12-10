     
subroutine write_fields3D

use vars
#ifdef CLUBB
use clubbvars
use sgs_params
#endif

use rad, only: qrad
use sgs, only: tke

#ifdef ATEX
use domain
use tracers, only: tracer, tracername
#elif DYCOMSRF01
use domain
use tracers, only: tracer, tracername
#elif BOMEX
use domain
use tracers, only: tracer, tracername
#elif HISCALE
use domain
use tracers, only: tracer, tracername
#elif LASSO_ENA
use domain
use tracers, only: tracer, tracername
#endif

use params
use microphysics, only: nmicro_fields, micro_field, flag_number, &
     flag_micro3Dout, mkname, mklongname, mkunits, mkoutputscale, &
     index_water_vapor, GET_reffc, Get_reffi
use rad, only: rel_rad, rei_rad
implicit none
character *120 filename
character *80 long_name
character *8 name
character *10 timechar
character *4 rankchar
character *5 sepchar
character *6 filetype
character *10 units
character *12 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
integer i,j,k,n,nfields,nfields1
real(4) tmp(nx,ny,nzm)
#ifdef CLUBB
real, external :: LIN_INT
#endif

nfields= 9 ! number of 3D fields to save
#ifdef CLUBB
if( .not.docloud .and. .not.(doclubb.or.doclubbnoninter ) ) nfields=nfields-1 ! dschanen UWM 19 June
#else
if(.not.docloud) nfields=nfields-1
#endif
if(.not.doprecip) nfields=nfields-1
#ifdef CLUBB
if( doclubb ) nfields=nfields+11 ! dschanen UWM 28 May 2008
#endif
!bloss: add 3D outputs for microphysical fields specified by flag_micro3Dout
!       except for water vapor (already output as a SAM default).
! turning off 3D micro prognostic variable output for HISCALE tracking run - HX
! if(docloud) nfields=nfields+SUM(flag_micro3Dout)-flag_micro3Dout(index_water_vapor)
if((dolongwave.or.doshortwave).and..not.doradhomo) nfields=nfields+1
if(compute_reffc.and.(dolongwave.or.doshortwave).and.rad3Dout) nfields=nfields+1
if(compute_reffi.and.(dolongwave.or.doshortwave).and.rad3Dout) nfields=nfields+1

#ifdef ATEX
nfields = nfields + ntracers
#elif DYCOMSRF01
nfields = nfields + ntracers
#elif BOMEX
nfields = nfields + ntracers
#elif HISCALE
nfields = nfields + ntracers
#elif LASSO_ENA
nfields = nfields + ntracers
#endif

nfields1=0


if(masterproc.or.output_sep) then

  if(output_sep) then
     write(rankchar,'(i4)') rank
     sepchar="_"//rankchar(5-lenstr(rankchar):4)
  else
     sepchar=""
  end if
  write(rankchar,'(i4)') nsubdomains
  write(timechar,'(i10)') nstep
  do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
  end do

  if(RUN3D) then
    if(save3Dbin) then
      filetype = '.bin3D'
    else
      filetype = '.com3D'
    end if
    filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
    open(46,file=filename,status='unknown',form='unformatted')

  else
    if(save3Dbin) then
     if(save3Dsep) then
       filetype = '.bin3D'
     else
       filetype = '.bin2D'
     end if
    else
     if(save3Dsep) then
       filetype = '.com3D'
     else
       filetype = '.com2D'
     end if
    end if
    if(save3Dsep) then
      filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
      open(46,file=filename,status='unknown',form='unformatted')	
    else
      filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//filetype//sepchar
      if(nrestart.eq.0.and.notopened3D) then
         open(46,file=filename,status='unknown',form='unformatted')	
      else
         open(46,file=filename,status='unknown', &
                              form='unformatted', position='append')
      end if
      notopened3D=.false.
    end if  

  end if

  if(masterproc) then

    if(save3Dbin) then

      write(46) nx,ny,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields
      do k=1,nzm
        write(46) z(k) 
      end do
      do k=1,nzm
        write(46) pres(k)
      end do
      write(46) dx
      write(46) dy
      write(46) nstep*dt/(3600.*24.)+day0

    else

      write(long_name,'(8i4)') nx,ny,nzm,nsubdomains, &
                                   nsubdomains_x,nsubdomains_y,nfields
      do k=1,nzm
         write(c_z(k),'(f12.3)') z(k)
      end do
      do k=1,nzm
         write(c_p(k),'(f12.3)') pres(k)
      end do
      write(c_dx,'(f12.0)') dx
      write(c_dy,'(f12.0)') dy
      write(c_time,'(f12.5)') nstep*dt/(3600.*24.)+day0
	
      write(46) long_name(1:32)
      write(46) c_time,c_dx,c_dy, (c_z(k),k=1,nzm),(c_p(k),k=1,nzm)

    end if ! save3Dbin

  end if ! masterproc
 
end if ! masterproc.or.output_sep

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=u(i,j,k) + ug
    end do
   end do
  end do
  name='U'
  long_name='X Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=v(i,j,k) + vg
    end do
   end do
  end do
  name='V'
  long_name='Y Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=w(i,j,k)
    end do
   end do
  end do
  name='W'
  long_name='Z Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=p(i,j,k)
    end do
   end do
  end do
  name='PP'
  long_name='Pressure Perturbation'
  units='Pa'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tke(i,j,k)
    end do
   end do
  end do
  name='TKES'
  long_name='SGS TKE'
  units='m^2/s^2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

if((dolongwave.or.doshortwave).and..not.doradhomo) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qrad(i,j,k)*86400.
    end do
   end do
  end do
  name='QRAD'
  long_name='Radiative heating rate'
  units='K/day'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if

if(compute_reffc.and.(dolongwave.or.doshortwave).and.rad3Dout) then
  nfields1=nfields1+1
  tmp(1:nx,1:ny,1:nzm)=Get_reffc()
  name='REL'
  long_name='Effective Radius for Cloud Liquid Water'
  units='mkm'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if

if(compute_reffi.and.(dolongwave.or.doshortwave).and.rad3Dout) then
  nfields1=nfields1+1
  tmp(1:nx,1:ny,1:nzm)=Get_reffi()
  name='REI'
  long_name='Effective Radius for Cloud Ice'
  units='mkm'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if


  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tabs(i,j,k)
    end do
   end do
  end do
  name='TABS'
  long_name='Absolute Temperature'
  units='K'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qv(i,j,k)*1.e3
    end do
   end do
  end do
  name='QV'
  long_name='Water Vapor'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

#ifdef CLUBB
if(docloud .or.doclubb.or.doclubbnoninter) then ! dschanen UWM 19 June 2007
#else
if(docloud) then
#endif  
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=(qcl(i,j,k)+qci(i,j,k))*1.e3
    end do
   end do
  end do
  name='QN'
  long_name='Non-precipitating Condensate (Water+Ice)'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if

#ifdef ATEX
  do n = 1, ntracers
    nfields1=nfields1+1
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,k)=tracer(i,j,k,n)
        end do
      end do
    end do
    name=trim(tracername(n))
    long_name=trim(tracername(n))
    units='kg/kg'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)
  end do
#elif DYCOMSRF01
  do n = 1, ntracers
    nfields1=nfields1+1
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,k)=tracer(i,j,k,n)
        end do
      end do
    end do
    name=trim(tracername(n))
    long_name=trim(tracername(n))
    units='kg/kg'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)
  end do
#elif BOMEX
  do n = 1, ntracers
    nfields1=nfields1+1
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,k)=tracer(i,j,k,n)
        end do
      end do
    end do
    name=trim(tracername(n))
    long_name=trim(tracername(n))
    units='kg/kg'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)
  end do
#elif HISCALE
  do n = 1, ntracers
    nfields1=nfields1+1
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,k)=tracer(i,j,k,n)
        end do
      end do
    end do
    name=trim(tracername(n))
    long_name=trim(tracername(n))
    units='kg/kg'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)
  end do
#elif LASSO_ENA
>>>>>>> 29eec2b (LASSO_ENA testrun setup on cirrus.)
  do n = 1, ntracers
    nfields1=nfields1+1
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,k)=tracer(i,j,k,n)
        end do
      end do
    end do
    name=trim(tracername(n))
    long_name=trim(tracername(n))
    units='kg/kg'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                  save3Dbin,dompi,rank,nsubdomains)
  end do
#endif

if(doprecip) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=(qpl(i,j,k)+qpi(i,j,k))*1.e3
    end do
   end do
  end do
  name='QP'
  long_name='Precipitating Water (Rain+Snow)'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if

! turning off 3D micro prognostic variable output for HISCALE tracking run - HX
! do n = 1,nmicro_fields
!    if(docloud.AND.flag_micro3Dout(n).gt.0.AND.n.ne.index_water_vapor) then
!       nfields1=nfields1+1
!       do k=1,nzm
!          do j=1,ny
!             do i=1,nx
!                tmp(i,j,k)=micro_field(i,j,k,n)*mkoutputscale(n)
!             end do
!          end do
!          ! remove factor of rho from number, if this field is a number concentration
!          if(flag_number(n).gt.0) tmp(:,:,k) = tmp(:,:,k)*rho(k)
!       end do
!       name=TRIM(mkname(n))
!       long_name=TRIM(mklongname(n))
!       units=TRIM(mkunits(n))
!       call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
!            save3Dbin,dompi,rank,nsubdomains)
!    end if
! end do
#ifdef CLUBB
if( doclubb ) then ! dschanen UWM 28 May 2008

  nfields1=nfields1+1
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=lin_int( real( upwp(i,j,k+1) ), real( upwp(i,j,k) ), zi(k+1), zi(k), z(k) )
     end do
    end do
   end do
   name='UPWP'
   long_name="U'W'"
   units='m^2/s^2'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                   save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=lin_int( real( vpwp(i,j,k+1) ), real( vpwp(i,j,k) ), zi(k+1), zi(k), z(k) )

     end do
    end do
   end do
   name='VPWP'
   long_name="U'W'"
   units='m^2/s^2'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                   save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=lin_int( real( up2(i,j,k+1) ), real( up2(i,j,k) ), zi(k+1), zi(k), z(k) )
     end do
    end do
   end do
   name='UP2'
   long_name="U'^2'"
   units='m^2/s^2'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                   save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=lin_int( real( vp2(i,j,k+1) ), real( vp2(i,j,k) ), zi(k+1), zi(k), z(k) )
     end do
    end do
   end do
   name='VP2'
   long_name="V'^2'"
   units='m^2/s^2'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                   save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=lin_int( real( wp2(i,j,k+1) ), real( wp2(i,j,k) ), zi(k+1), zi(k), z(k) )
     end do
    end do
   end do
   name='WP2'
   long_name="W'^2'"
   units='m^2/s^2'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                   save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=lin_int( real( wpthlp(i,j,k+1) ), real( wpthlp(i,j,k) ), zi(k+1), zi(k), z(k) )
     end do
    end do
   end do
   name='WPTHLP'
   long_name="w'theta_l'"
   units='m K/s'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                   save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=lin_int( real( wprtp(i,j,k+1) ), real( wprtp(i,j,k) ), zi(k+1), zi(k), z(k) )
     end do
    end do
   end do
   name='WPRTP'
   long_name="w'r_t'"
   units='kg m/kg s'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                   save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=real( cloud_frac(i,j,k+1) )
     end do
    end do
   end do
   name='CLD_FRAC'
   long_name="Cloud Fraction"
   units='-'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                   save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=real( wp3(i,j,k+1) )
     end do
    end do
   end do
   name='WP3'
   long_name="w'^3"
   units='m^3/s^3'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                   save3Dbin,dompi,rank,nsubdomains)


  nfields1=nfields1+1
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=real( um(i,j,k+1) )
     end do
    end do
   end do
   name='UM'
   long_name='CLUBB u wind'
   units='m/s'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                   save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       tmp(i,j,k)=real( vm(i,j,k+1) )
     end do
    end do
   end do
   name='VM'
   long_name='CLUBB v wind'
   units='m/s'
   call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                   save3Dbin,dompi,rank,nsubdomains)
end if 
#endif

  call task_barrier()

  if(nfields.ne.nfields1) then
    if(masterproc) print*,'write_fields3D error: nfields=',nfields,'  nfields1=',nfields1
    call task_abort()
  end if
  if(masterproc) then
    close (46)
    if(RUN3D.or.save3Dsep) then
       if(dogzip3D) call systemf('gzip -f '//filename)
       print*, 'Writting 3D data. file:'//filename
    else
       print*, 'Appending 3D data. file:'//filename
    end if
  endif

 
end
#ifdef CLUBB
!-------------------------------------------------------------------------------
FUNCTION LIN_INT( var_high, var_low, height_high, height_low, height_int )

! This function computes a linear interpolation of the value of variable.
! Given two known values of a variable at two height values, the value
! of that variable at a height between those two height levels (rather 
! than a height outside of those two height levels) is computed.
!
! Here is a diagram:
!
!  ################################ Height high, know variable value
!
!
!
!  -------------------------------- Height to be interpolated to; linear interpolation
!
!
!
!
!
!  ################################ Height low, know variable value
!
!
! FORMULA:
!
! variable(@ Height interpolation) =
!
! [ (variable(@ Height high) - variable(@ Height low)) / (Height high - Height low) ]
! * (Height interpolation - Height low)  +  variable(@ Height low)

! Author: Brian Griffin, UW-Milwaukee

! References: None

IMPLICIT NONE

! Input Variables
REAL, INTENT(IN):: var_high
REAL, INTENT(IN):: var_low
REAL, INTENT(IN):: height_high
REAL, INTENT(IN):: height_low
REAL, INTENT(IN):: height_int

! Output Variable
REAL:: lin_int

lin_int = ( var_high - var_low ) / ( height_high - height_low ) &
         * ( height_int - height_low ) + var_low


END FUNCTION LIN_INT
#endif /*CLUBB*/
