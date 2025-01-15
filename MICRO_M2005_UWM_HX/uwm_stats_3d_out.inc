!--------------------------------------------------------------------------------------------------
subroutine write_3d_micro_fields()

! Description:
!   This subroutine was created in response to sam_clubb:ticket:82.
!   We would like to look at 3D fields of various microphysical budgets, terms,
!   etc. However, the variables are local to Morrison. Therefore, we write a
!   seperate 3D field for the micro variables.
!
!   This is largely 'write_3D_fields.F90', but reused here. 
!--------------------------------------------------------------------------------------------------

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
character *12 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
integer i,j,k,n,nfields,nfields1
real tmp(nx,ny,nzm)

real, dimension(nx,ny,nzm) :: &
  thl, & ! Liquid water potential temperature                 [K]
  rt,  & ! Total water mixing ratio                           [kg/kg]
  chi, & ! Extended liquid water mixing ratio                 [kg/kg]
  eta    ! Coordinate orthogonal to chi in PDF transformation [kg/kg]

call t_startf('3D_microrain')

#ifdef LASSO_ENA_3D_MICRO_OUTPUT
nfields = 37-8
#else
nfields=37 ! number of 3D fields to save
#endif
nfields1=0 ! assertion check

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
      filetype = '_micro.bin3D'
    else
      filetype = '_micro.com3D'
    end if
    
    filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
    
    open(46,file=filename,status='unknown',form='unformatted')

  else
    
    if(save3Dbin) then
      if(save3Dsep) then
        filetype = '_micro.bin3D'
      else
        filetype = '_micro.bin2D'
      end if
    else
      if(save3Dsep) then
        filetype = '_micro.com3D'
      else
        filetype = '_micro.com2D'
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
    end if!save3Dsep
  end if!RUN3D

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


!--------------------------------------
! Micro fields
!--------------------------------------
#ifndef LASSO_ENA_3D_MICRO_OUTPUT
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
         tmp(i,j,k) = qcl(i,j,k)
      enddo
   enddo
enddo
name = 'RC'
long_name = 'Cloud water mixing ratio'
units = 'kg/kg'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )



if (doprecip) then

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
           tmp(i,j,k) = micro_field(i,j,k,inr)
        enddo
     enddo
  enddo
  name = 'NR'
  long_name = 'Rain drop concentration'
  units = 'num/kg'
  call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )
endif ! doprecip
# endif
  !--------------------------------------
  ! RR Budget Terms
  !--------------------------------------
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRE_3D(i,j,k)
    end do
   end do
  end do
  name='PRE'
  long_name='Evaporation of rain'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRA_3D(i,j,k)
    end do
   end do
  end do
  name='PRA'
  long_name='Accretion of cloud droplets by rain'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRC_3D(i,j,k)
    end do
   end do
  end do
  name='PRC'
  long_name='Autoconversion of cloud droplets to rain'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRACS_3D(i,j,k)
    end do
   end do
  end do
  name='PRACS'
  long_name='Collection of rain by snow to form snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mMNUCCR_3D(i,j,k)
    end do
   end do
  end do
  name='MNUCCR'
  long_name='Contact freezing of rain drops'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mQMULTR_3D(i,j,k)
    end do
   end do
  end do
  name='QMULTR'
  long_name='Splintering  from rain drops accreted onto snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mQMULTRG_3D(i,j,k)
    end do
   end do
  end do
  name='QMULTRG'
  long_name='Splintering from rain drops accreted onto graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPIACR_3D(i,j,k)
    end do
   end do
  end do
  name='PIACR'
  long_name='Collection of cloud ice by rain to form graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPIACRS_3D(i,j,k)
    end do
   end do
  end do
  name='PIACRS'
  long_name='Collection of cloud ice by rain to form snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRACG_3D(i,j,k)
    end do
   end do
  end do
  name='PRACG'
  long_name='Collection of rain by graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPGRACS_3D(i,j,k)
    end do
   end do
  end do
  name='PGRACS'
  long_name='Collection of rain by snow to form snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPSMLT_3D(i,j,k)
    end do
   end do
  end do
  name='PSMLT'
  long_name='Freezing of rain to form snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPGMLT_3D(i,j,k)
    end do
   end do
  end do
  name='PGMLT'
  long_name='Freezing of rain to form graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)


!--------------------------------------
! NR Budget Terms
!--------------------------------------
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSIZEFIX_NR_3D(i,j,k)
    end do
   end do
  end do
  name='SIZEFIX_NR'
  long_name='Adjust rain number when size is too large/small'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNEGFIX_NR_3D(i,j,k)
    end do
   end do
  end do
  name='NEGFIX_NR'
  long_name='Removal of negative rain drop number concentration'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNSUBR_3D(i,j,k)
    end do
   end do
  end do
  name='NSUBR'
  long_name='Evaporation of rain'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNSMLTR_3D(i,j,k)
    end do
   end do
  end do
  name='NSMLTR'
  long_name='Melting of snow to form rain'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNGMLTR_3D(i,j,k)
    end do
   end do
  end do
  name='NGMLTR'
  long_name='Melting of graupel to form rain'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNPRC1_3D(i,j,k)
    end do
   end do
  end do
  name='NPRC1'
  long_name='Change in rain due to autoconversion'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNPRACS_3D(i,j,k)
    end do
   end do
  end do
  name='NPRACS'
  long_name='Collection of rain by snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNNUCCR_3D(i,j,k)
    end do
   end do
  end do
  name='NNUCCR'
  long_name='Contact freezing of rain'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNRAGG_3D(i,j,k)
    end do
   end do
  end do
  name='NRAGG'
  long_name='Self collection of rain drops'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNIACR_3D(i,j,k)
    end do
   end do
  end do
  name='NIACR'
  long_name='Collection of ice by rain to form graupel'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNIACRS_3D(i,j,k)
    end do
   end do
  end do
  name='NIACRS'
  long_name='Collection of ice by rain to form snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNPRACG_3D(i,j,k)
    end do
   end do
  end do
  name='NPRACG'
  long_name='Collection of rain by graupel'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNGRACS_3D(i,j,k)
    end do
   end do
  end do
  name='NGRACS'
  long_name='Collection of rain by snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNR_INST_3D(i,j,k)
    end do
   end do
  end do
  name='NR_INST'
  long_name='Instantaneous processes'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=rain_vel_3D(i,j,k)
    end do
   end do
  end do
  name='rain_vel'
  long_name='Rain drop velocity'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=EFFR_3D(i,j,k)
    end do
   end do
  end do
  name='EFFR'
  long_name='Effective radius of rain'
  units='micron'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

call task_barrier()

if(nfields.ne.nfields1) then
  if( masterproc ) print*,'write_fields3D error: nfields=',nfields,'nfields1=',nfields1
  call task_abort()
endif
if ( masterproc ) then
   close (46)
   if(RUN3D .or. save3Dsep) then
     if(dogzip3D) call systemf('gzip -f '//filename)
     print*, 'Writting 3D data. file:'//filename
   else
     print*, 'Appending 3D data. file:'//filename
  end if
endif


call t_stopf('3D_microrain')

end subroutine write_3d_micro_fields


!--------------------------------------------------------------------------------------------------
subroutine write_3d_micro_fields_frzmr()
!   MWSWONG/Apr2015
!   SIMILAR TO WRITE_3D_MICRO_FIELDS BUT FOR FROZEN MIXING RATIO BUDGETS.
!
! Description:
!   This subroutine was created in response to sam_clubb:ticket:82.
!   We would like to look at 3D fields of various microphysical budgets, terms,
!   etc. However, the variables are local to Morrison. Therefore, we write a
!   seperate 3D field for the micro variables.
!
!   This is largely 'write_3D_fields.F90', but reused here. 
!--------------------------------------------------------------------------------------------------

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
character *12 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
integer i,j,k,n,nfields,nfields1
real tmp(nx,ny,nzm)

real, dimension(nx,ny,nzm) :: &
  thl, & ! Liquid water potential temperature                 [K]
  rt,  & ! Total water mixing ratio                           [kg/kg]
  chi, & ! Extended liquid water mixing ratio                 [kg/kg]
  eta    ! Coordinate orthogonal to chi in PDF transformation [kg/kg]

call t_startf('3D_frzmr')

nfields=24 ! number of 3D fields to save (frozen condensates mixing ratio budgets, MWSWONG)
nfields1=0 ! assertion check

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
      filetype = '_frzmr.bin3D'
    else
      filetype = '_frzmr.com3D'
    end if
    
    filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
    
    open(46,file=filename,status='unknown',form='unformatted')

  else
    
    if(save3Dbin) then
      if(save3Dsep) then
        filetype = '_frzmr.bin3D'
      else
        filetype = '_frzmr.bin2D'
      end if
    else
      if(save3Dsep) then
        filetype = '_frzmr.com3D'
      else
        filetype = '_frzmr.com2D'
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
    end if!save3Dsep
  end if!RUN3D

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


!--------------------------------------
! QG Budget Terms
!--------------------------------------
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPSACWG_3D(i,j,k)
    end do
   end do
  end do
  name='PSACWG'
  long_name='Collection of cloud water by graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPGSACW_3D(i,j,k)
    end do
   end do
  end do
  name='PGSACWG'
  long_name='Riming of cloud water by snow to form graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRDG_3D(i,j,k)
    end do
   end do
  end do
  name='PRDG'
  long_name='Depositional growth of graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mEPRDG_3D(i,j,k)
    end do
   end do
  end do
  name='EPRDG'
  long_name='Sublimation of graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRACI_3D(i,j,k)
    end do
   end do
  end do
  name='PRACI'
  long_name='Collection of cloud ice by rain to form graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPSACR_3D(i,j,k)
    end do
   end do
  end do
  name='PSACR'
  long_name='Collection of snow by rain to form graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mEVPMG_3D(i,j,k)
    end do
   end do
  end do
  name='EVPMG'
  long_name='Evaporation of metled graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mQG_INST_3D(i,j,k)
    end do
   end do
  end do
  name='QG_INST'
  long_name='Changes due to instantaneous processes'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)


!--------------------------------------
! QS Budget Terms
!--------------------------------------
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRAI_3D(i,j,k)
    end do
   end do
  end do
  name='PRAI'
  long_name='Collection of cloud ice by snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPSACWS_3D(i,j,k)
    end do
   end do
  end do
  name='PSACWS'
  long_name='Collection of cloud water by snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRDS_3D(i,j,k)
    end do
   end do
  end do
  name='PRDS'
  long_name='Depositional growth of snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRCI_3D(i,j,k)
    end do
   end do
  end do
  name='PRCI'
  long_name='Autoconversion of cloud ice to snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mEPRDS_3D(i,j,k)
    end do
   end do
  end do
  name='EPRDS'
  long_name='Sublimation of snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRACIS_3D(i,j,k)
    end do
   end do
  end do
  name='PRACIS'
  long_name='Collection of cloud ice by rain to form snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mEVPMS_3D(i,j,k)
    end do
   end do
  end do
  name='EVPMS'
  long_name='Evaporation of melted snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mQS_INST_3D(i,j,k)
    end do
   end do
  end do
  name='QS_INST'
  long_name='Changes due to instantaneous processes'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

!--------------------------------------
! QI Budget Terms
!--------------------------------------
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRD_3D(i,j,k)
    end do
   end do
  end do
  name='PRD'
  long_name='Depositional growth of cloud ice'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mEPRD_3D(i,j,k)
    end do
   end do
  end do
  name='EPRD'
  long_name='Sublimation of cloud ice'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPSACWI_3D(i,j,k)
    end do
   end do
  end do
  name='PSACWI'
  long_name='Collection of cloud water by cloud ice'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mMNUCCC_3D(i,j,k)
    end do
   end do
  end do
  name='MNUCCC'
  long_name='Contact freezing of cloud droplets'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mQMULTS_3D(i,j,k)
    end do
   end do
  end do
  name='QMULTS'
  long_name='Splintering from cloud droplets accreted onto snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mQMULTG_3D(i,j,k)
    end do
   end do
  end do
  name='QMULTG'
  long_name='Splintering from droplets accreted onto graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mMNUCCD_3D(i,j,k)
    end do
   end do
  end do
  name='MNUCCD'
  long_name='Freezing of aerosol'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mQI_INST_3D(i,j,k)
    end do
   end do
  end do
  name='QI_INST'
  long_name='Changes due to instantaneous processes'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)


call task_barrier()

if(nfields.ne.nfields1) then
  if( masterproc ) print*,'write_fields3D error: nfields=',nfields,'nfields1=',nfields1
  call task_abort()
endif
if ( masterproc ) then
   close (46)
   if(RUN3D .or. save3Dsep) then
     if(dogzip3D) call systemf('gzip -f '//filename)
     print*, 'Writting 3D data. file:'//filename
   else
     print*, 'Appending 3D data. file:'//filename
  end if
endif


call t_stopf('3D_frzmr')

end subroutine write_3d_micro_fields_frzmr


!--------------------------------------------------------------------------------------------------
subroutine write_3d_micro_fields_frznc()
!   MWSWONG/Apr2015
!   SIMILAR TO WRITE_3D_MICRO_FIELDS BUT FOR FROZEN MIXING RATIO BUDGETS.
!
! Description:
!   This subroutine was created in response to sam_clubb:ticket:82.
!   We would like to look at 3D fields of various microphysical budgets, terms,
!   etc. However, the variables are local to Morrison. Therefore, we write a
!   seperate 3D field for the micro variables.
!
!   This is largely 'write_3D_fields.F90', but reused here. 
!--------------------------------------------------------------------------------------------------

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
character *12 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
integer i,j,k,n,nfields,nfields1
real tmp(nx,ny,nzm)

real, dimension(nx,ny,nzm) :: &
  thl, & ! Liquid water potential temperature                 [K]
  rt,  & ! Total water mixing ratio                           [kg/kg]
  chi, & ! Extended liquid water mixing ratio                 [kg/kg]
  eta    ! Coordinate orthogonal to chi in PDF transformation [kg/kg]

call t_startf('3D_frznc')

nfields=25 ! number of 3D fields to save (frozen condensates mixing ratio budgets, MWSWONG)
nfields1=0 ! assertion check

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
      filetype = '_frznc.bin3D'
    else
      filetype = '_frznc.com3D'
    end if
    
    filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
    
    open(46,file=filename,status='unknown',form='unformatted')

  else
    
    if(save3Dbin) then
      if(save3Dsep) then
        filetype = '_frznc.bin3D'
      else
        filetype = '_frznc.bin2D'
      end if
    else
      if(save3Dsep) then
        filetype = '_frznc.com3D'
      else
        filetype = '_frznc.com2D'
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
    end if!save3Dsep
  end if!RUN3D

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

!--------------------------------------
! NG Budget Terms
!--------------------------------------
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSIZEFIX_NG_3D(i,j,k)
    end do
   end do
  end do
  name='SZFIX_NG'
  long_name='Adjust graupel number when size is too large/small'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNEGFIX_NG_3D(i,j,k)
    end do
   end do
  end do
  name='NGFIX_NG'
  long_name='Removal of negative graupel concentration'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNGMLTG_3D(i,j,k)
    end do
   end do
  end do
  name='NGMLTG'
  long_name='Melting of graupel'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNSCNG_3D(i,j,k)
    end do
   end do
  end do
  name='NSCNG'
  long_name='Conversion of snow to graupel'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNSUBG_3D(i,j,k)
    end do
   end do
  end do
  name='NSUBG'
  long_name='Sublimation of graupel'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNG_INST_3D(i,j,k)
    end do
   end do
  end do
  name='NG_INST'
  long_name='Changes due to instantaneous processes'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

!--------------------------------------
! NS Budget Terms
!--------------------------------------
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSIZEFIX_NS_3D(i,j,k)
    end do
   end do
  end do
  name='SZFIX_NS'
  long_name='Adjust snow number when size is too large/small'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNEGFIX_NS_3D(i,j,k)
    end do
   end do
  end do
  name='NGFIX_NS'
  long_name='Removal of negative snow concentration'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNSMLTS_3D(i,j,k)
    end do
   end do
  end do
  name='NSMLTS'
  long_name='Melting of snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNSAGG_3D(i,j,k)
    end do
   end do
  end do
  name='NSAGG'
  long_name='Aggregation of snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNPRCI_3D(i,j,k)
    end do
   end do
  end do
  name='NPRCI'
  long_name='Autoconversion of cloud ice to snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNSUBS_3D(i,j,k)
    end do
   end do
  end do
  name='NSUBS'
  long_name='Sublimination of snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNS_INST_3D(i,j,k)
    end do
   end do
  end do
  name='NS_INST'
  long_name='Changes due to instantaneous processes'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
!--------------------------------------
! NI Budget Terms
!--------------------------------------
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSIZEFIX_NI_3D(i,j,k)
    end do
   end do
  end do
  name='SZFIX_NI'
  long_name='Adjust cloud ice number when size is too large/small'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNEGFIX_NI_3D(i,j,k)
    end do
   end do
  end do
  name='NGFIX_NI'
  long_name='Removal of negative ice concentration'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNNUCCC_3D(i,j,k)
    end do
   end do
  end do
  name='NNUCCC'
  long_name='Contact freezing'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNPRAI_3D(i,j,k)
    end do
   end do
  end do
  name='NPRAI'
  long_name='Accretion of cloud ice by snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNMULTS_3D(i,j,k)
    end do
   end do
  end do
  name='NMULTS'
  long_name='Splintering due to accretion of droplets by snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNMULTG_3D(i,j,k)
    end do
   end do
  end do
  name='NMULTG'
  long_name='Splintering due to accretion of droplets by graupel'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNMULTR_3D(i,j,k)
    end do
   end do
  end do
  name='NMULTR'
  long_name='Splintering due to rain riming snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNMULTRG_3D(i,j,k)
    end do
   end do
  end do
  name='NMULTRG'
  long_name='Splintering due to rain riming graupel'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNNUCCD_3D(i,j,k)
    end do
   end do
  end do
  name='NNUCCD'
  long_name='Ice nucleation by freezing aerosol'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNSUBI_3D(i,j,k)
    end do
   end do
  end do
  name='NSUBI'
  long_name='Sublimation of cloud ice'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNIM_MORR_CL_3D(i,j,k)
    end do
   end do
  end do
  name='NIM_MORR_CL'
  long_name='Clipping of large ice concentration'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNI_INST_3D(i,j,k)
    end do
   end do
  end do
  name='NI_INST'
  long_name='Changes due to instantaneous processes'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

call task_barrier()

if(nfields.ne.nfields1) then
  if( masterproc ) print*,'write_fields3D error: nfields=',nfields,'nfields1=',nfields1
  call task_abort()
endif
if ( masterproc ) then
   close (46)
   if(RUN3D .or. save3Dsep) then
     if(dogzip3D) call systemf('gzip -f '//filename)
     print*, 'Writting 3D data. file:'//filename
   else
     print*, 'Appending 3D data. file:'//filename
  end if
endif


call t_stopf('3D_frznc')

end subroutine write_3d_micro_fields_frznc


subroutine write_3d_micro_fields_sed()
!   MWSWONG/Apr2015
!   SIMILAR TO WRITE_3D_MICRO_FIELDS BUT FOR SEDIMENTATION TERMS. 
!
! Description:
!   This subroutine was created in response to sam_clubb:ticket:82.
!   We would like to look at 3D fields of various microphysical budgets, terms,
!   etc. However, the variables are local to Morrison. Therefore, we write a
!   seperate 3D field for the micro variables.
!
!   This is largely 'write_3D_fields.F90', but reused here. 
!--------------------------------------------------------------------------------------------------

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
character *12 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
integer i,j,k,n,nfields,nfields1
real tmp(nx,ny,nzm)

real, dimension(nx,ny,nzm) :: &
  thl, & ! Liquid water potential temperature                 [K]
  rt,  & ! Total water mixing ratio                           [kg/kg]
  chi, & ! Extended liquid water mixing ratio                 [kg/kg]
  eta    ! Coordinate orthogonal to chi in PDF transformation [kg/kg]

call t_startf('3D_sed')

nfields=8 ! number of 3D fields to save (frozen condensates mixing ratio budgets, MWSWONG)
nfields1=0 ! assertion check

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
      filetype = '_sedim.bin3D'
    else
      filetype = '_sedim.com3D'
    end if
    
    filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
    
    open(46,file=filename,status='unknown',form='unformatted')

  else
    
    if(save3Dbin) then
      if(save3Dsep) then
        filetype = '_sedim.bin3D'
      else
        filetype = '_sedim.bin2D'
      end if
    else
      if(save3Dsep) then
        filetype = '_sedim.com3D'
      else
        filetype = '_sedim.com2D'
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
    end if!save3Dsep
  end if!RUN3D

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


!--------------------------------------
! Sedimentation terms
!--------------------------------------
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSTENDQR_3D(i,j,k)
    end do
   end do
  end do
  name='STENDQR'
  long_name='Tendency of rain due to sedimentation'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSTENDQCI_3D(i,j,k)
    end do
   end do
  end do
  name='STENDQCI'
  long_name='Tendency of cloud ice due to sedimentation'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSTENDQG_3D(i,j,k)
    end do
   end do
  end do
  name='STENDQG'
  long_name='Tendency of graupel due to sedimentation'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSTENDQS_3D(i,j,k)
    end do
   end do
  end do
  name='STENDQS'
  long_name='Tendency of snow due to sedimentation'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSTENDNR_3D(i,j,k)
    end do
   end do
  end do
  name='STENDNR'
  long_name='Tendency of rain number conc due to sedimentation'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSTENDNCI_3D(i,j,k)
    end do
   end do
  end do
  name='STENDNCI'
  long_name='Tendency of cloud ice number conc due to sedimentation'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSTENDNG_3D(i,j,k)
    end do
   end do
  end do
  name='STENDNG'
  long_name='Tendency of graupel number conc due to sedimentation'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSTENDNS_3D(i,j,k)
    end do
   end do
  end do
  name='STENDNS'
  long_name='Tendency of snow number conc due to sedimentation'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

call task_barrier()

if(nfields.ne.nfields1) then
  if( masterproc ) print*,'write_fields3D error: nfields=',nfields,'nfields1=',nfields1
  call task_abort()
endif
if ( masterproc ) then
   close (46)
   if(RUN3D .or. save3Dsep) then
     if(dogzip3D) call systemf('gzip -f '//filename)
     print*, 'Writting 3D data. file:'//filename
   else
     print*, 'Appending 3D data. file:'//filename
  end if
endif


call t_stopf('3D_sed')

end subroutine write_3d_micro_fields_sed
