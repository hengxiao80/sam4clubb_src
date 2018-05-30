     
subroutine write_fields3D_micro
	
use vars
use microphysics
 
implicit none
character *80 filename,long_name
character *8 name
character *10 timechar
character *4 rankchar
character *6 filetype
character *10 units
character *10 c_z(nzm),c_p(nzm),c_rho(nzm),c_dx, c_dy, c_time
integer i,j,k,nfields,nfields1
real tmp(nx,ny,nzm)

character *3 binumber
integer m


nfields= 13         ! number of 3D fields to save

! if(.not.docloud) nfields=nfields-1
! if(.not.doprecip) nfields=nfields-1
nfields1=0

if(masterproc) then

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
    filename='./DATA3D/'//trim(case)//'_'//trim(caseid)//'_micro_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype
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
      filename='./DATA3D/'//trim(case)//'_'//trim(caseid)//'_micro_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype
      open(46,file=filename,status='unknown',form='unformatted')	
    else
      filename='./DATA3D/'//trim(case)//'_'//trim(caseid)//'_micro_'// &
        rankchar(5-lenstr(rankchar):4)//filetype
      if(nrestart.eq.0.and.notopened3D) then
         open(46,file=filename,status='unknown',form='unformatted')	
      else
         open(46,file=filename,status='unknown', &
                              form='unformatted', position='append')
      end if
      notopened3D=.false.
    end if  

  end if

  if(save3Dbin) then

    write(46) nx,ny,nzm,   &
            nsubdomains,nsubdomains_x,nsubdomains_y,nfields
    do k=1,nzm
      write(46) z(k) 
    end do
    do k=1,nzm
      write(46) pres(k)
    end do
    do k=1,nzm
      write(46) rho(k)
    end do
    write(46) dx
    write(46) dy
    write(46) nstep*dt/(3600.*24.)+day0

  else

    write(long_name,'(10i4)') nx,ny,nzm, &
             nsubdomains,nsubdomains_x,nsubdomains_y,nfields
    do k=1,nzm
       write(c_z(k),'(f10.3)') z(k)
    end do
    do k=1,nzm
       write(c_p(k),'(f10.3)') pres(k)
    end do
    do k=1,nzm
       write(c_rho(k),'(f10.3)') rho(k)
    end do
    write(c_dx,'(f10.0)') dx
    write(c_dy,'(f10.0)') dy
    write(c_time,'(f10.5)') nstep*dt/(3600.*24.)+day0
	
    write(46) long_name(1:40)
    write(46) c_time,c_dx,c_dy,                                    &
         (c_z(k),k=1,nzm),(c_p(k),k=1,nzm),(c_rho(k),k=1,nzm)
  end if ! save3Dbin

end if ! masterproc

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qq(i,j,k)*1000.  ! kg/kg to g/kg
    end do
   end do
  end do
  name='QV'
  long_name='Water Vapor Mixing Ratio'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write qv:',tmp,1,nx,1,ny,nzm)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=super_WRF(i,j,k)
    end do
   end do
  end do
  name='SSW'
  long_name='Saturation ratio (water)'
  units=' '
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write super_WRF:',tmp,1,nx,1,ny,nzm)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
!      tmp(i,j,k)=qc(i,j,k)*1000.*rho(k)
      tmp(i,j,k)=qc(i,j,k)*1000.  ! kg/kg to g/kg
    end do
   end do
  end do
  name='QL'
  long_name='Cloud Water Mixing Ratio'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write qc:',tmp,1,nx,1,ny,nzm)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qr(i,j,k)*1000.  ! kg/kg to g/kg
    end do
   end do
  end do
  name='QR'
  long_name='Rain Water Mixing Ratio'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write qr:',tmp,1,nx,1,ny,nzm)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qi(i,j,k)*1000.  ! kg/kg to g/kg
    end do
   end do
  end do
  name='QI'
  long_name='Cloud Ice Mixing Ratio'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write qi:',tmp,1,nx,1,ny,nzm)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qs(i,j,k)*1000.  ! kg/kg to g/kg
    end do
   end do
  end do
  name='QS'
  long_name='Snow Mixing Ratio'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write qs:',tmp,1,nx,1,ny,nzm)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qg(i,j,k)*1000.  ! kg/kg to g/kg
    end do
   end do
  end do
  name='QG'
  long_name='Graupel Mixing Ratio'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write qg:',tmp,1,nx,1,ny,nzm)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=ncn(i,j,k)*rho(k)*1.e-6  ! #/kg to #/cm3
    end do
   end do
  end do
  name='NCN'
  long_name='CCN Number Concentration'
  units='#/cm3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write ncn:',tmp,1,nx,1,ny,nzm)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=nin(i,j,k)*rho(k)*1.e-6  ! #/kg to #/cm3
    end do
   end do
  end do
  name='NIN'
  long_name='IN Number Concentration'
  units='#/cm3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write nin:',tmp,1,nx,1,ny,nzm)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=nw(i,j,k)*rho(k)*1.e-6  ! #/kg to #/cm3
    end do
   end do
  end do
  name='NCD'
  long_name='Droplet Number Concentration'
  units='#/cm3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write ncd:',tmp,1,nx,1,ny,nzm)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=ni(i,j,k)*rho(k)*1.e-6  ! #/kg to #/cm3
    end do
   end do
  end do
  name='NIC'
  long_name='Ice Crystal Number Concentration'
  units='#/cm3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write nic:',tmp,1,nx,1,ny,nzm)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=ncn_a(i,j,k)*rho(k)*1.e-6  ! #/kg to #/cm3
    end do
   end do
  end do
  name='NCN_a'
  long_name='Activated CCN Number Concentration'
  units='#/cm3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write ncn_a:',tmp,1,nx,1,ny,nzm)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=nin_a(i,j,k)*rho(k)*1.e-6  ! #/kg to #/cm3
    end do
   end do
  end do
  name='NIN_a'
  long_name='Activated IN Number Concentration'
  units='#/cm3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  call fminmax_print(' write nin_a:',tmp,1,nx,1,ny,nzm)

  call task_barrier()
  print*, nfields, nfields1

  if(nfields.ne.nfields1) then
    if(masterproc) print*,'write_fields3D_micro error: nfields'
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

 
end subroutine write_fields3D_micro
