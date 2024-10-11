f( doclubb ) then ! dschanen UWM 28 May 2008

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