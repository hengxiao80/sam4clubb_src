subroutine homogenize_cld_env

  ! This subroutine homogenizes the environment in the cloud layer
  ! and collects statistics for the environment and cloud layer.
  ! It is called at the beginning of each time step.

  use vars
  use domain
  use params
  use tracers, only: tracer
  use microphysics, only: micro_field

  implicit none

  ! mean cloud water limit in kg/kg for cloud layer determination
  real, parameter :: qn0_limit = 1.0e-6
  ! cloud water limit in kg/kg for cloudy grid points 
  real, parameter :: qc_limit = 1.0e-18
  ! parameter controlling the minimum tracer concentration anomaly
  ! in a plume grid point, in terms of fraction of
  ! the horizontal standard deviation 
  real, parameter :: tr_frac = 0.5
  ! the relaxation time scale in number of time steps
  ! for our standard BOMEX run, dt = 0.3 s
  real, parameter :: relax_steps = 1.0 ! instantaneous homogenization
  ! real, parameter :: relax_steps = 10.0 ! 3 s homogenization
  ! real, parameter :: relax_steps = 100.0 ! 30 s homogenization
  ! real, parameter :: relax_steps = 300.0 ! 90 s homogenization
  ! real, parameter :: relax_steps = 1000.0 ! 300 s homogenization
  ! real, parameter :: relax_steps = 6000.0 ! 1800 s homogenization

  ! env grid point flag
  logical :: l_env(nx,ny,nzm)

  ! fields before homogenization for output
  real :: qt_before(nx,ny,nzm), t_before(nx,ny,nzm)

  ! horizontal mean and standard deviation of tracer concentration
  real :: tr0(nzm), tr_sd(nzm)
  ! minimum tracer concentration required for a plume grid point
  real :: tr_min(nzm)

  ! env counts and means
  real :: env_counts(nzm), mqt_env(nzm), mtabs_env(nzm)
  ! column cloud base counts
  real :: ccb_counts(nzm)

  ! kcb is the domain-wide cloud-base level,
  ! based on the column-wise cloud base occurrence frequency
  integer :: kcb
  ! top and bottom of the lowest layer of cloud (in terms of qn0)
  ! all cases we examine will be single-cloud-layer cases
  integer :: cl_top, cl_base
  ! top and bottom of the homogenization layer
  integer :: hl_top, hl_base

  ! variables for collecting statistics
  real(8) coef, coef1, buffer(nzm,3), buffer1(nzm,3)

  ! local logicals
  logical :: l_c
  integer :: i, j, k

  ! for 3d output
  character(len=120) :: filename
  character(len=80) :: long_name
  character(len=8) :: name
  character(len=10) :: timechar
  character(len=4) :: rankchar
  character(len=5) :: sepchar
  character(len=6) :: filetype
  character(len=10) :: units
  character(len=12) :: c_z(nzm), c_p(nzm), c_dx, c_dy, c_time
  integer, parameter :: nfields = 5 ! number of 3d output fields
  real(4) :: tmp(nx,ny,nzm) ! temporary array for 3d output

  ! initialize local variables
  env_counts(:) = 0.0
  mqt_env(:) = 0.0
  mtabs_env(:) = 0.0

  if (nstep .gt. nstep_homo1 .and. nstep .le. nstep_homo2) then

    ! calculate mean and variance of tracer concentration in the subdomain
    coef = 1./float(nx*ny)
    do k = 1, nzm
      tr0(k) = 0.0
      tr_sd(k) = 0.0
      do i = 1, nx
        do j = 1, ny
          tr0(k) = tr0(k) +  tracer(i,j,k,1)
          tr_sd(k) = tr_sd(k) + tracer(i,j,k,1)**2
        enddo 
      enddo
      tr0(k) = tr0(k)*coef
      tr_sd(k) = tr_sd(k)*coef
    enddo
    ! calculate the mean and variance over the entire horizontal domain
    if(dompi) then
      coef1 = 1./float(nsubdomains)
      do k = 1, nzm
        buffer(k,1) = tr0(k)
        buffer(k,2) = tr_sd(k)
      enddo
      call task_sum_real8(buffer,buffer1,nzm*3)
      do k = 1, nzm
        tr0(k) = buffer1(k, 1)*coef1
        tr_sd(k) = buffer1(k, 2)*coef1
      enddo
    endif
    tr_min(:) = 0.0
    do k = 1, nzm
      tr_sd(k) = tr_sd(k) - tr0(k)**2
      if (tr_sd(k) .gt. 0.0) then
        tr_sd(k) = sqrt(tr_sd(k))
      else
        tr_sd(k) = 0.0
      endif ! tr_sd(k) .gt. 0.0
      ! calculate tr_min following Dawe and Austin (2012, ACP)
      ! and the cloud_tracker implementation
      if (k .eq. 1) then
        tr_min(k) = 0.05*tr_sd(k)
      else
        tr_min(k) = (tr_min(k-1)*float(k-1) + 0.05*tr_sd(k))/float(k)
      endif
    enddo

    ! Determine cl_top and cl_base using qn0
    ! qn0 has just been updated in diagnose()
    ! at the end of the last time step.
    cl_top = 0
    cl_base = 0
    l_c = .False.
    do k = 1, nzm
      if (l_c) then
        if (qn0(k) .gt. qn0_limit) then
          cl_top = k
        else
          exit
        endif
      else
        if (qn0(k) .gt. qn0_limit) then
          l_c = .True.
          cl_base = k
          cl_top = k
        endif
      endif
    enddo

    ! set hl_top to a fixed level to avoid variations of cloud top height
    ! hl_top = 100 ! ~2.5 km with dz = 25 m
    ! hl_top = 50 ! ~ 1.25 km
    ! hl_top = 38 ! ~ (25 (cl_base mean) + 50) / 2
    hl_top = 32 
    ! hl_base = 33
    ! hl_base = 29 
    hl_base = cl_base
    ! set hl_base to be the mid of the cloud layer
    ! if ((hl_top - cl_base) .ge. 2) hl_base = floor(cl_base + (hl_top - cl_base)/2.0)

    ! Find, kcb, the level where the column-wise cloud base occuring frequency
    ! maximizes, i.e., the level where the column-wise cloud base occurs most often. 
    ! First sample through the subdomain
    ccb_counts(:) = 0.0
    do i = 1, nx
      do j = 1, ny
        kcb = 0
        do k = nzm, 1, -1
          if ((qcl(i,j,k)+qci(i,j,k)) .gt. qc_limit) kcb = k
        enddo
        if (kcb .gt. 0) ccb_counts(kcb) = ccb_counts(kcb) + 1.0
      enddo 
    enddo 
    ! Then for the entire horizontal domain
    if(dompi) then
      do k = 1, nzm
        buffer(k,1) = ccb_counts(k)
      enddo
      call task_sum_real8(buffer,buffer1,nzm*3)
      do k = 1, nzm
        ccb_counts(k) = buffer1(k, 1)
      enddo
    endif ! dompi
    kcb = 0 
    do k = 1, nzm
      if (kcb .eq. 0) then
        if (ccb_counts(k) .gt. 0.5) kcb = k
      else
        if (ccb_counts(k) .gt. ccb_counts(kcb)) kcb = k 
      endif ! kcb .eq. 0
    enddo

    ! store the before homogenization fields
    do k = 1, nzm
      do j = 1, ny
        do i = 1, nx
          qt_before(i,j,k) = micro_field(i,j,k,1)
          t_before(i,j,k) = t(i,j,k)
        enddo
      enddo 
    end do

    ! only homogenize the levels at and above kcb
    ! so that we are only homogenizing the environment
    ! in the layers where most clouds have already formed.
    ! (actually most of the time kcb is lower than cl_base)
    if (kcb .gt. 1) then
      if (hl_base .lt. kcb) hl_base = kcb 
    endif ! kcb .gt. 1
    ! if (hl_base .lt. 50) hl_base = 50

    ! calculate env counts and mean qt and tabs within the subdomain
    do k = hl_base, hl_top
      do i = 1, nx
        do j = 1, ny
          l_env(i,j,k) = .False.
          ! a grid point is in the environment if it is not cloudy (qcl+qci .lt. 1.0e-18)
          ! and also not in the plume.
          if (((qcl(i,j,k)+qci(i,j,k)) .lt. qc_limit) .and. &
              ! following Dawe and Austin (2012, ACP)
              (tracer(i,j,k,1) .lt. max((tr0(k)+tr_frac*tr_sd(k)), tr_min(k)))) then
            l_env(i,j,k) = .True.
            env_counts(k) = env_counts(k) + 1.0
            ! qt in micro_field(:,:,:,1) in M2005
            ! all the qts here should be just qv because qcl+qci < qc_limit
            mqt_env(k) = mqt_env(k) + micro_field(i,j,k,1)
            ! tabs, just diagnosed in diagnose()
            mtabs_env(k) = mtabs_env(k) + tabs(i,j,k)
          endif
        enddo 
      enddo
    enddo
    ! calculate the env means over the entire horizontal domain
    if(dompi) then
      do k = 1, nzm
        buffer(k,1) = env_counts(k)
        buffer(k,2) = mqt_env(k)
        buffer(k,3) = mtabs_env(k)
      enddo
      call task_sum_real8(buffer,buffer1,nzm*3)
      do k = hl_base, hl_top
        env_counts(k) = buffer1(k, 1)
        mqt_env(k) = buffer1(k, 2)
        mtabs_env(k) = buffer1(k, 3)
      enddo
    endif ! dompi

    do k = hl_base, hl_top
        mqt_env(k) = mqt_env(k)/env_counts(k)
        mtabs_env(k) = mtabs_env(k)/env_counts(k)
    enddo

    ! smooth out the prognostic variables in the environment
    do k = hl_base, hl_top
      do i = 1, nx
        do j = 1, ny
          if (l_env(i,j,k)) then
            micro_field(i,j,k,1) = (mqt_env(k) + &
              micro_field(i,j,k,1) * (relax_steps - 1.0))/relax_steps
            ! we only homogenize actual temperature
            ! the part of 't' associated with potential energy
            ! and latent heat are not touched
            ! t(i,j,k) = t(i,j,k) - (tabs(i,j,k) - mtabs_env(k))/relax_steps
          endif ! l_env(i,j,k)
        enddo
      enddo
    enddo

    ! output stats on masterproc
    if (masterproc) then
      if (no_ehe_file) then
        open(168, file='./OUT_STAT/ehe_stats.ascii', status='unknown', &
             form='formatted')
        no_ehe_file = .False.
      else
        open(168, file='./OUT_STAT/ehe_stats.ascii', status='unknown', &
            form='formatted', position='append')
      endif ! no_ehe_file
      write(168, '(7i10)') nstep, nzm, kcb, cl_base, cl_top, hl_base, hl_top
      do k = 1, nzm 
        write(168, '(7e20.12)') &
             tr0(k), tr_sd(k), tr_min(k), &
             env_counts(k), ccb_counts(k), &
             mtabs_env(k), mqt_env(k)
      enddo
      close(168)
    endif ! masterproc

    ! output before and after homogenization fields
    if(mod(nstep,nsave3D).eq.0.and.nstep.ge.nstep_homo1.and.nstep.le.nstep_homo2 ) then
      ! create the file name 
      ! and open the file for writing
      if(masterproc.or.output_sep) then
        if(output_sep) then
          write(rankchar,'(i4)') rank
          sepchar="_"//rankchar(5-lenstr(rankchar):4)
        else
          sepchar=""
        end if ! output_sep
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
          filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_homo_'// &
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
          end if ! save3Dbin
          if(save3Dsep) then
            filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_homo_'// &
              rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
            open(46,file=filename,status='unknown',form='unformatted')	
          else
            filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_homo_'// &
              rankchar(5-lenstr(rankchar):4)//filetype//sepchar
            if(nrestart.eq.0.and.notopened3D_homo) then
              open(46,file=filename,status='unknown',form='unformatted')	
            else
              open(46,file=filename,status='unknown', &
                      form='unformatted', position='append')
            end if
            notopened3D_homo =.false.
          end if ! save3Dsep
        end if ! RUN3D
        ! write the header information on masterproc
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

      ! write the fields
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            tmp(i,j,k)=qt_before(i,j,k)*1.e3
          end do
        end do
      end do
      name='QT_BEFOR'
      long_name='Total water mixing ratio before homogenization'
      units='g/kg'
      call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                      save3Dbin,dompi,rank,nsubdomains)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            tmp(i,j,k)=t_before(i,j,k)
          end do
        end do
      end do
      name='TL_BEFOR'
      long_name='TL (prognostic) before homogenization'
      units='K'
      call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                      save3Dbin,dompi,rank,nsubdomains) 
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            tmp(i,j,k)=micro_field(i,j,k, 1)*1.e3
          end do
        end do
      end do
      name='QT_AFTER'
      long_name='Total water mixing ratio after homogenization'
      units='g/kg'
      call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                      save3Dbin,dompi,rank,nsubdomains)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            tmp(i,j,k)=t(i,j,k)
          end do
        end do
      end do
      name='TL_AFTER'
      long_name='TL (prognostic) after homogenization'
      units='K'
      call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                      save3Dbin,dompi,rank,nsubdomains) 
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if (l_env(i,j,k)) then
              tmp(i,j,k)=1.0
            else
              tmp(i,j,k)=0.0
            end if
          end do
        end do
      end do
      name='L_ENV'
      long_name='Environment grid point flag (1.0 for env)'
      units=''
      call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                      save3Dbin,dompi,rank,nsubdomains) 

      ! close out
      call task_barrier()
      if (masterproc) close(46)
    end if ! mod(nstep,nsave3D).eq.0.and.nstep.ge.nstep_homo1.and.nstep.le.nstep_homo2
  endif ! nstep .gt. nstep_homo1 .and. nstep .le. nstep_homo2

  return
end subroutine homogenize_cld_env