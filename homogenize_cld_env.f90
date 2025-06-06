subroutine homogenize_cld_env

  use vars
  use tracers, only: tracer
  use microphysics, only: micro_field

  implicit none

  ! cloud water limit in kg/kg 
  real, parameter :: qn_limit = 1.0e-6

  ! flags for homogenization
  logical :: l_env(nx, ny, nzm)

  ! horizontal mean and standard deviation of tracer concentration
  real :: tr0(nzm), tr_sd(nzm)

  ! env counts and means
  real :: env_counts(nzm), mqt_env(nzm), mtabs_env(nzm)
  ! column cloud base counts
  real :: ccb_counts(nzm)

  ! kcb is the domain-wide cloud-base level,
  ! also the lowest level for homogenization.
  integer :: kcb
  ! top and bottom of the lowest layer of cloud (in terms of qn0)
  ! all cases we examine will be single-cloud-layer cases
  integer :: cl_top, cl_base
  ! top and bottom of the homogenization layer
  integer :: hl_top, hl_base, hl_top_new

  ! variables for collecting statistics
  real(8) coef, coef1, buffer(nzm,3), buffer1(nzm,3)

  ! local logicals
  logical :: l_c
  integer :: i, j, k

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
    do k = 1, nzm
      tr_sd(k) = tr_sd(k) - tr0(k)**2
      if (tr_sd(k) .gt. 0.0) then
        tr_sd(k) = sqrt(tr_sd(k))
      else
        tr_sd(k) = 0.0
      endif ! tr_sd(k) .gt. 0.0
    enddo

    ! Determine cl_top and cl_base using qn0
    cl_top = 0
    cl_base = 0
    l_c = .False.
    do k = 1, nzm
      if (l_c) then
        if (qn0(k) .gt. qn_limit) then
          cl_top = k
        else
          exit
        endif
      else
        if (qn0(k) .gt. qn_limit) then
          l_c = .True.
          cl_base = k
          cl_top = k
        endif
      endif
    enddo

    ! set hl_top to a fixed level to avoid variations of cloud top height
    hl_top = 100 ! ~2.5 km with dz = 25 m
    hl_base = cl_base
    ! set hl_base to be the mid of the cloud layer
    ! if ((hl_top - cl_base) .ge. 2) hl_base = floor(cl_base + (hl_top - cl_base)/2.0)

    ! Find the level where the column-wise cloud base occuring frequency
    ! maximizes, i.e., the level where the column-wise cloud base occurs most often. 
    ! First sample through the subdomain
    ccb_counts(:) = 0.0
    do i = 1, nx
      do j = 1, ny
        kcb = 0
        do k = nzm, 1, -1
          if ((qcl(i,j,k)+qci(i,j,k)) .gt. 1.0e-18) kcb = k
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

    ! only homogenize the levels at and above kcb
    if (kcb .gt. 1) then
      if (hl_base .lt. kcb) hl_base = kcb 
    endif ! kcb .gt. 1

    ! calculate env counts and mean qt and tabs within the subdomain
    do k = hl_base, hl_top
      do i = 1, nx
        do j = 1, ny
          l_env(i,j,k) = .True.
          ! not env ?
          if ((qcl(i,j,k)+qci(i,j,k) .gt. 0.0) .or. &
              ((tracer(i,j,k,1)-tr0(k)) .gt. tr_sd(k))) then
                l_env(i,j,k) = .False.
          endif
          if (l_env(i,j,k)) then
            env_counts(k) = env_counts(k) + 1.0
            ! qt in micro_field(:,:,:,1) in M2005
            ! all the qts here should be just qv because qn .le. 0.0
            mqt_env(k) = mqt_env(k) + micro_field(i,j,k,1)
            ! tabs, just diagnosed in diagnose()
            mtabs_env(k) = mtabs_env(k) + tabs(i,j,k)
          endif ! l_env(i,j,k)
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
    ! adjust hl_top in case it is too high and there are no more env/cloud pts
    do k = hl_base, hl_top
      if (env_counts(k) .gt. 0.5) then
        hl_top_new = k
      else
        exit
      endif ! env_counts(k) .le. 0.5
    enddo

    do k = hl_base, hl_top_new
        mqt_env(k) = mqt_env(k)/env_counts(k)
        mtabs_env(k) = mtabs_env(k)/env_counts(k)
    enddo

    ! smooth out the prognostic variables in the environment
    do k = hl_base, hl_top_new
      do i = 1, nx
        do j = 1, ny
          if (l_env(i,j,k)) then
            micro_field(i,j,k,1) = mqt_env(k)
            ! we only homogenize actual temperature or potential temperature
            ! the part of TL associated with latent heat is untouched
            t(i,j,k) = t(i,j,k) - tabs(i,j,k) + mtabs_env(k)
          endif ! l_env(i,j,k)
        enddo
      enddo
    enddo

    ! output on masterproc
    if (masterproc) then
      if (no_ehe_file) then
        open(168, file='./OUT_STAT/ehe_stats.ascii', status='unknown', &
             form='formatted')
        no_ehe_file = .False.
      else
        open(168, file='./OUT_STAT/ehe_stats.ascii', status='unknown', &
            form='formatted', position='append')
      endif ! no_ehe_file
      write(168, '(8i10)') nstep, nzm, kcb, cl_base, cl_top, hl_base, hl_top, hl_top_new
      do k = 1, nzm 
        write(168, '(6e20.12)') &
             tr0(k), tr_sd(k), &
             env_counts(k), ccb_counts(k), &
             mtabs_env(k), mqt_env(k)
      enddo
      close(168)
    endif ! masterproc

  endif ! nstep .gt. nstep_homo1 .and. nstep .le. nstep_homo2

  return
end subroutine homogenize_cld_env