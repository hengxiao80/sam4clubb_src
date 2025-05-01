subroutine homogenize_cld_env

  use vars
  use tracers, only: tracer
  use microphysics, only: micro_field

  implicit none
  ! cloud water limit in kg/kg 
  real, parameter :: qn_limit = 1.0e-6

  ! flags for homogenization
  logical :: l_homo(nzm)
  logical :: l_env(nx, ny, nzm)

  ! horizontal mean and standard deviation of tracer concentration
  real :: tr0(nzm), tr_sd(nzm)

  ! env counts and means
  real :: env_counts(nzm), mqt_env(nzm), mtabs_env(nzm)

  ! variables for collecting statistics
  real(8) coef, coef1, buffer(nzm,3), buffer1(nzm,3)

  integer :: i, j, k

  if (nstep .ge. nstep_homo) then

    ! calculate mean and variance of tracer concentration in the subdomain
    coef = 1./float(nx*ny)
    do k = 1, nzm
      tr0(k) = 0.0
      tr_sd(k) = 0.0
      do i = 1, nx
        do j = 1, ny
          tr0(k) = tr0(k) +  tracer(i,j,k,1)
          tr_sd(k) = tr_sd(k) + tracer(i,j,k,1)**2
        end do 
      end do
      tr0(k) = tr0(k)*coef
      tr_sd(k) = tr_sd(k)*coef
    end do
   
    ! calculate the mean and variance over the entire horizontal domain
    if(dompi) then
      coef1 = 1./float(nsubdomains)
      do k = 1, nzm
        buffer(k,1) = tr0(k)
        buffer(k,2) = tr_sd(k)
      end do
      call task_sum_real8(buffer,buffer1,nzm*nb)
      do k = 1, nzm
        tr0(k) = buffer1(k, 1)*coef1
        tr_sd(k) = buffer1(k, 2)*coef1
        tr_sd(k) = tr_sd(k) - tr0(k)**2
        if (tr_sd(k) .gt. 0.0) then
          tr_sd(k) = sqrt(tr_sd(k))
        else
          tr_sd(k) = 0.0
        end if ! tr_sd(k) .gt. 0.0
      end do
    end if ! dompi

    do k = 1, nzm
      ! First, decide here whether we need to homogenize at this level.

      ! Option #1: Homogenize whenever/wherever mean cloud water qn0 
      ! (cloud liquid+ice for M2005) exceeds the limit.
      ! We should be able to use qn0 directly because it has just been updated
      ! in diagnose() right above the call to this subroutine in main().
      if (qn0(k) .gt. qn_limit) l_homo(k) = .True.

      ! Option #2: Homogenize only where <w'tv'> is negative.
      ! The design here tries to homogenize only the so-called "transition 
      ! layer" (Albright et al. 2023, JAS).

      env_counts(k) = 0.0
      mqt_env(k) = 0.0
      mtabs_env(k) = 0.0
      if (l_homo(k)) then
        do i = 1, nx
          do j= 1, ny
            l_env(i,j,k) = .True.
            ! cloud?
            if ((qcl(i,j,k)+qci(i,j,k) .gt. 0.0) .or. &
                ((tracer(i,j,k,1)-tr0(k)) .gt. tr_sd(k))) then
                  l_env(i,j,k) = .False.
            end if
            if (l_env(i,j,k)) then
              env_counts(k) = env_counts(k) + 1.0
              ! qt in micro_field(:,:,:,1) in M2005
              ! all the qts here should be just qv because qn .le. 0.0
              mqt_env(k) = mqt_env(k) + micro_field(i,j,k,1)
              ! tabs, just diagnosed in diagnose()
              mtabs_env(k) = mtabs_env(k) + tabs(i,j,k)
            end if ! l_env(i,j,k)
          end do 
        end do
      end if ! l_homo(k)
    enddo

    ! calculate the env means over the entire horizontal domain
    if(dompi) then
      do k = 1, nzm
        buffer(k,1) = env_counts(k)
        buffer(k,2) = mqt_env(k)
        buffer(k,3) = mtabs_env(k)
      end do
      call task_sum_real8(buffer,buffer1,nzm*nb)
      do k = 1, nzm
        if (l_homo(k)) then
          env_counts(k) = buffer1(k, 1)
          mqt_env(k) = buffer1(k, 2)
          mtabs_env(k) = buffer1(k, 3)
          ! very unlikely to have zero env point, but ...
          if (env_counts(k) .gt. 0.5) then
            mqt_env(k) = mqt_env(k)/env_counts(k)
            mtabs_env(k) = mtabs_env(k)/env_counts(k)
          else
            ! no need to homogenize
            l_homo(k) = .False.
          end if ! env_counts(k) .gt. 0.5
        end if ! l_homo(k)
      end do
    end if ! dompi

    ! update the prognostic variables for the smoothing
    do k = 1, nzm
      if (l_homo(k)) then
        do i = 1, nx
          do j = 1, ny
            if (l_env(i,j,k)) then
              micro_field(i,j,k,1) = mqt_env(k)
              ! we only homogenize actual temperature or potential temperature
              ! the part due to latent heat is untouched
              t(i,j,k) = t(i,j,k) - tabs(i,j,k) + mtabs_env(k)
            end if ! l_env(i,j,k)
          end do
        end do
      end if ! l_homo(k)
    end do

  endif ! nstep .ge. nstep_homo

end subroutine homogenize_cld_env