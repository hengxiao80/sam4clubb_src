MODULE hoppel_transfer



contains  

REAL FUNCTION prob_lognorm_lt_Dmax(Dg, Sg, Dmax)
  implicit none
  real, intent(in) :: Dg ! geometric mean diameter, m
  real, intent(in) :: Sg ! standard deviation of lognormal size distribution = log(sigmag)
  real, intent(in) :: Dmax ! upper limit (diameter, m) for lognormal integration

  prob_lognorm_lt_Dmax = &
       0.5*(1. + erf( log(Dmax/Dg) / ( sqrt(2.)*Sg ) ) )
END FUNCTION prob_lognorm_lt_Dmax

REAL FUNCTION prob_lognorm_gt_Dmin(Dg, Sg, Dmin)
  implicit none
  real, intent(in) :: Dg ! geometric mean diameter, m
  real, intent(in) :: Sg ! standard deviation of lognormal size distribution = log(sigmag)
  real, intent(in) :: Dmin ! upper limit (diameter, m) for lognormal integration

  prob_lognorm_gt_Dmin = &
       0.5*(1. - erf( log(Dmin/Dg) / ( sqrt(2.)*Sg ) ) )
END FUNCTION prob_lognorm_gt_Dmin

REAL FUNCTION lognorm_conc_at_diam(N, Dg, Sg, Dtarget)
  implicit none
  ! returns the pdf value at Dtarget
  real, intent(in) :: N ! total number mixing ratio, #/kg
  real, intent(in) :: Dg ! geometric mean diameter, m
  real, intent(in) :: Sg ! standard deviation of lognormal size distribution = log(sigmag)
  real, intent(in) :: Dtarget ! diameter at which lognormal pdf is evaluated, m
  real ::  pi = 3.1415926535

  lognorm_conc_at_diam = &
       N * exp( -0.5 * (log(Dg/Dtarget)/Sg)**2)/(SQRT(2.*pi)*Sg)

END FUNCTION lognorm_conc_at_diam

REAL FUNCTION modal_diameter(N, M, Sg, rho_aerosol)
  implicit none
  ! returns the modal diameter
  real, intent(in) :: N, M, Sg, rho_aerosol
  real ::  pi = 3.1415926535
  modal_diameter = ((1/rho_aerosol) * (M/N) * (6/pi) * exp(-4.5 * Sg**2))**(1./3.)

END FUNCTION modal_diameter  

REAL FUNCTION mass_mixing_ratio(N,D,sg, rho_aero)
  implicit none
  real, intent(in) :: N,D,sg,rho_aero
  real :: pi =  3.1415926535
  mass_mixing_ratio = (pi/6)*rho_aero * (D**3) * exp(4.5*sg**2) * N

END FUNCTION mass_mixing_ratio  

REAL FUNCTION mass_fraction(nfrac, sigma)
  !    Compute mass fraction of lognormal distribution smaller than some size
  !   given a number fraction smaller than the same size
  implicit none
  double precision, external :: derfi
  double precision :: dnfrac, tmp1, xx, arg_erf, dmfrac
  real, intent(in) :: nfrac, sigma
  real :: sqrt2 = 2.**0.5

  ! use double precision throughout to preserve more of the tails
  !   of the error function and inverse error function as the argument
  !   goes to -1 and 1.
  dnfrac = nfrac ! derfi needs double precision input
  tmp1 = 2.*( dnfrac - 0.5 ) ! remap to [-1,1]

  ! bound the values away from -1 and 1:  -1 + 1e-8 < tmp1 < 1 - 1e-8
  tmp1 = MAX(-1.+SQRT(EPSILON(tmp1)), MIN( 1.-SQRT(EPSILON(tmp1)), tmp1 ) )
  xx = sqrt(2.) * derfi( tmp1 )
  arg_erf = (-3.*sigma+ xx)/sqrt2 
  dmfrac = 0.5 + 0.5 * erf( arg_erf ) ! truncated third moment
  mass_fraction = dmfrac

!!$  if (nfrac <= 0.) then
!!$     mass_fraction = 0.
!!$  else if (nfrac < 0.999999) then
!!$     xx = sqrt2 * derfi(2. * (dnfrac - 0.5))
!!$     mass_fraction = 0.5 + 0.5 * erf((-3.*sigma+ xx)/sqrt2)
!!$  else
!!$     mass_fraction = 1.
!!$  end if   

END FUNCTION mass_fraction

SUBROUTINE hoppel_aitken_accum_transfer(N0, N1, M0, M1, sg0, sg1, Dc0, Dc1, &
     rho_aero0, rho_aero1, Ntrans, Mtrans)
  ! Compute number and mass of transfer from aitken to accumulation
  ! if aitken number at critical radius is larger than accumulation number at
  !   critical radius
  ! All particles transfered from Aitken have diameter Dc0

  ! Does not check that inputs are nonzero, N0 and N1. (diameter is
  ! computed and n in denomintor).

  implicit none

  
  real, intent(in) :: N0, N1  !number in /kg   0 is 'Aitken' 1 is 'accumulation'
  real, intent(in) :: M0, M1  ! mass in kg/kg
  real, intent(in) :: sg0, sg1  !natural log of sigma shape parameter
  real, intent(in) :: Dc0, Dc1  !critical diameter of activation (from ARG)
  real, intent(in) :: rho_aero0 ! rho_aero1  aerosol density
  real, intent(in) :: rho_aero1
  real, intent(out) :: Ntrans
  real, intent(out) :: Mtrans  ! in /kg and kg/kg

  real :: Dg0, Dg1
  real :: N0n, N1n, M0n, M1n, NT, MT, kT 
  real :: T, Tnew  ! Amoutn transferred scaled by N0
  real :: eta0, eta1 ! number distribution functions of modes
  real :: deta0_DT, deta1_DT
  
  real ::  pi = 3.1415926535
  integer, parameter :: maxiter = 15
  real, parameter :: deltaN_conv = 50. ! /m^3

  integer ::  count = 0
  
  ! N0, MO, N1, M1 are unchanged -> N0 = N00 in writeup

  
  
  N0n = N0
  N1n = N1
  M0n = M0
  M1n = M1
  
  kT = (pi/6.) * rho_aero0*Dc0**3 
  NT = 0.
  MT = 0.
  
  T = 0.
  count = 0
  DO
     count = count+1
     Dg0 = modal_diameter(N0n, M0n, sg0, rho_aero0)
     Dg1 = modal_diameter(N1n, M1n, sg1, rho_aero1)
!     print*, 'Dg0, Dg1', Dg0*1.e6, Dg1*1.e6
     eta0 = N0n*exp(-0.5 * (log(Dg0/Dc0)/sg0)**2)/(SQRT(2.*pi)*sg0)
     eta1 = N1n*exp(-0.5 * (log(Dg1/Dc1)/sg1)**2)/(SQRT(2.*pi)*sg1)
!     print*, 'eta0, eta1', eta0/1.e6, eta1/1.e6
     
     IF ((eta0-eta1).LE.0.AND.count.EQ.1) THEN
        EXIT
     END IF   
     IF (abs(eta0-eta1).LT.deltaN_conv) THEN
        EXIT
     END IF

     ! N0n nad N1n should not be identically zero
     deta0_DT = -eta0 * (N0/N0n) * (1. - 1./3. * sg0**(-2) * log(Dg0/Dc0)*(exp(-4.5 * sg0**2)*(Dc0/Dg0)**3 - 1.))
     deta1_DT = eta1 * (N0/N1n) * (1. - 1./3. * sg1**(-2) * log(Dg1/Dc1)*(exp(-4.5 * sg1**2)*(Dc1/Dg1)**3 - 1.))   

     Tnew = T - (eta0-eta1)/(deta0_dT - deta1_dT)

     IF (Tnew < 0) THEN
!!$        !bloss(2023-11-11): Take a bisection step
!!$        !  with NT = 0 and NT from last iteration
!!$        Tnew = 0.5*T
!!$
        print*, 'Aitken transfer negative. Count= ', count
        print*, 'N0,N1,M0,M1,Dc0,Dc1=', N0, N1, M0, M1, Dc0, Dc1
        NT = 0.
        MT = 0.
        EXIT
     END IF
     
     NT = Tnew * N0

     !put limiters here
     IF (NT.GE.N0) THEN
        DO WHILE (NT.GE.N0) 
           NT = NT/2.
        END DO
        Tnew = NT/N0
     END IF
     MT = NT*kT

     IF (MT.GE.M0) THEN 
        DO WHILE (MT.GE.M0)
           MT = MT/2
        END DO
        NT = MT/kT
    END IF

     IF (maxiter.EQ.count) THEN
        NT = 0.
        MT = 0.
        print*, 'Aitken Transfer Convergence Failure'
        print*, 'N0,N1,M0,M1,Dc0,Dc1=', N0, N1, M0, M1, Dc0, Dc1
        EXIT
     END IF


     
     N0n = N0 - NT
     N1n = N1 + NT
     M0n = M0 - MT 
     M1n = M1 + MT

     T = Tnew
    
!     print*,  '      N0          M0              N1           M1           F         T        NT'
!     print*, N0n/1.e6, M0n*1.e9, N1n/1.e6, M1n*1.e9, (eta0-eta1)/1.e6, Tnew, NT/1.e6
     
  END DO   

  Ntrans = NT
  Mtrans = MT

  
END SUBROUTINE hoppel_aitken_accum_transfer
  
SUBROUTINE hoppel_aitken_accum_transfer_v2(N0, N1, M0, M1, sg0, sg1, Dc0, Dc1, &
     rho_aero0, rho_aero1, Ntrans, Mtrans)
  ! Compute number and mass of transfer from aitken to accumulation
  ! if aitken number at critical radius is larger than accumulation number at
  !   critical radius
  ! All particles transfered from Aitken have diameter Dc0

  ! Does not check that inputs are nonzero, N0 and N1. (diameter is
  ! computed and n in denomintor).

  ! Duplicated/Modified by Peter Blossey, 2023-10-13
  !   - Mass tranferred is mass of the largest Ntrans Aitken aerosols
  !       in the lognormal size distribution, computed using mass_fraction().
  !   - Iterate with secant method to find value of Ntrans/Mtrans that
  !       leads to the accumulation and Aitken mode size distributions
  !       having the same concentration at their respective Dc's.

  implicit none


  real, intent(in) :: N0, N1  !number in /kg   0 is 'Aitken' 1 is 'accumulation'
  real, intent(in) :: M0, M1  ! mass in kg/kg
  real, intent(in) :: sg0, sg1  !natural log of sigma shape parameter
  real, intent(in) :: Dc0, Dc1  !critical diameter of activation (from ARG)
  real, intent(in) :: rho_aero0 ! rho_aero1  aerosol density
  real, intent(in) :: rho_aero1
  real, intent(out) :: Ntrans
  real, intent(out) :: Mtrans  ! in /kg and kg/kg
  
  real :: Dg0, Dg1
  real :: Ntrans_old, Ntrans_new
  real :: eta0, eta1 ! number distribution functions of modes at critical diameters
  real :: deta, deta_old ! eta0 - eta1 (Aitken conc at Dc0 minus Accum. conc at Dc1)
  real :: Nlo, Nhi ! evolving lower/upper bounds for Ntrans
  real :: minM0perN0 ! mass mixing ratio of 1 nanometer diameter aerosol with N0=1.
  
  real ::  pi = 3.1415926535

  ! converged if |eta0-eta1| < MAX(deltaN_conv, eta_frac_conv * max(eta0,eta1) )
  real, parameter :: eta_frac_conv = 0.01 
  real, parameter :: deltaN_conv = 50. ! /m^3
  integer, parameter :: maxiter = 25

  integer ::  count

  logical :: debug = .false.

  ! Return if there is no mass in Aitken mode
  if(M0.le.0.) return

  Ntrans = 0.
  Mtrans = 0.

  ! Check initial gap between N0 and N1 at Dcrit
  Dg0 = modal_diameter(N0, M0, sg0, rho_aero0)
  eta0 = lognorm_conc_at_diam(N0, Dg0, sg0, Dc0)

  Dg1 = modal_diameter(N1, M1, sg1, rho_aero1)
  eta1 = lognorm_conc_at_diam(N1, Dg1, sg1, Dc1)

!!$  if(debug) then
!!$    write(*,*) 'count, Dg0, Dg1, Ntrans, Mtrans, eta0, eta1 = '
!!$    write(*,920) 0, Dg0, Dg1, 0., 0., eta0, eta1
920 format(I8,6E12.4)
!!$  end if

  deta = eta0 - eta1
  if ( (deta < 0) .OR. &
       (abs(deta) < MAX(deltaN_conv, eta_frac_conv*max(eta0,eta1) ) ) ) then
    ! return with no transfer if Naitken < Naccum at Dcrit
    !   of if convergence criteria is satisfied.
    Ntrans = 0
    Mtrans = 0
    return
  end if

  ! save initial values for eta0-eta1, Ntrans and Mtrans
  !  These will be used in the secant method
  deta_old = deta
  Ntrans_old = 0.

  ! make initial guess for Ntrans: equal to twice the Aitken concentration in excess of Dcrit
  Ntrans_new = MIN(N0-1., &
       N0*prob_lognorm_gt_Dmin(Dg0,sg0,Dc0), & ! twice Aitken concentration with D/Dcrit
       0.67*deta) ! (2/3)*(Nait(Dc0) - Naccum(Dc0))

  ! track the lower and upper bounds for Ntrans
  !   lower bound only includes Ntrans for which eta0 > eta1, opposite for upper bound
  Nlo = 0.
  Nhi = N0 - 1. ! Extreme transfers should preserve small Aitken number for numerical robustness
  
  minM0perN0 = mass_mixing_ratio(1.,1.e-9,sg0,rho_aero0)

  ! loop over iterations of secant method.
  do count = 1,maxiter

    ! Secant will predict Ntrans, seeking a value that gives eta0=eta1.
    Ntrans = Ntrans_new

    ! Assume Mtrans is the mass of largest Ntrans Aitken aerosols
    Mtrans = M0*( 1. - mass_fraction(1.-Ntrans/N0, sg0) )

    if(Mtrans.eq.0.) then
       ! The mass_fraction function has difficulty
       !    when Ntrans << N0 and sometimes produces zero
       !    Mtrans when Ntrans is non-zero.  Since this is
       !    non-physical, assume the transferred aerosols
       !    have the critical diameter to produce a
       !    reasonable estimate of Mtrans.
       Mtrans = mass_mixing_ratio(Ntrans,Dc0,sg0,rho_aero0)
!!$       write(*,*) 'Assumed Dc0 size for Mtrans'
    end if

    ! ensure that remaining Aitken mode has a modal diameter >= 1 nanometer
    Mtrans = MIN(M0 - minM0perN0*(N0-Ntrans), Mtrans)
    
    ! Check gap between N0 and N1 at Dcrit
    !   after accounting for predicted transfer
    Dg0 = modal_diameter(N0-Ntrans, M0-Mtrans, sg0, rho_aero0)
    eta0 = lognorm_conc_at_diam(N0-Ntrans, Dg0, sg0, Dc0)
    Dg1 = modal_diameter(N1+Ntrans, M1+Mtrans, sg1, rho_aero1)
    eta1 = lognorm_conc_at_diam(N1+Ntrans, Dg1, sg1, Dc1)

    deta = eta0 - eta1
    if(abs(deta) < MAX( deltaN_conv, eta_frac_conv*max(eta0,eta1) ) ) then
!!$      if(debug) then
!!$         write(*,*) 'count, Dg0, Dg1, Ntrans, Mtrans, eta0, eta1 = '
!!$         write(*,920) 0, Dg0, Dg1, Ntrans, Mtrans, eta0, eta1
!!$      end if

      ! reached convergence, exit loop
      exit
    end if

!!$      if(debug) then
!!$         write(*,*) 'count, Dg0, Dg1, Ntrans, Mtrans, eta0, eta1 = '
!!$         write(*,920) 0, Dg0, Dg1, Ntrans, Mtrans, eta0, eta1
!!$      end if

    ! If not converged, updated guess for Ntrans based on secant method
    !x0 = Ntrans_old
    !f0 = deta_old
    !x1 = Ntrans
    !f1 = deta
    !x2 = x1 - f1*(x1 - x0)/(f1 -f0)
    !Ntrans_new = x2

    ! Increase lower bound if deta is still positive but closer to zero.
    if(deta.gt.0.) Nlo = MAX(Nlo,Ntrans)
    if(deta.lt.0.) Nhi = MIN(Nhi,Ntrans)
    
    Ntrans_new = Ntrans - deta*(Ntrans - Ntrans_old)/(deta - deta_old)
    Ntrans_new = max(Nlo, min(Nhi, Ntrans_new) )

    if(Ntrans_new.eq.Ntrans) then
       ! No change from the last prediction of Ntrans
       !   despite the solution not being converged.
       !   Take a bisection step between Nlo and Nhi,
       !   noting that this could be a huge step if
       !   Nhi == N0.
       Ntrans_new = 0.5*(Nlo + Nhi)
!!$       write(*,*) 'Modified Ntrans_new because Ntrans_new == Ntrans'
    end if
    
    ! update old guesses with values at current step
    Ntrans_old = Ntrans
    deta_old = deta
  end do




  if(count >= maxiter) then
    write(*,*) 'Aitken Transfer Failed to Converge'  
    write(*,*) 'N0,N1=', N0, N1
    write(*,*) 'M0,M1=', M0, M1
    write(*,*) 'Dc0,Dc1=', Dc0, Dc1
    write(*,*) 'Nlo,Ntrans,Nhi=', Nlo, Ntrans, Nhi
    write(*,*) 'count, Dg0, Dg1, Ntrans, Mtrans, eta0, eta1 = '
    write(*,920) count, Dg0, Dg1, Ntrans, Mtrans, eta0, eta1

    ! Only stop the simulation if Ntrans or Mtrans are NaN.
    if(isnan(Ntrans).OR.isnan(Mtrans)) then
       STOP 'in hoppel_transfer_v2'
    end if
 end if


END SUBROUTINE hoppel_aitken_accum_transfer_v2

SUBROUTINE hoppel_aitken_accum_transfer_W2022b(N0, N1, M0, M1, sg0, sg1, Dc0, Dc1, &
     rho_aero0, rho_aero1, Ntrans, Mtrans)
  ! Compute number and mass of transfer from aitken to accumulation
  ! if aitken number at critical radius is larger than accumulation number at
  !   critical radius
  ! All particles transfered from Aitken have diameter MAX(Dg0,Dc0)

  ! Does not check that inputs are nonzero, N0 and N1. (diameter is
  ! computed and n in denomintor).

  ! Duplicated/Modified by Peter Blossey, 2023-10-13
  !   - Choose diameter of transferred aerosols to be MAX(Dv0,Dc0)
  !       rather than Dc0
  !   - Iterate with secant method to find value of Ntrans/Mtrans that
  !       leads to the accumulation and Aitken mode size distributions
  !       having the same concentration at their respective Dc's.

  ! Duplicated/Modified by Peter Blossey, 2023-11-29
  !   - The mass of transferred aerosols is the average of the mass
  !       of the largest Ntrans Aitken aerosols and Ntrans
  !       aerosols with diameter Dc0

  implicit none


  real, intent(in) :: N0, N1  !number in /kg   0 is 'Aitken' 1 is 'accumulation'
  real, intent(in) :: M0, M1  ! mass in kg/kg
  real, intent(in) :: sg0, sg1  !natural log of sigma shape parameter
  real, intent(in) :: Dc0, Dc1  !critical diameter of activation (from ARG)
  real, intent(in) :: rho_aero0 ! rho_aero1  aerosol density
  real, intent(in) :: rho_aero1
  real, intent(out) :: Ntrans
  real, intent(out) :: Mtrans  ! in /kg and kg/kg

  real :: Dg0, Dg1
  real :: Ntrans_old, Ntrans_new
  real :: eta0, eta1 ! number distribution functions of modes at critical diameters
  real :: deta, deta_old ! eta0 - eta1 (Aitken conc at Dc0 minus Accum. conc at Dc1)
  real :: Nlo, Nhi ! evolving lower/upper bounds for Ntrans
  real :: kT ! mass mixing ratio of transferred Aitken aerosol is kT*Ntrans
  real :: minM0perN0 ! mass mixing ratio of 3 nanometer diameter aerosol with N0=1.
  
  real ::  pi = 3.1415926535

  ! converged if |eta0-eta1| < MAX(deltaN_conv, eta_frac_conv * max(eta0,eta1) )
  real, parameter :: eta_frac_conv = 0.01 
  real, parameter :: deltaN_conv = 50. ! /m^3
  integer, parameter :: maxiter = 25

  integer ::  count

  logical :: debug = .false.

  Ntrans = 0.
  Mtrans = 0.

  ! Check initial gap between N0 and N1 at Dcrit
  Dg0 = modal_diameter(N0, M0, sg0, rho_aero0)
  eta0 = lognorm_conc_at_diam(N0, Dg0, sg0, Dc0)

  Dg1 = modal_diameter(N1, M1, sg1, rho_aero1)
  eta1 = lognorm_conc_at_diam(N1, Dg1, sg1, Dc1)

  if(debug) then
    write(*,*) 'count, Dg0, Dg1, Ntrans, M/N, eta0, eta1 = '
    write(*,920) 0, Dg0, Dg1, 0., M0/N0, eta0, eta1
920 format(I8,7E12.4)
  end if

  deta = eta0 - eta1
  if ( (deta < 0) .OR. &
       (abs(deta) < MAX(deltaN_conv, eta_frac_conv*max(eta0,eta1) ) ) ) then
    ! return with no transfer if Naitken < Naccum at Dcrit
    !   of if convergence criteria is satisfied.
    Ntrans = 0
    Mtrans = 0
    return
  end if

  ! save initial values for eta0-eta1, Ntrans and Mtrans
  !  These will be used in the secant method
  deta_old = deta
  Ntrans_old = 0.

  ! make initial guess for Ntrans: equal to Aitken concentration in excess of Dcrit
  Ntrans_new = MIN(N0-1., &
       N0*prob_lognorm_gt_Dmin(Dg0,sg0,Dc0) )

  ! track the lower and upper bounds for Ntrans
  !   lower bound only includes Ntrans for which eta0 > eta1, opposite for upper bound
  Nlo = 0.
  Nhi = N0 - 1. ! Extreme transfers should preserve small Aitken number for numerical robustness
  
  ! Mass of transferred aerosols is the average of Ntrans aerosols with the
  !   critical diameter, Dc0, and the mass of the largest Ntrans Aitken aerosols.
  !   Here, compute the multiplier to compute the mass of aerosols with diameter Dc0.
  kT = (pi/6.) * rho_aero0*Dc0**3

  ! Minimum per-particle Aitken mass after transfer
  minM0perN0 = mass_mixing_ratio(1.,3.e-9,sg0,rho_aero0)

  ! loop over iterations of secant method.
  do count = 1,maxiter

    ! Secant will predict Ntrans, seeking a value that gives eta0=eta1.
    Ntrans = Ntrans_new

    ! Mtrans is the average of the mass of Ntrans aerosols with diameter Dc0
    !   and of the largest Ntrans aerosols in the Aitken size distribution.
    Mtrans = 0.5*( kT*Ntrans + M0*( 1. - mass_fraction(1.-Ntrans/N0, sg0) ) )
    
    ! Prevent transfer from exhausting all mass in Aitken mode.
    !   and ensure that the mean mass of the transferred particles
    !   is at least as big as that of the whole Aitken mode
    Mtrans = MAX(Ntrans*M0/N0, MIN( Mtrans, M0 - (N0-Ntrans)*minM0perN0 ) )

    if(debug.AND.(count.eq.1)) then
       write(*,920) count, Dg0, Dg1, Ntrans, Mtrans/Ntrans
    end if

    ! Check gap between N0 and N1 at Dcrit
    !   after accounting for predicted transfer
    Dg0 = modal_diameter(N0-Ntrans, M0-Mtrans, sg0, rho_aero0)
    eta0 = lognorm_conc_at_diam(N0-Ntrans, Dg0, sg0, Dc0)
    Dg1 = modal_diameter(N1+Ntrans, M1+Mtrans, sg1, rho_aero1)
    eta1 = lognorm_conc_at_diam(N1+Ntrans, Dg1, sg1, Dc1)

    deta = eta0 - eta1
    if(abs(deta) < MAX( deltaN_conv, eta_frac_conv*max(eta0,eta1) ) ) then
      if(debug) then
         write(*,920) count, Dg0, Dg1, Ntrans, Mtrans/Ntrans, eta0, eta1, Nlo
      end if

      ! reached convergence, exit loop
      exit
    end if

    ! If not converged, updated guess for Ntrans based on secant method
    !x0 = Ntrans_old
    !f0 = deta_old
    !x1 = Ntrans
    !f1 = deta
    !x2 = x1 - f1*(x1 - x0)/(f1 -f0)
    !Ntrans_new = x2

    ! Increase lower bound if deta is still positive but closer to zero.
    if(deta.gt.0.) Nlo = MAX(Nlo,Ntrans)
    if(deta.lt.0.) Nhi = MIN(Nhi,Ntrans)
    
    if(debug) then
       write(*,920) count, Dg0, Dg1, Ntrans, Mtrans/Ntrans, eta0, eta1, Nlo
    end if

    Ntrans_new = Ntrans - deta*(Ntrans - Ntrans_old)/(deta - deta_old)

    if( (Ntrans_new.eq.Ntrans) &
         .OR. (Ntrans_new.lt.Nlo) &
         .OR.  (Ntrans_new.gt.Nhi) )then
       ! No change from the last prediction of Ntrans
       !   despite the solution not being converged.
       !   Take a bisection step between Nlo and Nhi,
       !   noting that this could be a huge step if
       !   Nhi == N0.
       Ntrans_new = 0.5*(Nlo + Nhi)
!!$       write(*,*) 'Modified Ntrans_new because Ntrans_new == Ntrans or outside [Nlo,Nhi]'
    end if
    
    ! update old guesses with values at current step
    Ntrans_old = Ntrans
    deta_old = deta
  end do

  if(count >= maxiter) then
    write(*,*) 'Aitken Transfer Convergence Failure'  
    write(*,*) 'N0,N1=', N0, N1
    write(*,*) 'M0,M1=', M0, M1
    write(*,*) 'Dc0,Dc1=', Dc0, Dc1
    write(*,*) 'Nlo,Ntrans,Nhi=', Nlo, Ntrans, Nhi
    write(*,*) 'count, Dg0, Dg1, Ntrans, Mtrans, eta0, eta1 = '
    write(*,920) count, Dg0, Dg1, Ntrans, Mtrans, eta0, eta1
!!$    STOP 'in hoppel_transfer_W2022b'
  end if


END SUBROUTINE hoppel_aitken_accum_transfer_W2022b

END MODULE hoppel_transfer
