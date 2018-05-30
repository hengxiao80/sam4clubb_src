!==================================================================================================
!==================================================================================================
! This module contains subroutines from the
!	DOUBLE-MOMENT VERSION OF LIN MICROPHYSICS SCHEME 
!		by Vaughan Phillips (GFDL, 2006)
!		(developed from single-moment version by Lord et al., 1984, JAS, vol 41, 2836, and
!		Krueger et al, 1995, JAM, vol 34, 281)
! adopted for SAM by Mikhail Ovtchinnikov (PNNL, 2008)
!==================================================================================================
!==================================================================================================
! Modifications for software compatibilty done by dschanen UWM 15 May 2008
! $Id: module_lin_2mom.f90,v 1.2 2008-05-15 21:54:16 dschanen Exp $

 MODULE module_lin_2mom
  

CONTAINS


REAL FUNCTION find_lambdadstar(p_w, fraction_frozen)
implicit none
REAL, INTENT(IN) :: p_w, fraction_frozen
integer :: iter 
 REAL :: fdiff, lam_dstar_zero, lam_dstar_one, lam_dstar_two, a, junk
! REAL :: gammq_frac
a = p_w + 1.;
iter = 0;
if(p_w == 0.) then
lam_dstar_zero = -log(fraction_frozen)
if(ABS(gammq_frac(a, lam_dstar_zero) - fraction_frozen) > 1.e-6) stop 711
find_lambdadstar = lam_dstar_zero
RETURN
else
lam_dstar_zero = p_w
endif
if(gammq_frac(a, lam_dstar_zero) - fraction_frozen < 0.) then 
	lam_dstar_one = lam_dstar_zero/2.
else
	lam_dstar_one = lam_dstar_zero*2.
endif
lam_dstar_two = lam_dstar_one + 1.
do while(ABS(lam_dstar_one - lam_dstar_zero) > 1.e-6) 
fdiff = gammq_frac(a, lam_dstar_one) - gammq_frac(a, lam_dstar_zero)
if(fdiff == 0.) then
	exit
endif
lam_dstar_two = lam_dstar_one - (gammq_frac(a, lam_dstar_one) - fraction_frozen) *( lam_dstar_one -  lam_dstar_zero)/fdiff
if(lam_dstar_two < 0.) lam_dstar_two = 0.
if(lam_dstar_two > 10.*p_w+1.) lam_dstar_two = 10.*p_w+1.

lam_dstar_zero = lam_dstar_one
lam_dstar_one = lam_dstar_two
! print *, gammq_frac(a, lam_dstar_one) - fraction_frozen, iter
iter = iter +1
if(iter > 20) exit 
enddo
if(ABS(gammq_frac(a, lam_dstar_two) - fraction_frozen) > 1.e-5) then
	junk = gammq_frac(a, lam_dstar_two)
	print *,' WARNING:  find_lambdadstar not accurate', ABS(gammq_frac(a, lam_dstar_two) - fraction_frozen)
!	stop 91111
endif
find_lambdadstar = lam_dstar_two
RETURN


END FUNCTION find_lambdadstar





SUBROUTINE reset_aerosols_in_environment(qlz_x, nwz_x, ncn_az_x, &
        dbar_cw, qiz_x, niz_x, nin_az_x, dbar_ci, C_by_rho_x, C_env_x, rho_x, &
	temp_x, N_w_crit_env, qv_x, XLV, XLS, CP)
implicit none
REAL :: T_FRZ_HOM_DEGC
PARAMETER(T_FRZ_HOM_DEGC = -36.)
REAL, INTENT(INOUT) :: qlz_x, nwz_x, ncn_az_x, qiz_x, niz_x, nin_az_x, &
           C_by_rho_x, dbar_cw, dbar_ci, temp_x, qv_x
REAL, INTENT(IN) :: C_env_x, rho_x, N_w_crit_env, XLV, XLS, CP
REAL :: XLVOCP, XLSOCP

XLVOCP = XLV/CP
XLSOCP = XLS/CP


! ENVIRONMENT::  reset # of activated CCN and IN to zero.
	if(qiz_x + qlz_x <= 0.) then
		ncn_az_x = 0.
		qlz_x = 0.
		nwz_x = 0.
                dbar_cw = 0.
		nin_az_x = 0.
  		qiz_x = 0.
		niz_x = 0.
                dbar_ci = 0.
		C_by_rho_x = C_env_x/rho_x;
		if(C_env_x <= 0.) stop 9112
	endif



! NO CLOUD-WATER (below anvil)::  reset # of activated CCN to zero.
	if( (qlz_x <= 0. .or. nwz_x*rho_x < N_w_crit_env*1.e6) .and. temp_x > T_FRZ_HOM_DEGC + 273.15) then
		
		qv_x = qv_x + qlz_x
		temp_x = temp_x - XLVOCP*qlz_x
		ncn_az_x = 0.
		qlz_x = 0.
		nwz_x = 0.
                dbar_cw = 0.
	endif

! NO ICE ::  reset # of activated IN to zero

	if( qiz_x <= 0. .or. niz_x*rho_x < 1. ) then
		qv_x = qv_x + qiz_x
		temp_x = temp_x - XLSOCP*qiz_x
		nin_az_x = 0.
		qiz_x = 0.
		niz_x = 0.
                dbar_ci = 0.
	endif

END SUBROUTINE reset_aerosols_in_environment



SUBROUTINE condense_and_sublime_driver( qv_x, T_x, pre_x, rho_x, qv_zero_x, T_zero_x, &
						p_zero_x, F_q_x, F_T_x,  &
						ql_x, nw_x, ncn_a_x, qi_x, ni_x, nin_a_x, qr_x, qs_x, qg_x, &
						rvapor, a_i, b_i, c_i, d_i, a_w, b_w, c_w, a_r, b_r, rhowater, &
				a_s, b_s, rhosnow, a_g, b_g, rhograupel, LAMBDA_FAC_w,  LAMBDA_FAC_i, DBAR_FAC_w, &
				DBAR_FAC_i, N0_FAC_w, N0_FAC_i,  p_w, p_i, GRAUPEL_FAC, RAIN_FAC, SNOW_FAC, CP, &
						PIE, EPS, XLS, XLV, Rair, ESIT, &
						DESIT, ESWT, DESWT, temp_K, dt, super_WRF, s_FD_test,  s_FD_diff, N_0S, N_0G, N_0R)
IMPLICIT NONE

REAL, INTENT(IN):: qv_zero_x, T_zero_x, pre_x, rho_x, p_zero_x, rvapor, a_i, b_i, c_i, d_i, a_w, b_w, c_w, a_r, b_r, rhowater, &
					a_s, b_s, rhosnow, a_g, b_g, rhograupel, LAMBDA_FAC_w,  LAMBDA_FAC_i,  p_w, p_i, CP, &
						PIE, EPS, XLS, XLV, Rair, dt, DBAR_FAC_w, DBAR_FAC_i, N0_FAC_w, N0_FAC_i, &
					GRAUPEL_FAC, RAIN_FAC, SNOW_FAC, N_0S, N_0G, N_0R

DOUBLE PRECISION, INTENT(IN)::  ESIT(150001), DESIT(150001), ESWT(150001), DESWT(150001), temp_K(150001)

REAL, INTENT(INOUT):: qv_x, T_x, F_q_x, F_T_x, ql_x, nw_x, ncn_a_x, qi_x, ni_x, nin_a_x, qr_x, qs_x, qg_x

REAL, INTENT(OUT)::	super_WRF, s_FD_test,  s_FD_diff

INTEGER :: ikj, NUM_CYCLES

DOUBLE PRECISION:: qv_zero_num, T_zero_num, qv_zero_hi, t_zero_hi, p_zero_hi, F_q_hi, F_T_hi, F_p_hi, qv_hi, T_hi, &
	 ql_hi, qr_hi, qi_hi, qs_hi, qg_hi, nw_hi, ni_hi, ncn_a_hi, nin_a_hi, dtfine

REAL :: super_WRF_fine, s_numerical

INTRINSIC DBLE, REAL

		F_q_x = ( qv_x - qv_zero_x)/dt
		F_T_x = ( t_x - T_zero_x)/dt

		if((ql_x > 0. .or. qi_x > 0. .or. qr_x > 0. .or. qs_x > 0. .or. qg_x > 0.) .and. qv_x > 0.) then 

				F_p_hi = (pre_x - p_zero_x)/dt
                	        if(qi_x > 1. .or. ql_x > 1.) stop 137

				qv_zero_hi = DBLE(qv_zero_x)
				t_zero_hi = DBLE(t_zero_x)
				p_zero_hi = DBLE(p_zero_x)
				qv_zero_num = qv_zero_hi
				T_zero_num = t_zero_hi
				F_q_hi = DBLE(F_q_x)
				F_T_hi = DBLE(F_T_x)


				num_cycles = NUMBER_SUBCYCLES(nw_x*rho_x, ni_x*rho_x, F_T_x)


				num_cycles = 1;
				if(num_cycles < 1) num_cycles = 1;
				super_WRF_fine = 0.; 
				dtfine = DBLE(dt)/DBLE(NUM_CYCLES)

				qv_hi = DBLE(qv_x)
 				T_hi = DBLE(t_x)
				ql_hi = DBLE(ql_x)
				qr_hi = DBLE(qr_x)
				qi_hi = DBLE(qi_x)
				qs_hi = DBLE(qs_x)
				qg_hi = DBLE(qg_x)
				nw_hi = DBLE(nw_x)
				ni_hi = DBLE(ni_x)
				ncn_a_hi = DBLE(ncn_a_x)
				nin_a_hi = DBLE(nin_a_x)
				s_numerical = 0.;
				do ikj = 1, NUM_CYCLES
					qv_hi = F_q_hi*dtfine + qv_zero_hi
					T_hi = F_T_hi*dtfine + t_zero_hi

					if( t_x < 170.) stop 5235

					call condensation_and_sublimation(qv_hi, T_hi, qv_zero_hi, t_zero_hi, &
						p_zero_hi, F_q_hi, F_T_hi,  &
						ql_hi, nw_hi, ncn_a_hi, qi_hi, ni_hi, nin_a_hi, qr_hi, qs_hi, qg_hi, rvapor, &
                        		        a_i, b_i, c_i, d_i, a_w, b_w, c_w, a_r, b_r, rhowater, a_s, b_s, rhosnow, &
						a_g, b_g, rhograupel, LAMBDA_FAC_w,  LAMBDA_FAC_i, DBAR_FAC_w, DBAR_FAC_i, &
						N0_FAC_w, N0_FAC_i, p_w, p_i, GRAUPEL_FAC, RAIN_FAC, SNOW_FAC, CP, &
						PIE, EPS, XLS, XLV, Rair, ESIT, &
						DESIT, ESWT, DESWT, temp_K, dtfine, super_WRF_fine,  s_numerical, &
						qv_zero_num, t_zero_num, N_0S, N_0G, N_0R) 

					qv_zero_hi = qv_hi
					t_zero_hi = T_hi
					p_zero_hi = p_zero_hi + F_p_hi*dtfine
				enddo
				super_WRF = super_WRF_fine
				s_FD_test = s_numerical
				s_FD_diff = super_WRF_fine - s_numerical

				qv_x = REAL(qv_hi)
				t_x = REAL(T_hi)
				ql_x = REAL(ql_hi)
				qr_x = REAL(qr_hi)
				qi_x = REAL(qi_hi)
				qs_x = REAL(qs_hi)
				qg_x = REAL(qg_hi) 
				nw_x = REAL(nw_hi)
				ni_x = REAL(ni_hi) 
				ncn_a_x = REAL(ncn_a_hi)
				nin_a_x = REAL(nin_a_hi)

			

                      		if(qi_x > 1. .or. ql_x > 1.) then 
                              		stop 138
                        	endif
			endif

END SUBROUTINE condense_and_sublime_driver


integer FUNCTION NUMBER_SUBCYCLES(N_w, N_i, F_t)
implicit none
REAL, INTENT(IN):: N_w, N_i, F_t
REAL :: w
INTEGER :: NUM_CYCLES_liq, NUM_CYCLES_ice
INTRINSIC  ABS

	w = -F_T*1004/9.8

	NUM_CYCLES_liq = 2
	if(N_w > 30.e6 .and. N_w <= 300.e6 ) then
 		if(ABS(w) >= 5.) then
			NUM_CYCLES_liq = 5
		endif
	endif
	if(N_w > 300.e6 .and. N_w <= 3000.e6 ) then
 		if(ABS(w) >= 10.) then
			NUM_CYCLES_liq = 5
		endif
 		if(ABS(w) >= 20.) then
			NUM_CYCLES_liq = 10
		endif
	endif
	if(N_w > 3000.e6 ) then
 		NUM_CYCLES_liq = 5
		if(ABS(w) >= 1.) then
			NUM_CYCLES_liq = 10
		endif
 		if(ABS(w) >= 10.) then
			NUM_CYCLES_liq = 20
		endif
	endif




	NUM_CYCLES_ice = 2
	if(N_i > 3.e6 .and.  N_i <= 30.e6) then
		if(ABS(w) >= 20) then
 			NUM_CYCLES_ice = 5
		endif
	endif

	if(N_i > 30.e6 .and. N_i <= 300.e6) then
		if(ABS(w) >= 1.) then
			NUM_CYCLES_ice = 5
		endif
		if(ABS(w) >= 20.) then
			NUM_CYCLES_ice = 10
		endif
	endif
	if(N_i > 300.e6 .and. N_i <= 3000.e6) then
		if(ABS(w) >= 1.) then
			NUM_CYCLES_ice = 10
		endif
		if(ABS(w) >= 10.) then
			NUM_CYCLES_ice = 20
		endif
		if(ABS(w) >= 20.) then
			NUM_CYCLES_ice = 40
		endif
        endif

	if(N_i > 3000.e6) then
		if(ABS(w) >= 0.1) then
			NUM_CYCLES_ice = 10
		endif
		if(ABS(w) >= 1.) then
			NUM_CYCLES_ice = 20
		endif
		if(ABS(w) >= 10.) then
			NUM_CYCLES_ice = 40
		endif
		if(ABS(w) >= 20.) then
			NUM_CYCLES_ice = 80
		endif
	endif

	NUMBER_SUBCYCLES = MAX(NUM_CYCLES_liq, NUM_CYCLES_ice)
       return
       END FUNCTION NUMBER_SUBCYCLES


subroutine fillz(km, q_lo, rho, dz8w)
implicit none

! !INPUT PARAMETERS:
   integer,  intent(in):: km                ! No. of levels

   real , intent(in)::  rho(km), dz8w(km)  
! !INPUT/OUTPUT PARAMETERS:
   real , intent(inout) :: q_lo(km)   ! tracer mixing ratio

!
! !LOCAL VARIABLES:
   integer i, k
   double precision ::  q(km), mass_air(km), qup, qdn,qly, dup, ddn, sum_mass1, sum_mass2, &
	totalmass_available
    intrinsic DMIN1, DBLE, DABS

q=DBLE(q_lo)

mass_air(:) = DBLE(dz8w(:)) * DBLE(rho(:)) 

    do k=1,km
	if(DABS(q(k)) < 1.e-20) q(k) = 0.
    enddo

	sum_mass1 = 0.
    do k=1,km
	sum_mass1 = sum_mass1 + q(k)*mass_air(k)
    enddo

if(sum_mass1 == 0.) return

! Top layer
          if( q(km) < 0.) then
             q(km-1) = q(km-1) + q(km)*mass_air(km)/mass_air(km-1)
             q(km) = 0.
          endif

! Interior
      do k=km-1,2,-1
         if( q(k) < 0. ) then
             qdn =  q(k-1)*mass_air(k-1)
! Borrow from above
             qup =  q(k+1)*mass_air(k+1)
	     totalmass_available = qup + qdn
	     if(totalmass_available <= 0.) continue

             qly = -q(k)*mass_air(k)*2.
		qly = DMIN1(qly, totalmass_available) 
             dup =  DMIN1(0.5*qly, qup )        !borrow no more than 50% from top
             q(k+1) = q(k+1) - dup/mass_air(k+1)
	     q(k) = q(k) +  dup/mass_air(k)

! Borrow from below: q(k) is still negative at this stage
              ddn =  DMIN1(qly-dup, qdn )      
             q(k-1) = q(k-1) - ddn/mass_air(k-1)
             q(k) = q(k) + ddn/mass_air(k)
          endif
       enddo
 
! Bottom layer
      k = 1
         if( q(k) < 0.) then
! Borrow from above
            qup =  q(k+1)*mass_air(k+1)
 	     if(qup > 0.) then 

            qly = -q(k)*mass_air(k)*2.
             dup =  DMIN1( qly, qup)
             q(k+1) = q(k+1) - dup/mass_air(k+1) 
             q(k) =  q(k) + dup/mass_air(k) 
          endif
          endif
 
	sum_mass2 = 0.
    do k=1,km
	sum_mass2 = sum_mass2 + q(k)*mass_air(k)
    enddo
    if(sum_mass1 > 0.) then
    if(DABS(sum_mass1 - sum_mass2)/sum_mass1 > 1.e-10) then
		print *, 'WARNING(fillz)', sum_mass1, sum_mass2
    endif
    endif
q_lo=REAL(q)
end subroutine fillz





SUBROUTINE condensation_and_sublimation(qvz_vap, tz_vap, qv_zero_vap, tz_zero_vap, p_zero_vap, F_q_vap, F_T_vap, &
	qlz_vap, nwz_vap, ncn_a_vap, qiz_vap, niz_vap, nin_a_vap, qrz_vap, qsz_vap, qgz_vap, rvapor, a_i, b_i, c_i, &
	d_i, a_w, b_w, c_w, a_r, b_r, RHOW, a_s, b_s, RHOS, a_g, b_g, RHOG, &
	LAMBDA_FAC_w,  LAMBDA_FAC_i, DBAR_FAC_w, DBAR_FAC_i, N0_FAC_w, N0_FAC_i, p_w, p_i, &
	GRAUPEL_FAC, RAIN_FAC, SNOW_FAC, CP, PIE, EPS, XLS, XLV, Rair, ESIT, DESIT, &
	ESWT, DESWT, TEMP_K, dt_vap, &
	super_final_WRF,  s_numerical, qv_zero_numerical, tz_zero_numerical, N_0S, N_0G, N_0R)
implicit none
REAL:: DCRIT_EVAP_UM_LIQ = 10., DCRIT_EVAP_UM_ICE = 40., FUKUTA_FAC
INTEGER :: THRESHOLD_LAMBDA, VENTILATION_ON, TESTING, KOGAN_PINTY_PARTIAL_EVAP = 1 
PARAMETER(THRESHOLD_LAMBDA = 1, VENTILATION_ON = 1, FUKUTA_FAC = 1., TESTING = 0)

DOUBLE PRECISION, INTENT(IN):: F_q_vap, F_T_vap, tz_zero_vap, qv_zero_vap, p_zero_vap, &
	 dt_vap

DOUBLE PRECISION, INTENT(IN)::  ESIT(150001), DESIT(150001), ESWT(150001), DESWT(150001), temp_K(150001)

REAL, INTENT(IN)::         XLS, XLV,  a_r, b_r, RHOW, a_s, b_s, RHOS, a_g, b_g, RHOG, rvapor, a_i, &
	b_i, c_i, d_i, a_w, b_w, c_w, LAMBDA_FAC_w,  LAMBDA_FAC_i, DBAR_FAC_w, DBAR_FAC_i, CP, PIE, &
	EPS, Rair, p_w, p_i, N0_FAC_w, N0_FAC_i, GRAUPEL_FAC, RAIN_FAC, SNOW_FAC, N_0S, N_0G, N_0R

REAL, INTENT(INOUT):: super_final_WRF, s_numerical

 
DOUBLE PRECISION, INTENT(INOUT):: qv_zero_numerical, tz_zero_numerical, tz_vap, qvz_vap, qlz_vap, &
	nwz_vap, qiz_vap, niz_vap, ncn_a_vap, &
       nin_a_vap, qrz_vap, qsz_vap, qgz_vap




DOUBLE PRECISION :: z_vap(2), x_1(2), x_2(2), X_vap(2,2), X_inverse(2,2), M_vap(2,2), h_vap(2), V_vap(2), &
        R_vap(2), norm, QSW_zero, QSI_zero, G_q, G_T, A_q(2), C_q, D_q(2), E_q(2), A_T(2), C_T, &
	D_T(2), E_T(2), Gamw, Gami, r_ice, r_liq, XLVOCP, XLSOCP, tau(2) 
DOUBLE PRECISION :: qvz_vap_save,  tz_vap_save, delta_qi, rwoqswz, rioqsiz, rgoqs2_w, rgoqs2_i, &
	delta_qv, delta_ql, delta_T, ex(2), delta_T_condense, delta_qv_diffusional_growth, qv_dash, T_dash, lambda(2), &	
	dqlovdqi, delta_qliq_approx, delta_qice_approx, exdum(2), &
        super, QSW, delta_qs, delta_qg, &
	delta_qr, r_w, r_r, r_i, r_s, r_g, delta_qliq, delta_qice, ESW, P_SAT, ESI, QSI, ESW_zero, ESI_zero

REAL :: lambda_i, n_0i, AB_ice, k_a, D_v, dbar, lambda_s, lambda_g,  lambda_r
REAL :: lambda_w, n_0w, AB_liq
REAL :: VENT_w, VENT_i, VENT_r, VENT_s, VENT_g, f_vent, delta_ni, delta_nw, rhoz_vap
DOUBLE PRECISION :: T_theory_bar, qv_theory_bar, T_theory(2), qv_theory(2), dts
INTEGER :: evap_flag, j, num_int, kk
! DOUBLE PRECISION :: esatw, esati, deswdt, desidt, ggesw
! REAL :: vent_x, vent_snow, vent_sphere
INTRINSIC  ABS, DEXP, DABS

if((qlz_vap <= 0. .and. qiz_vap <= 0. .and. qrz_vap <= 0. .and. qsz_vap <= 0. .and. qgz_vap <= 0.) .or. qvz_vap <= 0.) return 

if(qlz_vap < 0.) stop 1376
if(qiz_vap < 0.) stop 1375
if(qsz_vap < 0.) stop 13767
if(qgz_vap < 0.) stop 13753
if(qrz_vap < 0.) stop 137

P_SAT = p_zero_vap/100.
rhoz_vap = p_zero_vap/(Rair*tz_zero_vap)
! UWM change -dschanen put numbers > 5 digits in quotes for Sun Fortran compat.
if(RHOW < 900) stop "165897"
if(RHOW > 1100) stop "165897"
if(RHOS > 920) stop 16587
if(RHOS < 10) stop 16587
if(RHOG > 920) stop 16587
if(RHOG < 10) stop 16587


XLVOCP = DBLE(XLV/CP)
XLSOCP = DBLE(XLS/CP)
k_a = (5.69 + 0.017 * (tz_zero_vap - 273.15)) * 1.0e-3 * 4.187
D_v = (tz_zero_vap / 273.15)**1.94
D_v = 0.211 * D_v * (101325. / p_zero_vap) / 1.e4


ESW = ESATW_int(tz_zero_vap, ESWT, TEMP_K)
ESW_zero = 100.*ESW
QSW_zero = DBLE(EPS)*ESW/(P_SAT - ESW)
ESI = ESATI_int(tz_zero_vap, ESIT, TEMP_K)
ESI_zero = 100.*ESI
QSI_zero = DBLE(EPS)*ESI/(P_SAT - ESI)
Gamw = DBLE(EPS)*P_SAT*DESWDT(tz_zero_vap, DESWT, temp_K)/((P_SAT - ESW)**2.)
Gami = DBLE(EPS)*P_SAT*DESIDT(tz_zero_vap, DESIT, temp_K)/((P_SAT - ESI)**2.)


qvz_vap_save = qvz_vap
tz_vap_save = tz_vap

! AB_ice = XLS*XLS/(k_a*rvapor*tz_zero_vap*tz_zero_vap) + 1./(rhoz_vap*QSI_zero*D_v)
AB_ice = (XLS/(rvapor*tz_zero_vap) -1.)*(XLS/(k_a*tz_zero_vap)) + (rvapor*tz_zero_vap/(ESI_zero*D_v))

if(qiz_vap > 0.) then
	lambda_i = LAMBDA_FAC_i*c_i*niz_vap/qiz_vap
	lambda_i = lambda_i**(1./d_i) ! CORRECT
!	n_0i = (niz_vap*rhoz_vap)*(lambda_i**(1.+p_i))*N0_FAC_i ! CORRECT
	n_0i = (niz_vap*rhoz_vap)*N0_FAC_i ! CORRECT
	if(VENTILATION_ON == 1 .and. qiz_vap > 1.e-20) then
        	VENT_i = VENT_x(lambda_i, p_i, a_i, b_i, REAL(tz_zero_vap), rhoz_vap, D_v, DBAR_FAC_i/N0_FAC_i, 1)
		f_vent = VENT_i*lambda_i/(DBAR_FAC_i/N0_FAC_i)
		if(f_vent > 5.) stop 1398
	else
        	VENT_i = DBAR_FAC_i/(N0_FAC_i*lambda_i)
	endif

	r_i = FUKUTA_FAC*DBLE(2.*PIE* n_0i * VENT_i/(AB_ice*rhoz_vap))
else
	r_i = 0.
endif


if(qsz_vap > 0.) then
	if(RHOS < 2.) stop "165986"
	lambda_s = PIE*RHOS*n_0s/(rhoz_vap*qsz_vap)
	lambda_s = lambda_s**(0.25) ! CORRECT
	if(VENTILATION_ON == 1 .and. qsz_vap > 1.e-20) then
        	VENT_s = VENT_snow(lambda_s, a_s, b_s, REAL(tz_zero_vap), rhoz_vap, SNOW_FAC, D_v)
		f_vent = VENT_s*(lambda_s**2.)
		if(f_vent > 50.) stop 1398
	else
        	VENT_s = 1./(lambda_s**2.)
	endif

	r_s = DBLE(2.*PIE* n_0s * VENT_s/(AB_ice*rhoz_vap))
else
	r_s = 0.
endif



if(qgz_vap > 0.) then
	if(RHOG < 2.) stop "165986"
	lambda_g = PIE*RHOG*n_0g/(rhoz_vap*qgz_vap)
	lambda_g = lambda_g**(0.25) ! CORRECT
	if(VENTILATION_ON == 1 .and. qgz_vap > 1.e-20) then
        	VENT_g = VENT_sphere(lambda_g, a_g, b_g, REAL(tz_zero_vap), rhoz_vap, GRAUPEL_FAC, D_v)
		f_vent = VENT_g*(lambda_g**2.)
		if(f_vent > 50.) stop 1398
	else
        	VENT_g = 1./(lambda_g**2.)
	endif

	r_g = DBLE(2.*PIE* n_0g * VENT_g/(AB_ice*rhoz_vap))
else
	r_g = 0.
endif



! AB_liq = XLV*XLV/(k_a*rvapor*tz_zero_vap*tz_zero_vap) + 1./(rhoz_vap*QSW_zero*D_v)
AB_liq = (XLV/(rvapor*tz_zero_vap) -1.)*(XLV/(k_a*tz_zero_vap)) + (rvapor*tz_zero_vap/(ESW_zero*D_v))
if(qlz_vap > 0.) then
	lambda_w = LAMBDA_FAC_w*c_w*nwz_vap/qlz_vap
	lambda_w = lambda_w**(1./3.) ! CORRECT
	n_0w = (nwz_vap*rhoz_vap)*N0_FAC_w ! CORRECT
!	if(VENTILATION_ON == 1) then
!        	VENT_w = VENT_x(lambda_w, p_w, a_w, b_w, REAL(tz_zero_vap), rhoz_vap, D_v, DBAR_FAC_w/N0_FAC_w, 0)
!		f_vent = VENT_w*lambda_w/( DBAR_FAC_w/N0_FAC_w)
!		if(f_vent > 1.5) stop 13998
!		if(f_vent < 1.) stop 13980
!	else
       		VENT_w = DBAR_FAC_w/(N0_FAC_w*lambda_w)
!	endif
	
	r_w = DBLE(2.*PIE* n_0w* VENT_w /(AB_liq*rhoz_vap)) ! neglect ventilation factor for the time being (see Rogers and Yau 1991)
else
	r_w = 0.
endif

if(qrz_vap > 0.) then
	if(RHOW < 2.) stop "165986"
	lambda_r = PIE*RHOW*n_0r/(rhoz_vap*qrz_vap)
	lambda_r = lambda_r**0.25 ! CORRECT

	if(VENTILATION_ON == 1 .and. qrz_vap > 1.e-20) then
        	VENT_r = VENT_sphere(lambda_r, a_r, b_r, REAL(tz_zero_vap), rhoz_vap, RAIN_FAC, D_v)
		f_vent = VENT_r*(lambda_r**2.)
		if(f_vent > 50.) stop 1398
	else
       		VENT_r = 1./(lambda_r**2.)
	endif
	
	r_r = DBLE(2.*PIE* n_0r* VENT_r /(AB_liq*rhoz_vap)) ! neglect ventilation factor for the time being (see Rogers and Yau 1991)
else
	r_r = 0.
endif

r_liq = r_w + r_r
r_ice = r_i + r_s + r_g

if(r_ice <= 0. .and. r_liq <= 0.) return

rwoqswz = r_liq/QSW_zero;  rioqsiz = r_ice/QSI_zero;  
rgoqs2_w = r_liq*Gamw/(QSW_zero**2.);  rgoqs2_i = r_ice*Gami/(QSI_zero**2.); 

M_vap(1,1) = -(rwoqswz+rioqsiz)
M_vap(1,2) = DBLE(qv_zero_vap)*(rgoqs2_w+rgoqs2_i)
M_vap(2,1) = XLVOCP*rwoqswz + XLSOCP*rioqsiz
M_vap(2,2) = -DBLE(qv_zero_vap)*(XLVOCP*rgoqs2_w + XLSOCP*rgoqs2_i)

CALL eigen(M_vap(1,1), M_vap(1,2), M_vap(2,1), M_vap(2,2), x_1, lambda(1), x_2, lambda(2))
ex(:) = DEXP(lambda(:) * DBLE(dt_vap))
do j=1,2
	if(ex(j) == 1.) lambda(j) = 0.
	if(THRESHOLD_LAMBDA == 1) then
		if(DABS(lambda(j)) < 1.e-12) lambda(j) = 0.;
	endif
enddo



X_vap(1,1) = x_1(1)
X_vap(2,1) = x_1(2)
X_vap(1,2) = x_2(1)
X_vap(2,2) = x_2(2)

G_q = DBLE(F_q_vap) - (rwoqswz + rioqsiz)*DBLE(qv_zero_vap) + (r_liq + r_ice )
G_T = DBLE(F_T_vap) + (XLVOCP*rwoqswz + XLSOCP*rioqsiz)*DBLE(qv_zero_vap)  - (XLVOCP*r_liq + XLSOCP*r_ice)

h_vap(1) = G_q
h_vap(2) = G_T

norm = 1./(X_vap(1,1)*X_vap(2,2) - X_vap(1,2)*X_vap(2,1))
X_inverse(1,1) = X_vap(2,2)*norm
X_inverse(1,2) = -X_vap(1,2)*norm
X_inverse(2,1) = -X_vap(2,1)*norm
X_inverse(2,2) = X_vap(1,1)*norm

R_vap(1) = X_inverse(1,1)*h_vap(1) + X_inverse(1,2)*h_vap(2)
R_vap(2) = X_inverse(2,1)*h_vap(1) + X_inverse(2,2)*h_vap(2)



do j=1,2
	if(lambda(j) /= 0.) then
		ex(j)= DEXP(lambda(j) * DBLE(dt_vap))
		tau(j) = 1./lambda(j)
		A_q(j)= X_vap(1,j)* R_vap(j)*tau(j)
		A_T(j)= X_vap(2,j)* R_vap(j)*tau(j)
		E_q(j)= A_q(j)*(ex(j)- 1.)
		E_T(j)= A_T(j)*(ex(j)- 1.)
		D_q(j)= A_q(j)*(ex(j)- 1.)*tau(j)- A_q(j)*DBLE(dt_vap)
		D_T(j)= A_T(j)*(ex(j)- 1.)*tau(j)- A_T(j)*DBLE(dt_vap)
	else
		E_q(j)= X_vap(1,j)*R_vap(j)*DBLE(dt_vap)
		E_T(j)= X_vap(2,j)*R_vap(j)*DBLE(dt_vap)
		D_q(j)= X_vap(1,j)*R_vap(j)*DBLE(dt_vap)*DBLE(dt_vap)/2.
		D_T(j)= X_vap(2,j)*R_vap(j)*DBLE(dt_vap)*DBLE(dt_vap)/2.
	endif
enddo

qv_dash = E_q(1) + E_q(2)
T_dash = E_T(1) + E_T(2)

delta_qv_diffusional_growth = qv_zero_vap + qv_dash - qvz_vap_save
delta_T_condense = tz_zero_vap + T_dash - tz_vap_save

if(delta_qv_diffusional_growth == 0. .and. delta_T_condense == 0.) return


!
! REPRESENT CHANGES IN MASS: ALL ICE AND LIQUID 
!
delta_ql = 0.
delta_qr = 0.
delta_qi = 0.
delta_qs = 0.
delta_qg = 0.
delta_qliq = 0.
delta_qice=  0.
delta_T = 0.


evap_flag = 0
if(r_liq > 0. .and. r_ice > 0.) then
	C_q = D_q(1) + D_q(2)
	C_T = D_T(1) + D_T(2)
	delta_qliq_approx = C_q - DBLE(dt_vap)*QSW_zero + DBLE(dt_vap)*DBLE(qv_zero_vap) - C_T*Gamw*DBLE(qv_zero_vap)/QSW_zero
	delta_qliq_approx = rwoqswz * delta_qliq_approx

	delta_qice_approx = C_q - DBLE(dt_vap)*QSI_zero + DBLE(dt_vap)*DBLE(qv_zero_vap) - C_T*Gami*DBLE(qv_zero_vap)/QSI_zero
	delta_qice_approx = rioqsiz * delta_qice_approx

	if(DABS(delta_qice_approx) < DABS(delta_qliq_approx)) then
		Delta_qice = delta_qice_approx
		Delta_qliq = -delta_qv_diffusional_growth - Delta_qice
	else
		Delta_qliq = delta_qliq_approx
		Delta_qice = -delta_qv_diffusional_growth - Delta_qliq
	endif

else
	if(r_liq > 0.) then
		delta_qliq = -delta_qv_diffusional_growth
	endif
	if(r_ice > 0.) then
		delta_qice = -delta_qv_diffusional_growth
	endif
endif

if( delta_qliq > qvz_vap/2. .and. delta_qliq > 0.) delta_qliq = qvz_vap/2.
if( delta_qice > qvz_vap/2. .and. delta_qice > 0.) delta_qice = qvz_vap/2.

if(r_ice > 0.) then
	delta_qi = delta_qice*r_i/r_ice
	delta_qs = delta_qice*r_s/r_ice
	delta_qg = delta_qice*r_g/r_ice
endif

if(r_liq > 0.) then
	delta_ql = delta_qliq*r_w/r_liq
	delta_qr = delta_qliq*r_r/r_liq
endif

!
! REPRESENT CHANGES IN MASS : CLOUD-LIQUID
!

if(qlz_vap > 0. .and. nwz_vap > 0.) then
	dbar = DBAR_FAC_w/lambda_w	
else
	dbar = 0.
endif

if((DBLE(qlz_vap) + delta_ql <= 1.e-10 .or. dbar < 1.e-6) .and. delta_ql < 0.) then

! TOTAL EVAPORATION
	evap_flag = 1
	delta_ql = -DBLE(qlz_vap)
	qlz_vap = 0.
	nwz_vap = 0.
!	call code_storage()
	
else
! PARTIAL EVAPORATION

	if(dbar < 1.e-6*DCRIT_EVAP_UM_LIQ .and. delta_ql < 0. .and. qlz_vap > 0. .and. nwz_vap > 0.) then
		delta_nw = delta_ql*nwz_vap/qlz_vap
		if( nwz_vap + delta_nw < 0.) then
			 delta_nw = -nwz_vap
		endif
		nwz_vap =  nwz_vap + delta_nw

! every droplet contains a CCN
		ncn_a_vap = ncn_a_vap + delta_nw
	endif
	qlz_vap = qlz_vap + delta_ql
endif

!
! REPRESENT CHANGES IN MASS : RAIN
!

if(qrz_vap > 0.) then
	dbar = 1./lambda_r	
else
	dbar = 0.
endif
if((DBLE(qrz_vap) + delta_qr <= 1.e-10 .or. dbar < 1.e-6) .and. delta_qr < 0.) then
	evap_flag = 1
	delta_qr = -DBLE(qrz_vap)
	qrz_vap = 0.
else
	qrz_vap = qrz_vap + delta_qr
endif

!
! REPRESENT CHANGES IN MASS : CLOUD-ICE
!


if(qiz_vap > 0. .and. niz_vap > 0.) then
	dbar = DBAR_FAC_i/lambda_i	
else
	dbar = 0.
endif

if( (DBLE(qiz_vap) + delta_qi <= 1.e-10 .or. dbar < 1.e-6) .and. delta_qi < 0.) then

! TOTAL SUBLIMATION

	evap_flag = 2
	delta_qi = -DBLE(qiz_vap)
	niz_vap = 0.
	qiz_vap = 0.
else

! PARTIAL SUBLIMATION

	if(dbar < 1.e-6*DCRIT_EVAP_UM_ICE .and. delta_qi < 0. .and. qiz_vap > 0. .and. niz_vap > 0. ) then
		delta_ni = delta_qi*niz_vap/qiz_vap
		if( niz_vap + delta_ni < 0.) then
			 delta_ni = -niz_vap
		endif
		niz_vap =  niz_vap + delta_ni
! not every crystal contains an IN - in fact most dont => do nothing here to nin_a_vap 
	endif
	qiz_vap = qiz_vap + delta_qi

endif
!
! REPRESENT CHANGES IN MASS : SNOW
!

if(qsz_vap > 0.) then
	dbar = 1./lambda_s	
else
	dbar = 0.
endif
if( (DBLE(qsz_vap) + delta_qs <= 1.e-10 .or. dbar < 1.e-6) .and. delta_qs < 0.) then
	evap_flag = 2
	delta_qs = -DBLE(qsz_vap)
	qsz_vap = 0.
else
	qsz_vap = qsz_vap + delta_qs
endif

!
! REPRESENT CHANGES IN MASS : GRAUPEL
!

if(qgz_vap > 0.) then
	dbar = 1./lambda_g	
else
	dbar = 0.
endif
if( (DBLE(qgz_vap) + delta_qg <= 1.e-10 .or. dbar < 1.e-6) .and. delta_qg < 0.) then
	evap_flag = 2
	delta_qg = -DBLE(qgz_vap)
	qgz_vap = 0.
else
	qgz_vap = qgz_vap + delta_qg

endif

!
! REPRESENT CHANGES IN VAPOUR AND TEMPERATURE
!


delta_T = (delta_ql+delta_qr)*XLVOCP +  (delta_qi+delta_qs+delta_qg)*XLSOCP
delta_qv = -( delta_ql+delta_qr + delta_qi+delta_qs+delta_qg )
qvz_vap = qvz_vap + delta_qv
tz_vap = tz_vap + delta_T




if(qvz_vap < 0.) stop 14368
if(tz_vap < 150.) stop 14369
if(tz_vap > 350.) stop 14370

if(evap_flag == 0) then
	if(ABS(delta_qv - delta_qv_diffusional_growth) > 1.E-10) then
!		stop 13655
	endif
	
!	if(ABS(delta_T)/tz_vap > 1.e-8) then
!		if(ABS((delta_T - delta_T_condense)/delta_T) > 0.01) then
!			stop 13656
!		endif
!	endif

endif

ESW = GGESW(tz_vap)
QSW = DBLE(EPS)*ESW/(P_SAT - ESW)
super_final_WRF = 100.*REAL(qvz_vap/QSW - 1.)


if(TESTING == 1) then

 
	num_int = 1000
	dts = dt_vap/DBLE(num_int)		

	T_theory(1) = tz_zero_numerical
	qv_theory(1) = qv_zero_numerical

	do kk = 1, num_int	

! First step

		ESW =  GGESW(T_theory(1))
		QSW = DBLE(EPS)*ESW/(P_SAT - ESW)
		qv_theory_bar = qv_theory(1) + dts*(F_q_vap - r_liq*(qv_theory(1)/QSW - 1.))
		T_theory_bar = T_theory(1) + dts*(F_T_vap + r_liq*XLVOCP*(qv_theory(1)/QSW - 1.))
		qv_theory_bar = (qv_theory_bar + qv_theory(1) )*0.5
		T_theory_bar = (T_theory_bar + T_theory(1) )*0.5

! Second step

		ESW =  GGESW(T_theory_bar)
		QSW = DBLE(EPS)*ESW/(P_SAT - ESW)
		qv_theory(2) = qv_theory(1) + dts*(F_q_vap - r_liq*(qv_theory_bar/QSW - 1.))
		T_theory(2) = T_theory(1) + dts*(F_T_vap + r_liq*XLVOCP*(qv_theory_bar/QSW - 1.))

		T_theory(1) = T_theory(2)
		qv_theory(1) = qv_theory(2)

	enddo

	tz_zero_numerical = T_theory(2) 
	qv_zero_numerical = qv_theory(2)


	ESW = GGESW(T_theory(2))
	QSW = DBLE(EPS)*ESW/(P_SAT - ESW)
	s_numerical = 100.*REAL(qv_theory(2)/QSW - 1.) 

endif



END SUBROUTINE condensation_and_sublimation

















SUBROUTINE zero_trace_cloud_in_subsat_env(qlz_vap, nwz_vap, qiz_vap, niz_vap,  qvz_vap, tz_vap, prez_vap, &
	CP, EPS, XLS, XLV, ESIT, DESIT, ESWT, DESWT, temp_K)
implicit none
REAL, INTENT(IN) :: CP, EPS, XLS, XLV
REAL, INTENT(INOUT) :: qlz_vap, nwz_vap, qiz_vap, niz_vap,  qvz_vap, tz_vap, prez_vap
REAL :: XLVOCP, XLSOCP, P_SAT, ESW, QSW, ESI, QSI
DOUBLE PRECISION, INTENT(IN) :: ESIT(150001), DESIT(150001), ESWT(150001), DESWT(150001), temp_K(150001)

XLVOCP = XLV/CP
XLSOCP = XLS/CP
P_SAT = prez_vap/100.
ESW = REAL(ESATW_int(DBLE(tz_vap), ESWT, temp_K))
QSW = DBLE(EPS*ESW/(P_SAT - ESW))
if(qvz_vap/QSW < 1. .and. qlz_vap > 0. .and. qlz_vap < 1.e-10)then
	tz_vap = tz_vap - XLVOCP * qlz_vap 
	qvz_vap = qvz_vap + qlz_vap
	qlz_vap = 0.
	nwz_vap = 0.
endif

ESI = REAL(ESATI_int(DBLE(tz_vap), ESIT, temp_K))
QSI = EPS*ESI/(P_SAT-ESI)
if(qvz_vap/QSI < 1. .and. qiz_vap > 0. .and. qiz_vap < 1.e-10)then
	tz_vap = tz_vap - XLSOCP * qiz_vap 
	qvz_vap = qvz_vap + qiz_vap
	qiz_vap = 0.
	niz_vap = 0.
endif

END SUBROUTINE zero_trace_cloud_in_subsat_env


REAL FUNCTION VENT_x( lambda_x, p_x, a_x, b_x, temp_x, rho_x, D_vx, GAM2PLUSp, ice_flag)
implicit none
REAL, INTENT(IN) :: a_x, b_x, lambda_x, temp_x, rho_x, D_vx, p_x, GAM2PLUSp
INTEGER, INTENT(IN) :: ice_flag
REAL :: te, eta, rhofac, Sc, D_xf, chi_cut, psi_x
REAL :: arg1, arg2, arg3, vent_min

	 vent_min = GAM2PLUSp/lambda_x
         te = temp_x - 273.15
         eta = (1.718 + 0.0049 * te - 1.2e-5 * te * te) * 1.e-5/rho_x;
         rhofac = ((1.2/rho_x)**0.5)
         Sc = eta/D_vx
	 psi_x = ((rhofac*a_x)**0.5)*(((eta**0.5) * D_vx)**(-1./3.))

         if(ice_flag == 0) then

		D_xf = (1.4/(Sc**(1./3.)))**2.
		D_xf = (eta * D_xf/(a_x*rhofac ))**(1./(b_x+1.))
	
 		arg1 = vent_min*(0.78 + (1.-0.78)*gammp(2.+p_x,lambda_x*D_xf))
	 	arg2 = (0.108*psi_x*psi_x/(lambda_x**(2.+b_x)))*gammp(3.+p_x+b_x,lambda_x*D_xf)
	 	arg3 = (0.308*psi_x/(lambda_x**(1.5+0.5*b_x))) &
			*(exp(gammln(2.5+p_x+0.5*b_x)) - gammp(2.5+p_x+0.5*b_x,lambda_x*D_xf))

        else

		D_xf = (1./(Sc**(1./3.)))**2.
		D_xf = (eta * D_xf/(a_x*rhofac ))**(1./(b_x+1.))
	
		arg1 = vent_min*(0.86 + (1.-0.86)*gammp(2.+p_x,lambda_x*D_xf))
		arg2 = (0.14*psi_x*psi_x/(lambda_x**(2.+b_x)))*gammp(3.+b_x+p_x,lambda_x*D_xf)
		arg3 = (0.28*psi_x/(lambda_x**(1.5+0.5*b_x))) &
			*(exp(gammln(2.5+p_x+0.5*b_x)) - gammp(2.5+p_x+0.5*b_x,lambda_x*D_xf))

        endif
	VENT_x = arg1 + arg2 + arg3
	if(VENT_x < vent_min) VENT_x = vent_min
RETURN
END FUNCTION VENT_x

REAL FUNCTION VENT_sphere( lambda_x, a_x, b_x, temp_x, rho_x, PRECIP_FAC, D_vx)
implicit none
REAL, INTENT(IN) :: a_x, b_x, lambda_x, temp_x, rho_x, D_vx, PRECIP_FAC
REAL :: te, eta, rhofac, Sc
REAL :: arg1, arg2, vent_min
! PK97, page 541

  	vent_min = (1./(lambda_x**2.))
       te = temp_x - 273.15
         eta = (1.718 + 0.0049 * te - 1.2e-5 * te * te) * 1.e-5/rho_x;
         rhofac = ((1.2/rho_x)**0.5)
         Sc = eta/D_vx

	
 	arg1 = vent_min*0.78
	arg2 = 0.308*(Sc**(1./3.))*((rhofac*a_x/eta)**0.5) * PRECIP_FAC/(lambda_x**(2.5+0.5*b_x))

  	VENT_sphere = arg1 + arg2
	if(VENT_sphere < vent_min) VENT_sphere = vent_min
RETURN
END FUNCTION VENT_sphere

REAL FUNCTION VENT_snow( lambda_x, a_x, b_x, temp_x, rho_x, SNOW_FAC, D_vx)
implicit none
REAL, INTENT(IN) :: a_x, b_x, lambda_x, temp_x, rho_x, D_vx, SNOW_FAC
REAL :: te, eta, rhofac, Sc
REAL :: arg1, arg2, vent_min
! PK97, page 552

	if(ABS(exp(gammln(2.)) - 1.) > 1.e-8) stop 65432
 
 	 vent_min = (1./(lambda_x**2.))
         te = temp_x - 273.15
         eta = (1.718 + 0.0049 * te - 1.2e-5 * te * te) * 1.e-5/rho_x;
         rhofac = ((1.2/rho_x)**0.5)
         Sc = eta/D_vx

	
 	arg1 = vent_min*0.86
	arg2 = 0.28*(Sc**(1./3.))*((rhofac*a_x/eta)**0.5) * SNOW_FAC/(lambda_x**(2.5+0.5*b_x))

  	VENT_snow = arg1 + arg2
	if(VENT_snow < vent_min) VENT_snow = vent_min
RETURN
END FUNCTION VENT_snow

SUBROUTINE CCN_activity_spectrum(tz, zz, z_pbl_top, s_pc_cut_off, kaero, C_by_rhoz, &
	rhoz, PIE, N_env_perm3, CCN_conc_env,  rccn_cm, kts, kte)
implicit none
INTEGER :: AMMONIUM_SULPHATE_IN_PBL, NKR
REAL :: COL
PARAMETER(COL = 0.23105, NKR=33, AMMONIUM_SULPHATE_IN_PBL = 1)
REAL, INTENT(IN) ::  s_pc_cut_off, kaero, PIE

INTEGER, INTENT(IN) :: kts, kte
REAL, DIMENSION (NKR), INTENT(OUT) :: rccn_cm 
REAL, DIMENSION (NKR, kts:kte), INTENT(OUT) ::  N_env_perm3
REAL, DIMENSION (kts:kte), INTENT(OUT) ::  CCN_conc_env
REAL, DIMENSION (kts:kte), INTENT(IN) ::  rhoz, zz, tz
REAL, DIMENSION (kts:kte), INTENT(INOUT) ::  C_by_rhoz
REAL, INTENT(IN) :: z_pbl_top

REAL, DIMENSION(NKR) :: roccn, fccnr_percm3
DOUBLE PRECISION ::  N_env_perm3_hi(NKR), sum_ncn, del_ncn
REAL :: CCN_conc_total, s_kr, ro_solute, x0ccn, &
aa, a,b, a1,a2, x0drop, x0, r0, rccnkr_cm, caero

INTEGER:: kr, k, ammonium_sulphate
 
! Bill says:  for tropical Pacific, winds are light and DMS production will prevail for PBL from plankton etc

  DO k = kts, kte
	if(zz(k) > z_pbl_top) then
		ammonium_sulphate = 1
	else
		if(AMMONIUM_SULPHATE_IN_PBL == 1) then
			ammonium_sulphate = 1
		else
			ammonium_sulphate = 0
		endif	
	endif

	if(ammonium_sulphate == 0)then 
        	ro_solute=2.16;
	else
        	ro_solute=1.77;
	endif
	x0drop= 3.3510000d-11 ;
	aa = REAL(NKR - 1);
	x0ccn =x0drop/(2.**aa);
	a=3.3d-05/tz(k); 

!	b=2.*4.3/(22.9+35.5)
!

	if(ammonium_sulphate == 0) then
		b=2.*4.3/(58.44)
	else
		b=3.*4.3/(132.14)
	endif
	b=b*(4./3.)*PIE*ro_solute;
	a1=2.*((a/3.)**1.5)/(b**0.5); 
	a2=a1*100.; 

	do kr = 1,NKR
		roccn(kr)= ro_solute;
        	x0=x0ccn*(2.** REAL(kr-1));
        	r0=(3.*x0/4./PIE/roccn(kr))** (1./3.);
         	rccn_cm(kr)=r0;
        enddo


  	if(C_by_rhoz(k) < 1.) then
		print *, 'WARNING :: THRESHOLDING C_by_rho'
		C_by_rhoz(k) = 1.
	endif
  	if(rhoz(k) < 0.) stop "6543653"
	caero = C_by_rhoz(k) *  rhoz(k)
 	N_env_perm3(:,k) = 0.
	do kr = 1,NKR
        	rccnkr_cm= rccn_cm(kr);
        	s_kr=a2/(rccnkr_cm** 1.5);
        	fccnr_percm3(kr)=1.5*(caero/1.e6)*kaero*((s_kr)**(kaero));
		if(s_kr > s_pc_cut_off) then
			fccnr_percm3(kr)= 0.
		endif
	enddo
	CCN_conc_env(k) = caero*( s_pc_cut_off**kaero)
	N_env_perm3(:,k) = fccnr_percm3(:)*COL*1.e6
	N_env_perm3_hi = DBLE(N_env_perm3(:,k))


	sum_ncn = 0.
	do kr = 1,NKR
		sum_ncn = sum_ncn + N_env_perm3_hi(kr)
	enddo
	del_ncn = DBLE(CCN_conc_env(k)) - sum_ncn

	do kr = 1,NKR
		if( N_env_perm3_hi(kr) > 0.) then
			if( N_env_perm3_hi(kr)+ del_ncn > 0.) then
				N_env_perm3_hi(kr) = N_env_perm3_hi(kr)+ del_ncn
				exit
			else
				del_ncn = del_ncn + N_env_perm3_hi(kr)
				N_env_perm3_hi(kr) = 0.
			endif
		endif
	enddo

	sum_ncn = 0.
	do kr = 1,NKR
		sum_ncn = sum_ncn + N_env_perm3_hi(kr)
	enddo
	if( ABS(sum_ncn - CCN_conc_env(k)) > 100.) then
		stop "136580"
	endif
	N_env_perm3(:,k) = REAL(N_env_perm3_hi)
   END DO

END SUBROUTINE CCN_activity_spectrum

SUBROUTINE eigen(a, b, c, d, x_1, lambda_1, x_2, lambda_2)
DOUBLE PRECISION, INTENT(IN):: a,b,c,d
DOUBLE PRECISION, INTENT(OUT):: x_1(2), lambda_1, x_2(2), lambda_2
DOUBLE PRECISION :: roots

if(a**2. + 4.*b*c - 2.*a*d + d**2. < 0.) stop 14376

if(c == 0.) stop 14398

roots = (a**2. + 4.*b*c - 2.*a*d + d**2.)**(1./2.)

x_1(1) = -(-a + d + roots)/(2.*c)
x_1(2) = 1.
x_2(1) = -(-a + d - roots)/(2.*c)
x_2(2) = 1.

lambda_1 = (a + d - roots)/2.
lambda_2 = (a + d + roots)/2.

if(lambda_1 == lambda_2) stop 14365



if(lambda_1 == 0. .and. lambda_2 == 0.) print *,'WARNING:  both eigenvalues are zero in SUBROUTINE eigen'


END SUBROUTINE eigen


SUBROUTINE droplet_nucleation(supersat, s_vw_saved, s_pc_start, s_pc_cut_off, temperature, pressure, rho_a, &
	 w_speed, qv_nuc, ql_nuc, nw_nuc, ncn_a_nuc, &
	caero, kaero, XLV, CP, EPS, Rair, rvapor, grav, &
	BETA_SUPER, PIE, ESIT, DESIT, ESWT, DESWT, temp_K, c_w, n_w_crit_env)
implicit none
REAL, INTENT(IN) :: s_vw_saved, s_pc_start, s_pc_cut_off, &
	 pressure, rho_a,w_speed, caero, kaero, XLV, n_w_crit_env, &
	CP, EPS, Rair, rvapor, grav, BETA_SUPER, PIE, &
          c_w
REAL, INTENT(INOUT) :: supersat, temperature, qv_nuc, ql_nuc, nw_nuc, ncn_a_nuc

DOUBLE PRECISION, INTENT(IN) ::ESIT(150001), DESIT(150001), ESWT(150001), DESWT(150001), temp_K(150001)
REAL :: caero_pinty, Q_1, Q_2, super_fac, r_mean, s_vw_eqm, s_vw_cb
REAL :: psi_1, psi_2,  G_factor, D_v, k_a, s_vw, D_crit, CVHENC, RVHENC
REAL ::  QSW, ESW, P_SAT, TT, supersaturation, ee_nuc

        caero_pinty = caero*(100.**kaero)
	TT = temperature
	P_SAT= pressure/100.
  	ESW = REAL(ESATW_int(DBLE(TT), ESWT, temp_K))
      	QSW = EPS*ESW/(P_SAT-ESW)
	ESW = 100.*ESW
	ee_nuc = qv_nuc*pressure/(qv_nuc + EPS)
	supersaturation = (ee_nuc/ESW - 1.)
	if(supersaturation < 0.) then
		supersat = 0.
		return
	endif	
	if( supersaturation > 0. &
			 .and. temperature > 273.15 - 100.) then
		if(w_speed > 0.) then

			k_a = (5.69 + 0.017 * (temperature - 273.15)) * 1.0e-3 * 4.187
			D_v = (temperature / 273.15)**1.94
			D_v = 0.211 * D_v * (101325. / pressure) / 1.e4


			if(nw_nuc > 0. .and. ql_nuc > 0.) then
				r_mean = (ql_nuc/nw_nuc)/c_w
				r_mean = r_mean**(1./3.)
				r_mean = r_mean*0.5

			else
				r_mean = 0.
			endif

			if(nw_nuc*rho_a > n_w_crit_env*1.e6 .and. r_mean > 1.e-6) then	
! IN-CLOUD
				s_vw = supersaturation
				
			else	
! OUT-OF-CLOUD
				psi_1 = (grav/(temperature*Rair)) * (EPS*XLV/(cp*temperature) - 1.)
				psi_2 = pressure/(EPS* ESW) + EPS* XLV * XLV/( Rair * temperature*temperature* cp) 
				G_factor = rvapor*temperature/ (ESW * D_v) + (XLV/(k_a*temperature)) *(XLV/(rvapor*temperature)  - 1.)
				G_factor = (1./1000.)*(1./G_factor )
				s_vw_cb = 2.* kaero*caero_pinty*PIE*1000.*psi_2*( G_factor**(3./2.))* BETA_SUPER
			        s_vw_cb = rho_a * ( ( psi_1 * w_speed)**(3./2.))/s_vw_cb
				s_vw_cb = s_vw_cb**(1./(kaero+2.))
				s_vw = s_vw_cb

				if(s_vw < supersaturation) s_vw = supersaturation			
				

			endif

!			print *, SUPERSAT (Pinty) ::, s_vw_cb*100., SUPERSAT (Twomey) ::, twomey_s_cb(w_speed, caero, kaero)
		
		else
			s_vw = supersaturation
			
		endif
		if(s_vw > s_pc_cut_off/100.) then
			s_vw = s_pc_cut_off/100.
		endif
		supersat = s_vw*100.
	else		
		supersat = 0.
	endif

  	if(supersat > s_pc_start) then
		CVHENC = Caero*(supersat**kaero) - ncn_a_nuc * rho_a
		if(CVHENC/rho_a < 0.) CVHENC = 0.
		D_crit = 7.5E-2/(rvapor*temperature*1000.)
		D_crit = 4.*D_crit/(3.*supersat/100.)
		if(D_crit > 10.e-6) then
			print *,'WARNING::  thresholding D_crit = ', D_crit 
			D_crit = 10.e-6
		endif
		if(D_crit < 1.e-6) D_crit = 1.e-6
		RVHENC = CVHENC*(c_w/rho_a) * (D_crit**3.)  ! RVHENC > 0.
		
		if(qv_nuc - RVHENC > 0.) then		
			if(temperature < -45.+ 273.15)then
				print *, 'secondary droplet nucleation at ', temperature - 273.15, 'degC'
			endif
			nw_nuc =  nw_nuc + CVHENC/rho_a
			ncn_a_nuc = ncn_a_nuc + CVHENC/rho_a
			ql_nuc = ql_nuc + RVHENC
			qv_nuc = qv_nuc - RVHENC
			if(qv_nuc < 0.) then
				print *, temperature, supersat, D_crit
				stop 3156
			endif
			temperature = temperature + XLV*RVHENC/cp	
		else
			CVHENC = 0.
			RVHENC = 0.
		endif	
	else
		CVHENC = 0.
		RVHENC = 0.
	endif


END SUBROUTINE droplet_nucleation


SUBROUTINE homogeneous_aerosol_freezing(temperature, pressure, rho_a,  w_prim, &
	qv_prim, qi_prim, ni_prim, ncn_a_prim, &
	 XLF, XLS, CP, PIE, EPS, ESIT, DESIT, ESWT, DESWT, temp_K, &
        N_env_perm3, CCN_conc_env, rccn_cm, num_aero)
implicit none
INTEGER :: NKR
REAL :: COL
PARAMETER(NKR=33, COL = 0.23105)
REAL, INTENT(INOUT) :: num_aero
REAL, INTENT(IN) :: XLF, XLS, CP, EPS, pressure, w_prim, rho_a, PIE, &
          N_env_perm3(NKR), rccn_cm(NKR), CCN_conc_env
DOUBLE PRECISION, INTENT(IN) :: ESIT(150001), &
	 DESIT(150001), ESWT(150001), DESWT(150001), temp_K(150001)
REAL, INTENT(INOUT) :: temperature, qv_prim, qi_prim, ni_prim, ncn_a_prim
REAL, DIMENSION(NKR) :: n_ccn
REAL :: CCN_conc_total
REAL :: DEL_NCN,  sum_N, S_i_crit, dmm, ESW, ESI, P_SAT, TT, QSI, QSW, Q_is, Q_ws, S_i
INTEGER :: kr
REAL :: m_i0, CIHENC, RIHENC, DEL_QV, S_i_crit_next, dmm_next, del_n_ccn 

	if(temperature > 235.) then
		return
	endif
	TT = temperature
	P_SAT= pressure/100.
 
  	ESI = REAL(ESATI_int(DBLE(TT), ESIT, temp_K))
      	QSI = EPS*ESI/(P_SAT-ESI)
	if(qv_prim <  QSI) then
		return
	endif
  	ESW = REAL(ESATW_int(DBLE(TT), ESWT, temp_K))
      	QSW = EPS*ESW/(P_SAT-ESW)
	Q_is = QSI
	Q_ws = QSW

 	CCN_conc_total =  CCN_conc_env - ncn_a_prim*rho_a
	if(CCN_conc_total < 0.) CCN_conc_total = 0.

	n_ccn(:) = 0.
	if(CCN_conc_total > 0.) then
		sum_N = 0.
		do kr = 1, NKR
			DEL_NCN = N_env_perm3(kr)
			if(sum_N + DEL_NCN >= CCN_conc_total) then
				n_ccn(kr) = CCN_conc_total - sum_N
				sum_N = sum_N + n_ccn(kr)
				exit
			endif
			n_ccn(kr) = DEL_NCN
			sum_N = sum_N + n_ccn(kr)
		enddo
!		if( ABS(CCN_conc_total-SUM(n_ccn(:))) > 1.)then
!			print *, sum_N, rho_a, rho_a*SUM(n_ccn(:)),  ABS(sum_N-rho_a*SUM(n_ccn(:)))
!			stop 13567	
!		endif

!		print *, n_ccn
!		print *, sum_N =, sum_N/1.e6, cm-3 
		n_ccn = n_ccn/rho_a	
	endif



	if(qv_prim > Q_is) then
		if(qv_prim < Q_ws) then
			S_i = qv_prim/Q_is 
		else
			S_i = Q_ws/Q_is 
		endif	
		do kr = NKR,1,-1

!going to smaller and smaller sizes, larger and larger S_i_{crit}
			TT = temperature
			P_SAT= pressure/100.
 
  			ESI = REAL(ESATI_int(DBLE(TT), ESIT, temp_K))
      			QSI = EPS*ESI/(P_SAT-ESI)
			if(qv_prim <  QSI) then
				return
			endif
  			ESW = REAL(ESATW_int(DBLE(TT), ESWT, temp_K))
      			QSW = EPS*ESW/(P_SAT-ESW)
			Q_is = QSI
			Q_ws = QSW

			if(qv_prim < Q_ws) then
				S_i = qv_prim/Q_is 
			else
				S_i = Q_ws/Q_is 
			endif	
			if(n_ccn(kr) > 0.) then
				dmm = rccn_cm(kr)*20.;
				S_i_crit =  S_i_homog_koop(dmm, temperature);
				if(kr > 1) then
					dmm_next = rccn_cm(kr-1)*20.;
					S_i_crit_next =  S_i_homog_koop(dmm_next, temperature);
				else
					S_i_crit_next = S_i_crit;
				endif
				if(S_i_crit_next < S_i_crit) stop 13579

				if(S_i > S_i_crit ) then
					m_i0 = 920. * (4./3.)*PIE*( (rccn_cm(kr)/100.) **3.)
					if(S_i > S_i_crit_next) then
						CIHENC =  n_ccn(kr)* rho_a
						n_ccn(kr) = 0.
						if(CIHENC/rho_a < 0.) CIHENC = 0.
					else
						if((S_i - S_i_crit)/(S_i_crit_next - S_i_crit) < 0.) stop 1356
						if((S_i - S_i_crit)/(S_i_crit_next - S_i_crit) > 1.) stop 1357
						if(S_i_crit_next - S_i_crit  == 0.) stop "134576"

						del_n_ccn = n_ccn(kr) * (S_i - S_i_crit)/(S_i_crit_next - S_i_crit)
						CIHENC =  del_n_ccn* rho_a
						n_ccn(kr) = n_ccn(kr) - CIHENC
					endif
					num_aero = num_aero + CIHENC*500.*2000.*2000.

					ni_prim =  ni_prim + CIHENC/rho_a
					ncn_a_prim = ncn_a_prim + CIHENC/rho_a
					RIHENC = CIHENC*m_i0/rho_a
					qi_prim = qi_prim + RIHENC
					qv_prim = qv_prim - RIHENC
					temperature = temperature + XLS*RIHENC/cp	
					print *,'HOMOGENEOUS FREEZING OF AEROSOLS at N = ', CIHENC/1.e6,'cm-3 at T =', &
						temperature - 273.15,' degC and w =', w_prim, ' m/s'


				endif



			endif
		enddo	 			
	endif



END SUBROUTINE homogeneous_aerosol_freezing

REAL FUNCTION S_i_homog_koop(dmmx, temx)
implicit none
REAL :: AMMONIUM_SULPHATE_CORRECTION
! PARAMETER(AMMONIUM_SULPHATE_CORRECTION = 0.3)  ! AIDA cloud chamber results presented at EGU 2005 by A. Mangold
PARAMETER(AMMONIUM_SULPHATE_CORRECTION = 0.)
REAL, INTENT(IN) :: dmmx, temx

	integer:: curve, ix;
	real :: rad_um, S_ix(7), rad_um_ix(7), ans;


	rad_um_ix(1) = 0.05;  
	rad_um_ix(2) = 0.1; 
	rad_um_ix(3) = 0.2; 
	rad_um_ix(4) = 0.5; 
	rad_um_ix(5) = 1.0;
	rad_um_ix(6) = 2.; 
	rad_um_ix(7) = 5.; 
	
	rad_um = 1.e3*dmmx/2.;

	S_ix(1) = 1.71 + (temx - 175.)*(1.485 - 1.71)/(230. - 175.);
	S_ix(2) = 1.70 + (temx - 175.)*(1.48 - 1.70)/(230. - 175.);
	S_ix(3) = 1.69 + (temx - 175.)*(1.475 - 1.69)/(230. - 175.);
	S_ix(4) = 1.675 + (temx - 175.)*(1.465 - 1.675)/(230. - 175.);
	S_ix(5) = 1.67 + (temx - 175.)*(1.46 - 1.67)/(230. - 175.);
	S_ix(6) = 1.66 + (temx - 175.)*(1.455 - 1.66)/(230. - 175.);
	S_ix(7) = 1.65 + (temx - 175.)*(1.445 - 1.65)/(230. - 175.);

	if(rad_um < rad_um_ix(1)) then
		S_i_homog_koop = S_ix(1)
		return;
	endif
	if(rad_um > rad_um_ix(7)) then
		S_i_homog_koop = S_ix(7)
		return;
	endif
	do ix = 1, 7
		if (rad_um_ix(ix) > rad_um) then
			exit
		endif
	enddo
	if(ix == 1) stop "135675"
	
	S_i_homog_koop = S_ix(ix-1) + (log10(rad_um) - log10(rad_um_ix(ix-1))) * (S_ix(ix) - S_ix(ix-1)) / &
		(log10(rad_um_ix(ix)) - log10(rad_um_ix(ix-1)));

	if(temx < 230.) then
	S_i_homog_koop = S_i_homog_koop - AMMONIUM_SULPHATE_CORRECTION
	endif

	return;

END FUNCTION S_i_homog_koop


REAL FUNCTION frozen_fraction(wx, number_mr_D, s_w_pc)
 IMPLICIT NONE
REAL, INTENT(IN) :: wx, number_mr_D, s_w_pc 
INTEGER :: i, j, ismax, is
REAL :: minfrac, maxfrac, froz_frac_1, froz_frac_2, s_norm, s_centre(5)
! REAL :: frozen_fraction_2D
INTRINSIC :: log10


if(wx /= 0.) then
s_norm = s_w_pc * ABS(number_mr_D/wx)
else
s_norm = 0.
endif

 minfrac = 0.01
 maxfrac = 1.

s_centre(1) = 0.
s_centre(2) = 150.
s_centre(3) = 250.
s_centre(4) = 350.
s_centre(5) = 700.
ismax = 5;

if(s_norm <= s_centre(1)) then
frozen_fraction = frozen_fraction_2D(wx, number_mr_D,1)
return
endif

if(s_norm >= s_centre(ismax)) then
frozen_fraction = frozen_fraction_2D(wx, number_mr_D,ismax)
return
endif

is = ismax
do j  = 2,ismax
if(s_centre(j) > s_norm)then
	is = j-1
	exit
endif
enddo
if(is == ismax) stop 4321

froz_frac_1 = frozen_fraction_2D(wx, number_mr_D,is)
froz_frac_2 = frozen_fraction_2D(wx, number_mr_D,is+1)

frozen_fraction = froz_frac_1 + (s_norm - s_centre(is))* &
		(froz_frac_2 - froz_frac_1)/(s_centre(is+1) - s_centre(is))

if(frozen_fraction > maxfrac)  frozen_fraction = maxfrac
if(frozen_fraction < minfrac)  frozen_fraction = minfrac

       return
       END FUNCTION frozen_fraction

   REAL FUNCTION frozen_fraction_2D(w_vel, n_mr_D, isx)
implicit none
INTEGER, INTENT(IN) :: isx
REAL, INTENT(IN) :: w_vel, n_mr_D
REAL :: w_min, nd_centre(6), m(6), c(6), froz_frac_1, froz_frac_2, maxfrac, minfrac, wx
INTEGER :: ndmax, inmin, inmax, in, j

 minfrac = 0.01
 maxfrac = 1.

wx = w_vel;
w_min = 0.1
if(wx < w_min) wx = w_min


nd_centre(1) = 20.*0.316
nd_centre(2) = 20.*3.16
nd_centre(3) = 20.*17.32
nd_centre(4) = 20.*54.77
nd_centre(5) = 20.*316.3
nd_centre(6) = 20.*3162.
ndmax = 6

if(isx == 1) then
	m(4) = 0.4628;
	m(5) = 0.3760;
	m(6) = 0.4435
	c(4) = -0.6281
	c(5) = -0.8107
	c(6) = -0.8609
	inmin = 4
	inmax = 6
endif

if(isx == 2) then
	m(3) = 0.6218
	m(4) = 0.4019
	m(5) = 0.4009 
	c(3) = -0.4199 
	c(4) = -0.6302
	c(5) = -0.7756
	inmin = 3
	inmax = 5
endif

if(isx == 3) then
	m(3) = 0.5334
	m(4) = 0.4541
	m(5) =  0.4250
	m(6) = 0.5678 
	c(3) = -0.4004 
	c(4) = -0.5549 
	c(5) = -0.7446
	c(6) = -0.9655
	inmin = 3
	inmax = 6

endif

if(isx == 4) then
	m(2) =  0.5683
	m(3) = 0.4299
	m(4) =  0.5200
	m(5) = 0.4627 
	c(2) = -0.2305
	c(3) = -0.3579
	c(4) = -0.5233
	c(5) = -0.7417
	inmin = 2
	inmax = 5

endif

if(isx == 5) then
	m(3) = 0.3314
	m(4) = 0.3156
	m(5) = 0.3235
	m(6) = 0.3670
	c(3) = -0.2759
	c(4) = -0.3647
	c(5) = -0.5378
	c(6) = -0.7189
	inmin = 3
	inmax = 6

endif

if(n_mr_D <= nd_centre(inmin)) then
	frozen_fraction_2D = 10.**(m(inmin)*log10(wx) + c(inmin))
	if(frozen_fraction_2D > maxfrac)  frozen_fraction_2D = maxfrac
	if(frozen_fraction_2D < minfrac)  frozen_fraction_2D = minfrac
	return
endif

if(n_mr_D >= nd_centre(inmax)) then
	frozen_fraction_2D = 10.**(m(inmax)*log10(wx) + c(inmax))
	if(frozen_fraction_2D > maxfrac)  frozen_fraction_2D = maxfrac
	if(frozen_fraction_2D < minfrac)  frozen_fraction_2D = minfrac
	return
endif

in = inmax
do j  = 2,inmax
if(nd_centre(j) > n_mr_D)then
	in = j-1
	exit
endif
enddo

froz_frac_1 = 10.**(m(in)*log10(wx) + c(in))
if(froz_frac_1 > maxfrac) froz_frac_1 = maxfrac
if(froz_frac_1 < minfrac) froz_frac_1 = minfrac

froz_frac_2 = 10.**(m(in+1)*log10(wx) + c(in+1))
if(froz_frac_2 > maxfrac) froz_frac_2 = maxfrac
if(froz_frac_2 < minfrac) froz_frac_2 = minfrac

frozen_fraction_2D = froz_frac_1 + (n_mr_D - nd_centre(in)) *  &
		(froz_frac_2 - froz_frac_1)/(nd_centre(in+1) - nd_centre(in))

if(frozen_fraction_2D > maxfrac)  frozen_fraction_2D = maxfrac
if(frozen_fraction_2D < minfrac)  frozen_fraction_2D = minfrac


      return
       END FUNCTION frozen_fraction_2D



SUBROUTINE primary_ice_nucleation(temperature, pressure, rho_a,  w_prim, qv_prim, ql_prim, &
	qi_prim, nw_prim, ni_prim, nin_a_prim, &
	 XLF, XLS, XLV, CP, PIE, EPS, ESIT, DESIT, &
	 ESWT, DESWT, temp_K, rvapor, dt_prim, p_w, &
	LAMBDA_FAC_W, a_w, b_w, c_w, N0_FAC_w)
implicit none
integer :: CONTACT, ULTRA_COLD_HET_ICE_NUCL_DEMOTT
real :: PSI_DUST, DEMOTT_CORRECTION_FACTOR
!PARAMETER(CONTACT = 1, DEMOTT_CORRECTION_FACTOR = 0.063, ULTRA_COLD_HET_ICE_NUCL_DEMOTT = 1, PARAMETER(CONTACT = 1, DEMOTT_CORRECTION_FACTOR = 1., ULTRA_COLD_HET_ICE_NUCL_DEMOTT = 1, PSI_dust = 1.)
REAL, INTENT(IN) :: XLF, XLS, XLV, CP, EPS, pressure, w_prim, rho_a, PIE, rvapor, dt_prim, N0_FAC_w
REAL, INTENT(IN) :: p_w, LAMBDA_FAC_W, a_w, b_w, c_w
DOUBLE PRECISION, INTENT(IN) :: ESIT(150001), &
	 DESIT(150001), ESWT(150001), DESWT(150001), temp_K(150001)
REAL, INTENT(INOUT) :: temperature, qv_prim, ql_prim, qi_prim, nw_prim, ni_prim, nin_a_prim

REAL ::  SS_i, QSI, QSW, ESW, P_SAT, TT, ESI, m_i0, N_in, Q_is, Q_ws, RIHENC, CIHENC, N_in_ultra

REAL :: DEL_QW, DEL_QV, DEL_NW, N_in_contact, supersat, delta_ni_contact, delta_qi_contact, &
	 d_delta_qi_contact, d_delta_ni_contact
 

	if(ni_prim > 0. .and. qi_prim <= 0.) stop 1366		
	if(qi_prim > 0. .and. ni_prim <= 0.) stop 1367		
	

	TT = temperature
	P_SAT= pressure/100.
 
  	ESI = REAL(ESATI_int(DBLE(TT), ESIT, temp_K))
      	QSI = EPS*ESI/(P_SAT-ESI)
  	ESW = REAL(ESATW_int(DBLE(TT), ESWT, temp_K))
      	QSW = EPS*ESW/(P_SAT-ESW)
	Q_is = QSI
	Q_ws = QSW
         supersat = 100.*(qv_prim/Q_ws - 1.)
	if(qv_prim >  QSI) then
		if(temperature < 273.15 .and. temperature > 273.15 - 80. ) then	
			if(qv_prim < Q_ws) then
				SS_i = qv_prim/Q_is - 1.
			else
				SS_i = Q_ws/Q_is - 1.
			endif

			if(w_prim > 0. .and. temperature < 273.15 - 5. .and. temperature >= 273.15 - 30. ) then
				N_in =PSI_dust* 1.E3* exp(12.96*SS_i - 0.639) * DEMOTT_CORRECTION_FACTOR
			else
				N_in =  0.
			endif
			
			if(w_prim > 0. .and. temperature < 273.15 - 30. .and. ULTRA_COLD_HET_ICE_NUCL_DEMOTT == 1) then
				N_in_ultra = PSI_dust*1000.*(exp(0.1296*(SS_i*100.-10.))**0.3)  ! email on 23 Apr 2005 from Paul DeMott
			else
				N_in_ultra =  0.
			endif
			N_in = N_in + N_in_ultra;

			delta_ni_contact = 0.; delta_qi_contact = 0.;
                        if(CONTACT == 1) then
! Meyers et al. (1992)
! For now, make no attempt to track the number of contact-IN activated (Ovtchinnikov and Kogan 2000) state that fraction of contact-IN lost in one evaporation cycle is tiny )

				if(ql_prim > 0. .and. nw_prim > 0. .and. temperature < 273.15 - 3. .and. temperature > 273.15 - 35. ) then
					N_in_contact = -2.80 + 0.262 * (273.15 - temperature);
 	                             	N_in_contact = PSI_dust*1000.*(2.718281**N_in_contact) * DEMOTT_CORRECTION_FACTOR
					call contact_driver(dt_prim, N_in_contact, nw_prim, ql_prim, N0_FAC_w,  LAMBDA_FAC_w, temperature, pressure, &
						rho_a, supersat, ESW, XLV, EPS, rvapor, PIE, a_w, b_w, c_w, p_w, delta_ni_contact, delta_qi_contact)					
				else
					N_in_contact = 0.
				endif
                       endif


			CIHENC = N_in - nin_a_prim * rho_a
			if(CIHENC/rho_a < 0.) CIHENC = 0.

			ni_prim =  ni_prim + CIHENC/rho_a
			nin_a_prim = nin_a_prim + CIHENC/rho_a
			ni_prim =  ni_prim + delta_ni_contact	

			m_i0 = 920. * (4./3.)*PIE*(5.E-6 **3.)
			RIHENC = CIHENC*m_i0/rho_a
			qi_prim = qi_prim + RIHENC + delta_qi_contact

			if(qi_prim > 0.1) then
				print *, temperature, CIHENC, ql_prim, qv_prim/Q_is, qv_prim/Q_ws, ni_prim, nin_a_prim 
				stop 3164
			endif



			if( ql_prim > 0.) then
				if(ql_prim >= RIHENC + delta_qi_contact) then
					DEL_QW =  -(RIHENC + delta_qi_contact)
					DEL_QV = 0.
					DEL_NW = -nw_prim*RIHENC/ql_prim - delta_ni_contact
					if(DEL_NW + nw_prim < 0.) DEL_NW = -nw_prim
				else
					DEL_QV = -(RIHENC + delta_qi_contact) + ql_prim 
					DEL_QW = -ql_prim	
					DEL_NW = -nw_prim
				endif
			else
				DEL_QW = 0.			
				DEL_QV = -RIHENC
				DEL_NW = 0.
			endif

			if(DEL_QW > 0.) stop 4156
			if(DEL_QV > 0.) stop 4157

			ql_prim = ql_prim + DEL_QW
			if(ql_prim < 0.) stop 667

			nw_prim = nw_prim + DEL_NW	! DEL_NW <= 0
			if(nw_prim < 0.) stop 668

			qv_prim = qv_prim + DEL_QV

			temperature = temperature - XLF*DEL_QW/cp	! DEL_QW <= 0
			temperature = temperature - XLS*DEL_QV/cp	! DEL_QV <= 0
			
		else
			CIHENC = 0.
			RIHENC = 0.
		endif
	endif

	if(ni_prim > 0. .and. qi_prim <= 0.) stop 1376		
	if(nw_prim > 0. .and. ql_prim <= 0.) stop 1375		
	if(qi_prim > 0. .and. ni_prim <= 0.) stop 1377		



  END SUBROUTINE primary_ice_nucleation

SUBROUTINE contact_driver(dt_prim, N_in_contact, nw_prim, ql_prim, N0_FAC_w, LAMBDA_FAC_w, temperature, pressure,rho_a, &
					supersat, ESW_mb, XLV, EPS, rvapor, PIE, a_w, b_w, c_w, p_w, delta_ni, delta_qi)
implicit none
INTEGER :: NKR
PARAMETER(NKR=33)
real, intent(in) :: dt_prim, nw_prim, ql_prim, N0_FAC_w, LAMBDA_FAC_w, temperature, pressure,rho_a, &
			supersat, ESW_mb, XLV, EPS, rvapor, PIE, a_w, b_w, c_w,  p_w, N_in_contact

real, intent(out) :: delta_ni, delta_qi
real :: m_water(NKR), r_water(NKR), n_water(NKR), D_water(NKR), n_0w_ov_lambdapower, &
	lambda_w, sum_q, dD_water(NKR), sum_n, NUB, NUC, NUD, d_delta_ni, d_delta_qi
INTEGER :: i


	delta_ni = 0.;  delta_qi = 0.;
	NUB = 0.; NUC = 0.; NUD = 0.;
	lambda_w = LAMBDA_FAC_w*c_w*nw_prim/ql_prim ! CORRECT
	lambda_w = lambda_w**(1./3.) ! CORRECT
	n_0w_ov_lambdapower = (nw_prim*rho_a)*N0_FAC_w ! CORRECT


	m_water(1) = 1000.*(4./3.)*PIE*(1.e-6**3.)
	do i=2,NKR
		m_water(i) = 1.5*m_water(i-1)
	enddo
	do i=1,NKR
		r_water(i) = (m_water(i)/1000.)/(4.*PIE/3.)
		r_water(i) = r_water(i)**(1./3.)
		D_water(i) = 2.*r_water(i)

!		print *, D_water(i) =, D_water(i), i

	enddo
        dD_water(1) =  D_water(2) - D_water(1) 
        dD_water(NKR) =  D_water(NKR) - D_water(NKR-1)
 
	do i=2,NKR-1
		dD_water(i) = 0.5*(D_water(i+1) - D_water(i-1))
	enddo
	sum_n = 0.
	do i=1,NKR
		n_water(i) = dD_water(i)*n_0w_ov_lambdapower*(D_water(i)**p_w)*exp(-lambda_w*D_water(i))/rho_a
		sum_n = sum_n + n_water(i)
	enddo
	if(sum_n > 0.) then
		do i=1,NKR
			n_water(i) = n_water(i) * nw_prim/sum_n
		enddo
		sum_q = 0.
		do i=1,NKR
			sum_q = sum_q + n_water(i)*m_water(i)
		enddo
		if(sum_q > 0.) then
			do i=1,NKR
				n_water(i) = n_water(i) * ql_prim/sum_q
			enddo
			sum_q = 0.
			do i=1,NKR
				sum_q = sum_q + n_water(i)*m_water(i)
			enddo
			if(ABS(sum_q - ql_prim)/ql_prim > 0.01) stop 32345

			do i=1,NKR
				if(n_water(i) > 0. .and. ql_prim > 0.) then
					call contact_sources(NUB, NUC, NUD, m_water(i), r_water(i), N_in_contact, &
						supersat, ESW_mb*100., n_water(i), temperature, pressure, rho_a, XLV, &
						EPS, rvapor, PIE, a_w, b_w)

					d_delta_qi = dt_prim*(NUB + NUC + NUD)
					if(d_delta_qi > n_water(i)*m_water(i)) d_delta_qi = n_water(i)*m_water(i)

					d_delta_ni = d_delta_qi/m_water(i)
					if(d_delta_ni > n_water(i))d_delta_ni = n_water(i)

					if(NUB + NUC + NUD > 0.) then 
						delta_ni = delta_ni + d_delta_ni
						delta_qi = delta_qi + d_delta_qi
					endif
				else
					NUB = 0.;  NUC = 0.; NUD = 0.; 
				endif
				if(delta_ni < 0.) stop 65434
 			enddo
		else
			print *, 'contact freezing: sum_q =', sum_q
		endif
	else
		print *, 'contact freezing: sum_n =', sum_n
	endif
	if(delta_qi > ql_prim) delta_qi = ql_prim;
	if(delta_ni > nw_prim) delta_ni = nw_prim;


END SUBROUTINE contact_driver


SUBROUTINE contact_sources(NUB, NUC, NUD, mi, r_water, N_in_contact, supersat_pc, es, nw_co, temp_co, &
	pres_co, rhoa_co, XLV, EPS, RV, PIE, a_w, b_w)
implicit none

REAL, INTENT(IN) :: N_in_contact, supersat_pc, es, XLV, EPS, RV, PIE, nw_co, temp_co, &
	pres_co, rhoa_co, mi, r_water, a_w, b_w
REAL, INTENT(OUT) :: NUB, NUC, NUD
REAL:: AERORAD, JOULES_IN_AN_ERG, JOULES_IN_A_CAL, K, dd, Fk, Fd, F1, F2, fact, mfpath 
REAL:: 	Kn, etaa, tdiff, ft, K_aerosol, dfa2, alph,  B, te, rad, vt_drop, N_Peclet

! NUB, NUC, NUD are all mixing ratio tendencies of ice mass, with each crystal weighing 

	NUB = 0.; NUC = 0.; NUD = 0.;

!	AERORAD = 0.075e-6; ! EMM
!	AERORAD = 0.1e-6; ! Meyers et al. (1992):  they criticize Young for choosing too large a diameter (~1 um)
!	AERORAD = 0.2e-6; ! Swann 
	AERORAD = 0.3e-6; ! Young (1974), Cotton et al. (1986), Ovtchinnikov and Kogan (2000);  Hess et al. (1998, BAMS)

	JOULES_IN_AN_ERG = 1.0e-7
	JOULES_IN_A_CAL = 4.187

	rad = r_water

	if (mi > 1.0e-17) then
		fact = 4. * PIE * rad * nw_co*rhoa_co
	else
		fact = 0.
	endif
	F1 = fact * N_in_contact;
	te = temp_co - 273.15;

	K = (5.69 + 0.017 * te) * 1.0e-3 * JOULES_IN_A_CAL  ;
	dd = 0.211 * ((temp_co / 273.15)**1.94) * (101325. / pres_co) / 1.e4;
	Fk =  (XLV / (RV * temp_co) - 1.) * XLV * 1000. / (K * temp_co);
	Fd = 1000. * RV * temp_co / (dd * es);

! work out the aerosol diffusivity - it should be O(1.e-10) 

	Kn = 7.37e-8 * temp_co * 100000.0 / (288. * pres_co * AERORAD);
	mfpath = 6.6e-8 * (101325. / pres_co) * (temp_co / 293.15); 
	Kn = mfpath / AERORAD;
	etaa = (1.718 + 0.0049 * te - 1.2e-5 * te * te) * 1.e-5;

	alph = 1.45 + 0.4 * exp(-1./Kn);  
	B = (1. + alph * Kn) / (6. * PIE * AERORAD * etaa);
	dfa2 = 1.3804e-16 * JOULES_IN_AN_ERG * temp_co * B;

! Brownian motion

	NUB = mi * dfa2 * F1 / rhoa_co;

! tdiff = Tc - T  :  this formula for the diffusion coefficient (dd) is from Eqn 13-3, of P & K, which agrees with Table 7.1 of R & Y 

	tdiff =  (1000. / (Fk + Fd)) * (supersat_pc / 100.) * XLV / K;

!					tdiff = (LV/K)*RHOW*(SUPERSAT/100.) * 70e-12; 
!					printf("In the updraft the average cloud droplet is %f K warmer than the ambient cloudy air\n",tdiff);

	F2 = -K * tdiff / pres_co;

!	K_aerosol = 5.39 * 1.0e-4 * 100. * JOULES_IN_AN_ERG;
	K_aerosol = 0.25 ! Ovtchinnikov et al. (2000)
! Daviess constants */
!  ft = 0.4*(1.+1.257*Kn+0.4*Kn*exp(-1.10/Kn))*(K+2.5*Kn*K_aerosol)/( (1.+3.*Kn)*(2.*K+5.*K_aerosol*Kn+K_aerosol) );  

	ft = 0.4*(1.+1.45*Kn+0.4*Kn*exp(-1.0/Kn))*(K+2.5*Kn*K_aerosol)/((1.+3.*Kn)*(2.*K+5.*K_aerosol*Kn+K_aerosol));
! Thermophoresis
	NUC = mi * F1 * F2 * ft / rhoa_co;
! Diffusiophoresis
	NUD = -NUC * (RV * temp_co / (ft * XLV)) * (4.8 / 4.) * EPS;

	vt_drop = a_w*((2.*rad)**b_w)
	N_Peclet = vt_drop*2.*rad/dfa2

!	print *,  NUB/NUC, 2.*rad*1.e6, supersat_pc:::::, NUB/NUC, 2.*rad*1.e6, supersat_pc
!	print *, Peclet #, Brownian (NUB), Thermophoresis (NUC), Diffusiophoresis (NUD), D(um) :::::::, N_Peclet, NUB, NUC, NUD, 2.*rad*1.e6
	if(N_Peclet > 1.25 .and. ABS(NUC) > 2.*ABS(NUB)) then
		NUB = 0.
	endif


  END SUBROUTINE contact_sources

SUBROUTINE melting_cloudice(tz_mlt, qlz_mlt, nwz_mlt, qiz_mlt, niz_mlt,  XLF, CP)
implicit none
REAL, INTENT(IN) :: XLF, CP
REAL, INTENT(INOUT) :: tz_mlt, qlz_mlt, nwz_mlt, qiz_mlt, niz_mlt

REAL :: dtmlt, del_qw, del_nw


	if(niz_mlt > 0. .and. qiz_mlt <= 0.) stop 1456		
	if(qiz_mlt > 0. .and. niz_mlt <= 0.) stop 1457		

	
	if(tz_mlt > 273.15 .and. qiz_mlt > 0.) then

		if(tz_mlt - XLF*qiz_mlt/CP > 273.15) then
		
			tz_mlt = tz_mlt - XLF*qiz_mlt/CP
			qlz_mlt = qlz_mlt + qiz_mlt
			nwz_mlt = nwz_mlt + niz_mlt 
			qiz_mlt = 0.
			niz_mlt = 0.
		else
			dtmlt = tz_mlt - 273.15 		! dtmlt > 0. 
			tz_mlt = 273.15
			
			del_qw = CP*dtmlt/XLF		! del_qw > 0.
			if(del_qw < 0.) stop 1547

			del_nw = niz_mlt*del_qw/qiz_mlt  !  del_nw > 0.


			qiz_mlt = qiz_mlt - del_qw
			if(qiz_mlt < 0.) then
				qiz_mlt = 0.
			endif
			if(qiz_mlt > 1.) stop 3161

			qlz_mlt =  qlz_mlt + del_qw
			
			niz_mlt = niz_mlt - del_nw 
			nwz_mlt = nwz_mlt + del_nw  
			
		endif

	endif

END SUBROUTINE melting_cloudice


SUBROUTINE condensation_of_cloudwater(s_pc_start, temperature, pressure, qv_vap, ql_vap, nw_vap, EPS, &
     ESIT, DESIT, ESWT, DESWT, temp_K, CPLC)
implicit none
REAL, INTENT(IN) :: s_pc_start, CPLC, pressure, EPS
DOUBLE PRECISION, INTENT(IN) :: ESIT(150001), DESIT(150001), ESWT(150001), DESWT(150001), temp_K(150001)
REAL, INTENT(INOUT) :: temperature, qv_vap, ql_vap, nw_vap 
REAL :: TQLK,  TT, P_SAT, QVK, QLK, ESW, QSW 
INTEGER :: ISTOP

	TQLK = 0.
	ISTOP = 0
	if(temperature > 273.15 - 100.) then
		TT = temperature
		P_SAT= pressure/100.
		QVK = qv_vap
		QLK = ql_vap

  		ESW = REAL(ESATW_int(DBLE(TT), ESWT, temp_K))
      		QSW = EPS*ESW/(P_SAT-ESW)
		if(ABS(QVK/QSW - 1.) > 1.E-8) then
			if(QVK > QSW) then 
				TQLK = 0.

!				print *, QVK, QLK BEFORE:: , QVK, QLK
! 				print *, QSW, T BEFORE::, QSW, TT
				ISTOP = 0
   	  			CALL CONDNS(TT,P_SAT,QVK,QSW,ESW,TQLK,ISTOP, &
				CPLC, ESWT, DESWT, temp_K, EPS)
				if(TQLK < 0.) then 
					TQLK = 0.
					TT = temperature
					QVK = qv_vap
				endif
      				QLK = QLK+TQLK
!				print *, QVK, QLK AFTER:: , QVK, QLK

  				ESW = REAL(ESATW_int(DBLE(TT), ESWT, temp_K))
      				QSW = EPS*ESW/(P_SAT-ESW)
!				print *, QSW, T AFTER::, QSW, TT
				IF( ABS(QVK/QSW - 1.) > 1.E-5) STOP 1578
			else
      				TT = TT-CPLC*QLK
      				QVK = QVK+QLK
      				QLK = 0.0

			       ESW = REAL(ESATW_int(DBLE(TT), ESWT, temp_K))
				ISTOP = 0
      			       IF(QVK/QSW-1.0 > 1.E-8) then
      					CALL CONDNS(TT,P_SAT,QVK,QSW,ESW,QLK,ISTOP, &
				CPLC, ESWT, DESWT, temp_K, EPS)
 				endif
  				if(QLK <= 0.) then
					nw_vap = 0.
					QLK = 0.
				endif
			endif

			temperature = TT
			qv_vap = QVK
			if(qv_vap < 0.) stop 3160

			ql_vap = QLK
		endif
	endif
	



!	CALL SAT(i, k, tz(k), pressure/100., qvz(k), qlz(k), &
!	nwz(k), ncn_az(k), qiz(k), niz(k),  nin_az(k),&
!       CPLC,CPLS,CPLF,TDIF,ESW00,ESI00,ESWT,ESIT, &
!	DESWT,DESIT, TT,TTD, HLTC,HLTF,TICE,EPS,TTFRZ, & 
!	PVI, PIL, PVL)

!	if(PVL .gt. 0.) then
!	G_COND(i,k,j) = G_COND(i,k,j) + PVL/dt
!	else
!	G_EVAP(i,k,j) = G_EVAP(i,k,j) - PVL/dt
!	endif
!	if(PVI .gt. 0.) then
!	G_DEP(i,k,j) = G_DEP(i,k,j) + PVI/dt
!	else
!	G_SUBL(i,k,j) = G_SUBL(i,k,j) - PVI/dt
!	endif
!	if(PIL .lt. 0.) then
!	G_FREEZE(i,k,j) = G_FREEZE(i,k,j) - PIL/dt
!	else
!	G_MELT(i,k,j) = G_MELT(i,k,j) + PIL/dt
!	endif





END SUBROUTINE condensation_of_cloudwater


SUBROUTINE threshold_number(d_bar_cx, qiz_x, niz_x, naero_x, D_CX_MAX,  d_x, c_x, LAMBDA_FAC, DBAR_FAC)
implicit none

real, intent(in) :: D_CX_MAX,  d_x, c_x, LAMBDA_FAC, DBAR_FAC
real, intent(inout) ::  d_bar_cx, qiz_x, niz_x, naero_x

real :: lambda_min, niz_min, lambda_x

	if(qiz_x > 1.e-20) then
		lambda_min = DBAR_FAC/D_CX_MAX

		niz_min = (lambda_min**d_x)*qiz_x/(LAMBDA_FAC*c_x)	
		if(niz_x < niz_min) then
			naero_x = naero_x + niz_min - niz_x
  			niz_x = niz_min
		endif

		lambda_x = LAMBDA_FAC*c_x*niz_x/qiz_x ! CORRECT
		lambda_x = lambda_x**(1./d_x) ! CORRECT
		if(lambda_x <= 0.) stop 13678
		d_bar_cx = DBAR_FAC/lambda_x
		if(d_bar_cx > D_CX_MAX) then
			print *, d_bar_cx, D_CX_MAX, qiz_x, niz_x, niz_min 
			stop 13487
		endif
                
	else
		qiz_x = 0.
		niz_x = 0.
		d_bar_cx = 0.
	endif

END SUBROUTINE threshold_number



REAL FUNCTION twomey_s_cb(U_t, C_t, k_t)
REAL, INTENT(IN):: U_t, C_t, k_t
REAL :: s_cbx


	s_cbx= 1.6e-3 * ((100.*U_t)**(1.5))/(C_t/1.e6);
	s_cbx =  3.6 * ((s_cbx)**(1./(k_t + 2.)));
	twomey_s_cb = s_cbx

	return

END FUNCTION twomey_s_cb


SUBROUTINE homogeneous_freezing(tz_hom, qiz_hom,  qlz_hom, niz_hom,  nwz_hom, ncn_az_hom, &
	qvz_hom, LAMBDA_FAC_w, p_w, XLF, XLV, CP, T_FRZ_HOM_DEGC, c_w, w_vel, s_hom, tz, kts, kte, k, num_droplets, rho )
implicit none
integer :: SIMPLE_FRAC
PARAMETER(SIMPLE_FRAC = 0)
real, intent(in) ::   XLF, XLV,CP, T_FRZ_HOM_DEGC, LAMBDA_FAC_w, c_w, p_w, w_vel, s_hom(kts:kte), tz(kts:kte), rho
 real, intent(inout) :: tz_hom, qiz_hom,  qlz_hom, niz_hom,  nwz_hom, ncn_az_hom, qvz_hom, num_droplets
integer, intent(in) :: kts, kte, k
real :: dtfr, del_ni, del_qi,del_ql,del_qv,lambda_w, dstar, fraction_mass_frozen, fraction_num_frozen, mass_air
integer :: kx, kfreeze
	if(nwz_hom < 0.) stop 13655
	if(qlz_hom < 0.) then 
		stop 13656
	endif
	if(qiz_hom < 0.) stop 13657
	if(niz_hom < 0.) stop "136581"

	mass_air = 500.*2000.*2000.*rho
	if(tz_hom < 273.15 + T_FRZ_HOM_DEGC .and. qlz_hom > 0. .and. k > 1) then

	
		lambda_w = LAMBDA_FAC_w*c_w*nwz_hom/qlz_hom ! CORRECT
		lambda_w = lambda_w**(1./3.) ! CORRECT
! assume there is a grid level between T_FRZ_HOM_DEGC and T_FRZ_HOM_DEGC - 6.
		if(tz_hom > 273.15 + T_FRZ_HOM_DEGC - 6.) then

			if(SIMPLE_FRAC == 1) then
				fraction_num_frozen = 0.7
			else
				kfreeze = 0
				do kx = k,1, -1
					if(tz(kx) > 273.15 + T_FRZ_HOM_DEGC) then
						kfreeze = kx
						exit
					endif
				enddo 
				if(kx == 1) stop "991111"
				if(kfreeze == 0) stop 91191

				fraction_num_frozen = frozen_fraction(w_vel, nwz_hom*(1.+p_w)/lambda_w, s_hom(kfreeze))
				if(p_w /= 3.5) print *, 'WARNING:  need to re-do lookup table for new p_w'
				print *, 'HOMOGEN FREEZE: fraction of droplets freezing = ', fraction_num_frozen, &
					'with assumed s(%) at -35 degC = ', s_hom(kfreeze) 

			endif

			dstar=find_lambdadstar(p_w, fraction_num_frozen)/lambda_w
			if(dstar < 0.) stop "135675"
			if(dstar > 1.) stop "136582"
			fraction_mass_frozen = gammq_frac(4.+p_w,dstar*lambda_w)

			if(fraction_mass_frozen < 0.) stop "135635"
			if(fraction_mass_frozen > 1.) stop "136583"

		else			
			fraction_num_frozen = 1.
			fraction_mass_frozen = 1.
		endif

		del_ql = qlz_hom*fraction_mass_frozen
		del_qv = qlz_hom - del_ql

		if(del_ql < 0.) stop 3154
		if(del_ql > qlz_hom) stop 3157
		if(del_qv < 0.) stop 3184

		if(tz_hom + XLF*del_ql/CP - XLV*del_qv/CP < 273.15 + T_FRZ_HOM_DEGC) then
! FREEZING
			tz_hom = tz_hom + XLF*del_ql/CP
			qiz_hom = qiz_hom + del_ql
			if(qiz_hom > 1.) stop 3162
			niz_hom = niz_hom + fraction_num_frozen*nwz_hom 
! EVAPORATION
			tz_hom = tz_hom - XLV*del_qv/CP
			qvz_hom = qvz_hom + del_qv
			ncn_az_hom = ncn_az_hom - (1. - fraction_num_frozen)*nwz_hom 
			if(ncn_az_hom < 0.) ncn_az_hom = 0.

			num_droplets = num_droplets + (fraction_num_frozen*nwz_hom)*mass_air
			qlz_hom = 0.
			nwz_hom = 0.

		else
			print *, 'HOMOGEN FREEZE: marginal temp :: coming here 123'
			dtfr = (273.15 + T_FRZ_HOM_DEGC) - tz_hom   ! dtfr > 0.
			del_qi = CP*dtfr/XLF		! del_qi > 0.
                        if(del_qi > qlz_hom) del_qi = qlz_hom 

			del_ni = nwz_hom* del_qi/qlz_hom	   ! del_ni > 0.
			if(del_ni < 0.) stop 1367

			qlz_hom = qlz_hom - del_qi
			if(qlz_hom < 0.) stop 531
			qiz_hom = qiz_hom + del_qi
			if(qiz_hom > 1.) stop 3163
			tz_hom = tz_hom + XLF*del_qi/CP
			num_droplets = num_droplets + del_ni*mass_air
			nwz_hom = nwz_hom - del_ni  
			niz_hom = niz_hom + del_ni 

		endif
	endif
	if(niz_hom > 0. .and. qiz_hom <= 0.) stop 1656		
	if(qiz_hom > 0. .and. niz_hom <= 0.) stop 1657		

  END SUBROUTINE homogeneous_freezing


  
     SUBROUTINE IMICRO(dt, ttt,PRESR,QMIXY,TSSY, QNW, QNI, TSSNW,  TSSNI,RHOK, &
	   TDT, RHOS, CRACS, CSACR,CGACR,CGACS,ACCOY,CSACW &
          ,CRACI,CIACR, CSACI,CGACW,CGACI,CRACW,CSSUB,CGSUB,CREVP  &
          ,CGFR,CSMLT,CGMLT,CGWET,CRAUT,QI0,QS0,QL0,ES0 &
          ,CES0,C1BRG,C2BRG,RMI50,RMI40,ESWT,ESIT,DESWT,DESIT &
          ,temp_K,HLTS,HLTC,HLTF,CH2O,CICE,TICE,CP,EPS, &
	   VCONR,VCONS,VCONG, GAMMA_K, d_i, RHO_CLOUDICE, c_i, &
	  c_w, LAMBDA_FAC_w, LAMBDA_FAC_i, p_w, p_i, DBAR_FAC_w, lambda_i0, DIS)

	IMPLICIT NONE
	REAL :: HUCM_TURB_SACI, HUCM_TURB_SACW, HUCM_TURB_GACW
	INTEGER :: ICE_NUMBER, ADVANCED_HM_FORMULA, ADVANCED_AUTOCONV_RAIN
	PARAMETER(ICE_NUMBER = 1, ADVANCED_HM_FORMULA = 0, ADVANCED_AUTOCONV_RAIN = 1, &
		HUCM_TURB_SACI = 1., HUCM_TURB_SACW = 1., HUCM_TURB_GACW = 1.)
!		HUCM_TURB_SACI = 5., HUCM_TURB_SACW = 10., HUCM_TURB_GACW = 2.)

      REAL             ::   DD=0.15,D_W_AUTO=20.E-6, D_XW=15.E-6, D_XI=40.E-6, D_W_MAX = 50.e-6
	REAL, INTENT(IN)::  ttt,PRESR,RHOK, DIS
	REAL, INTENT(INOUT):: QMIXY(6), TSSY(7), ACCOY(3,4), GAMMA_K(6), TSSNW, TSSNI
 
	REAL, INTENT(IN) :: TDT, RHOS, dt, QNW, QNI, d_i, RHO_CLOUDICE, c_i, c_w, &
                 LAMBDA_FAC_w, LAMBDA_FAC_i, lambda_i0, p_w,  p_i, DBAR_FAC_w

	REAL, INTENT(IN) :: CRACS, CSACR,CGACR,CGACS,CSACW &
          ,CRACI,CIACR, CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5),CREVP(5)  &
          ,CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,ES0 &
          ,CES0,C1BRG,C2BRG,RMI50,RMI40
	REAL, INTENT(IN) :: HLTS,HLTC &
          ,HLTF,CH2O,CICE,TICE,CP,EPS
	DOUBLE PRECISION, INTENT(IN) :: ESIT(150001), DESIT(150001), ESWT(150001), DESWT(150001), temp_K(150001)
	REAL, INTENT(IN) ::  VCONR,VCONS,VCONG


	REAL :: RHOFAC,RHO, PSS(26), QMIX(6), TSS(7), ACCO(3,4)

	REAL :: epsi, NRAUT, NRACW, NSACW, NGACW, NIDW, NSFW, VAR_D_w, D_BAR_w
	REAL :: NSAUT, NSACI, NSFI, NRACI, NGACI, NWLONW_eff, NIHMS, PIHMS, NIHMG, PIHMG
!
!   Changes by Martin Koehler for 
!   I41: Ice Cloud Decay Sensitivity Experiment (ICDS-41)
!   * all thresholds neglecting small mixing ratio cases: * 10^-2
!   - threshold in any condensate to run microphysics: 10^-8 (normally 10^-6)
!   - threshold in each condensate to define its existance: 5*10^-8 (5*10^-6)
!   - other thresholds (CRIT, QVK/QS...-1) * 10^-2  (increased accuracy)
!
! *************************************************************************
! MODIFICATIONS BELOW ARE MADE TO LORDS ORIGINAL CODE (Q. Fu & S. Krueger)
!**************************************************************************
!       Updates (Jun. 8, 1993, fu comments):
!         1. RMINUC = QI * RHO / NC;  this  will give the depositional
!            growth of cloud ice at expense of cloud water  properly (
!            PIDW).
!         2. Replacing R50 by R100; then DT1 is the growth time needed
!            for ice crystal to grow from 40 to 100um (PSFI and PSFW).
!         3. Setting Psfi = 0.0 when QL is less than 5e-6 (PSFI).
!         4. A1 and A2 are obtained by using interpolation; A1T(14) is
!            0.1725e-9 for function realidw ( PSFI, PSFW and PIDW).
!         5. Setting QI0 = .6e-3 (threshold for cloud ice aggregation)
!            and QS0 = 1.0e-3 (threshold for snow aggregation ) (PSAUT
!            and PGAUT) (private communication with S. Chin).
!         6. GAM625 = 184.860962 (PWACS)
!         7. Esi and Egs ( T < T0 ) are set to be  0.1 following  RH84 
!            ( PSACI and PGACS ).
!**********************************************************************
!       Updates (Jun. 9, 1993, rh comments):
!         1. Modifying intercept parameters and densities.
!         2. Replacing Lins hail fall speed relation with an average
!            graupel relation (Ug=40.74*Dg**0.5) based on Locatelli &
!            Hobbs (1974).
!**********************************************************************
!       Updates (Jun. 10, 1993, lin comments):
!         1. Correcting C2BRG.
!         2. Egs = 1.0 for T > T0 and wet growth. Note:  modifications
!            here should be consistent to (7) under the comment cfu.
!**********************************************************************
!       Updates (Jun. 10, 1993, sk comments):
!         1. Lin et al. (1983) version is  installed for the following
!            processes: PSMLT, PGWET, PGMLT (constants only).
!            Change DIMENSION CSMLT(3),CGWET(2) => CSMLT(5), CGWET(4).
!         2. PWACS = 0.0 following Lin et al. (1983).
!         3. And more.
! ************************************************************************
!

 
 
!
      REAL :: C3RACS(3) &
       ,C3SACR(3),C3GACR(3),C3GACS(3)
!
       REAL :: QRHO(6),THRSH(6),SMAX(6)
!
!      REAL IDW
!
      LOGICAL QLOGIC(6),VAPR,LIQ,ICE,SNO,GRAUP,RAIN

      REAL :: QVK, QLK, QIK, QSK, QGK, QRK, QVR, QLR, QIR, &
	QSR, QGR, QRR, TSSV, TSSL, TSSI, TSSS, TSSG, TSSR, TTMP
      REAL :: SVMAX, SLMAX, SIMAX, SSMAX, SGMAX, SRMAX, &
	 PRAUT, PRACW, PRACS, PSACW, PWACS, PGACW, PGMLT
      REAL :: PSMLT, PREVP, PIACR, PSACR, PGACR, PGFR, PGACS, &
	 PGAUT, PGACI, PGWORD, PRACI, PSAUT, PSACI
      REAL :: PSSUB, PGSUB, PSFW, PSFI, PIDW, PGWET, QL0RHO, TC, &
	TSQ, EXPT1, EXPT2, ESW, QSW, DQS, DQS0, VTS, VTG, VTR, QEXRHO
      REAL :: PGDRY, SI,  QSI, ESI, C, QIKT, TEM,  SW, CTGACS, PSDEP
      INTEGER :: I, K
!
!	REAL :: VTRG, VTRS, VTRR, ACR1, ACR2, ACR3, GMLT, REVP, &
!	SMLT, REALIDW, ESATW, RAUT, AUT, GFR, ESATI, SSUB, GSUB, GWET

	REAL :: delt1, delt2, sources, sinks, factor
	REAL :: lambda_w, lambda_i
        REAL :: frac_num_collide_w, frac_mass_collide_w, QLK_eff, QNW_eff
        REAL :: frac_num_collide_i, frac_mass_collide_i, QIK_eff, QNI_eff
	LOGICAL :: rescale 
	INTRINSIC SQRT


      EQUIVALENCE (QMIX(1),QVK),(QMIX(2),QLK),(QMIX(3),QIK), &
       (QMIX(4),QSK),(QMIX(5),QGK),(QMIX(6),QRK),(QRHO(1),QVR), &
       (QRHO(2),QLR),(QRHO(3),QIR),(QRHO(4),QSR),(QRHO(5),QGR), &
       (QRHO(6),QRR),(QLOGIC(1),VAPR),(QLOGIC(2),LIQ), &
       (QLOGIC(3),ICE),(QLOGIC(4),SNO),(QLOGIC(5),GRAUP), &
       (QLOGIC(6),RAIN),(TSSV,TSS(1)),(TSSL,TSS(2)),(TSSI,TSS(3)), &
       (TSSS,TSS(4)),(TSSG,TSS(5)),(TSSR,TSS(6)),(TTMP,TSS(7)), &
       (SVMAX,SMAX(1)),(SLMAX,SMAX(2)),(SIMAX,SMAX(3)), &
       (SSMAX,SMAX(4)),(SGMAX,SMAX(5)),(SRMAX,SMAX(6))
!
 
      EQUIVALENCE (PRAUT,PSS(1)),(PRACW,PSS(2)),(PRACS,PSS(3)) &
       ,(PSACW,PSS(4)),(PWACS,PSS(5)),(PGACW,PSS(6)),(PGMLT,PSS(7)) &
       ,(PSMLT,PSS(8)),(PREVP,PSS(9)),(PIACR,PSS(10)),(PSACR,PSS(11)) &
       ,(PGACR,PSS(12)),(PGFR,PSS(13)),(PGACS,PSS(14)),(PGAUT,PSS(15)) &
       ,(PGACI,PSS(16)),(PGWORD,PSS(17)),(PRACI,PSS(18)),(PSAUT,PSS(19)) &
       ,(PSACI,PSS(20)),(PSSUB,PSS(21)),(PGSUB,PSS(22)),(PSFW,PSS(23)) &
       ,(PSFI,PSS(24)),(PIDW,PSS(25)),(PGWET,PSS(26))
!
  


      EQUIVALENCE (C3RACS(1),ACCO(1,1)),(C3SACR(1),ACCO(1,2)) &
       ,(C3GACR(1),ACCO(1,3)),(C3GACS(1),ACCO(1,4))


	rescale = .true.
!	rescale = .false.
	epsi = 1.e-5
!
!xmk  DATA THRSH/0.0,5*1.E-6/
      DATA THRSH/0.0,5*1.E-8/          !threshold / 100

      DO K=1,6
		QMIX(K) = QMIXY(K)
		TSS(K) = TSSY(K)
      ENDDO
      TSS(7) = TSSY(7)

	DO I = 1, 3
        DO K = 1,4
		ACCO(I,K) = ACCOY(I,K)
	ENDDO
	ENDDO

!
!     ***************************
!     ***************************
!     ****                   ****
!     ****   PRELIMINARIES   ****
!     ****                   ****
!     ***************************
!     ***************************
!
!
!     LOCAL VARIABLES -- TEMP, DENSITY, TEMP. DEPENDENT EFFICIENCIES
!     AND SATURATION VARIABLES
!
	if(RHOS > 2. .or. RHOS < 0.5) stop 9997

      RHO = RHOK
      RHOFAC = SQRT(RHOS/RHO)
!      QL0RHO = QL0 * RHO
      QL0RHO = c_w*(D_W_AUTO**3.)*QNW*RHO
	if(QLK > 0.) then
	      lambda_w = LAMBDA_FAC_w*c_w*QNW/QLK ! CORRECT
	      lambda_w = lambda_w**(1./3.) ! CORRECT
      		frac_mass_collide_w = gammq_frac(4.+p_w, lambda_w*D_XW)
      		if(frac_mass_collide_w > 1.) stop "1438799"
      		frac_num_collide_w = gammq_frac(1.+p_w, lambda_w*D_XW)
      		if(frac_num_collide_w > 1.) stop "1438798"
      		if(frac_num_collide_w > frac_mass_collide_w ) stop "1458498"
      		if( frac_mass_collide_w < 0.) stop "4438798"
      		if( frac_num_collide_w < 0.) stop "4434798"
	else
	      lambda_w = 0.
     		frac_mass_collide_w = 0.
      		frac_num_collide_w = 0.
	endif

	if(QIK > 0.) then
		lambda_i = LAMBDA_FAC_i*c_i*QNI/QIK
		lambda_i = lambda_i**(1./d_i)
       		frac_mass_collide_i = gammq_frac(4.+p_i, lambda_i*D_XI)
      		if(frac_mass_collide_i > 1.) stop "1438799"
      		frac_num_collide_i = gammq_frac(1.+p_i, lambda_i*D_XI)
      		if(frac_num_collide_i > 1.) stop "1438798"
      		if(frac_num_collide_i > frac_mass_collide_i ) stop "1458498"
      		if( frac_mass_collide_i < 0.) stop "4438798"
      		if( frac_num_collide_i < 0.) stop "4434798"
	else
	        lambda_i = 0.
     		frac_mass_collide_i = 0.
      		frac_num_collide_i = 0.
	endif

      TC = ttt-TICE
      TSQ = ttt**2
!
      EXPT1 = AMIN1(EXP(0.09*TC),1.0)
      EXPT2 = AMIN1(EXP(0.025*TC),1.0)
!
!      print *,HERE 12
 
      ESW   = REAL(ESATW_int(DBLE(ttt), ESWT, temp_K))
      QSW   = EPS*ESW/(PRESR-ESW)
      DQS   = QVK-QSW
      DQS0  = CES0/(PRESR-ES0)-QVK
!
!     ZERO SOURCE AND SINK TERMS
!
      DO K=1,26
  1   PSS(K) = 0.0
	ENDDO
!
!     DEFINE MIXING RATIOS GREATER THAN THRESHOLD FOR
!     ACCRETION PROCESSES AND MASS DENSITIES
!
!           1:  WATER VAPOR
!           2:  CLOUD (SUSPENDED LIQUID) WATER
!           3:  CLOUD ICE CRYSTALS
!           4:  SNOW
!           5:  GRAUPEL
!           6:  RAIN
!

	TSSNW = 0
	TSSNI = 0
      DO K=1,6
      TSS(K) = 0.0
      QLOGIC(K) = QMIX(K) .GT. THRSH(K)
      SMAX(K) = QMIX(K)/TDT
   10 QRHO(K) = QMIX(K) * RHO
	ENDDO
      TSS(7)  = 0.0
!
!     TERMINAL VELOCITIES
!
	VTS = 0.
	VTR = 0.
	VTG = 0.
 
      IF(SNO)    VTS = VTRS(QSR,VCONS,RHO,RHOFAC)
      IF(GRAUP)  VTG = VTRG(QGR,VCONG,RHO,RHOFAC)
      IF(RAIN)   VTR = VTRR(QRR,VCONR,RHO,RHOFAC)
     


!       IF(RAIN) print *, RAIN ALLOWED, QMIX(6)
!
!
!     *********************************************
!     ****                                     ****
!     ****   PROCESSES CALLED INDEPENDENT OF   ****
!     ****   TEMPERATURE ARE:                  ****
!     ****                                     ****
!     ****        GACW    SACW    RAUT         ****
!     ****        RACW    RACS    GACS         ****
!     ****                REVP                 ****
!     ****                                     ****
!     *********************************************
!     *********************************************
!
      IF(.NOT. LIQ)  GO TO 150
      IF(.NOT. GRAUP)  GO TO 110
!
!     PGACW
!
      QLK_eff = QLK * frac_mass_collide_w
      PGACW = HUCM_TURB_GACW*AMIN1(ACR2(QLK_eff,QGR,CGACW,0.875,RHOFAC,RHO),SLMAX)
!
  110 IF(.NOT. SNO)  GO TO 120
!
!     PSACW
!
      QLK_eff = QLK * frac_mass_collide_w
      PSACW = HUCM_TURB_SACW*AMIN1(ACR1(QLK_eff,QSR,CSACW,0.8125,RHOFAC,RHO),SLMAX)
!

!  120 IF(QLK .LE. QL0)  GO TO 130

120 IF(QLR - QL0RHO .LE. 0.)  GO TO 130

!     PRAUT

      QEXRHO = QLR - QL0RHO
	if(ADVANCED_AUTOCONV_RAIN == 0) then
      PRAUT = AMIN1(RAUT(CRAUT,QEXRHO, RHO,QNW*RHO/1.E6, DD),SLMAX)
	else
      PRAUT = AMIN1(RAUT(CRAUT,QEXRHO, RHO,QNW*RHO/1.E6, DIS),SLMAX)
	endif
!       PRAUT = AMIN1(RAUT(CRAUT,QEXRHO, RHO,N1, DD),SLMAX)


!      if(  PRAUT .GT. 0.) then
!	print *, PRAUT = ,   PRAUT , CRAUT,QEXRHO, RHO, SLMAX
!      endif

  130 IF(.NOT. RAIN)  GO TO 200

!     PRACW
      QLK_eff = QLK * frac_mass_collide_w
      PRACW = AMIN1(ACR1(QLK_eff,QRR,CRACW,0.95,RHOFAC,RHO),SLMAX)

  150 IF(.NOT. RAIN)  GO TO 200

  170 IF(.NOT. SNO)  GO TO 290

!     PRACS

      PRACS = AMIN1(ACR3(VTR,VTS,QSR,QRR,CRACS,C3RACS,RHOFAC,RHO),SSMAX)

      GO TO 210

  200 IF(.NOT. SNO)  GO TO 290

  210 IF(.NOT. GRAUP)  GO TO 290

!     PGACS
!
! fu ********************************************************************
! fu   CTGACS = CGACS*EXPT1
      CTGACS = CGACS*0.1
! fu ********************************************************************
      PGACS = AMIN1(ACR3(VTG,VTS,QSR,QGR,CTGACS,C3GACS,RHOFAC,RHO),SSMAX)
!
  290 IF(QRK .EQ. 0.0 .OR. ICE .OR. DQS .GE. 0.0)  GO TO 300
!
!     PREVP
!
      SW = QVK/QSW
      PREVP = 0.
!
  300 IF(TC .LT. 0.0)  GO TO 400
!
!     ***********************************
!     ***********************************
!     ****                           ****
!     ****     TC >= 0 PROCESSES     ****
!     ****                           ****
!     ****    SMLT   WACS   GMLT     ****
!     ****                           ****
!     ***********************************
!     ***********************************
!
!sk ********************************************************************
      IF(QSK .EQ. 0.0) GO TO 305

!     PSMLT

     IF ( SNO .AND. RAIN ) THEN
!sk
!sk   PSACR CALLED FOR SMLT HERE    (T < 0 PROCESS)
!sk
      PSACR = AMIN1(ACR3(VTS,VTR,QRR,QSR,CSACR,C3SACR,RHOFAC,RHO),SRMAX)
      END IF
!sk
!sk   PSMLT (Follow Lin et al. 1983)
!sk
      PSMLT = AMIN1(SMLT(TC,DQS0,QSR,PSACW,PSACR,CSMLT,RHO, RHOFAC),SSMAX)
      PSMLT = AMAX1(PSMLT,0.0)
      PSACR = 0. ! Not used except in PSMLT
!sk ********************************************************************
!
  305 IF(.NOT. LIQ .OR. .NOT. SNO)  GO TO 310
!
!     PWACS
!
!sk ******************************************************************
      PWACS = 0.0
!sk   PWACS = AMIN1(ACR1(QLK,QSR,CWACS,1.5625,RHOFAC,RHO),SSMAX)
!sk ******************************************************************
!
!
  310 IF(.NOT. GRAUP .OR. .NOT. RAIN)  GO TO 320
!
!     GACR CALLED FOR GMLT HERE
!
      PGACR = AMIN1(ACR3(VTG,VTR,QRR,QGR,CGACR,C3GACR,RHOFAC,RHO),SRMAX)
!
!     GMLT:  GACW HAS ALREADY BEEN CALLED IF APPROPRIATE
!            GUARD AGAINST NEGATIVE VALUES AT TEMP CLOSE TO 0 DEG C
!
  320 IF(QGK .EQ. 0.0)  GO TO 330
!
      PGMLT = AMIN1(GMLT(TC,DQS0,QGR,PGACW,PGACR,CGMLT, RHO),SGMAX)
      PGMLT = AMAX1(PGMLT,0.0)
!sk ********************************************************************
      PGACR = 0. ! Not used except in PGMLT
!sk ********************************************************************
!
!lin *******************************************************************
      PGACS = AMIN1(10.*PGACS,SSMAX)
!lin *******************************************************************
!     **************************************
!     ****                              ****
!     ****     ADD SOURCES AND SINKS    ****
!     ****             T>=0             ****
!     ****                              ****
!     **************************************
!
  330 CONTINUE

!
 IF(rescale) THEN
!
! combined cloud water depletions
	sinks = PRAUT + PSACW + PRACW + PGACW
      sources = 0.
	if(sources .lt. 0.) stop 701
	if(sinks .lt. 0.) stop 702
      IF (sinks .GT. (1. - epsi)*QLK/dt + sources) THEN
         factor =  ((1. - epsi)*QLK/dt + sources)/sinks
         PRAUT = factor * PRAUT
         PSACW = factor * PSACW
         PRACW = factor * PRACW
         PGACW = factor * PGACW 
      ENDIF

! combined snow processes
     sinks = PSMLT + PGACS + PRACS   ! PWACS is zero!!
      sources = 0.
	if(sources .lt. 0.) stop 703
	if(sinks .lt. 0.) stop 704
	if(PWACS .NE. 0.) stop 700
      IF (sinks .GT. (1. - epsi)*QSK/dt + sources) THEN       
         factor =  ((1. - epsi)*QSK/dt + sources)/sinks
         PSMLT  = factor * PSMLT 
         PGACS  = factor * PGACS
         PRACS  = factor * PRACS
      ENDIF

! combined graupel processes
      sources = PGACS
      sinks = PGMLT
	if(sources .lt. 0.) stop 705
	if(sinks .lt. 0.) stop 706
     IF ( sinks .GT. (1. - epsi)*QGK/dt + sources) THEN             
         factor =  ((1. - epsi)*QGK/dt + sources)/sinks
         PGMLT = factor * PGMLT
      ENDIF

! rain processes
      sinks = PREVP 
      sources = PRAUT+PRACW+PRACS+PSACW+PGACW+PSMLT+PGMLT 
	if(sources .lt. 0.) stop 707
	if(sinks .lt. 0.) stop 708
    IF ( sinks .GT. (1. - epsi)*QRK/dt + sources) THEN             
         factor =  ((1. - epsi)*QRK/dt + sources)/sinks
          PREVP = factor * PREVP
      ENDIF
  ENDIF
 

      TSSV = PREVP
      TSSL = -(PRAUT+PRACW+PSACW+PGACW)
	if(QLK > 0.) then
		NRAUT = PRAUT/(c_w*(D_W_MAX**3.))
!		NRAUT = 0.                 
                QLK_eff = QLK * frac_mass_collide_w
                QNW_eff = QNW * frac_num_collide_w
		if(QLK_eff > 0.) then
			NRACW = QNW_eff * (PRACW/QLK_eff)
			NSACW = QNW_eff * (PSACW/QLK_eff)
			NGACW = QNW_eff * (PGACW/QLK_eff)
		else
			NRACW = 0.
			NSACW = 0.
			NGACW = 0.
		endif
		TSSNW = -(NRAUT+NRACW+NSACW+NGACW)
	else
		TSSNW = 0.
	endif

      TSSS = -(PGACS+PRACS+PSMLT)
      TSSG = PGACS-PGMLT
      TSSR = PRAUT+PRACW+PRACS+PSACW+PGACW+PSMLT+PGMLT-PREVP
      TTMP = (-HLTF*(PRACS+PGMLT+PSMLT)-HLTC*PREVP)/CP

! PREVP > 0. for evaporation

	GAMMA_K(2) = GAMMA_K(2) + PREVP
	GAMMA_K(6) = GAMMA_K(6) + PRACS
	GAMMA_K(6) = GAMMA_K(6) + PGMLT
	GAMMA_K(6) = GAMMA_K(6) + PSMLT
	
	if(PWACS .ne. 0.) stop 705
	if(PREVP .lt. 0.) stop 710
	if(PRACS .lt. 0.) stop 711
	if(PGMLT .lt. 0.) stop 712
	if(PSMLT .lt. 0.) stop 713

!
! HLTF = Latent heat of fusion J / kg 
! PRACS = mixing ratio tendency (/sec)
! CP = specific heat at constant pressure (J/ kg/ K) = 1004.5
! =>   TTMP ~ K / sec ie.  TTMP is the tendency of potential temperature
!

!
!     WRITE(6,111) TSSV,TSSL,TSSI,TSSS,TSSG,TSSR,TTMP
  111 FORMAT('0',30X,'T>=0  '//10X,'TSSV=',E16.8/ &
     ,10X,'TSSL=',E16.8/10X,'TSSI=',E16.8/10X,'TSSS=',E16.8, &
     /10X,'TSSG=',E16.8/10X,'TSSR=',E16.8/10X,'TTMP=',E16.8) 
      TEM = 0.0
      DO K=1,6
 1112 TEM = TEM+TSS(K)
      ENDDO
!     WRITE(6,1113)  TEM
 1113 FORMAT('0',10X,'TOTAL SOURCE AND SINK=',E16.8)
!
      GO TO 1000
!
!     ************************************
!     ************************************
!     ****                            ****
!     ****      TC < 0 PROCESSES      ****
!     ****                            ****
!     ****   IDW    RIME_CI_TO_SNOW    RACI    ****
!     ****   GACI      SACI    SAUT   ****
!     ****   GFR       SACR    GACR   ****
!     ****   SSUB      GAUT    GSUB   ****
!     ****   WORD                     ****
!     ****   NOTE: IACR IS NOT USED   ****
!     ****                            ****
!     ************************************
!     ************************************
!
  400 IF(.NOT. LIQ)  GO TO 410
!
!     PIDW
!
! fu *******************************************************************
!      PIDW = AMIN1(REALIDW(TC,qik,RHO),SLMAX)

!
!	VTJP(2004) :  Follow Ferriers (1994) advice on scrapping B-F paramtrisation, and using
!	explicit vapour deposition instead
!

	PIDW = 0.

! fu   PIDW = AMIN1(IDW(TC),SLMAX)
! fu *******************************************************************
!
  410 IF(.NOT. ICE)  GO TO 450
!
!     BERGERON PROCESSES -- PSFW AND PSFI
!
! fu *******************************************************************

!
! VTJP (2004)  Follow Ferrier (1994) who recommends scrapping the B-F parametrisation when 
!	vapour growth is explicitly represented
!
	
      CALL RIME_CI_TO_SNOW(TC,QLK,QIK,QLR,PSFW, &
        C1BRG,C2BRG,RMI50,RMI40)
	PSFI = 0.


! fu   CALL RIME_CI_TO_SNOW(TC,QIK,QLR,PSFW,PSFI)
! fu *******************************************************************
      PSFW = AMIN1(PSFW,SLMAX)
!     PSFI = AMIN1(PSFI,SIMAX)


!
      IF(.NOT. RAIN)  GO TO 420
!
!     PRACI
!

      	QIK_eff = QIK * frac_mass_collide_i
      	QNI_eff = QNI * frac_num_collide_i
    PRACI = AMIN1(ACR1(QIK_eff,QRR,CRACI,0.95,RHOFAC,RHO),SIMAX)

      PIACR = AMIN1(ACR1(QNI_eff,QRR,CIACR,1.7,RHOFAC,RHO),SRMAX)

!
  420 IF(.NOT. GRAUP)  GO TO 430
!
!     PGACI
!
     	QIK_eff = QIK * frac_mass_collide_i
      PGACI = AMIN1(ACR2(QIK_eff,QGR,CGACI,0.875,RHOFAC,RHO),SIMAX)
!
  430 IF(.NOT. SNO)  GO TO 440
!
!     PSACI
!
! fu ******************************************************************
! fu   QIKT = QIK*EXPT2
       QIK_eff = QIK * frac_mass_collide_i
     QIKT = QIK_eff*0.1
! fu ******************************************************************
      PSACI = HUCM_TURB_SACI*AMIN1(ACR1(QIKT,QSR,CSACI,0.8125,RHOFAC,RHO),SIMAX)
!
  440 IF(QIK .LE. QI0)  GO TO 450
!
!     PSAUT
!
      C = 0.1*EXPT2
!      PSAUT = AMIN1(AUT(C,QIK,QI0),SIMAX)

! Follow Ferrier (1994) and Swann



	lambda_i = LAMBDA_FAC_i*c_i*QNI/QIK
	lambda_i = lambda_i**(1./d_i)


! Takes account of vapour deposition AND aggregation/autoconversion 

	if(lambda_i < lambda_i0) then
		PSAUT = 1. - ((lambda_i/lambda_i0)**3.)
		PSAUT = PSAUT * QIK/dt
	else
		PSAUT = 0.
	endif


!
  450 IF(.NOT. RAIN)  GO TO 470
!
!     PGFR
!
      PGFR = AMIN1(GFR(TC,QRR,CGFR, RHO),SRMAX)
! VTJP (30/3/05)		start       
      IF(TC < -35.) PGFR = SRMAX
! VTJP (30/3/05)		end


!
      IF(.NOT. SNO)  GO TO 460
!
!     PSACR
!
      PSACR = AMIN1(ACR3(VTS,VTR,QRR,QSR,CSACR,C3SACR,RHOFAC,RHO),SRMAX)
!
  460 IF(.NOT. GRAUP)  GO TO 470
!
!     PGACR
!
      PGACR = AMIN1(ACR3(VTG,VTR,QRR,QGR,CGACR,C3GACR,RHOFAC,RHO),SRMAX)
!
  470 ESI = REAL(ESATI_int(DBLE(ttt), ESIT, temp_K))
      QSI = EPS*ESI/(PRESR-ESI)
      SI = QVK/QSI
!
      IF(QSK .EQ. 0.0)  GO TO 480
!
!
	PSSUB = 0.
!
      IF(QSK .LE. 0.)  GO TO 480
!
!     PGAUT
!
      C = 1.E-3*EXPT1
!      PGAUT = AUT(C,QSK,QS0)

	if(QSR > 5.e-4 .and. TC < 0. .and. SI > 1. .and. PSACW > 0.) then
		PSDEP = AMIN1(SSUB(SI,TSQ,QSI,QSR,CSSUB, RHO, RHOFAC),SVMAX)
        	PGAUT = 0.5*AMAX1(0., PSACW - PSDEP - PSACI)
		if(PGAUT > SSMAX) PGAUT = SSMAX

		if(PSDEP < 0.) stop 43241
	else
	       PGAUT = 0.
	endif

  480 IF(QGK .EQ. 0.0)  GO TO 520
      IF(QIK .NE. 0.0 .OR. QLK .NE. 0.0 .OR. QVK .GE. QSI)  GO TO 500
!
!
	PGSUB = 0.
!!
!     PGDRY OR PGWET
!
  500 PGDRY = PGACW + PGACI + PGACR + PGACS
      IF(.NOT.LIQ .AND. .NOT.RAIN)  GO TO 510
!
!lin *******************************************************************
      PGWET = GWET(DQS0,TC,QGR,PGACI,PGACS,CGWET,SIMAX,SSMAX, &
	HLTF,CH2O,CICE,RHO)

	
!lin *******************************************************************
!
!
      IF(PGDRY .LE. PGWET)  GO TO 510
!
!lin  PGWET INVOLVES REDEFINITION OF PGACI, PGACS (INCREASE OF
!     COLLECTION EFFICIENCY BY FACTOR OF 10), AND RECALCULATION
!     OF PGACR (WHICH CAN NOW BECOME NEGATIVE DUE TO SHEDDING
!     OF CLOUD LIQUID WATER WHICH IS CONVERTED INTO RAINWATER)
!
!lin *******************************************************************
!lin  PGACI = PGACI * 10.0
      PGACI = AMIN1(PGACI * 10.0,SIMAX)
      PGACS = AMIN1(10.*PGACS,SSMAX)
!lin  PGACR = PGWET - PGACW - PGACI - PGACS
      PGACR = AMIN1(PGWET - PGACW - PGACI - PGACS,SRMAX)
!lin *******************************************************************
!
      PGWORD = PGWET
      GO TO 520
!
  510 PGWORD = PGDRY
!
!
!     **************************************
!     ****                              ****
!     ****     ADD SOURCES AND SINKS    ****
!     ****              T<0             ****
!     ****                              ****
!     **************************************
!
  520 CONTINUE



	if(QLK > 0.) then
		NIDW = 0.
		NSFW = QNW *PSFW/QLK
		NRAUT = PRAUT/(c_w*(D_W_MAX**3.))
!		NRAUT = 0. 
                QLK_eff = QLK * frac_mass_collide_w
                QNW_eff = QNW * frac_num_collide_w
		if(QLK_eff > 0.) then
			NRACW = QNW_eff * (PRACW/QLK_eff)
			NSACW = QNW_eff * (PSACW/QLK_eff)
			NGACW = QNW_eff * (PGACW/QLK_eff)
		else
			NRACW = 0.
			NSACW = 0.
			NGACW = 0.
		endif

		TSSNW = -(NSACW+NSFW+NRAUT+NRACW+NGACW+NIDW )
	else
		NIDW = 0.
		NSFW = 0.
		NRAUT = 0.
		NRACW = 0.
		NSACW = 0.
 		NGACW = 0.
		TSSNW = 0.
	endif

	if(QLK > 0. .and. ICE_NUMBER == 1) then

		if(ADVANCED_HM_FORMULA == 1) then
			NWLONW_eff = gammq_frac(1.+p_w, 25.E-6*lambda_w)/frac_num_collide_w  
			if(NWLONW_eff > 1.) then
				print *, 'D_XW > 25 um ?', NWLONW_eff, QNW, QLK
				stop 31456
			endif 
			if(NWLONW_eff < 0.) stop 31459

			NIHMS = F_HM(TC) * NSACW * NWLONW_eff/200. 
			NIHMG = F_HM(TC) * NGACW * NWLONW_eff/200. 
		else
			NIHMS = F_HM(TC) * PSACW * 350.e6 
			NIHMG = F_HM(TC) * PGACW * 350.e6 
		endif
		PIHMS = NIHMS * 920.* (4./3.)*3.14159265* (2.5e-6**3.)
		PIHMG = NIHMG * 920.* (4./3.)*3.14159265* (2.5e-6**3.)
	else
		NIHMS = 0.
		PIHMS = 0.
		NIHMG = 0.
		PIHMG = 0.
	endif

 	if(QIK > 0.) then
		NIDW = 0.
		NSFI = 0.
		NSAUT = 0.
     		QIK_eff = QIK * frac_mass_collide_i
      		QNI_eff = QNI * frac_num_collide_i
		if(QIK_eff > 0.) then
			NSACI = QNI_eff * PSACI/QIK_eff
			NRACI = QNI_eff * PRACI/QIK_eff
			NGACI = QNI_eff * PGACI/QIK_eff
		else
			NSACI = 0.
			NRACI = 0.
			NGACI = 0.
		endif
		TSSNI = -(NSAUT+NSACI+NSFI+NRACI+NGACI)+NIDW + NIHMS + NIHMG
	else

		TSSNI = 0. + NIHMS + NIHMG
	endif


 IF(rescale) THEN

! rescale conversion rates in order to conserve water mass as done 
! in Shuhua Chens version
!

delt1 = 1.
delt2 = 1.

IF(QRK .LT. 1.E-4 .AND. QSK .LT. 1.E-4) delt1 = 0. 
IF(QRK .GE. 1.E-4) delt2 = 0.

! combined water vapour depletions
!PGSUB is always -ve => SUBLIMATION CAN ONLY BE A SOURCE OF VAPOUR
 
      IF (PSSUB .GT. 0) THEN    
       sinks = PSSUB
	sources = -PGSUB + PREVP
	if(sources .lt. 0.) stop 800

       IF (sinks .GT. (1. - epsi)*QVK/dt + sources ) THEN     
           factor = ((1. - epsi)*QVK/dt + sources)/sinks
            PSSUB =  factor * PSSUB
       ENDIF
      ENDIF

! combined cloud water depletions 
	sources = 0.
       sinks = PRAUT + PSACW + PSFW + PRACW + PGACW + PIDW
       If (sinks .GT. (1. - epsi)*QLK/dt + sources) THEN
          factor = ((1. - epsi)*QLK/dt + sources)/sinks
!!$          IF(-PRAUT .GT. 0) WRITE(6,*) "!! Prod PRAUT"
!!$          IF(-PSACW .GT. 0) WRITE(6,*) "!! Prod PSACW"
!!$          IF(-PSFW .GT. 0) WRITE(6,*) "!! Prod PSFW"  
!!$          IF(-PRACW .GT. 0) WRITE(6,*) "!! Prod PRACW"
!!$          IF(-PGACW .GT. 0) WRITE(6,*) "!! Prod PGACW"
!!$          IF(-PIDW .GT. 0) WRITE(6,*) "!! Prod PIDW"
          PRAUT = factor * PRAUT
          PSACW = factor * PSACW
          PSFW = factor * PSFW
          PRACW = factor * PRACW 
          PGACW = factor * PGACW
          PIDW = factor * PIDW
       ENDIF


!     TSSI = -(PSAUT+PSACI+PSFI+PRACI+PGACI)+PIDW + PIHMS + PIHMG
! combined cloud ice depletions
	if( PIDW /= 0.) stop 54321

	sources = PIDW
       sinks = PSAUT + PSACI + PSFI + PRACI + PGACI
	if(sources .lt. 0.) stop 801
	if(sinks .lt. 0.) stop 802

       IF (sinks .GT. (1. - epsi)*QIK/dt + sources) THEN
          factor = ((1. - epsi)*QIK/dt + sources)/sinks
!!$          IF(-PSAUT .GT. 0) WRITE(6,*) "!! Prod PSAUT"
!!$          IF(-PSACI .GT. 0) WRITE(6,*) "!! Prod PSACI"
!!$          IF(-PSFI .GT. 0) WRITE(6,*) "!! Prod PSFI"
!!$          IF(-PRACI .GT. 0) WRITE(6,*) "!! Prod PRACI"
!!$          IF(-PGACI .GT. 0) WRITE(6,*) "!! Prod PGACI"
          PSAUT = factor * PSAUT
          PSACI = factor * PSACI
          PSFI  = factor * PSFI
          PRACI = factor * PRACI
          PGACI = factor * PGACI 
       ENDIF


! combined rain processes
	if(PGACR .ge. 0.) then
		sources = PRAUT+PRACW
		sinks = PSACR + PREVP + PGFR + PGACR + PIACR
	else
		sources = PRAUT+PRACW - PGACR 
		sinks = PSACR + PREVP + PGFR + PIACR

	endif
	if(sources .lt. 0.) stop 803
	if(sinks .lt. 0.) then
		print *, PSACR, PREVP, PGFR, PGACR 

		stop 804
	endif
       IF (sinks .GT. (1. - epsi)*QRK/dt + sources) THEN
          factor = ((1. - epsi)*QRK/dt + sources)/sinks
!!$          IF(-PSACR .GT. 0) WRITE(6,*) "!! Prod PSACR"
!!$          IF(-PREVP .GT. 0) WRITE(6,*) "!! Prod PREVP"
!!$          IF(-PGFR .GT. 0) WRITE(6,*) "!! Prod PGFR"
!!$          IF(-PGACR .GT. 0) WRITE(6,*) "!! Prod PGACR"
          PSACR = factor * PSACR          
          PREVP = factor * PREVP          
          PGFR = factor * PGFR     
          PIACR = factor * PIACR     
  	  if(PGACR .ge. 0.) then
          	PGACR = factor * PGACR
	  endif
        ENDIF


! combined snow processes
  
	if(PSSUB .ge. 0.) then
		sources = PSAUT+PSACI+PSACW+PSFW+PSFI+PSSUB + (1. - delt1) * PSACR + delt2 * (PRACI+PIACR)
		sinks = PRACS*delt1 + PGACS + PGAUT + PIHMS
	else
		sources = PSAUT+PSACI+PSACW+PSFW+PSFI + (1. - delt1) * PSACR + delt2 * (PRACI+PIACR)
		sinks = PRACS*delt1 + PGACS + PGAUT - PSSUB + PIHMS
	endif
	if(sources .lt. 0.) stop 805
	if(sinks .lt. 0.) stop 806
	IF (sinks .GT. (1. - epsi)*QSK/dt + sources) THEN
		factor = ((1. - epsi)*QSK/dt + sources)/sinks
        	if(PSSUB .LT. 0.) THEN
	        	PSSUB = factor * PSSUB 
          	endif          
		PRACS = factor * PRACS        
		PGAUT = factor * PGAUT           
		PGACS = factor * PGACS 
		PIHMS = factor * PIHMS 
	ENDIF


  ! adjust PGWORD
      IF( PGWORD .NE. pgacw+pgaci+pgacr+pgacs ) THEN
           PGWORD =  pgacw+pgaci+pgacr+pgacs
      ENDIF


! combined graupel processes

	if(PGWORD .ge. 0.) then
		sources = PGAUT + PGFR + delt1 * PSACR + delt1 * PRACS + (1. - delt2) * (PRACI+PIACR) + PGWORD
        	sinks =  -PGSUB  + PIHMG
	else
		sources = PGAUT + PGFR + delt1 * PSACR + delt1 * PRACS + (1. - delt2) * (PRACI+PIACR)
        	sinks =  -PGSUB - PGWORD + PIHMG
	endif

  	if(sources .lt. 0.) stop 807
	if(sinks .lt. 0.) stop 808
         
       IF(sinks .GT. (1. - epsi)*QGK/dt + sources) THEN
          factor = ((1. - epsi)*QGK/dt + sources)/sinks
          PGSUB = factor * PGSUB
	if(PGWORD .lt. 0.) then
		PGWORD = PGWORD *factor
		pgacw = pgacw*factor
		pgaci = pgaci*factor
		pgacr = pgacr*factor
		pgacs = pgacs*factor
		PIHMG = factor * PIHMG 
	ENDIF

       ENDIF 


  ENDIF



      TSSV = -(PSSUB+PGSUB)+PREVP
      TSSL = -(PSACW+PSFW+PRAUT+PRACW+PGACW+PIDW)
      TSSI = -(PSAUT+PSACI+PSFI+PRACI+PGACI)+PIDW + PIHMS + PIHMG

 

      TSSS = PSAUT+PSACI+PSACW+PSFW+PSFI+PSSUB-(PGACS+PGAUT + PIHMS)

      IF(QRK .LT. 1.E-4 .AND. QSK .LT. 1.E-4)  GO TO 530
!	delt1 = 1
      TSSS = TSSS-PRACS
      TSSG = PSACR+PRACS
      GO TO 540
  530 TSSS = TSSS+PSACR
! delt1 = 0

  540 IF(QRK .GE. 1.E-4)  GO TO 550
      TSSS = TSSS+PRACI+PIACR
! delt2 = 1
      GO TO 560
  550 TSSG = PRACI+PIACR+TSSG
!delt2 = 0
  560 TSSG = TSSG+PGSUB+PGAUT+PGFR+PGWORD - PIHMG
      TSSR = PRAUT+PRACW-(PSACR+PGACR+PREVP+PGFR+PIACR)

!       if(PRAUT .gt. 0.) then
!      print *, (2) TSSR = , TSSR, TSS(6), TSS(6)
!       endif
      TTMP = (HLTF*(PSACW+PGACW+PSACR+PGACR+PGFR+PSFW+PIDW+PIACR) &
      -HLTC*PREVP+HLTS*(PSSUB+PGSUB))/CP

! PREVP > 0. for evaporation
	if(PREVP .ge. 0.) then
		GAMMA_K(2) = GAMMA_K(2) + ABS(PREVP)
	else
		stop 610
	endif
!PSSUB +ve => vapour deposition
	if(PSSUB .ge. 0.)  then
		GAMMA_K(3) = GAMMA_K(3) + ABS(PSSUB)
	else 
		GAMMA_K(4) = GAMMA_K(4) + ABS(PSSUB)
	endif

!PGSUB +ve => vapour deposition
	if(PGSUB .ge. 0.)  then
		GAMMA_K(3) = GAMMA_K(3) + ABS(PGSUB)
	else 
		GAMMA_K(4) = GAMMA_K(4) + ABS(PGSUB)
	endif

! VTJP:
!Of the following processes, only PGACR is sometimes -ve,
! (in which case it is representing an inefficient form of freezing 
! rather than melting)
!
	GAMMA_K(5) = GAMMA_K(5) + PSACW
	GAMMA_K(5) = GAMMA_K(5) + PGACW
	GAMMA_K(5) = GAMMA_K(5) + PSACR
	GAMMA_K(5) = GAMMA_K(5) + PGACR
	GAMMA_K(5) = GAMMA_K(5) + PGFR
	GAMMA_K(5) = GAMMA_K(5) + PSFW
	GAMMA_K(5) = GAMMA_K(5) + PIDW

	if(PSACW .lt. 0.) stop 700	
	if(PGACW .lt. 0.) stop 701	
	if(PSACR .lt. 0.) stop 702	
	if(PGFR .lt. 0.) stop 703	
	if(PSFW .lt. 0.) stop 704	
	if(PIDW .lt. 0.) stop 705	
	
!
!     WRITE(6,666) TSSV,TSSL,TSSI,TSSS,TSSG,TSSR,TTMP
  666 FORMAT('1',30X,'T<0  '//10X,'TSSV=',E16.8/ &
       ,10X,'TSSL=',E16.8/10X,'TSSI=',E16.8/10X,'TSSS=',E16.8, &
       /10X,'TSSG=',E16.8/10X,'TSSR=',E16.8/10X,'TTMP=',E16.8)
      TEM = 0.0
      DO K=1,6
 6661 TEM = TEM+TSS(K)
      ENDDO
!     WRITE(6,1113)  TEM
!
!
 1000 CONTINUE
!
!     PRINT VALUES FOR EACH PROCESS
!
!     WRITE(6,2000)ALL TEMP,PGACW,PGACW,PSACW,PSACW,PRAUT,PRAUT,
!    1   PRACW,PRACW,PRACS,PRACS,PGACS,PGACS,PREVP,PREVP
!     WRITE(6,2000)TC > 0,PSMLT,PSMLT,PWACS,PWACS,PGMLT,PGMLT
!     WRITE(6,2000) TC < 0,PSFW,PSFW,PSFI,PSFI,PRACI,PRACI,
!    1   PIACR,PIACR,PGACI,PGACI,PSACI,PSACI,PSAUT,PSAUT,
!    2   PGFR,PGFR,PSACR,PSACR,PGACR,PGACR,PSSUB,PSSUB,
!    3   PGAUT,PGAUT,PGSUB,PGSUB,PIDW,PIDW
!     WRITE(6,2000)WET-DRY,PGWET,PGWET,PGDRY,PGDRY,PGWORD,PGWORD
 
     DO K=1,6
		QMIXY(K) = QMIX(K)
		TSSY(K) = TSS(K)
      ENDDO
      TSSY(7) = TSS(7)
  
	DO I = 1,3
        DO K = 1,4
		ACCOY(I,K) = ACCO(I,K)
	ENDDO
	ENDDO
     
       RETURN   
      END SUBROUTINE IMICRO




      REAL FUNCTION REALIDW(TC,qik,RHO )
       IMPLICIT NONE
      REAL, INTENT(IN):: TC, qik,  RHO

      REAL :: fnc, fmi, A1, A2, TC1, A1T, A2T, RMINUC

! fu   REAL FUNCTION REALIDW(TC,)
! fu *******************************************************************
  

! fu *******************************************************************
      DIMENSION A1T(0:31),A2T(0:31)
! fu   DIMENSION A1T(31),A2T(31)
! fu *******************************************************************
!
!     NOTE -- A2T IS IDENTICAL TO A2T IN SUBROUTINE RIME_CI_TO_SNOW, BUT A1T
!     IS MULTIPLIED BY A FACTOR OF 1.E-5.
!
! fu *******************************************************************
      DATA A1T/.0,.7939E-12,.7841E-11,.3369E-10,0.4336E-10,0.5285E-10, &
! fu   DATA A1T/0.7939E-12,0.7841E-11,0.3369E-10,0.4336E-10,0.5285E-10,&
! fu *******************************************************************
              0.3728E-10,0.1852E-10,0.2991E-11,0.4248E-11,0.7434E-11,&
! fu *******************************************************************
              0.1812E-10,0.4394E-10,0.9145E-10,0.1725E-9 ,0.3348E-9 ,&
! fu  2         0.1812E-10,0.4394E-10,0.9145E-10,0.1725E-11,0.3348E-9 ,&
! fu *******************************************************************
              0.1725E-9 ,0.9175E-10,0.4412E-10,0.2252E-10,0.9115E-11,&
              0.4876E-11,0.3473E-11,0.4758E-11,0.6306E-11,0.8573E-11,&
              0.7868E-11,0.7192E-11,0.6513E-11,0.5956E-11,0.5333E-11,&
              0.4834E-11/
! fu *******************************************************************
      DATA A2T/0.0,0.4006,0.4831,0.5320,0.5307,0.5319, &
! fu   DATA A2T/0.4006,0.4831,0.5320,0.5307,0.5319, &
! fu *******************************************************************
              0.5249,0.4888,0.3894,0.4047,0.4318, &
              0.4771,0.5183,0.5463,0.5651,0.5813, &
              0.5655,0.5478,0.5203,0.4906,0.4447, &
              0.4126,0.3960,0.4149,0.4320,0.4506, &
              0.4483,0.4460,0.4433,0.4413,0.4382, &
              0.4361/
      DATA RMINUC/1.05E-15/
!
      TC1 = AMAX1(TC,-30.0)
! fu ********************************************************************
      A1  = (A1T(-INT(TC1))-A1T(-INT(TC1)+1))*(TC1-INT(TC1)+1.0)+ &
           A1T(-INT(TC1)+1)
      A2  = (A2T(-INT(TC1))-A2T(-INT(TC1)+1))*(TC1-INT(TC1)+1.0)+ &
           A2T(-INT(TC1)+1)
! fu   A1  = A1T(-INT(TC1)+1)
! fu   A2  = A2T(-INT(TC1)+1)
! fu ********************************************************************
!
! fu ********************************************************************
      fnc = 0.01 * exp ( - 0.6 * tc )
      fmi = rho * qik / fnc * 1000.0
      REALIDW = EXP(-0.6*TC)*A1*fmi**A2/RHO
! fu   REALIDW = EXP(-0.6*TC)*A1*RMINUC**A2/RHO
! fu ********************************************************************
      RETURN
      END FUNCTION REALIDW

REAL FUNCTION F_HM(temperature)
    IMPLICIT NONE

REAL, INTENT(IN):: temperature
       

        if (temperature <= -8. .or. temperature > -3.) F_HM = 0.

        if (temperature <= -5.5 .and. temperature > -8.) F_HM =  ( 8. - abs(temperature)) / 2.5 

        if (temperature > -5.5 .and. temperature <= -3.) F_HM = ( abs(temperature) - 3.) / 2.5
    

	if(    F_HM > 1.) stop "154986"
	if(    F_HM < 0.) stop "144986"
      RETURN
      END FUNCTION F_HM



      SUBROUTINE BOMB(I,T)
       IMPLICIT NONE
        INTEGER, INTENT(IN) :: I
        REAL, INTENT(IN) :: T

!      COMMON/ERRDAT/NERR,IERR,JERR,PRESS,TO,QVO,QLO,QIO
      WRITE(6,*) 'BOMB: I, T',I,T
!      WRITE(6,*) BOMB: J,K,PRESS,T0,QV0,QL0,QI0,IERR,JERR,PRESS,TO,QVO,QLO,QIO
!      CALL OUTPUT ( 0 )
!      CALL TSPUT
      STOP 911
      END SUBROUTINE BOMB




      SUBROUTINE SETUPM ( DT, SFCRHO, grav_wrf, cp_wrf, Rair, rvapor, &
		XLS, XLV, XLF, rhowater, rhosnow, rhograupel, &
	CRACS,CSACR,CGACR,CGACS,ACCO,CSACW,CWACS, &
           CIACR,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB,CGSUB, &
!sk        CREVP,CGFR,CSMLT,CGMLT,CGWET,CRAUT, &
           CREVP,CGFR,CSMLT,CGMLT,CGWET,CRAUT, &
           QI0,QS0,QL0,ES0,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50, &
           CPLC,CPLS,CPLF,TDIF,ESW00,ESI00, &
 	PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC &
	  ,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ, &
	VCONR,VCONS,VCONG,VCONI, DTL,TDT, RHOSFC, BETA_SUPER, kaero, &
	s_pc_start, D_v, k_a, N_0S_OUT, N_0G_OUT, N_0R_OUT)

	IMPLICIT NONE

      REAL, INTENT(IN) :: DT, SFCRHO, kaero
      REAL, INTENT(IN)::  grav_wrf, cp_wrf, Rair, rvapor, XLS, XLV, XLF
       INTRINSIC ABS, SUM
	
!
!     CONSTANTS FOR MICROPHYSICS SOURCE AND SINK TERMS
!     AND TERMINAL VELOCITIES
!
      REAL, INTENT(OUT) :: CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS, &
           CIACR,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5), &
!sk        CREVP(5),CGFR(2),CSMLT(3),CGMLT(5),CGWET(2),CRAUT(2), &
           CREVP(5),CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2), &
           QI0,QS0,QL0,ES0,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50, &
           CPLC,CPLS,CPLF,TDIF,ESW00,ESI00,BETA_SUPER, D_v,k_a, s_pc_start
!
      REAL, INTENT (OUT) :: PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK, &
	   HLTS,HLTC,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ
!
      REAL, INTENT(OUT) :: VCONR,VCONS,VCONG,VCONI
!
      REAL, INTENT(OUT) :: DTL,TDT,RHOSFC, rhowater, rhosnow, rhograupel, &
		N_0S_OUT, N_0G_OUT, N_0R_OUT
!


      REAL :: GCON, CD, aa22, SCM3, PISQ, ACT, ACC, RNZR, & 
	RNZS, RNZG, RHOR, RHOS, RHOG, &
	ALIN, CLIN, GAM680, GAM625,GAM480, GAM450, GAM425, & 
	GAM380,GAM350, GAM325, GAM290, GAM275, GAM263, &
	GAMX, GAMY, GAMZ

	INTEGER :: ITC,  K, I, kr

      DIMENSION ACT(8),ACC(3)

!	DOUBLE PRECISION :: GGESW, GGESI, gammln

!
!     PHYSICAL CONSTANTS (MKS)
!
!*    DATA PIE/3.14159265/GRAV/9.81/VDIFU/2.11E-5/
!*    DATA TCOND/2.36E-2/RVAPR/4.615E2/RDRYA/2.87E2/VISK/1.259E-5/
!*    DATA HLTS/2.8336E6/HLTC/2.5E6/HLTF/3.336E5/CH2O/4.1855E3/
!*    DATA CICE/2.093E3/TICE/273.16/EPS/0.62197/TTFRZ/233.1/
!
!     LINS CONSTANTS(MKS) EXCEPT RMI50,RMI40 (CGS)
!
!rh ********************************************************************
!rh   DATA RNZR/20.E6/RNZS/3.E6/RNZG/4.E4/
!rh   DATA RHOR/1.E3/RHOS/1.E2/RHOG/0.3E3/
      DATA RNZR/8.E6/  ! Lin83
!      DATA RNZS/3.E6/  ! Lin83
       DATA RNZS/3.e7/  ! VTJP:   KWAJEX value corresponding to mean size for snow of about 500 microns (see Fig 19, Heymsfield et al 2002)
			! same value attained when McFarquhar et al. (1999) size distribution is fitted with an exponential at 0.6 - 1 mm:
			! (we choose 0.6 mm as lower threshold because this is needed for <D_s> to be at the autoconversion threshold (200 um)) 	
      DATA RNZG/4.E6/  ! RH84
!      DATA RNZG/4.E8/  ! VTJP:  CRYSTAL-FACE (11 km MSL, 18th July, 17:49:10-30 GMT, Citation)
      DATA RHOR/1.0E3/ ! Lin83
      DATA RHOS/0.1E3/ ! Lin83
      DATA RHOG/0.4E3/ ! RH84
!rh ********************************************************************
      DATA ALIN/841.99667/clin/4.836071224/
!*    DATA RMI50/4.80E-7/RMI40/2.46E-7/RI50/5.E-5/
!*    DATA QI0/1.E-3/QS0/0.6E-3/QL0/0.5E-3/
!*    DATA CRAUT/1.2E-1,1.064E-2/

      DATA ACC/5.0,2.0,0.5/
      DATA GAM263/1.456943/GAM275/1.608355/GAM290/1.827363/ &
          GAM325/2.54925/GAM350/3.323363/GAM380/4.694155/ &
          GAM425/8.285063/GAM450/11.631769/GAM480/17.837789/ &
          GAM625/184.860962/GAM680/496.604067/
!
      PIE=3.14159265

! MARITIME SETTINGS:  Covert (1992) for Equatorial Pacific Ocean
       rhowater = RHOR
       rhosnow = RHOS
       rhograupel = RHOG

! CLEAN MARITIME
!

        s_pc_start = 0.01 ! this corresponds to D_crit ~ 10 um, see Fig. 9-3, PK97

!      Caero = 200.e6
!      kaero = 0.6
!      s_pc_cut_off = 2.**(1./kaero)
!      s_pc_start = 0.03 ! Fig. 9-3, PK97


! CONTINENTAL SETTINGS
!      Caero = 3500.e6
!      kaero = 0.9
!      s_pc_cut_off = 2.





	N_0S_OUT = RNZS; N_0G_OUT = RNZG; N_0R_OUT = RNZR;







      D_v = 2.21E-5
      k_a = 2.40E-2

      if(  ABS( EXP(gammln(0.5)) - SQRT(PIE)) > 0.00001) STOP 4321

       GAMX = EXP(gammln(kaero/2.))
       GAMY = EXP(gammln(3./2.))
       GAMZ = EXP(gammln(kaero/2. + 3./2.))

      BETA_SUPER = GAMX*GAMY/GAMZ

 
      GRAV=9.81
!      GRAV = grav_wrf

      VDIFU=2.11E-5
      TCOND=2.36E-2
 
      RVAPR=4.615E2
!     RVAPR=Rvapor

      RDRYA=2.87E2
!	RDRYA = Rair

      VISK=1.259E-5
      HLTS=2.8336E6
!      HLTS=XLS
      HLTC=2.5E6
!      HLTC=XLV
      HLTF=3.336E5
!      HLTF=XLF

      CH2O=4.1855E3
      CICE=2.093E3
      TICE=273.16
      EPS=0.62197
      TTFRZ=233.1
! fu ******************************************************************
      RMI50=3.84E-6
! fu   RMI50=4.80E-7
! fu ******************************************************************
      RMI40=2.46E-7
! fu ******************************************************************
      RI50=1.E-4
! fu   RI50=5.E-5
! fu ******************************************************************
! fu ******************************************************************
! fu   QI0=1.E-3
! fu   QS0=0.6E-3
      QS0=1.E-3
      QI0=0.6E-3
! fu ******************************************************************
!      QL0=0.5E-3
      QL0=2.E-3
      CRAUT(1)=1.2E-1
      CRAUT(2)=1.064E-2
!
!     PARAMETERS PASSED FROM MODEL WHICH USES ADAMS-BASHFORTH SCHEME
!
      TDT = DT
      DTL = DT
      RHOSFC = SFCRHO
!
      CP = 3.5*RDRYA
!     CP = cp_wrf
!       print *, cp = , CP
 

      CPLC = HLTC/CP
      CPLS = HLTS/CP
      CPLF = HLTF/CP
      PISQ = PIE*PIE
      SCM3 = (VISK/VDIFU)**(1./3.)
!
!     ACR3:  FOUR LEAD CONSTANTS REQUIRED, THREE FOR EACH SUMMATION
!            FOUR SEPARATE PROCESSES:  RACS,SACR,GACR,GACS
!
      CRACS = PISQ*RNZR*RNZS*RHOS
      CSACR = PISQ*RNZR*RNZS*RHOR
      CGACR = PISQ*RNZR*RNZG*RHOR
      CGACS = PISQ*RNZG*RNZS*RHOS
!
!     ACT:  1-2:RACS(S-R); 3-4:SACR(R-S);
!           5-6:GACR(R-G); 7-8:GACS(S-G)
!
      ACT(1) = PIE * RNZS * RHOS
      ACT(2) = PIE * RNZR * RHOR
      ACT(6) = PIE * RNZG * RHOG
      ACT(3) = ACT(2)
      ACT(4) = ACT(1)
      ACT(5) = ACT(2)
      ACT(7) = ACT(1)
      ACT(8) = ACT(6)


!
      DO I=1,3
      DO K=1,4
      ACCO(I,K) = ACC(I)/(ACT(2*K-1)**((7-I)*0.25)*ACT(2*K)**(I*0.25))
      ENDDO
      ENDDO
!
!     TERMINAL VELOCITY CONSTANTS
!
      VCONR = ALIN*GAM480/(6.*ACT(2)**0.20)
      VCONS = clin*GAM425/(6.*ACT(1)**.0625)
!rh ********************************************************************
      aa22 = 40.74 * 40.74
      CD = 4. * GRAV * RHOG / 3.0 / RHOSFC / aa22
      GCON  = SQRT(4.*GRAV*RHOG/3.0/CD)
!rh   GCON  = SQRT(4.*GRAV*RHOG/1.8)
!rh ********************************************************************
      VCONG = GAM450*GCON/(6.*ACT(6)**0.125)
      VCONI = 3.29
!
!     ACR1:  SINGLE CONSTANT REQUIRED
!     FIVE SEPARATE PROCESSES:  SACW,WACS,IACR,RACI,SACI
!
      CSACW = PIE*RNZS*clin*GAM325/(4.*ACT(1)**0.8125)
      CWACS = PISQ*RHOS*RNZS*clin*GAM625/(1.0056E-10*ACT(1)**1.5625)
      CIACR = PISQ*RHOR*RNZR*ALIN*GAM680/(24.*ACT(2)**1.7)
      CRACI = PIE*RNZR*ALIN*GAM380/(4.*ACT(2)**0.95)
      CSACI = CSACW
!
!     ACR2:  SINGLE CONSTANT REQUIRED
!     CAGCI IS FOR DRY GROWTH
!
      CGACW = PIE*RNZG*GAM350*GCON/(4.*ACT(6)**0.875)
      CGACI = CGACW*0.1
!
!     RACW
!
      CRACW = CRACI
!
!     SUBL AND REVP:  FIVE CONSTANTS FOR THREE SEPARATE PROCESSES
!
      CSSUB(1) = 2.*PIE*VDIFU*TCOND*RVAPR*RNZS
      CGSUB(1) = 2.*PIE*VDIFU*TCOND*RVAPR*RNZG
      CREVP(1) = 2.*PIE*VDIFU*TCOND*RVAPR*RNZR
      CSSUB(2) = 0.78/SQRT(ACT(1))
      CGSUB(2) = 0.78/SQRT(ACT(6))
      CREVP(2) = 0.78/SQRT(ACT(2))
      CSSUB(3) = 0.31*SCM3*GAM263*SQRT(clin/VISK)/ACT(1)**0.65625
      CGSUB(3) = 0.31*SCM3*GAM275*SQRT(GCON/VISK)/ACT(6)**0.6875
      CREVP(3) = 0.31*SCM3*GAM290*SQRT(ALIN/VISK)/ACT(2)**0.725
      CSSUB(4) = TCOND*RVAPR
      CSSUB(5) = HLTS**2*VDIFU
      CGSUB(4) = CSSUB(4)
      CREVP(4) = CSSUB(4)
      CGSUB(5) = CSSUB(5)
      CREVP(5) = HLTC**2*VDIFU
!
!     GFR:  TWO CONSTANTS
!
      CGFR(1) = 20.E2*PISQ*RNZR*RHOR/ACT(2)**1.75
      CGFR(2) = 0.66
!
!sk   SMLT:  FIVE CONSTANTS ( Lin et al. 1983 )
!
!sk ********************************************************************
!sk   CSMLT(1) = 2.*PIE*TCOND*RNZS/HLTF
!sk   CSMLT(2) = CSSUB(2)
!sk   CSMLT(3) = CSSUB(3)
      CSMLT(1) = 2.*PIE*TCOND*RNZS/HLTF
      CSMLT(2) = 2.*PIE*VDIFU*RNZS*HLTC/HLTF
      CSMLT(3) = CSSUB(2)
      CSMLT(4) = CSSUB(3)
      CSMLT(5) = CH2O/HLTF
!sk ********************************************************************
!
!     GMLT:  FIVE CONSTANTS
!
      CGMLT(1) = 2.*PIE*TCOND*RNZG/HLTF
      CGMLT(2) = 2.*PIE*VDIFU*RNZG*HLTC/HLTF
!sk ********************************************************************
      CGMLT(3) = CGSUB(2)
      CGMLT(4) = CGSUB(3)
!sk   CGMLT(3) = 1.6/SQRT(ACT(6))
!sk   CGMLT(4) = CGSUB(3)/SCM3
!sk ********************************************************************
      CGMLT(5) = CH2O/HLTF
!
!     GWET:  TWO CONSTANTS PLUS CALC. OF SAT. VAPOR PRESSURE AT TC=0.0
!     VENTILATION COEFFICIENT OF 10 INCLUDED
!
!sk ********************************************************************
!sk   CGWET(1) = 20.*PIE*RNZG/SQRT(ACT(6))*HLTC*VDIFU
!sk   CGWET(2) = 20.*PIE*RNZG/SQRT(ACT(6))*TCOND
      CGWET(1) = CGMLT(1)
      CGWET(2) = CGMLT(2)
      CGWET(3) = CGMLT(3)
      CGWET(4) = CGMLT(4)
!sk ********************************************************************
      ES0 = 6.107799961
      CES0 = EPS*ES0
!
!     RIME_CI_TO_SNOW: TWO CONSTANTS
!     C2BRG HAS CONVERSION FACTOR OF 10**3
!
      C1BRG = DTL/RMI50
!lin ********************************************************************
!lin  C2BRG = RI50**2*1.E3 ! Error
      C2BRG = PIE*RI50**2*1.E3
!lin ********************************************************************
!
!     SATURATION VAPOR PRESSURE VALUES ARE TABULATED EACH WHOLE
!     DEG. (C) FOR -100 < TC < 50 FOR WATER AND -100 < TC < 0 FOR ICE.
!     DERIVATIVES ARE CENTERED AT EACH HALF DEG. (C).
!

!
!     SATURATION ADJUSTMENT:  RANGE OF MIXED ICE-WATER SATURATION
!     AND SATURATION VAPOR PRESSURE OVER WATER AND ICE AT T = -40 C.
!
      TDIF  = TICE-TTFRZ
      ESW00 = GGESW(DBLE(TICE-TDIF))
      ESI00 = GGESI(DBLE(TICE-TDIF))
!
      RETURN
      END SUBROUTINE SETUPM

      DOUBLE PRECISION FUNCTION ESATW_int(tt, ESWT, temp_K)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: ESWT(150001), temp_K(150001), tt
      DOUBLE PRECISION :: TTT
      INTEGER :: I

!	ESATW_int = 2.53e11*EXP(-5420./tt)/100.
!  RETURN

      TTT = tt-temp_K(1)
      I = INT(TTT/0.001)+1

!      IF(I.LT.1 .OR. I.GT.151) THEN
!	CALL BOMB(I,T)
!      ENDIF

	if(I < 1) then
		I = 1	
		ESATW_int = ESWT(I)
	else
      ESATW_int = ESWT(I)+ (tt - temp_K(I))*(ESWT(I+1)-ESWT(I))/0.001
	endif
	if(tt - temp_K(I) > 0.0015) stop "654321"
      RETURN
      END FUNCTION ESATW_int


      DOUBLE PRECISION FUNCTION ESATI_int (tt, ESIT,temp_K )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: ESIT(150001), temp_K(150001), tt
      DOUBLE PRECISION :: TTT
      INTEGER:: I
!     ESATI = 3.41e12*EXP(-6130./tt)/100.
!     RETURN


      TTT = tt-temp_K(1)
      I = INT(TTT/0.001)+1
!      IF(I.LT.1 .OR. I.GT.151) THEN
!	CALL BOMB(I,tt)
!      ENDIF
	if(I < 1) then
		I = 1	
     		ESATI_int = ESIT(I)
	else
      ESATI_int = ESIT(I)+(tt - temp_K(I))*(ESIT(I+1)-ESIT(I))/0.001
	endif
      RETURN
      END FUNCTION ESATI_int

      DOUBLE PRECISION FUNCTION DESWDT ( tt, DESWT, temp_K )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: DESWT(150001), temp_K(150001), tt
      DOUBLE PRECISION :: TTT
      INTEGER:: I

!       DESWDT = 2.5e6*(2.53e11*EXP(-5420./tt)/100.)/(461.5*tt*tt)
!    return

       TTT = tt-temp_K(1)
       I = INT(TTT/0.001)+1
!      IF(I.LT.1 .OR. I.GT.150) THEN
!	CALL BOMB(I,T)
!      ENDIF
	if(I < 1) then
		I = 1	
		DESWDT = DESWT(I)
	else
      DESWDT = DESWT(I)+(tt - temp_K(I))*(DESWT(I+1)-DESWT(I))/0.001
	endif
      RETURN
      END FUNCTION DESWDT
!
      DOUBLE PRECISION FUNCTION DESIDT ( tt, DESIT,temp_K )
 	IMPLICIT NONE
 	DOUBLE PRECISION, INTENT(IN) :: tt, DESIT(150001),temp_K(150001)
	INTEGER :: I
 	DOUBLE PRECISION :: TTT

!      DESIDT = 2.83e6*(3.41e12*EXP(-6130./tt)/100.)/(461.5*tt*tt)
! return

       TTT = tt-temp_K(1)
       I = INT(TTT/0.001)+1
!      IF(I.LT.1 .OR. I.GT.150) THEN
!	CALL BOMB(I,tt)
!      ENDIF
	if(I < 1) then
		I = 1
 	      DESIDT = DESIT(I)
	else
      DESIDT = DESIT(I)+(tt - temp_K(I))*(DESIT(I+1)-DESIT(I))/0.001
 	endif
       RETURN
      END FUNCTION DESIDT




      DOUBLE PRECISION FUNCTION GGESI(T)
      IMPLICIT NONE
	DOUBLE PRECISION :: T
      DOUBLE PRECISION :: A,B, C1, C2, C3, C4
!
!     SATURATION VAPOR PRESSURE OVER ICE
!      (GOFF AND GRATCH)
!
!     ESI     SATURATION VAPOR PRESSURE  (MB)
!     T       TEMP  (KELVIN)
!
      DATA C1/-9.09718/C2/-3.56654/C3/0.876793/C4/0.78583503/
!
      A = 273.16/T
      B = C1*(A-1.0)+C2*DLOG10(A)+C3*(1.0-1.0/A)+C4
      GGESI = 10.0**B
      RETURN
      END FUNCTION GGESI



      DOUBLE PRECISION FUNCTION GGESW (TA)
	IMPLICIT NONE 
	DOUBLE PRECISION, INTENT(IN) :: TA
	DOUBLE PRECISION :: TS, F5, E1, E2, F1, F2, F3, F4, F
!
!     SATURATION VAPOR PRESSURE OVER WATER
!          (GOFF AND GRATCH)
!
!     TA IS TEMPERATURE IN DEG KELVIN
!     ES IS SATURATION VAPOR PRESSURE IN MB
!
      DATA TS/373.16/F5/3.0057149/

! 	GGESW = 2.53e11*EXP(-5420./TA)/100.
! RETURN

      E1 = 11.344*(1.0 - TA/TS)
      E2 = -3.49149*(TS/TA - 1.0)
      F1 = -7.90298*(TS/TA - 1.0)
      F2 = 5.02808*DLOG10(TS/TA)
      F3 = -1.3816E-7*(10.0**E1 -1.0)
      F4 = 8.1328E-3*(10.0**E2 - 1.0)
      F = F1 + F2 + F3 + F4 + F5
      GGESW = 10.0**F
      RETURN
      END FUNCTION GGESW


      SUBROUTINE CONDNS(T,P,QVK,QSW,ESW,QLK,ISTOP, &
    	  CPLC, ESWT,DESWT, temp_K, EPS)

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: ESWT(150001), DESWT(150001), temp_K(150001)
      REAL, INTENT(IN) :: EPS, P, CPLC
      INTEGER, INTENT(INOUT) :: ISTOP
      REAL, INTENT(INOUT) :: T,QVK,QSW,ESW,QLK

      REAL:: EX, GAMMAW, SS, GAMFAC, CRIT
      INTEGER :: N
!
!xmk  DATA CRIT/1.E-6/
      DATA CRIT/1.E-8/    !thresholds / 100
!
!      print *, INPUT FOR CONDNS:: T,p,q_v, q_v,s, e_s, q_l, ISTOP , &
!		T,P,QVK,QSW,ESW,QLK,ISTOP

      N = 0
      GAMFAC = EPS*CPLC*P
!     WRITE(6,*)   10  ADJUST TO SAT---NO ICE   CRIT=,CRIT
   10 N = N+1
      IF(N .GT. 10)  GO TO 1000
      SS = QVK-QSW
      GAMMAW = GAMFAC*DESWDT(DBLE(T), DESWT, temp_K)/(P-ESW)**2
      EX = SS/(1.+GAMMAW)
      IF(N .EQ.  1)  GO TO  20
      IF(ABS(EX/QLK) .LT. CRIT)  RETURN
!       print *, T, CPLC, EX

   20 T = T+CPLC*EX
!       print *, T, CPLC, EX
!      print *,HERE 11.1
 
      ESW = REAL(ESATW_int(DBLE(T), ESWT, temp_K))
!      print *,HERE 11.2

      QSW = EPS*ESW/(P-ESW)
      QVK = QVK-EX
      QLK = QLK+EX
!     WRITE(6,111) N,T,SS,EX,QSW,QVK,QLK,GAMMAW
  111 FORMAT('0',I5,F14.7,6E16.8)
      GO TO  10
!
 1000 ISTOP = 1111
      RETURN
      END SUBROUTINE CONDNS


      SUBROUTINE SUBVAP(T,P,QVK,QSI,ESI,QIK,ISTOP, &
	CPLS,ESIT,DESIT,temp_K,EPS)

      IMPLICIT NONE
      REAL, INTENT(INOUT) :: T, QVK,QSI,ESI,QIK
      REAL, INTENT(IN) :: P
      INTEGER, INTENT(INOUT) :: ISTOP 
      REAL, INTENT(IN) :: CPLS,EPS
 	DOUBLE PRECISION, INTENT(IN) :: ESIT(150001), DESIT(150001), temp_K(150001)    
      REAL :: EX, GAMMAI, SS , GAMFAC, CRIT 
      INTEGER :: N

!    
!xmk  DATA CRIT/0.5E-7/
      DATA CRIT/0.5E-9/    !thresholds / 100
!
      N = 0
      GAMFAC = EPS*CPLS*P
!     WRITE(6,*)       ADJUST TO SAT---ICE ONLY   CRIT=,CRIT
   10 N = N+1
      IF(N .GT. 10)  GO TO 1000
      SS = QVK-QSI
      GAMMAI = GAMFAC*DESIDT(DBLE(T), DESIT, temp_K)/(P-ESI)**2
      EX = SS/(1.+GAMMAI)
      IF(N .EQ.  1)  GO TO 20
      IF(ABS(EX/QIK) .LT. CRIT)  RETURN
   20 T = T+CPLS*EX
      ESI = REAL(ESATI_int(DBLE(T), ESIT, temp_K))
      QSI = EPS*ESI/(P-ESI)
      QVK = QVK-EX
      QIK = QIK+EX
!     WRITE(6,111) N,T,SS,EX,QSI,QVK,QIK,GAMMAI
  111 FORMAT('0',I5,F14.7,6E16.8)
      GO TO 10
!
 1000 ISTOP = 2222
      RETURN
      END SUBROUTINE SUBVAP

      REAL FUNCTION gammln(xx)
       IMPLICIT NONE
       REAL, INTENT(IN):: xx
   
       INTEGER ::  j
       DOUBLE PRECISION :: ser, stp, tmp, x,y,cof(6)
       SAVE cof, stp
       DATA cof, stp/76.18009172947146d0,-86.50532032941677d0, &
	24.01409824083091d0, -1.231739572450155d0, &
        0.1208650973866179d-2, -0.5395239384953d-5, &
        2.5066282746310005d0/

       x = xx
       y = x
       tmp = x+5.5d0
       tmp=(x+0.5d0)*log(tmp)-tmp
       ser=1.000000000190015d0
       do j=1,6
          y=y+1.d0
          ser=ser+cof(j)/y
       enddo
       gammln = tmp+log(stp*ser/x)
       return
      END FUNCTION gammln

      REAL FUNCTION gammp(a,x)
        IMPLICIT NONE
	REAL, INTENT(IN):: a, x
	REAL gammcf, gamser,gln
	if(x.lt.0. .or. a .le. 0.) pause 'bad args in gammp'
	if(x .lt. a+1.)then
		call gser( gamser,a,x,gln)
		gammp = gamser
	else
		call gcf(gammcf,a,x,gln)
		gammp = 1. - gammcf
	endif
	gammp = gammp*exp(gln)
	return
      END FUNCTION gammp

       REAL FUNCTION gammq_frac(a,x)
       IMPLICIT NONE
	REAL, INTENT(IN):: a, x
	REAL gammcf, gamser,gln
	if(x.lt.0. .or. a .le. 0.) pause 'bad args in gammq_frac'
	if(x .lt. a+1.)then
		call gser( gamser,a,x,gln)
		gammq_frac = 1.-gamser
	else
		call gcf(gammcf,a,x,gln)
		gammq_frac = gammcf
	endif
	return
	END FUNCTION gammq_frac

	SUBROUTINE gser(gamser,a,x,gln)
       IMPLICIT NONE
	INTEGER :: ITMAX
	REAL, INTENT(IN):: a,x
	REAL, INTENT(INOUT):: gamser,gln
	REAL EPS
	PARAMETER(ITMAX=100,EPS=3.E-7)
	INTEGER n
	REAL :: ap,del,sum
!       real :: gammln
	gln=gammln(a)
	if(x .le. 0.) then
		if(x .lt. 0.)pause 'x < 0 in gser'
		gamser = 0.
		return
	endif
	ap=a
	sum = 1./a
	del = sum
	do n=1,ITMAX
		ap=ap+1.
		del=del*x/ap
		sum=sum+del
		if(abs(del) .lt. abs(sum)*EPS) goto 1
	enddo
	pause 'a too large, ITMAX too small in gser'
  1     gamser=sum*exp(-x+a*log(x)-gln)
	return
	END SUBROUTINE gser

	SUBROUTINE gcf(gammcf,a,x,gln)
       IMPLICIT NONE
	INTEGER :: ITMAX
	REAL, INTENT(IN):: a,x
	REAL, INTENT(INOUT)::gammcf,gln
	REAL EPS, FPMIN
	PARAMETER(ITMAX=100,EPS=3.E-7,FPMIN=1.E-30)
	INTEGER :: i
	REAL :: an,b,c,d,del,h
!        real :: gammln
	gln=gammln(a)
	b=x+1.-a
	c=1./FPMIN
	d=1./b
	h=d
	do i=1,ITMAX
		an=-i*(i-a)
		b=b+2.
		d=an*d+b
		if(abs(d) .lt. FPMIN)d=FPMIN
		c=b+an/c
		if(abs(c) .lt. FPMIN)c=FPMIN
		d=1./d
		del=d*c
		h=h*del
		if(abs(del-1.) .lt. EPS) goto 1
	enddo
	pause 'a too large, ITMAX too small in gcf'
  1     gammcf=exp(-x+a*log(x)-gln)*h
	return
	END SUBROUTINE gcf









! fu *******************************************************************
      SUBROUTINE RIME_CI_TO_SNOW(TC,QLL,QI,QLRHO,PSFW, &
        C1BRG,C2BRG,RMI50,RMI40)
      IMPLICIT NONE
	REAL :: T_FRZ_HOM_DEGC
	PARAMETER(T_FRZ_HOM_DEGC = -36.)
! fu ******************************************************************
!
       REAL, INTENT(INOUT) :: PSFW
      REAL, INTENT(IN):: QLRHO, QI, QLL, TC
      REAL, INTENT(IN)  :: C1BRG,C2BRG,RMI50,RMI40
      REAL :: DT1, A21, a2,  a1, TC1, A1T, A2T

! fu ******************************************************************
      DIMENSION A1T(0:31),A2T(0:31)
! fu   DIMENSION A1T(31),A2T(31)
! fu ******************************************************************
      DATA A1T/0.,0.7939E-7,0.7841E-6,0.3369E-5,0.4336E-5,0.5285E-5, &
              0.3728E-5,0.1852E-5,0.2991E-6,0.4248E-6,0.7434E-6, &
              0.1812E-5,0.4394E-5,0.9145E-5,0.1725E-4,0.3348E-4, &
              0.1725E-4,0.9175E-5,0.4412E-5,0.2252E-5,0.9115E-6, &
              0.4876E-6,0.3473E-6,0.4758E-6,0.6306E-6,0.8573E-6, &
              0.7868E-6,0.7192E-6,0.6513E-6,0.5956E-6,0.5333E-6, &
              0.4834E-6/
      DATA A2T/0.,0.4006,0.4831,0.5320,0.5307,0.5319, &
              0.5249,0.4888,0.3894,0.4047,0.4318, &
              0.4771,0.5183,0.5463,0.5651,0.5813, &
              0.5655,0.5478,0.5203,0.4906,0.4447, &
              0.4126,0.3960,0.4149,0.4320,0.4506, &
              0.4483,0.4460,0.4433,0.4413,0.4382, &
              0.4361/
!
      TC1 = AMAX1(TC,-30.0)
! fu ****************************************************************
      if ( tc1 .gt. -1.0 ) then
      a1 = a1t(1)
      a2 = a2t(1)
      else
      A1  = (A1T(-INT(TC1))-A1T(-INT(TC1)+1))*(TC1-INT(TC1)+1.0)+ &
           A1T(-INT(TC1)+1)
      A2  = (A2T(-INT(TC1))-A2T(-INT(TC1)+1))*(TC1-INT(TC1)+1.0)+ &
           A2T(-INT(TC1)+1)
      endif
! fu   A1  = A1T(-INT(TC1)+1)
! fu   A2  = A2T(-INT(TC1)+1)
! fu ****************************************************************
      A21 = 1.0 - A2
!
      DT1 = (RMI50**A21 - RMI40**A21) / (A1*A21)
!
!     NOTE:  MKS UNITS, UI50=1.0 M/SEC, EIW=1.0
!
      PSFW = C1BRG*QI/DT1*(QLRHO*C2BRG)
      if(TC < T_FRZ_HOM_DEGC) PSFW = 0.

! fu *******************************************************************
! fu ******************************************************************
      RETURN
      END SUBROUTINE RIME_CI_TO_SNOW
! fu *******************************************************************


      REAL FUNCTION ACR1(Q,QRHO,C,PWR,RHOFAC,RHO)
      IMPLICIT NONE
	REAL, INTENT(IN):: Q,QRHO,C,PWR,RHOFAC,RHO

!
!     RHO FACTOR FOR SOURCE TERMS:  PSACW,PWACS,PRACI,PSACI
!
      ACR1 = C * Q * QRHO**PWR * RHOFAC
      RETURN
      END FUNCTION ACR1



      REAL FUNCTION ACR2(Q,QRHO,C,PWR,RHOFAC,RHO)
      IMPLICIT NONE
      REAL, INTENT(IN) :: Q,QRHO,C,PWR,RHOFAC,RHO

      ACR2 = C * Q * QRHO**PWR /SQRT(RHO)
      RETURN
      END FUNCTION ACR2





      REAL FUNCTION ACR3(V1,V2,QRHO1,QRHO2,C,CAC,RHOFAC,RHO)
      IMPLICIT NONE
	REAL, INTENT(IN) :: V1,V2,QRHO1,QRHO2,C,CAC,RHOFAC,RHO
	REAL :: A
	INTEGER :: K
!
!      SOURCE TERMS:  PRACS,PSACR,PGACR,PGACS
!
      DIMENSION CAC(3)
      A=0.0
      DO 10 K=1,3
      A = A + CAC(K) * (QRHO1**((7-K)*0.25) * QRHO2**(K*0.25))
   10 CONTINUE
      ACR3 = C * ABS(V1-V2) * A / RHO
      RETURN
      END FUNCTION ACR3


      REAL FUNCTION REVP(S,TSQ,QS,QRHO,C, RHO, RHOFAC)
      IMPLICIT NONE
      REAL, INTENT(IN)  ::  S,TSQ,QS,QRHO,C(5), RHO, RHOFAC    
!
!     EVAPORATION
!
!      REVP = C(1)*TSQ*QS*(1.0-S) * (C(2)*SQRT(QRHO)+C(3)*QRHO**0.725) &
!        / (C(4)*TSQ+C(5)*QS*RHO)
! VTJP:  ssl says RHOFAC**1/2 is missing 

      REVP = C(1)*TSQ*QS*(1.0-S) * (C(2)*SQRT(QRHO)+C(3)*(QRHO**0.725)*SQRT(RHOFAC)) &
        / (C(4)*TSQ+C(5)*QS*RHO)


      RETURN
      END FUNCTION REVP



      REAL FUNCTION AUT(C,Q,Q0)
       IMPLICIT NONE
	REAL, INTENT(IN)  :: C,Q,Q0
!
!      SOURCE TERMS:  PGAUT,PSAUT
!
      AUT = AMAX1(C*(Q-Q0),0.0)
      RETURN
      END FUNCTION AUT


      REAL FUNCTION RAUT(C,QEXRHO, RHO, N1X, DDX)
      IMPLICIT NONE

      REAL, INTENT(IN) ::  C(2), QEXRHO, RHO, N1X, DDX

!
!     SEE ORVILLE AND KOPP (1977)
!
!      RAUT = QEXRHO**3 / ((C(1)*QEXRHO+C(2))*RHO)
!
!   VTJP:: ssl finds error in autoconversion constants (see Lin et al  (1983) and Orville and Kopp (1977)

      RAUT = QEXRHO**3 / ((1.2E-4*QEXRHO+1.596E-12*N1X/DDX)*RHO)       



      RETURN
      END FUNCTION RAUT

      REAL FUNCTION GFR(TC,QRR,C, RHO)
       IMPLICIT NONE
      REAL, INTENT(IN) ::  TC, QRR, C(2), RHO
    
      GFR = C(1)*(EXP(-C(2)*TC)-1.0)*QRR**1.75/RHO
      RETURN
      END FUNCTION GFR
!sk ********************************************************************


      REAL FUNCTION SMLT(TC,DQS,QSRHO,PSACW,PSACR,C,RHO, RHOFAC)
      IMPLICIT NONE
      REAL, INTENT(IN) ::  TC,DQS,QSRHO,PSACW,PSACR,C(5),RHO, RHOFAC
     
      SMLT = (C(1)*TC/RHO-C(2)*DQS) * (C(3)*SQRT(QSRHO)+ &
        C(4)*QSRHO**0.65625*SQRT(RHOFAC)) + C(5)*TC*(PSACW+PSACR)
!sk ********************************************************************
      RETURN
      END FUNCTION SMLT
 

     REAL FUNCTION GMLT(TC,DQS,QGRHO,PGACW,PGACR,C, RHO)
      IMPLICIT NONE
      REAL, INTENT(IN) ::  TC,DQS,QGRHO,PGACW,PGACR,C(5),RHO
     
!
!     NOTE:  PGACW AND PGACR MUST BE CALC BEFORE GMLT IS CALLED
!
      GMLT = (C(1)*TC/RHO-C(2)*DQS) * (C(3)*SQRT(QGRHO)+ &
        C(4)*QGRHO**0.6875/RHO**0.25) + C(5)*TC*(PGACW+PGACR)
      RETURN
      END FUNCTION GMLT
!lin *******************************************************************

      REAL FUNCTION GWET(DQS,TC,QGRHO,PGACI,PGACS,C,SIMAX,SSMAX, &
	HLTF,CH2O,CICE,RHO)
	IMPLICIT NONE
	REAL, INTENT(IN) :: DQS,TC,QGRHO,PGACI,PGACS,SIMAX,SSMAX

!lin *******************************************************************
!csk
!
      REAL, INTENT(IN) ::  HLTF,CH2O,CICE,C(4)
!
      REAL, INTENT(IN) :: RHO
      REAL :: x, y

!
!     NOTE:  PGACI AND PGACS MUST BE CALC BEFORE GWET IS CALLED
!      FACTOR OF 10 CONVERTS PGACI TO PGACI PRIME
!lin   FACTOR OF 10 CONVERTS PGACS TO PGACS PRIME
!
!lin *******************************************************************
      x = AMIN1(PGACI * 10.0,SIMAX)
      y = AMIN1(10.*PGACS,SSMAX)
      GWET = AMAX1((C(2)*DQS-C(1)*TC/RHO)*HLTF/(HLTF+CH2O*TC)* & 
        (C(3)*SQRT(QGRHO)+C(4)*QGRHO**0.6875/RHO**0.25) +  & 
        (x+y)*(1.0-CICE*TC/(HLTF+CH2O*TC)),0.0) 
!lin    (10.*PGACI+PGACS)*(1.0-CICE*TC/(HLTF+CH2O*TC)),0.0)
!lin *******************************************************************
      RETURN
      END FUNCTION GWET

      REAL FUNCTION GSUB(S,TSQ,QS,QRHO,C, RHO, RHOFAC)
      IMPLICIT NONE
      REAL, INTENT(IN) :: S,TSQ,QS,QRHO,C(5), RHO, RHOFAC
 
!
!     GRAUPEL SUBLIMATION
!
      GSUB = C(1)*TSQ*QS*(S-1.0) * (C(2)*SQRT(QRHO) + &
            C(3)*QRHO**0.6875/RHO**0.25) / (C(4)*TSQ+C(5)*QS*RHO)
      RETURN
      END FUNCTION GSUB


      REAL FUNCTION SSUB(S,TSQ,QS,QRHO,C, RHO, RHOFAC)
      IMPLICIT NONE
      REAL, INTENT(IN) :: S,TSQ,QS,QRHO, C(5), RHO, RHOFAC 
!
!     SNOW SUBLIMATION
!
      SSUB = C(1)*TSQ*QS*(S-1.0) * (C(2)*SQRT(QRHO) + &
            C(3)*QRHO**0.65625*SQRT(RHOFAC)) / (C(4)*TSQ+C(5)*QS*RHO)
      RETURN
      END FUNCTION SSUB



  
  
!
!==================================================================================
!		PRECIPITATION FLUX (VTJP)
!			based on algorithm by Shuhua Chen at NCAR
!==================================================================================
!


   SUBROUTINE PRECIP_FALLOUT(qrz, qgz, qsz, qiz, qlz, niz, nwz, dt, VCONR, VCONS, &
	VCONG, VCONI, a_i, b_i, c_i, d_i, a_w, b_w, c_w, LAMBDA_FAC_w, LAMBDA_FAC_i, &
        GAMMA_RAT_w_mass, GAMMA_RAT_i_mass, &
	GAMMA_RAT_w_number, GAMMA_RAT_i_number, z_sfc, RHOS, rhoz, &
        dzw, zz, pptrain, pptsnow, pptgraul, kts, kte)

   IMPLICIT NONE
	INTEGER ICE_NUMBER, HEYMSFIELD_DONNER_VT
	PARAMETER(ICE_NUMBER = 1, HEYMSFIELD_DONNER_VT = 0)

   INTEGER, INTENT(IN) ::  kts, kte
   REAL, DIMENSION(kts:kte), INTENT(INOUT)  :: qrz, qgz, qsz, qiz, niz, nwz, qlz  
   REAL, DIMENSION(kts:kte), INTENT(IN)  :: rhoz, dzw, zz
   REAL, INTENT(INOUT) :: pptrain, pptsnow, pptgraul
   REAL, INTENT(IN) :: dt, VCONR, VCONS, VCONG, VCONI, z_sfc, RHOS, a_i, &
            b_i, c_i, d_i, a_w, b_w, c_w,  GAMMA_RAT_i_mass,GAMMA_RAT_w_mass, &
	    GAMMA_RAT_i_number, LAMBDA_FAC_i, GAMMA_RAT_w_number, LAMBDA_FAC_w

   INTEGER                       :: min_q, max_q, k
   REAL				:: t_del_tv, del_tv, flux, fluxin, fluxout, &
                                      lambda_i, lambda_w
   LOGICAL			:: notlast
   REAL, DIMENSION(kts:kte) :: vtrold, vtsold, vtgold, vtiold, vtlold
   

    t_del_tv=0.
    del_tv=dt
    notlast=.true.
    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
!        if (qrzold(k) .gt. 1.0e-12) then
!!       if (qrzold(k) .gt. 1.0e-10) then
         if (qrz(k) .gt. 1.0e-10) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)

            vtrold(k)= VTRR(qrz(k)*rhoz(k),VCONR,rhoz(k),SQRT(RHOS/rhoz(k)))
	   
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(zz(k)-z_sfc)/vtrold(k))
            else
               del_tv=amin1(del_tv,0.9*(zz(k)-zz(k-1))/vtrold(k))
            endif
         else
            vtrold(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dt)

         t_del_tv=t_del_tv+del_tv
!
         if ( t_del_tv .ge. dt ) then
              notlast=.false.
              del_tv=dt+del_tv-t_del_tv
         endif

! use small delta t to calculate the q?z flux 
! termi is the qrzold flux pass in the grid box through the upper boundary
! termo is the qrzold flux pass out the grid box through the lower boundary
!
!
!
! fluxin, fluxout are:-  kg/m2/sec
!	del_tv is in sec
!
	if(1./2./3. .ne. 1./(2.*3.)) stop 4321 

         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rhoz(k)*vtrold(k)*qrz(k)
            flux=(fluxin-fluxout)/rhoz(k)/dzw(k)
            if(qrz(k)+del_tv*flux < 0.) then
		flux = -qrz(k)/del_tv
		fluxout = fluxin-flux*rhoz(k)*dzw(k)
            	qrz(k)=0.
	    else
		qrz(k)=qrz(k)+del_tv*flux
	    endif
            fluxin=fluxout
         enddo 
         if (min_q .eq. 1) then
            pptrain=pptrain+fluxin*del_tv
         else
            qrz(min_q-1)=qrz(min_q-1)+del_tv*  &
	                  fluxin/rhoz(min_q-1)/dzw(min_q-1)
         endif
!
      else
         notlast=.false.
      endif
    ENDDO

!
!-- snow
!
    t_del_tv=0.
    del_tv=dt
    notlast=.true.

    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
!        if (qszold(k) .gt. 1.0e-12) then
!!       if (qszold(k) .gt. 1.0e-10) then
         if (qsz(k) .gt. 1.0e-10) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)

            vtsold(k)= VTRS(qsz(k)*rhoz(k),VCONS,rhoz(k),SQRT(RHOS/rhoz(k)))
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(zz(k)-z_sfc)/vtsold(k))
            else
               del_tv=amin1(del_tv,0.9*(zz(k)-zz(k-1))/vtsold(k))
            endif
         else
            vtsold(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dt)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dt ) then
              notlast=.false.
              del_tv=dt+del_tv-t_del_tv
         endif

! use small delta t to calculate the qszold flux 
! termi is the qszold flux pass in the grid box through the upper boundary
! termo is the qszold flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rhoz(k)*vtsold(k)*qsz(k)
            flux=(fluxin-fluxout)/rhoz(k)/dzw(k)
            if(qsz(k)+del_tv*flux < 0.) then
		flux = -qsz(k)/del_tv
		fluxout = fluxin-flux*rhoz(k)*dzw(k)
            	qsz(k)=0.
	    else
		qsz(k)=qsz(k)+del_tv*flux
	    endif
            fluxin=fluxout
         enddo 
         if (min_q .eq. 1) then
            pptsnow=pptsnow+fluxin*del_tv
         else
            qsz(min_q-1)=qsz(min_q-1)+del_tv*  &
	                 fluxin/rhoz(min_q-1)/dzw(min_q-1)
         endif
!
      else
         notlast=.false.
      endif

    ENDDO
!
!-- graupel
!
    t_del_tv=0.
    del_tv=dt
    notlast=.true.
!
    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
!        if (qgzold(k) .gt. 1.0e-12) then
!!       if (qgzold(k) .gt. 1.0e-10) then
         if (qgz(k) .gt. 1.0e-10) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)
  
            vtgold(k)= VTRG(qgz(k)*rhoz(k),VCONG,rhoz(k),SQRT(RHOS/rhoz(k)))
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(zz(k)-z_sfc)/vtgold(k))
            else
               del_tv=amin1(del_tv,0.9*(zz(k)-zz(k-1))/vtgold(k))
            endif
         else
            vtgold(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dt)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dt ) then
              notlast=.false.
              del_tv=dt+del_tv-t_del_tv
         endif

! use small delta t to calculate the qgzold flux 
! termi is the qgzold flux pass in the grid box through the upper boundary
! termo is the qgzold flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rhoz(k)*vtgold(k)*qgz(k)
            flux=(fluxin-fluxout)/rhoz(k)/dzw(k)
            if(qgz(k)+del_tv*flux < 0.) then
		flux = -qgz(k)/del_tv
		fluxout = fluxin-flux*rhoz(k)*dzw(k)
            	qgz(k)=0.
	    else
		qgz(k)=qgz(k)+del_tv*flux
	    endif

            fluxin=fluxout
         enddo 
         if (min_q .eq. 1) then
            pptgraul=pptgraul+fluxin*del_tv
         else
            qgz(min_q-1)=qgz(min_q-1)+del_tv*  &
	                 fluxin/rhoz(min_q-1)/dzw(min_q-1)
         endif
!
      else
         notlast=.false.
      endif
!
   ENDDO


!=========================================================================
!	Precipitation flux of cloud-ice mass - VTJP (Jan/2002)
! 
!=========================================================================

!
!-- cloud-ice mass
!

    t_del_tv=0.
    del_tv=dt
    notlast=.true.
!
    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
         if (qiz(k) .gt. 1.0e-10) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)
  	    lambda_i = LAMBDA_FAC_i*c_i*niz(k)/qiz(k)
	    lambda_i = lambda_i**(1./d_i) ! CORRECT
            vtiold(k)= VTRI(qiz(k), VCONI, rhoz(k), a_i, b_i, GAMMA_RAT_i_mass, lambda_i, HEYMSFIELD_DONNER_VT)
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(zz(k)-z_sfc)/vtiold(k))
            else
               del_tv=amin1(del_tv,0.9*(zz(k)-zz(k-1))/vtiold(k))
            endif
         else
            vtiold(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dt)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dt ) then
              notlast=.false.
              del_tv=dt+del_tv-t_del_tv
         endif

! use small delta t to calculate the qgzold flux 
! termi is the qgzold flux pass in the grid box through the upper boundary
! termo is the qgzold flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rhoz(k)*vtiold(k)*qiz(k)
            flux=(fluxin-fluxout)/rhoz(k)/dzw(k)
            if(qiz(k)+del_tv*flux < 0.) then
		flux = -qiz(k)/del_tv
		fluxout = fluxin-flux*rhoz(k)*dzw(k)
            	qiz(k)=0.
	    else
		qiz(k)=qiz(k)+del_tv*flux
	    endif
            fluxin=fluxout
         enddo 
         if (min_q .eq. 1) then
            pptsnow=pptsnow+fluxin*del_tv
         else
            qiz(min_q-1)=qiz(min_q-1)+del_tv*  &
	                 fluxin/rhoz(min_q-1)/dzw(min_q-1)
         endif
  
!
      else
         notlast=.false.
      endif
!


   ENDDO

!
!-- cloud-ice number
!

    t_del_tv=0.
    del_tv=dt
    notlast=.true.
!
    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
         if (qiz(k) .gt. 1.0e-10 .and. niz(k) .gt. 1.0e-10) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)
  	    lambda_i = LAMBDA_FAC_i*c_i*niz(k)/qiz(k)
	    lambda_i = lambda_i**(1./d_i) ! CORRECT
            vtiold(k)= VTRI(qiz(k), VCONI, rhoz(k), a_i, b_i, GAMMA_RAT_i_number, lambda_i, HEYMSFIELD_DONNER_VT)
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(zz(k)-z_sfc)/vtiold(k))
            else
               del_tv=amin1(del_tv,0.9*(zz(k)-zz(k-1))/vtiold(k))
            endif
         else
            vtiold(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dt)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dt ) then
              notlast=.false.
              del_tv=dt+del_tv-t_del_tv
         endif

! use small delta t to calculate the qgzold flux 
! termi is the qgzold flux pass in the grid box through the upper boundary
! termo is the qgzold flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rhoz(k)*vtiold(k)*niz(k)
            flux=(fluxin-fluxout)/rhoz(k)/dzw(k)
            niz(k)=niz(k)+del_tv*flux
            niz(k)=amax1(0.,niz(k))
            fluxin=fluxout
         enddo 
         if (min_q .ne. 1) then
            niz(min_q-1)=niz(min_q-1)+del_tv*  &
	                 fluxin/rhoz(min_q-1)/dzw(min_q-1)
         endif
  
!
      else
         notlast=.false.
      endif
!


   ENDDO


!
!-- cloud-water mass
!

    t_del_tv=0.
    del_tv=dt
    notlast=.true.
!
    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
         if (qlz(k) .gt. 1.0e-10) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)
   	    lambda_w = LAMBDA_FAC_w*c_w*nwz(k)/qlz(k)
	    lambda_w = lambda_w**(1./3.) ! CORRECT
            vtlold(k)= VTRL(rhoz(k),a_w, b_w, GAMMA_RAT_w_mass, lambda_w)
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(zz(k)-z_sfc)/vtlold(k))
            else
               del_tv=amin1(del_tv,0.9*(zz(k)-zz(k-1))/vtlold(k))
            endif
         else
            vtlold(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dt)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dt ) then
              notlast=.false.
              del_tv=dt+del_tv-t_del_tv
         endif

! use small delta t to calculate the qgzold flux 
! termi is the qgzold flux pass in the grid box through the upper boundary
! termo is the qgzold flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rhoz(k)*vtlold(k)*qlz(k)
            flux=(fluxin-fluxout)/rhoz(k)/dzw(k)
            if(qlz(k)+del_tv*flux < 0.) then
		flux = -qlz(k)/del_tv
		fluxout = fluxin-flux*rhoz(k)*dzw(k)
            	qlz(k)=0.
	    else
		qlz(k)=qlz(k)+del_tv*flux
	    endif
            fluxin=fluxout
         enddo 
         if (min_q .eq. 1) then
            pptsnow=pptsnow+fluxin*del_tv
         else
            qlz(min_q-1)=qlz(min_q-1)+del_tv*  &
	                 fluxin/rhoz(min_q-1)/dzw(min_q-1)
         endif

   
!
      else
         notlast=.false.
      endif
!


   ENDDO


!
!-- cloud-water number
!

    t_del_tv=0.
    del_tv=dt
    notlast=.true.
!
    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
         if (qlz(k) .gt. 1.0e-10 .and. nwz(k) .gt. 1.0e-10) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)
   	    lambda_w = LAMBDA_FAC_w*c_w*nwz(k)/qlz(k)
	    lambda_w = lambda_w**(1./3.) ! CORRECT
            vtlold(k)= VTRL(rhoz(k),a_w, b_w, GAMMA_RAT_w_number, lambda_w)
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(zz(k)-z_sfc)/vtlold(k))
            else
               del_tv=amin1(del_tv,0.9*(zz(k)-zz(k-1))/vtlold(k))
            endif
         else
            vtlold(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dt)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dt ) then
              notlast=.false.
              del_tv=dt+del_tv-t_del_tv
         endif

! use small delta t to calculate the qgzold flux 
! termi is the qgzold flux pass in the grid box through the upper boundary
! termo is the qgzold flux pass out the grid box through the lower boundary
!
 
         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rhoz(k)*vtlold(k)*nwz(k)
            flux=(fluxin-fluxout)/rhoz(k)/dzw(k)
            nwz(k)=nwz(k)+del_tv*flux
            nwz(k)=amax1(0.,nwz(k))
            fluxin=fluxout
         enddo 
         if (min_q .ne. 1) then
            nwz(min_q-1)=nwz(min_q-1)+del_tv*  &
	                 fluxin/rhoz(min_q-1)/dzw(min_q-1)
         endif


!
      else
         notlast=.false.
      endif
!


   ENDDO



   END SUBROUTINE PRECIP_FALLOUT




    REAL FUNCTION VTRR(QRHO,VCONR,RHO,RHOFAC)
	IMPLICIT NONE
	REAL, INTENT(IN) :: QRHO,VCONR,RHO,RHOFAC
!
!     CALC MASS-WEIGHTED MEAN TERMINAL VELOCITIES
!
      VTRR = AMIN1(VCONR * QRHO**0.2 * RHOFAC,10.0)
      RETURN
      END FUNCTION VTRR


      REAL FUNCTION VTRS(QRHO,VCONS,RHO,RHOFAC)
	IMPLICIT NONE
	REAL, INTENT(IN) :: QRHO,VCONS,RHO,RHOFAC
!
!    C SNOW
!
      VTRS = VCONS * QRHO**0.0625 * RHOFAC
      RETURN
      END FUNCTION VTRS


      REAL FUNCTION VTRG(QRHO,VCONG,RHO,RHOFAC)
	IMPLICIT NONE
	REAL, INTENT(IN) :: QRHO,VCONG,RHO,RHOFAC
!
!     GRAUPEL
!
      VTRG = AMIN1(VCONG * QRHO**0.125 / SQRT(RHO),20.0)
      RETURN
      END FUNCTION VTRG
  
      REAL FUNCTION VTRI(QICE,VCONI,RHO, a_i, b_i, GAMMA_RAT_i_mass, lambda_i, heymsfield_flag)
	IMPLICIT NONE
	REAL, INTENT(IN) :: QICE,VCONI,RHO, a_i, b_i, GAMMA_RAT_i_mass, lambda_i
        INTEGER, INTENT(IN) :: heymsfield_flag
!
!     GRAUPEL
!
        if(heymsfield_flag == 1) then
 		VTRI = VCONI * (RHO * QICE)** 0.16
         else
		VTRI =  SQRT(1.2/RHO)*a_i*GAMMA_RAT_i_mass/(lambda_i**b_i)

        endif
      RETURN
      END FUNCTION VTRI

      REAL FUNCTION VTRL(RHO, a_w, b_w, GAMMA_RAT_w_mass, lambda_w)
	IMPLICIT NONE
	REAL, INTENT(IN) :: RHO, a_w, b_w, GAMMA_RAT_w_mass, lambda_w
  
!
!     GRAUPEL
!
		VTRL =  SQRT(1.2/RHO)*a_w*GAMMA_RAT_w_mass/(lambda_w**b_w)

       RETURN
      END FUNCTION VTRL

   SUBROUTINE code_storage


!	if(r_liq > 0. .and. r_ice > 0. .and. 0 == 1) then
!		if(lambda(1) /= 0. .and. lambda(2) /= 0.) then
!			Gr = Gamw*DBLE(qv_zero_vap)/QSW_zero;  ql_star = DBLE(qlz_vap)/rwoqswz;
!			dt_cut = dt_vap * ABS(delta_ql/delta_qliq_approx)
!			support_chelsea = 0
!			x_diff = 1.;!
!
!			arg1 = lambda(1)*lambda(2)*ql_star - A_q(2)*lambda(1) - A_q(1)*lambda(2) &
!					+A_T(2)*Gr*lambda(1) + A_T(1)*Gr*lambda(2)
!			arg2 = lambda(1)*lambda(2)*(A_q(1) + A_q(2) - A_T(1)*Gr - A_T(2)*Gr + QSW_zero - qv_zero_vap)
!			arg3 = (A_q(1) - A_T(1)*Gr)*lambda(2)
!			arg4 = (A_q(2) - A_T(2)*Gr)*lambda(1)!
!
!			dt_cut = (-arg1 - arg3 - arg4)/(-arg2 + arg3*lambda(1) + arg4*lambda(2))
!
!			root = (-2.*arg2 + 2.*arg3*lambda(1) + 2.*arg4*lambda(2))**2. - 4.*(2.*arg1 + 2.*arg3 + 2.*arg4)* &
!				(arg3*(lambda(1)**2.) + arg4*(lambda(2)**2.))
!			if(root >= 0.) then
!				root = DSQRT(root)
!				dt_cut_second_1 = (2.*arg2 - 2.*arg3*lambda(1) - 2*arg4*lambda(2) + root)&
!					/(2.*(arg3*(lambda(1)**2.) + arg4*(lambda(2)**2.)))
!				dt_cut_second_2 = (2.*arg2 - 2.*arg3*lambda(1) - 2*arg4*lambda(2) - root)&
!					/(2.*(arg3*(lambda(1)**2.) + arg4*(lambda(2)**2.)))
!				if(dt_cut_second_1 > 0. .and. dt_cut_second_1 < dt_vap) dt_cut = dt_cut_second_1
!				if(dt_cut_second_2 > 0. .and. dt_cut_second_2 < dt_vap) dt_cut = dt_cut_second_2
!			endif
!			print *, "DT_CUT = ", 	dt_cut		
!
!			do while(x_diff > 1.e-6)  
!				exdum(:)= DEXP(lambda(:) * DBLE(dt_cut))
!				arg1 = lambda(1)*lambda(2)*ql_star - A_q(2)*lambda(1) - A_q(1)*lambda(2) &
!					+A_T(2)*Gr*lambda(1) + A_T(1)*Gr*lambda(2)
!				arg2 = lambda(1)*lambda(2)*(A_q(1) + A_q(2) - A_T(1)*Gr - A_T(2)*Gr + QSW_zero - qv_zero_vap)
!				arg3 = (A_q(1) - A_T(1)*Gr)*lambda(2)
!				arg4 = (A_q(2) - A_T(2)*Gr)*lambda(1)
!				f_x = arg1 - arg2*dt_cut + arg3*exdum(1) + arg4*exdum(2)
!				fdash = -arg2 +arg3*lambda(1)*exdum(1) + arg4*lambda(2)*exdum(2)
!				dt_cut_save = dt_cut
!				if(fdash /= 0.) then
!					dt_cut = dt_cut - f_x/fdash
!				endif
!				x_diff = DABS(dt_cut - dt_cut_save)
!				support_chelsea = support_chelsea + 1
!				print *, ITER::, x_diff, support_chelsea, dt_cut
!				if(support_chelsea > 1000 .or. DABS(dt_cut) > 1e6 .or. fdash == 0. .or. dt_cut <= 0. .or. dt_cut > dt_vap) then
!					print *,  FAILURE TO CONVERGE
!					exit
!				endif
!			enddo
!
!		endif
!	endif


  END SUBROUTINE code_storage


END MODULE module_lin_2mom
