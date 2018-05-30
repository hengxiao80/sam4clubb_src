module microphysics
!
!==================================================================================================
!==================================================================================================
!==================================================================================================
! This module is an adaptation of the 
!
!	DOUBLE-MOMENT VERSION OF LIN MICROPHYSICS SCHEME 
!		by Vaughan Phillips (GFDL, 2006)
!
!		(developed from single-moment version by Lord et al., 1984, JAS, vol 41, 2836, and
!		Krueger et al, 1995, JAM, vol 34, 281)
! to SAM by Mikhail Ovtchinnikov (PNNL, Richland, WA, USA 2008, mikhail@pnl.gov)

!
! It should be located in the SRC/MICRO_lin_2mom directory.
! The following flags must be set
! docloud 	= .true.,
! doprecip 	= .true.,  
!==================================================================================================
!==================================================================================================
!==================================================================================================
! UWM changes are for software compatibility, some compilers follow for
! standard that stop cannot output more than a 5 digit integer 
! -dschanen 15 May 2008
! $Id: microphysics.F90,v 1.6 2008-08-11 15:26:43 dschanen Exp $ 

! grid is SAM module which contains the required grid information

use grid, only: nx, ny, nz, nzm, & ! grid dimensions; nzm=nz-1 - # of levels for all scalars
                dimx1_s,dimx2_s,dimy1_s,dimy2_s ! actual scalar-array dimensions 
use module_lin_2mom
use micro_params
implicit none

! It is probably easier to understand the procedure of microphysics implementation
! by following a specific example. Let's assume that we want to implement a
! 2-moment bulk microphysics package that has two prognostic variables for aerosol number 
! concentration (say CCN and IN), one prognostic water-vapor variable, and 
! 5 species of liquid/ice water each having mass and concentration characteristics; 
! therefore, our new microphysics scheme will have 10 prognostic liquid/ice water variables.  

! By prognostic variable I mean the one that needis to be explicitly advected and mixed
! by the dynamical core. In contrast, the diagnostic variables are those that can
! simply be computed from prognostic variables along; therefore, the diagnostic
! variable arrays can have simple (nx,ny,nzm) dimensions.

integer, parameter :: nmicro_fields = 12   ! total number of prognostic water vars

! Allocate the required memory for all the prognostic microphysics arrays:

real micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields)

! We assume that our prognostic variables are positioned in micro_field array as follows:

!  1 - CCN
!  2 - IN
!  3 - CCN_a number mixing ratio of CCN that have been activated as droplets 
!  4 - IN_a number mixing ratio of IN that have been activated as crystals 
!  5 - water vapor mixing ratio (kg/kg)
!  6 - cloud water mixing ratio (kg/kg)
!  7 - cloud ice mixing ratio (kg/kg)
!  8 - rain water mixing ratio (kg/kg)
!  9 - snow water mixing ratio (kg/kg)
! 10 - graupel/hail water mixing ratio (kg/kg)
! 11 - cloud water number mixing ratio (#/kg[air])
! 12 - cloud ice number mixing ratio (#/kg[air])

! For many reasons, for example, to be able to compute water budget, we may need to
! know which variables among the prognostic ones represent water mixing ratio, regardless
! of water species. We use a simple array of flags with 1 marking the water mass
! variable:
!                              variable             1 2 3 4 5 6 7 8 9 10 11 12
integer, parameter :: flag_wmass(nmicro_fields) = (/0,0,0,0,1,1,1,1,1, 1, 0, 0/) 

! To implement large-scale forcing, surface fluxes, etc, SAM needs to know
! which variable has a water vapor information. 

integer, parameter :: index_water_vapor = 5 ! index for variable that contains water vapor

! Now, we need to specify which variables describe precipitation. This is needed because
! SAM has two logical flags to deal with microphysics proceses - docloud and doprecip.
! docloud set to true means that condensation/sublimation processes are allowed to
! form clouds. However, the possibility of rain, snow, heil, etc., is controled by
! a second flag: doprecip. If doprecip=.false. than no precipitation is allowed, hence 
! no advection, diffusion, and fallout of corresponding variables should be done; 
! therefore, SAM needs an array of flags that mark the prognostic variables which
! only make sense when doprecip=.true. :

integer, parameter :: flag_precip(nmicro_fields) = (/0,0,0,0,0,1,1,1,1,1,1,1/)

! Sometimes, a cloud ice (or even cloud water) is allowed to be a subject of
! gravitational sedimentation, usually quite slow compared to the precipitation
! drops. SAM calls a special routine, ice_fall() that computes sedimentation of cloud ice.
! However, it is a rudiment from SAM's original single-moment microphysics.
! Instead, you may want to handle sedimentation of cloud water/ice yourself similarly
! to precipitation variables. In this case, set the index for falling cloud ice to -1, which
! means that no default ice mixing ration sedimentation is done. 

integer, parameter :: index_cloud_ice = -1   ! index for cloud ice (sedimentation)

#ifdef PNNL_STATS
integer :: &
  indx_qr = 8,  & ! Index for rain water mixing ratio (theta-l/qtog adv. budgets)
  indx_qi = 7,  & ! Index for ice mixing ratio (theta-l/qtog adv. budgets)
  indx_qs = 9,  & ! Index for snow mixing ratio (theta-l/qtog adv. budgets)
  indx_qg = 10    ! Index for graupel mixing ratio (theta-l/qtog adv. budgets)
#endif /* PNNL_STATS */

! The following arrays are needed to hold the turbulent surface and domain-top fluxes 
! for the microphysics prognostic variables:

real fluxbmk (nx, ny, 1:nmicro_fields) ! surface fluxes 
real fluxtmk (nx, ny, 1:nmicro_fields) ! top boundary fluxes 

! these arrays are needed for output statistics from advection and diffusion routines:

real mkwle(nz,1:nmicro_fields)  ! resolved vertical flux
real mkwsb(nz,1:nmicro_fields)  ! SGS vertical flux
real mkadv(nz,1:nmicro_fields)  ! tendency due to vertical advection
real mkdiff(nz,1:nmicro_fields)  ! tendency due to vertical diffusion

!------------------------------------------------------------------

! It would be quite inconvenient to work with the micro_field array itself. Besides,
! your original microphysics routines use some specific names for the prognostic variables
! that you don't wanna change. Therefore, you need to make aliases for prognostic variables. 
! You can use the aliases as you would ordinary arrays.

real ncn(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! CCN                            [#/kg] 
real nin(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! IN                             [#/kg]
real ncn_a(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! CCN activated                [#/kg]
real nin_a(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! IN activated                 [#/kg]
real qq(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! water vapor mixing ratio       [kg/kg]
real qc(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! cloud water mixing ratio       [kg/kg]
real qi(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! cloud ice mixing ratio         [kg/kg]
real qr(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! rain mixing ratio              [kg/kg]
real qs(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! snow mixing ratio              [kg/kg]
real qg(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! graupel mixin ratio            [kg/kg]
real nw(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! cloud water drop concentration  [#/kg]
real ni(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! cloud ice number concentration  [#/kg]
equivalence (ncn(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,1))
equivalence (nin(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,2))
equivalence (ncn_a(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,3))
equivalence (nin_a(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,4))
equivalence (qq(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,5))
equivalence (qc(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,6))
equivalence (qi(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,7))
equivalence (qr(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,8))
equivalence (qs(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,9))
equivalence (qg(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,10))
equivalence (nw(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,11))
equivalence (ni(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,12))

! Several of the above array names are different from original Phillips code in WRF
! due to name convention in SAM
!    this SAM module    Phillips WRF       
!           qq              qv     water vapor mixing ratio
!           qc              ql     cloud-liquid water mixing ratio
!

! You may also want to have some additional, diagnostic, arrays; 
! for example, total nonprecipitating cloud water, etc:

real qn(nx,ny,nzm)  ! cloud condensate (liquid + ice)

real qpsrc(nz)  ! source of precipitation microphysical processes
real qpfall(nz) ! source of precipitating water due to fall out in a given level
real qpevp(nz)  ! sink of precipitating water due to evaporation

real vrain, vsnow, vgrau, crain, csnow, cgrau  ! precomputed coefs for precip terminal velocity

!===============================================================vvvvvvvvvvvvv
	real ::                      &
         T_FRZ_HOM_DEGC,             &  !
         D_CI_MAX,                   &  !
         D_CW_MAX,                   &  !
         DT_VAP,                     &  !
         N_W_CM3_ENVIRONMENT,        &  !
         D_BAR_CI_TO_SNOW               !

	integer ::                   &
         ICE_NUMBER,                 &
         NKR,                        &
         FILL_HOLES_CW_WITH_SATADJ,  &
         BORROWING_SCHEME

	PARAMETER(                   &
         ICE_NUMBER = 1,             &
         T_FRZ_HOM_DEGC = -36.,      &
         D_CI_MAX = 500.E-6,         &
         D_BAR_CI_TO_SNOW = 200.e-6, &
         D_CW_MAX = 50.E-6,          &
         DT_VAP = 2.,                &
         NKR=33,                     & 
	 FILL_HOLES_CW_WITH_SATADJ=1,&
         BORROWING_SCHEME = 1,       &
         N_W_CM3_ENVIRONMENT = 0.1)

!  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , & 
!                                      ims,ime, jms,jme, kms,kme , &
!                                      its,ite, jts,jte, kts,kte

  INTEGER ::        &
          P_QI,     &
          P_QS,     &
          P_QG,     &
          P_FIRST_SCALAR

  REAL    ::  num_droplets, num_aero

  real, dimension(nx, ny, nzm) ::                 &
         th,        &  ! potential temperature (K)
         C_by_rho      ! number mixing ratio (#/kg[air]) of CN that would be activated 
                       !   at supersaturation of 1% (= normalised CN concentration 
                       !                                     divided by air density)

! Diagnistic variables
  real, dimension(nx, ny, nzm) ::                 &
         relhum,    &  ! output diagnostic for relative humidity (%)
         ni100um,   &  ! output diagnostic for number mixing ratio  (#/kg)
                       !   of cloud-ice particles > 70 microns
	 super_WRF, &  ! supersaturation (%) from linearised scheme
         s_FD_test, &  ! output diagnostic for supersaturation (%) from ultra-high
                       !   resolution finite difference scheme (switched off normally)
         s_FD_diff, &  ! output diagnostic for difference between
                       !   both supersaturation predictions (switched off normally)
         F_q,       &  ! output diagnostics for forcing terms in 
         F_T           !                        linearised supersaturation scheme

! look-up tables, set up at start of  simulation
  double precision,  dimension(150001) ::       &
         ESWT,      &  ! for saturated vapour pressures (mb) with respect to water
	 ESIT,      &  ! for saturated vapour pressures (mb) with respect to ice
         DESWT,     &  ! derivative of ESWT with respect to temperature
         DESIT,     &  ! derivative of ESIT with respect to temperature
         temp_K        ! physical temperature (K) of look-up table for ESWT, ESIT, etc.
 
! unused arrays for diagnostic output
! These are now stats output -dschanen 6 Aug 2008
  real, dimension(nx, ny, nzm) ::               &
         g_cond, g_evap, g_dep, g_subl, g_freeze, g_melt 

  real, dimension(nx, ny, nzm) ::                &
         s_vw_max,  & ! diagnostic output for supersaturation (%) applied in 
                      !    droplet nucleation (zero whenever there is no nucleation)
         d_bar_cw,  & ! output diagnostic for number-weighted mean droplet size (m)
         d_bar_ci     ! output diagnostic for number-weighted mean cloud-ice size (m)

  real, dimension(nx, ny, nzm) ::                &
         qrold,     &
         qsold,     &
         qgold,     &
         rho_wrf,   &  ! air density (kg/m3)
         pii,       &  ! ratio of physical temperature to potential temperature
!         p,         &  ! pressure (Pa)
         dz8w          ! spacing of adjacent vertical levels (m)

  real, dimension(nx, ny, nzm) ::                &
!         z,         &  !
         rw_vel,    &  !
         t_zero,    &  !
         qv_zero,   &  !
         p_zero        !

  real, dimension(nx, ny) ::  RAINNC

  real, dimension(nx, ny) ::  ht
  real, dimension(nzm)    ::  C_env_z
  real, dimension(nzm)    ::  C_by_rho_1d, super_WRF_1d    ! temporary arrays for subroutine arguments
 
  REAL :: 	    &
!	 dt_in,     &
 	 grav_wrf,  &
 	 Rair,      &
 	 rvapor,    &
 	 cp_wrf,    &
     	 XLS,       &
 	 XLV,       &
 	 XLF,       &
	 kaero,     &
         s_pc_cut_off, &
         z_pbl_top, &
         p_w,       &
         p_i

  REAL :: EP2,SVP1,SVP2,SVP3,SVPT0


!===============================================================^^^^^^^^^^^^^

CONTAINS

! Below are the required subroutines and functions that you need to fill in.

!----------------------------------------------------------------------
!!! Initialize microphysics:
!
! this one is guaranteed to be called by SAM at the beginning of each run, initial or restart:

subroutine micro_init()

  use vars
  use params

  double precision :: T_degC
  integer :: i,j,k, itc, itc_stop
!  double precision :: esatw_dble,esati_dble
!  external esatw_dble,esati_dble
  
  if(nrestart.eq.0) then

     micro_field = 0.
     do k=1,nzm
      qq(:,:,k) = q0(k)
     end do
     qn = 0.
     fluxbmk = 0.
     fluxtmk = 0.
  
  end if
  

  mkwle = 0.
  mkwsb = 0.
  mkadv = 0.
  mkdiff = 0.
 
  qpsrc = 0.
  qpevp = 0.

! your initialization calls are here. What is known at this point are the following 
! vertical profiles:

! temperature tabs0(1:nzm), 
! air density rho(1:nzm), 
! pressure pres(1:nzm), 
! water vapor (possibly supersaturated, so handle with caution) q0(1:nzm).
! Height of scalar levels is given by z(1:nzm), 
! height of leyer interfaces zi(1:nz).
! Thickness of each layer k is computed as dz*adz(k)
! So, the mass of each layer k per unit area is dz*adz(k)*rho(k)

! All the arrays above are available here through the module vars (vars.f90).

! Your additional initialization calls are placed below. 
! Remember that all your new files that contain actual microphysics subroutines 
! and functions should be added only to the microphysics directory.

! Initial CCN profile

  kaero    = 0.5   ! dimensionless exponent for supersaturation in power-law
                   !           activity spectrum (N_w = C s^k) for droplets
  Rair     = rgas  ! gas constant of air (R')
  grav_wrf = ggr   ! gravitational constant (~ 9.81 m/sec2)
  rvapor   = rv    ! gas constant of vapour (R_v)
  cp_wrf   = cp    ! specific heat capacity of air at constant pressure (~1005 J/kg/K)
  XLS=2.8336E6
!      HLTS=XLS
  XLV=2.5E6
!      HLTC=XLV
  XLF=3.336E5
!      HLTF=XLF
  p_w      = 3.5   ! shape parameter in droplet size distribution
  p_i      = 1.    ! shape parameter in cloud-ice size distribution

  ht(:,:)      = 0.    ! altitude of ground (m)
  z_pbl_top    = 0.
  s_pc_cut_off = 5.

  P_QI = 2    !  ??
  P_QS = 3
  P_QG = 4 

  DO k = 1,nzm
!    C_env_z(k) = 140.     ! 2mom high CCN
    C_env_z(k) = 80.     ! 2mom low CCN
    C_env_z(k) = C_env_z(k)*1.e6  ! #/cm3 to #/m3 
    print *, 'C = ', C_env_z(k)/1.e6, 'at ', z(k)/1000., ' km'

    C_by_rho(:,:,k) = C_env_z(k) / rho(k)

     do j=1,ny
       do i=1,nx
         tabs(i,j,k) = t(i,j,k) - gamaz(k) + fac_cond * (qc(i,j,k)+ qr(i,j,k)) &
   &                  + fac_sub * (qi(i,j,k)+ qs(i,j,k)+qg(i,j,k))

         t_zero(i,j,k)  = tabs(i,j,k)
         qv_zero(i,j,k) = qq(i,j,k) 
         p_zero(i,j,k)  = pres(k) * 100.   ! mb to Pa

       end do    ! i
     end do      ! j

  ENDDO

! At the start of the simulation, this code is needed to initialise lookup tables:-

	 T_degC = -100.

	do itc = 1, 150001
	 TEMP_K(itc) = 273.16 + T_degC

	 ESWT(itc) =  esatw_dble(TEMP_K(itc))
!mo	 ESWT(itc) = GGESW_DOUBLE(TEMP_K(itc))
	 if(T_degC < 0.) then
	   itc_stop = itc
           ESIT(itc) =  esati_dble(TEMP_K(itc))
!mo         ESIT(itc) = GGESI_DOUBLE(TEMP_K(itc))
	  else
           ESIT(itc) =  ESWT(itc)
	  endif

	  T_degC = T_degC + 0.001
	ENDDO
	DESWT(1) =  (ESWT(2)-ESWT(1))/0.001
 	DESWT(150001) =  (ESWT(150001)-ESWT(150000))/0.001
	DO itc = 2, 150000
		DESWT(itc) = (ESWT(itc+1)-ESWT(itc-1))/0.002
	ENDDO
	DESIT(1) =  (ESIT(2)-ESIT(1))/0.001
	DESIT(itc_stop) =  (ESIT(itc_stop)-ESIT(itc_stop - 1))/0.001
	DO itc = 2, itc_stop-1
		DESIT(itc) = (ESIT(itc+1)-ESIT(itc-1))/0.002
	ENDDO
	DO itc = itc_stop+1, 150001
		DESIT(itc) = DESWT(itc)
	ENDDO


!Preprocessor directives added to ensure pristine code -nielsenb UWM 30 July
!2008 
#ifdef CLUBB  
 if(docloud.or.doclubbnoninter.or.doclubb) call micro_diagnose()   ! leave this call here ! modified for CLUBB -dschanen UWM 16 May 2008
#else
 if(docloud) call micro_diagnose()   ! leave this call here
#endif

end subroutine micro_init

!----------------------------------------------------------------------
!!! fill-in the surface and top boundary fluxes here:
!
! Obviously, for liquid/ice water variables those fluxes are zero. They are not zero
! only for water vapor variable and, possibly, for CCN and IN if you have those.

subroutine micro_flux()

  use vars, only: fluxbq, fluxtq
  ! Added by dschanen UWM
#ifdef CLUBB
  use grid, only: doclubb
  if ( .not. doclubb ) then
    fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
    fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)
  else
    ! Add this in later
    fluxbmk(:,:,index_water_vapor) = 0.0
    fluxtmk(:,:,index_water_vapor) = 0.0
  end if
#else
  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
  fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)
#endif

end subroutine micro_flux

!----------------------------------------------------------------------
!!! compute local microphysics processes (beyond advection and SGS diffusion):
!
!  This is the place where the condensation/sublimation, accretion, coagulation, freezing,
!  melting, etc., are done, that is  all the microphysics processes except for the spatial 
!  transport and mixing.

! IMPORTANT: For consistancy, you need to use thermodynamic constants like specific heat,
! specific heat of condensation, gas constant, etc, the same as in file params.f90.
! Also, you should assume that the conservative thermodynamic variable during these
! proceses is the liquid/ice water static energy: t = tabs + gz - Lc (qc+qr) - Ls (qi+qs+qg) 
! It should not change during all of your point microphysical processes!

subroutine micro_proc()

! Preprocessor directives added -nielsenb UWM 30 July 2008
#ifdef CLUBB
      use vars, only: &
      w,dudt,dvdt,dwdt,tabs,t,pres,rho,nrestart,docloud,qv,gamaz, &
      doclubb, doclubbnoninter ! for CLUBB -dschanen UWM 16 May 2008
#else
      use vars, only: w,dudt,dvdt,dwdt,tabs,t,pres,rho,nrestart,docloud,qv,gamaz
#endif
      use grid, only: nc,dx,dy,dz,dt,nstep,z,adz,adzw, &
                       dimx1_w,dimx2_w, dimy1_w, dimy2_w 
      use params, only: fac_cond, fac_sub

      INTEGER, PARAMETER :: its=1, ite=nx, jts=1,     &
                            jte=ny, kts=1, kte=nzm    

! LOCAL VARIABLES

  real, dimension(nx, ny) ::  rain_wrf, snow_wrf, graupel_wrf
 
  real, dimension(nzm)    ::  qvz, qlz, niz, nin_az, nwz,  ncn_az,qrz, &
                                                   qiz, qsz, qgz, &
                                                  qrzold, qszold, &
                                                  qgzold,thz, tz, &
                                                     tothz, rhoz, &
                                                   orhoz, sqrhoz, &
                                                        prez, zz, &
                                         dzw, temp_tend, ql_tend, & 
					qi_tend, qs_tend, 	  & 
	qg_tend, qr_tend,qv_tend, nw_tend, ni_tend, s_vw_save, wz

 
  REAL    ::         pptrain, pptsnow, pptgraul, rhoe_s, DIS
  REAL                          :: qvmin=1.e-20,  N_0S, N_0G, N_0R
  LOGICAL			:: precip_flux
  INTEGER ::               i,j,k


	REAL :: ttt,QMIX(6),TSS(7)
	REAL :: DTL,TDT, RHOS
	REAL :: CRACS, CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,CIACR &
          ,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB(5),CGSUB(5),CREVP(5)  &
          ,CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,ES0 &
          ,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50,CPLC,CPLS,CPLF &
          ,TDIF,ESW00,ESI00
	REAL ::  PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC &
          ,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ, RHO_CLOUDICE, &
		GAMMA_RAT_i_mass, GAMMA_RAT_i_number, DBAR_FAC_i, &
		GAMMA_RAT_w_mass, GAMMA_RAT_w_number, DBAR_FAC_w, &  
           	N0_FAC_w, LAMBDA_FAC_w, &
 		N0_FAC_i, LAMBDA_FAC_i, GRAUPEL_FAC, RAIN_FAC, SNOW_FAC

	REAL ::  VCONR,VCONS,VCONG, VCONI, lambda_r, lambda_w, niz_min, nwz_min
	REAL :: PVI, PIL, PVL, GAMMA_K(6)
        REAL :: BETA_SUPER, Caero,  psi_1, psi_2 
        REAL :: ESW, G_factor, D_v, k_a, TSSNW, TSSNI
        REAL :: ESI, m_i0, N_in, Q_is, Q_ws, RIHENC, CIHENC, QIDEP
        REAL :: lambda_i, c_i, d_i, del_nw, P_SAT, QSW, &
		QSI, QVK, QLK, QIK, dtfr, TQIK, TQLK, TSV, DEL_QI, Q_ws_2, &
		s_pc_start, dtmlt, del_ni, &
		del_qw, del_qv, dtfine, lambda_i0, &
		rccn_cm(NKR), a_i, b_i, a_w, b_w, c_w, &
		a_r, b_r, a_s, b_s, a_g, b_g, rhograupel, rhosnow, rhowater

	REAL :: N_env_perm3(NKR,kts:kte), CCN_conc_env(kts:kte)
	REAL :: sum_mass_before, sum_mass_after, violation_pc, mass_rained_out
 	INTEGER :: ISTOP, i_fine, tot_fine
	INTRINSIC ABS, EXP

!===================================================================
!
!mo   dt=dt_in
! This is the surface air density - to be prescribed from observations :
   rhoe_s=1.20

! This is the switch governing code doing the precipitation fall-out
!   precip_flux = .true.
   precip_flux = .false.  ! precip is done by SAM

!
!    IF (P_QI .lt. P_FIRST_SCALAR .or. P_QS .lt. P_FIRST_SCALAR) THEN
!      CALL wrf_error_fatal ( 'module_mp_lin: Improper use of Lin et al scheme; no ice phase. Please chose another one.')
!   ENDIF

! only estimate rho_sfc for now

!=======================================================================
!		INITIALISATION 	
!=======================================================================
   
!mo   if((its .le. ite) .and. (jts .le. jte)) then
   	CALL SETUPM(dt, rhoe_s, grav_wrf, cp_wrf, &
	Rair, rvapor, XLS, XLV, XLF, rhowater, rhosnow, rhograupel, &
	CRACS,CSACR,CGACR,CGACS,ACCO,CSACW,CWACS, &
           CIACR,CRACI,CSACI,CGACW,CGACI,CRACW,CSSUB,CGSUB, &
!sk        CREVP,CGFR,CSMLT,CGMLT,CGWET,CRAUT, &
           CREVP,CGFR,CSMLT,CGMLT,CGWET,CRAUT, &
           QI0,QS0,QL0,ES0,CES0,C1BRG,C2BRG,RMI50,RMI40,RI50, &
           CPLC,CPLS,CPLF,TDIF,ESW00,ESI00, &
	PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC &
	  ,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ, &
	VCONR,VCONS,VCONG,VCONI, DTL,TDT, RHOS, BETA_SUPER, kaero, &
	 s_pc_start, D_v, k_a, N_0S, N_0G, N_0R)
!	   RHO_CLOUDICE = 200.

	RHO_CLOUDICE = 900  ! McFarquars recommendation (private comm. 2004; see MH96, JAS, page 2417)
   	c_i = RHO_CLOUDICE * PIE/6.    

	a_w = 1.19E6*100./4.; b_w = 2.; c_w = rhowater*PIE/6.;
! Ferriers values
	b_i = 0.6635;	
	a_i = (336./100.)*(100.**b_i) ;   ! = 70;  V_cm = "a_i" D_cm^"b_i"   =>  V*100 = "a_i" (D*100)^"b_i"  =>   V = ("a_i"/100) * 100^"b_i"  D^"b_i"

	d_i = 3.
! For terminal velocity of cloud-ice, assume C1f columns (density is about 900 kg/m3 at cloud-ice sizes) from PK97, Table 10.3b, (ultimately from Kajikawa 1976, J Met Soc Japan, Fig 10)

	b_i = 0.271 * 3.
	a_i = ((1.112/RHOS)**0.5)*2.53*(( c_i*1.e6)**0.271)  ! = 559

	b_r = 0.8
	a_r = (2115./100.)*(100.**b_r) 

	b_s = 0.25
	a_s = (152.93/100.)*(100.**b_s) 

	b_g = 0.5
	a_g = (4.*9.81*rhograupel/(3.*0.6*RHOS))**0.5
        if(RHOS > 2.) stop 17694
        if(RHOS < 0.5) stop 18694
	if(rhograupel < 100. .or. rhograupel > 950.) stop 65879

 ! see MH97: p_i = 1.
	N0_FAC_i = 1./exp(gammln(1.+p_i))
    	LAMBDA_FAC_i = exp(gammln(1.+d_i+p_i))*N0_FAC_i
    	DBAR_FAC_i = exp(gammln(2.+p_i))*N0_FAC_i
        GAMMA_RAT_i_mass = exp(gammln(1.+ p_i + d_i + b_i))/exp(gammln(1.+ p_i + d_i))
        GAMMA_RAT_i_number = exp(gammln(1.+b_i+p_i))*N0_FAC_i

! see PK97: p_w = 2.
 	N0_FAC_w = 1./exp(gammln(1.+p_w))
 	LAMBDA_FAC_w = exp(gammln(4.+p_w))*N0_FAC_w
   	DBAR_FAC_w = exp(gammln(2.+p_w))*N0_FAC_w
        GAMMA_RAT_w_mass = exp(gammln(4.+p_w+ b_w))/exp(gammln(4.+p_w))
        GAMMA_RAT_w_number = exp(gammln(1.+b_w+p_w))*N0_FAC_w

	GRAUPEL_FAC = exp(gammln(2.5+0.5*b_g))
	RAIN_FAC = exp(gammln(2.5+0.5*b_r))
	SNOW_FAC = exp(gammln(2.5+0.5*b_s))
	lambda_i0 = DBAR_FAC_i/D_BAR_CI_TO_SNOW	! when d_bar > D_bar_0 = 200 um we convert to snow
	if(1./lambda_i0 > D_CI_MAX) stop "1439876" ! UWM

	DIS = SQRT(1./(p_w+1.))
        if(DIS < 0.) stop "543657" ! UWM
	if(DIS > 10.) stop 65789

!mo   endif
 
   j_loop:  DO j = jts, jte
   i_loop:  DO i = its, ite

!- write data from 3-D to 1-D
!
   DO k = kts, kte
      qvz(k)=qq(i,j,k)
      qlz(k)=qc(i,j,k)
      nwz(k)=nw(i,j,k)
      ncn_az(k) = ncn_a(i,j,k)
      qrz(k)=qr(i,j,k)
!mo      qrzold(k)=qrold(i,j,k)
!mo      thz(k)=th(i,j,k)
!
      rhoz(k)=rho(k)
      orhoz(k)=1./rhoz(k)
      prez(k)=pres(K)*100.              ! mb to Pa
      sqrhoz(k)=sqrt(rhoe_s*orhoz(k))
!mo      tothz(k)=pii(i,j,k)
      zz(k)=z(k)
      dzw(k)=dz*adz(k)
	wz(k) = (w(i,j,k) + w(i,j,k+1))/2.

     tz(k) = t(i,j,k) - gamaz(k) + fac_cond * (qc(i,j,k)+ qr(i,j,k)) &
   &                  + fac_sub * (qi(i,j,k)+ qs(i,j,k)+qg(i,j,k))
   END DO

   C_by_rho_1d(:) = C_by_rho(i,j,:)
 CALL CCN_activity_spectrum(tz, zz, z_pbl_top, s_pc_cut_off, kaero, C_by_rho_1d(:), rhoz, PIE, N_env_perm3, &
	CCN_conc_env,  rccn_cm, kts, kte)
   C_by_rho(i,j,:) = C_by_rho_1d(:)

!   print *,rho at sfc::, rhoz(1) - zz(1) * (rhoz(2) - rhoz(1))/(zz(2) - zz(1)) 
!    print *,rho at 1::, rhoz(1),  zz(1)
   IF (P_QI .ge. P_FIRST_SCALAR) THEN
      DO k = kts, kte
         qiz(k)=qi(i,j,k)
         niz(k)=ni(i,j,k)
     	 nin_az(k) = nin_a(i,j,k)
       ENDDO
   ELSE
      DO k = kts, kte
         qiz(k)=0.
         niz(k)=0.
      ENDDO
   ENDIF

   IF (P_QS .ge. P_FIRST_SCALAR) THEN
      DO k = kts, kte
         qsz(k)=qs(i,j,k)
         qszold(k)=qsold(i,j,k)
      ENDDO
   ELSE
      DO k = kts, kte
         qsz(k)=0.
         qszold(k)=0.
      ENDDO
   ENDIF

   IF (P_QG .ge. P_FIRST_SCALAR) THEN
      DO k = kts, kte
         qgz(k)=qg(i,j,k)
         qgzold(k)=qgold(i,j,k)
      ENDDO
   ELSE
      DO k = kts, kte
         qgz(k)=0.
         qgzold(k)=0.
      ENDDO
   ENDIF
!
   pptrain=0.
   pptsnow=0.
   pptgraul=0.

!===============================================================================
!		DEAL WITH ANY NEGATIVE MIXING RATIOS FROM DYNAMICAL CORE
!===============================================================================

!mo   tz(:) = thz(:) * tothz(:)
   if(BORROWING_SCHEME == 1) then
	if(kts /= 1) stop 1234
	call fillz(kte, qlz, rhoz, dzw)
	call fillz(kte, nwz, rhoz, dzw)
	call fillz(kte, ncn_az, rhoz, dzw)
	call fillz(kte, qiz, rhoz, dzw)
	call fillz(kte, niz, rhoz, dzw)
	call fillz(kte, nin_az, rhoz, dzw)
	call fillz(kte, qrz, rhoz, dzw)
	call fillz(kte, qsz, rhoz, dzw)
	call fillz(kte, qgz, rhoz, dzw)
	call fillz(kte, qvz, rhoz, dzw)
	do k = kts, kte
		if(C_by_rho(i,j,k) < 0.) stop "154879"
 	enddo  
   endif


   if(FILL_HOLES_CW_WITH_SATADJ == 1) then 
   	do k = kts, kte
		if(C_by_rho(i,j,k) < 0.) stop "154879"
		ttt = tz(k)
		P_SAT= prez(k)/100.
		QVK = qvz(k)
  		ESW = REAL(ESATW_int(DBLE(ttt), ESWT, temp_K))
  	    	QSW = EPS*ESW/(P_SAT-ESW)

! is this a grid-point of cloud that has been rendered negative by the advection scheme ? If so,
!	fill the hole with the saturation adjustment scheme.

		if(QVK/QSW - 1. > s_pc_start/100. .and.  qlz(k) < -1.e-20 ) then 
			if(nwz(k) < 0.) nwz(k) = 0.
			nwz(k) = nwz(k) + C_by_rho(i,j,k) * (((QVK/QSW - 1.)*100.)**kaero)
			ncn_az(k) = ncn_az(k) + C_by_rho(i,j,k) * (((QVK/QSW - 1.)*100.)**kaero)
			call condensation_of_cloudwater(s_pc_start, tz(k), & 
				prez(k), qvz(k), qlz(k), nwz(k), &
				EPS, ESIT, DESIT, ESWT, DESWT, temp_K, CPLC)
		endif
	   enddo
!mo	   thz(:) = tz(:)/tothz(:) 
  endif
	
 
  do k = kts, kte
         qlz(k)=amax1( 0.0,qlz(k) )
         qiz(k)=amax1( 0.0,qiz(k) )
         niz(k)=amax1( 0.0,niz(k) )
         nin_az(k)=amax1( 0.0,nin_az(k) )
         nwz(k)=amax1( 0.,nwz(k) )
         ncn_az(k)=amax1( 0.0,ncn_az(k) )
	 qvz(k)=amax1( qvmin,qvz(k) )
         qsz(k)=amax1( 0.0,qsz(k) )
         qrz(k)=amax1( 0.0,qrz(k) )
         qgz(k)=amax1( 0.0,qgz(k) )
   enddo


!
!===============================================================================
!		THRESHOLDING OF NUMBERS
!===============================================================================
!



    do k = kts, kte
	call threshold_number(d_bar_ci(i,j,k), qiz(k), niz(k), nin_az(k), D_CI_MAX,  d_i, c_i, LAMBDA_FAC_i, DBAR_FAC_i)
	call threshold_number(d_bar_cw(i,j,k), qlz(k), nwz(k), ncn_az(k), D_CW_MAX,  3., c_w, LAMBDA_FAC_w, DBAR_FAC_w)
    enddo

!========================================================================================
!		PRECIPITATION FLUX (from Shuhua Chen, adapted by VTJP) 	
!========================================================================================
    
     if(precip_flux) then
   	CALL PRECIP_FALLOUT(qrz, qgz, qsz, qiz, qlz, niz, nwz, dt, VCONR, VCONS, &
		VCONG, VCONI, a_i, b_i, c_i, d_i, a_w, b_w, c_w, LAMBDA_FAC_w, &
                LAMBDA_FAC_i, GAMMA_RAT_w_mass, GAMMA_RAT_i_mass, &
		GAMMA_RAT_w_number, GAMMA_RAT_i_number, ht(i,j), RHOS, &
                rhoz, dzw, zz, pptrain, pptsnow, pptgraul, kts, kte)
		mass_rained_out = pptrain + pptsnow + pptgraul
     endif
!
!===============================================================================
!		THRESHOLDING OF NUMBERS
!===============================================================================
!



    do k = kts, kte
	call threshold_number(d_bar_ci(i,j,k), qiz(k), niz(k), nin_az(k), D_CI_MAX,  d_i, c_i, LAMBDA_FAC_i, DBAR_FAC_i)
	call threshold_number(d_bar_cw(i,j,k), qlz(k), nwz(k), ncn_az(k), D_CW_MAX,  3., c_w, LAMBDA_FAC_w, DBAR_FAC_w)
    enddo

!===============================================================================
!		MICROPHYSICS:: INTERACTIONS BETWEEN HYDROMETEOR SPECIES (IMICRO)
!===============================================================================

	temp_tend = 0.     
	qv_tend = 0.
	ql_tend = 0.
	qs_tend = 0.
	qg_tend = 0.
	qr_tend = 0.
	qi_tend = 0.
	nw_tend = 0.
	ni_tend = 0.

      do k = kts, kte

	QMIX(1) = qvz(k)
	QMIX(2) = qlz(k)
	QMIX(3) = qiz(k)
 	QMIX(4) = qsz(k)
	QMIX(5) = qgz(k)
	QMIX(6) = qrz(k)
  	TSSNW = 0.
	TSSNI = 0.
	ttt = tz(k)
	
	GAMMA_K = 0.
	TSS = 0.
	CALL IMICRO(dt, ttt,prez(k)/100.,QMIX,TSS, nwz(k), niz(k), TSSNW, TSSNI, rhoz(k), &
	   TDT, RHOS, CRACS, CSACR,CGACR,CGACS,ACCO,CSACW &
          ,CRACI,CIACR,CSACI,CGACW,CGACI,CRACW,CSSUB,CGSUB,CREVP  &
          ,CGFR,CSMLT,CGMLT,CGWET,CRAUT,QI0,QS0,QL0,ES0 &
          ,CES0,C1BRG,C2BRG,RMI50,RMI40,ESWT,ESIT,DESWT,DESIT &
          ,temp_K,HLTS,HLTC,HLTF,CH2O,CICE,TICE,CP,EPS, &
	   VCONR,VCONS,VCONG,GAMMA_K, d_i, RHO_CLOUDICE, c_i, c_w, &
	   LAMBDA_FAC_w, LAMBDA_FAC_i, p_w, p_i, DBAR_FAC_w, lambda_i0, &
	 DIS)
	

	G_COND(i,j,k) = GAMMA_K(1)
	G_EVAP(i,j,k) = GAMMA_K(2)
	G_DEP(i,j,k) = GAMMA_K(3)
 	G_SUBL(i,j,k) = GAMMA_K(4)
 	G_FREEZE(i,j,k) = GAMMA_K(5)
  	G_MELT(i,j,k) = GAMMA_K(6)


	temp_tend(k) = TSS(7)
	qv_tend(k) = TSS(1)
	ql_tend(k) = TSS(2)
	qi_tend(k) = TSS(3)
	qs_tend(k) = TSS(4)
	qg_tend(k) = TSS(5)
	qr_tend(k) = TSS(6) 
	nw_tend(k) = TSSNW
	ni_tend(k) = TSSNI

   enddo
  

	
     do k = kts, kte
	if(qvz(k) + dt * qv_tend(k) < 0.) print *,'WARNING - AFTER IMICRO: qv = ', qvz(k) + dt * qv_tend(k)
	qvz(k) = amax1(0.0, qvz(k) + dt * qv_tend(k))
	if(qlz(k) + dt * ql_tend(k) < 0.) print *,'WARNING - AFTER IMICRO: ql = ', qlz(k) + dt * ql_tend(k)
	qlz(k) = amax1(0.0, qlz(k) + dt * ql_tend(k))
	if(qlz(k) > 0.) then
		nwz(k) = amax1(0.0, nwz(k) + dt*nw_tend(k))
		if(nwz(k) == 0.) then
			tz(k) = tz(k) - (XLV/CP)*qlz(k)
			qvz(k) = qvz(k) + qlz(k)
			qlz(k) = 0.
		endif
	else
		nwz(k) = 0.
	endif
	if(qiz(k) + dt * qi_tend(k) < 0.) print *,'WARNING - AFTER IMICRO: qi = ', qiz(k) + dt * qi_tend(k)
	qiz(k) = amax1(0.0, qiz(k) + dt * qi_tend(k))
	if(qiz(k) > 0.) then
		niz(k) = amax1(0.0, niz(k) + dt*ni_tend(k))
		if(niz(k) == 0.) then
			tz(k) = tz(k) - (XLS/CP)*qiz(k)
			qvz(k) = qvz(k) + qiz(k)
			qiz(k) = 0.
		endif
	else
		niz(k) = 0.
	endif
	if(qsz(k) + dt * qs_tend(k) < 0.) print *,'WARNING - AFTER IMICRO: qs = ', qsz(k) + dt * qs_tend(k)
	qsz(k) = amax1(0.0, qsz(k) + dt * qs_tend(k))
	if(qgz(k) + dt * qg_tend(k) < 0.) print *,'WARNING - AFTER IMICRO: qg = ', qgz(k) + dt * qg_tend(k)
	qgz(k) = amax1(0.0, qgz(k) + dt * qg_tend(k))
	if(qrz(k) + dt * qr_tend(k) < 0.) print *,'WARNING - AFTER IMICRO: qr = ', qrz(k) + dt * qr_tend(k)
	qrz(k) = amax1(0.0, qrz(k) + dt * qr_tend(k))
	tz(k) = tz(k) + dt * temp_tend(k)
!mo        thz(k) = tz(k)/tothz(k)
     enddo

!===============================================================================
!		MICROPHYSICS::  NUCLEATION AND GROWTH OF CLOUD PARTICLES 
!                               (replaces saturation adjustment)
!===============================================================================

!===============================================================================
!		THRESHOLDING OF NUMBERS
!===============================================================================

    do k = kts, kte
	call threshold_number(d_bar_ci(i,j,k), qiz(k), niz(k), nin_az(k), D_CI_MAX,  d_i, c_i, LAMBDA_FAC_i, DBAR_FAC_i)
	call threshold_number(d_bar_cw(i,j,k), qlz(k), nwz(k), ncn_az(k), D_CW_MAX,  3., c_w, LAMBDA_FAC_w, DBAR_FAC_w)
                       if(qiz(k) > 1. .or. qlz(k) > 1.) stop 125
    enddo

!===========================================
!		RESET ACTIVATED CCN AND IN OUTSIDE CLOUD
!===========================================

     do k = kts, kte
	call reset_aerosols_in_environment(qlz(k), nwz(k), ncn_az(k),d_bar_cw(i,j,k), qiz(k), niz(k), &
                nin_az(k), d_bar_ci(i,j,k), C_by_rho(i,j,k), C_env_z(k), &
		rhoz(k), tz(k), N_W_CM3_ENVIRONMENT, qvz(k), XLV, XLS, CP)
     enddo


     do k = kts, kte
		ttt = tz(k)
		P_SAT= prez(k)/100.
		QVK = qvz(k)
  		ESW = REAL(esatw_int(DBLE(ttt), ESWT, temp_K))
      		QSW = EPS*ESW/(P_SAT-ESW)
		s_vw_save(k) = QVK/QSW - 1.      
      enddo


!
!===========================================
!		HOMOGENEOUS FREEZING/MELTING 
!===========================================
!
 
    do k = kts, kte
      super_WRF_1d(:) = super_WRF(i,j,:)
	call homogeneous_freezing( tz(k), qiz(k),  qlz(k), niz(k),  nwz(k),  ncn_az(k), qvz(k), &
		LAMBDA_FAC_w, p_w, XLF, XLV, CP, T_FRZ_HOM_DEGC, PIE, wz(k), super_WRF_1d(:), tz(:), kts, kte, k, num_droplets, rhoz(k))
 	call melting_cloudice(tz(k), qlz(k), nwz(k), qiz(k), niz(k),  XLF, CP)
      super_WRF(i,j,:) =  super_WRF_1d(:)
!mo       thz(k) = tz(k)/tothz(k)
                       if(qiz(k) > 1. .or. qlz(k) > 1.) stop 125

    enddo


!
!===============================================================================
!		THRESHOLDING OF NUMBERS
!===============================================================================
!

    do k = kts, kte
	call threshold_number(d_bar_ci(i,j,k), qiz(k), niz(k), nin_az(k), D_CI_MAX,  d_i, c_i, LAMBDA_FAC_i, DBAR_FAC_i)
	call threshold_number(d_bar_cw(i,j,k), qlz(k), nwz(k), ncn_az(k), D_CW_MAX,  3., c_w, LAMBDA_FAC_w, DBAR_FAC_w)
                       if(qiz(k) > 1. .or. qlz(k) > 1.) stop 125
    enddo


! Neglect ventilation coefficient for now

    tot_fine = INT(dt/dt_vap)
    if(tot_fine < 1) tot_fine = 1


 
     		do k = kts, kte


			super_WRF(i,j,k) = 0.
                        if(qiz(k) > 1. .or. qlz(k) > 1.) then 
                              stop 137
                        endif
			call condense_and_sublime_driver(qvz(k), Tz(k), prez(k), rhoz(k), qv_zero(i,j,k), T_zero(i,j,k), &
					p_zero(i,j,k), F_q(i,j,k), F_T(i,j,k),  &
					qlz(k), nwz(k), ncn_az(k), qiz(k), niz(k), nin_az(k), qrz(k), qsz(k), qgz(k), &
					rvapor, a_i, b_i, c_i, d_i, a_w, b_w, c_w, a_r, b_r, rhowater, &
					a_s, b_s, rhosnow, a_g, b_g, rhograupel, LAMBDA_FAC_w, LAMBDA_FAC_i, DBAR_FAC_w, &
					DBAR_FAC_i, N0_FAC_w, N0_FAC_i, p_w, p_i, GRAUPEL_FAC, RAIN_FAC, SNOW_FAC, CP, &
					PIE, EPS, XLS, XLV, Rair, ESIT, DESIT, ESWT, DESWT, temp_K, &
					dt, super_WRF(i,j,k), s_FD_test(i,j,k),  s_FD_diff(i,j,k), N_0S, N_0G, N_0R)

                        if(qiz(k) > 1. .or. qlz(k) > 1.) then 
                              stop 138
                        endif
			call homogeneous_aerosol_freezing(tz(k), prez(k), rhoz(k),  wz(k), &
				qvz(k), qiz(k), niz(k), ncn_az(k), &
				XLF, XLS, CP, PIE, EPS, ESIT, DESIT, ESWT, DESWT, temp_K,  &
				N_env_perm3(:,k), CCN_conc_env(k), rccn_cm, num_aero)
                      if(qiz(k) > 1. .or. qlz(k) > 1.) stop 139
!mo			thz(k) = tz(k)/tothz(k) 
		enddo

 
  
!
!===========================================
!		DROPLET and CRYSTAL NUCLEATION (supersaturation predicted at end of timestep)
!
!===========================================
!

      do k = kts, kte
	 caero = C_by_rho(i,j,k) *  rhoz(k);
	if(caero == 0.) stop 911
 	 call droplet_nucleation(s_vw_max(i,j,k), s_vw_save(k), s_pc_start, s_pc_cut_off, &
		tz(k), prez(k), rhoz(k), &
		 wz(k), qvz(k), qlz(k), nwz(k), ncn_az(k), &
		caero, kaero, XLV, CP, EPS, Rair, rvapor, grav, BETA_SUPER, &
		PIE, ESIT, DESIT, ESWT, DESWT, temp_K, c_w, N_W_CM3_ENVIRONMENT)
                if(qiz(k) > 1. .or. qlz(k) > 1.) stop 140

     enddo
     do k = kts, kte
              if(qiz(k) > 1. .or. qlz(k) > 1.) stop 141
			call primary_ice_nucleation(tz(k), prez(k), rhoz(k),  wz(k), & 
				qvz(k),qlz(k), qiz(k), nwz(k), niz(k), nin_az(k), &
				XLF, XLS, XLV, CP, PIE, EPS, ESIT, DESIT, ESWT, DESWT, &
				temp_K, rvapor, dt, p_w, LAMBDA_FAC_W, a_w, b_w, c_w, N0_FAC_w)
              if(qiz(k) > 1. .or. qlz(k) > 1.) stop 142
     enddo
 

!
!===============================================================================
!		THRESHOLDING OF NUMBERS
!===============================================================================
!

    do k = kts, kte
	call threshold_number(d_bar_ci(i,j,k), qiz(k), niz(k), nin_az(k), D_CI_MAX,  d_i, c_i, LAMBDA_FAC_i, DBAR_FAC_i)
	call threshold_number(d_bar_cw(i,j,k), qlz(k), nwz(k), ncn_az(k), D_CW_MAX,  3., c_w, LAMBDA_FAC_w, DBAR_FAC_w)
                   if(qiz(k) > 1. .or. qlz(k) > 1.) stop 143
    enddo


!
!===========================================
!		HOMOGENEOUS FREEZING/MELTING 
!===========================================
!

 
    do k = kts, kte
      super_WRF_1d(:) = super_WRF(i,j,:)
	call homogeneous_freezing( tz(k), qiz(k),  qlz(k), niz(k),  nwz(k),  ncn_az(k), qvz(k), &
		LAMBDA_FAC_w, p_w,   XLF, XLV, CP, T_FRZ_HOM_DEGC, c_w, wz(k), super_WRF_1d(:), &
		tz(:), kts, kte, k, num_droplets, rhoz(k) )
 	call melting_cloudice(tz(k), qlz(k), nwz(k), qiz(k), niz(k),  XLF, CP)
	 call zero_trace_cloud_in_subsat_env(qlz(k), nwz(k), qiz(k), niz(k),  qvz(k), tz(k), prez(k), &
			CP, EPS, XLS, XLV, ESIT, DESIT, ESWT, DESWT, temp_K)
      super_WRF(i,j,:) = super_WRF_1d(:) 

!mo        thz(k) = tz(k)/tothz(k)
                   if(qiz(k) > 1. .or. qlz(k) > 1.) stop 144
    enddo

!
!===============================================================================
!		THRESHOLDING OF NUMBERS
!===============================================================================
!


    do k = kts, kte
	call threshold_number(d_bar_ci(i,j,k), qiz(k), niz(k), nin_az(k), D_CI_MAX,  d_i, c_i, LAMBDA_FAC_i, DBAR_FAC_i)
	call threshold_number(d_bar_cw(i,j,k), qlz(k), nwz(k), ncn_az(k), D_CW_MAX,  3., c_w, LAMBDA_FAC_w, DBAR_FAC_w)
                  if(qiz(k) > 1. .or. qlz(k) > 1.) stop 145
    enddo

 
!
!===========================================
!		RESET ACTIVATED CCN AND IN OUTSIDE CLOUD
!===========================================
!

   do k = kts, kte
	call reset_aerosols_in_environment(qlz(k), nwz(k), ncn_az(k),d_bar_cw(i,j,k), qiz(k), niz(k), &
                nin_az(k), d_bar_ci(i,j,k), C_by_rho(i,j,k), C_env_z(k), &
		rhoz(k), tz(k), N_W_CM3_ENVIRONMENT, qvz(k), XLV, XLS, CP)
     enddo

  do k = kts, kte
	if(qlz(k) == 0. .and. d_bar_cw(i,j,k) > 0.) stop "147645"
  enddo

!=========================================================================================
! 		OUTPUT 
!=========================================================================================
 
     do k = kts, kte
		ttt = tz(k)
		P_SAT= prez(k)/100.
		QVK = qvz(k)
  		ESW = REAL(esatw_int(DBLE(ttt), ESWT, temp_K))
      		QSW = EPS*ESW/(P_SAT-ESW)
		relhum(i,j,k) = 100.*QVK/QSW   
		if(qiz(k) > 0. .and. niz(k) > 0.) then
			lambda_i = LAMBDA_FAC_i*c_i*niz(k)/qiz(k) ! CORRECT
			lambda_i = lambda_i**(1./d_i) ! CORRECT
			ni100um(i,j,k) = niz(k)*gammq_frac(1.+p_i,70.e-6*lambda_i)
		else
			ni100um(i,j,k) = 0.
		endif
		if(C_by_rho(i,j,k) < 0.) stop 53221
      enddo

!
! Precipitation from cloud microphysics -- only for one time step
!
! unit is transferred from m to mm

   rain_wrf(i,j)=pptrain
   snow_wrf(i,j)=pptsnow
   graupel_wrf(i,j)=pptgraul

   RAINNC(i,j)=RAINNC(i,j) + pptrain + pptsnow + pptgraul
 
!
!- update data from 1-D back to 3-D
!
!
   DO k = kts, kte
      qq(i,j,k)=qvz(k)
      qc(i,j,k)=qlz(k)
      nw(i,j,k)=nwz(k)
      ncn_a(i,j,k)=ncn_az(k)
      qr(i,j,k)=qrz(k)
!mo      th(i,j,k)=thz(k)
   END DO
!
   IF (P_QI .ge. P_FIRST_SCALAR) THEN
      DO k = kts, kte
         qi(i,j,k)=qiz(k)
         ni(i,j,k)=niz(k)
         nin_a(i,j,k)=nin_az(k)
      ENDDO
   ENDIF

   IF (P_QS .ge. P_FIRST_SCALAR) THEN
      DO k = kts, kte
         qs(i,j,k)=qsz(k)
      ENDDO
   ENDIF

   IF (P_QG .ge. P_FIRST_SCALAR) THEN
      DO k = kts, kte
         qg(i,j,k)=qgz(k)
      ENDDO
   ENDIF

   ENDDO i_loop
   ENDDO j_loop



!===================================================================
! Added preprocessor directives - nielsenb UWM 30 July 2008 
! modified for CLUBB -dschanen UWM 16 May 2008
#ifdef CLUBB  
   if(docloud.or.doclubbnoninter.or.doclubb) then
       call micro_diagnose()  
   end if
#else
   if (docloud)  call micro_diagnose()   ! leave this line here
#endif

! Save values for the next time step

      DO k = kts, kte
      do j = jts, jte
      do i = its, ite
         tabs(i,j,k) = t(i,j,k) - gamaz(k) + fac_cond * (qc(i,j,k)+ qr(i,j,k)) &
   &                  + fac_sub * (qi(i,j,k)+ qs(i,j,k)+qg(i,j,k))
         t_zero(i,j,k)  = tabs(i,j,k) 
         qv_zero(i,j,k) = qq(i,j,k) 
         p_zero(i,j,k)  = pres(k) * 100.   ! mb to Pa
      end do
      end do
      END DO

end subroutine micro_proc

!----------------------------------------------------------------------
!!! Diagnose arrays necessary for dynamical core and radiation:
!
! This is the pace where the microphysics field that SAM actually cares about
! are diagnosed. You need to compute all the arrays on the left-hand-side in the loop below
! for SAM dynamical core to see your microphysics (that is to see the cloud and precipitation).

subroutine micro_diagnose( )
 
   use vars, &
       only: qv, qcl, qci, qpl, qpi ! -dschanen UWM made this ``only:'' 20 May 08

   real omn, omp
   integer i,j,k

   do k=1,nzm
    do j=1,ny
     do i=1,nx
       qv(i,j,k) = qq(i,j,k)
!       qv(i,j,k) = qt(i,j,k) - (qc(i,j,k)+qr(i,j,k)+qi(i,j,k)+qs(i,j,k)+qg(i,j,k))
       qcl(i,j,k) = qc(i,j,k)
       qci(i,j,k) = qi(i,j,k)
       qpl(i,j,k) = qr(i,j,k)
       qpi(i,j,k) = qs(i,j,k) + qg(i,j,k)

     end do
    end do
   end do
       


end subroutine micro_diagnose

#ifdef CLUBB
!---------------------------------------------------------------------
subroutine micro_update( )
! Description: Added by dschanen UWM

  implicit none

  call micro_diagnose()
  return
end subroutine micro_update

!---------------------------------------------------------------------
subroutine micro_adjust( new_qv, new_qc )
! Description:
! Adjust vapor and liquid water.
! Microphysical variables are stored separately in
!    SAM's dynamics + CLUBB ( e.g. qv, qcl, qci) and
!    SAM's microphysics. (e.g. qq and qc).
! This subroutine stores values of qv, qcl updated by CLUBB
!   in the double-moment microphysical variables qq and qc.
!
! dschanen UWM 20 May 2008
!---------------------------------------------------------------------

  implicit none

  real, dimension(nx,ny,nzm), intent(in) :: &
  new_qv, & ! Water vapor mixing ratio that has been adjusted by CLUBB [kg/kg]
  new_qc    ! Cloud water mixing ratio that has been adjusted by CLUBB [kg/kg]

  qq(1:nx,1:ny,1:nzm) = new_qv ! Vapor
  qc(1:nx,1:ny,1:nzm) = new_qc ! Liquid

  return
end subroutine micro_adjust
#endif

!----------------------------------------------------------------------
!!! functions to compute terminal velocity for precipitating variables:
!
! you need to supply your own functions functions to compute terminal velocity 
! for all of your precipitating prognostic variables. Note that all functions should
! compute vertical velocity given two microphysics parameters var1, var2, 
! and temperature, and water vapor (single values, not arrays). Also, for
! bin-microphysics implementation, there is a fifth variable with the type of
! integer that can be used for bin index. Var1 and var2 
! are some microphysics variables like water content and concentration.
! IMPORTANT: Don't change the number of arguments or their meaning!

real function term_vel_qr(qrr,dummy,tabs,rho,ind)
  real, intent(in) ::  qrr, dummy, tabs, rho
  integer ind ! dummy variable
  term_vel_qr = vrain*(rho*qrr)**crain
 end function term_vel_qr

real function term_vel_qs(qss,dummy,tabs,rho,ind)
  real, intent(in) ::  qss, dummy, tabs, rho
  integer ind ! dummy variable
  term_vel_qs = vsnow*(rho*qss)**csnow
end function term_vel_qs

real function term_vel_qg(qgg,dummy,tabs,rho,ind)
  real, intent(in) ::  qgg, dummy, tabs, rho
  integer ind ! dummy variable
  term_vel_qg = vgrau*(rho*qgg)**cgrau 
end function term_vel_qg

!----------------------------------------------------------------------
!!! compute sedimentation 
!
!  The purpose of this subroutine is to prepare variables needed to call
! the precip_all() for each of the falling hydrometeor varibles

subroutine micro_precip_fall()
  
use vars
use params, only : pi

! before calling precip_fall() for each of falling prognostic variables,
! you need to set hydro_type and omega(:,:,:) variables.
! hydro_type can have four values:
! 0 - variable is liquid water mixing ratio
! 1 - hydrometeor is ice mixing ratio
! 2 - hydrometeor is mixture-of-liquid-and-ice mixing ratio. (As in original SAM microphysics).
! 3 - variable is not mixing ratio, but, for example, rain drop concentration
! OMEGA(:,:,:) is used only for hydro_type=2, and is the fraction of liquid phase (0-1).
! for our hypothetical case, there is no mixed hydrometeor, so omega is not actually used.
! In default SAM microphysics, omega is a mass partition between liquid and ice phases.

  integer hydro_type
  real omega(nx,ny,nzm) 
  real dummy  !
  integer ind ! variable that is reserved for bin-microphysics use (bin index).
  integer i,j,k
!  real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
!  real f0(nzm),df0(nzm)

  crain = b_rain / 4.
  csnow = b_snow / 4.
  cgrau = b_grau / 4.
  vrain = a_rain * gamr3 / 6. / (pi * rhor * nzeror) ** crain
  vsnow = a_snow * gams3 / 6. / (pi * rhos * nzeros) ** csnow
  vgrau = a_grau * gamg3 / 6. / (pi * rhog * nzerog) ** cgrau

! Initialize arrays that accumulate surface precipitation flux

 if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
   do j=1,ny
    do i=1,nx
     precsfc(i,j)=0.
    end do
   end do
   do k=1,nzm
    precflux(k) = 0.
   end do
 end if

 do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
    qpfall(k)=0.
    tlat(k) = 0.
 end do
   
! Compute sedimentation of falling variables:

 hydro_type=0
 call precip_fall(qr,dummy, term_vel_qr, hydro_type, omega, ind)
 hydro_type=1
 call precip_fall(qs,dummy, term_vel_qs, hydro_type, omega, ind)
 hydro_type=1
 call precip_fall(qg,dummy, term_vel_qg, hydro_type, omega, ind)

! hydro_type=3
! call precip_fall(Nr, term_vel_Nr, hydro_type, omega, ind)
! hydro_type=3
! call precip_fall(Ns, term_vel_Ns, hydro_type, omega, ind)
! hydro_type=3
! call precip_fall(Ng, term_vel_Ng, hydro_type, omega, ind)

! compute surface precipitation area fraction statistics
 do j=1,ny
   do i=1,nx
     if((qr(i,j,1)+qs(i,j,1)+qg(i,j,1).gt.1.e-6)) s_ar=s_ar+dtfactor
   end do
 end do

!=================================vvvvvvvvvvvvvvvvvvvvvvvvv
! if(dostatis) then
!   do k=1,nzm
!     do j=dimy1_s,dimy2_s
!       do i=dimx1_s,dimx2_s
!          df(i,j,k) = t(i,j,k)
!       end do
!     end do
!   end do
!   call stat_varscalar(t,df,f0,df0,t2leprec)
!   call setvalue(twleprec,nzm,0.)
!   call stat_sw2(t,df,twleprec)
! endif

! do j=1,ny
!   do i=1,nx
!     if(qp(i,j,1).gt.1.e-6) s_ar=s_ar+dtfactor
!   end do
! end do
!=================================^^^^^^^^^^^^^^^^^^^^^^^^^

end subroutine micro_precip_fall

!----------------------------------------------------------------------
!!! Initialize the list of microphysics statistics that will be outputted
!!  to *.stat statistics file

subroutine micro_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)


   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count,trcount
   integer ntr

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QVFLUX'
   deflist(count) = 'Water vapor flux (Resolved+SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QVFLUXS'
   deflist(count) = 'Water Vapor flux flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

! Commented out to prevent and issue with multiple definitions of QCFLUX when
!   running with the SAM1MOM microphysics scheme, in which QCFLUX occurs only in
!   the lst file.
! -dschanen UWM 11 Aug 2008
!
!  count = count + 1
!  trcount = trcount + 1
!  namelist(count) = 'QCFLUX'
!  deflist(count) = 'Cloud-water turbulent flux (Resolved+SGS)'
!  unitlist(count) = 'W/m2'
!  status(count) = 1    
!  average_type(count) = 0

!  count = count + 1
!  trcount = trcount + 1
!  namelist(count) = 'QCFLUXS'
!  deflist(count) = 'Cloud-water turbulent flux (SGS)'
!  unitlist(count) = 'W/m2'
!  status(count) = 1    
!  average_type(count) = 0
!
! End comment out

! Variables added by dschanen UWM

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QRFLUX'
   deflist(count) = 'Precipitating liquid water turbulent flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QRFLUXS'
   deflist(count) = 'Precipitating liquid water turbulent flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QSFLUX'
   deflist(count) = 'Precipitating snow water turbulent flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QSFLUXS'
   deflist(count) = 'Precipitating snow water turbulent flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QGFLUX'
   deflist(count) = 'Precipitating graupel water turbulent flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QGFLUXS'
   deflist(count) = 'Precipitating graupel water turbulent flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QC'
   deflist(count) = 'Liquid cloud water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QI'
   deflist(count) = 'Icy cloud water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QR'
   deflist(count) = 'Rain water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QS'
   deflist(count) = 'Snow water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QG'
   deflist(count) = 'Graupel water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

! Not sure about the units on these -dschanen UWM
   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QT_COND'
   deflist(count) = 'Change in total water due to condensation?'
   unitlist(count) = 'g/kg/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QT_EVAP'
   deflist(count) = 'Change in total water due to evaporation?'
   unitlist(count) = 'g/kg/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QT_DEP'
   deflist(count) = 'Change in total water due to deposition?'
   unitlist(count) = 'g/kg/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QT_SUBL'
   deflist(count) = 'Change in total water due to sublimation?'
   unitlist(count) = 'g/kg/s'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QT_FREEZ'
   deflist(count) = 'Change in total water due to freezing of liquid?'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QT_MELT'
   deflist(count) = 'Change in total water due to melting of ice?'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0

! end dschanen UWM

! ...
! etc.

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NCD'
   deflist(count) = 'CD number mix ratio'
   unitlist(count) = '#/mg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NIC'
   deflist(count) = 'IC number mix ratio'
   unitlist(count) = '#/mg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NCN_a'
   deflist(count) = 'activated CCN number mix ratio'
   unitlist(count) = '#/mg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NIN_a'
   deflist(count) = 'activated IN number mix ratio'
   unitlist(count) = '#/mg'
   status(count) = 1    
   average_type(count) = 0

  return
end subroutine micro_hbuf_init

!----------------------------------------------------------------------
!! Collect microphysics history statistics (vertical profiles)
!! Note that only the fields declared in micro_hbuf_init() are allowed to
!! be collected

subroutine micro_statistics()
  
  use vars
  use hbuffer, only: hbuf_put, hbuf_avg_put
  use params, only : lcond

  real tmp(2), factor_xy 
  real qcz(nzm), qiz(nzm), qrz(nzm), qsz(nzm), qgz(nzm), omg
  real qqz(nzm), nwz(nzm), niz(nzm), ncn_az(nzm), nin_az(nzm)
  integer i,j,k,m

  factor_xy = 1./float(nx*ny)

  do k=1,nzm
      tmp(1) = dz/rhow(k)
      tmp(2) = tmp(1) / dtn
      mkwsb(k,1) = mkwsb(k,1) * tmp(1) * rhow(k) * lcond
      mkwle(k,1) = mkwle(k,1)*tmp(2)*rhow(k)*lcond + mkwsb(k,1)

#ifdef CLUBB
      if( (docloud.or.doclubb.or.doclubbnoninter) .and.doprecip) then ! for CLUBB -dschanen UWM 16 May 2008
        mkwsb(k,2) = mkwsb(k,2) * tmp(1) * rhow(k) * lcond
        mkwle(k,2) = mkwle(k,2)*tmp(2)*rhow(k)*lcond + mkwsb(k,2)
      endif
#else
      if(docloud.and.doprecip) then
        mkwsb(k,2) = mkwsb(k,2) * tmp(1) * rhow(k) * lcond
        mkwle(k,2) = mkwle(k,2)*tmp(2)*rhow(k)*lcond + mkwsb(k,2)
      endif
#endif
  end do

  do k=1,nzm
      tmp(1) = dz/rhow(k)
      tmp(2) = tmp(1) / dtn
      mkwsb(k,:) = mkwsb(k,:) * tmp(1) * rhow(k) * lcond
      mkwle(k,:) = mkwle(k,:)*tmp(2)*rhow(k)*lcond + mkwsb(k,:)
  end do
  call hbuf_put('QVFLUX',mkwle(:,3),factor_xy)
  call hbuf_put('QVFLUXS',mkwsb(:,3),factor_xy)
  call hbuf_put('QCFLUX',mkwle(:,4),factor_xy)
  call hbuf_put('QCFLUXS',mkwsb(:,4),factor_xy)
! Added by dschanen UWM
  ! Rain Flux
  call hbuf_put('QRFLUX',mkwle(:,8),factor_xy)
  call hbuf_put('QRFLUXS',mkwsb(:,8),factor_xy)
  ! Snow Flux
  call hbuf_put('QSFLUX',mkwle(:,9),factor_xy)
  call hbuf_put('QSFLUXS',mkwsb(:,9),factor_xy)
  ! Graupel Flux
  call hbuf_put('QGFLUX',mkwle(:,10),factor_xy)
  call hbuf_put('QGFLUXS',mkwsb(:,10),factor_xy)

  ! Diagnostics from Phillips scheme for which we're not sure about the units
  ! Probably setup for kg/kg/s?
  call hbuf_avg_put('QT_COND',g_cond,1,nx,1,ny,nzm,1.e3)
  call hbuf_avg_put('QT_EVAP',g_evap,1,nx,1,ny,nzm,1.e3)
  call hbuf_avg_put('QT_DEP',g_dep,1,nx,1,ny,nzm,1.e3)
  call hbuf_avg_put('QT_SUBL',g_subl,1,nx,1,ny,nzm,1.e3)
  call hbuf_avg_put('QT_FREEZ',g_freeze,1,nx,1,ny,nzm,1.e3)
  call hbuf_avg_put('QT_MELT',g_melt,1,nx,1,ny,nzm,1.e3)

!G_COND(i,j,k) = GAMMA_K(1)
!G_EVAP(i,j,k) = GAMMA_K(2)
!G_DEP(i,j,k) = GAMMA_K(3)
!G_SUBL(i,j,k) = GAMMA_K(4)
!G_FREEZE(i,j,k) = GAMMA_K(5)
!G_MELT(i,j,k) = GAMMA_K(6)

! End added by dschanen UWM

! ... etc

  do k=1,nzm
    qcz(k) = 0.
    qiz(k) = 0.
    qrz(k) = 0.
    qsz(k) = 0.
    qgz(k) = 0.
    nwz(k) = 0. 
    niz(k) = 0. 
    ncn_az(k) = 0.
    nin_az(k) = 0.
    do j=1,ny
    do i=1,nx
      qcz(k)=qcz(k)+qc(i,j,k)
      qiz(k)=qiz(k)+qi(i,j,k)
      qrz(k)=qrz(k)+qr(i,j,k)
      qsz(k)=qsz(k)+qs(i,j,k)
      qgz(k)=qgz(k)+qg(i,j,k)
      nwz(k)=nwz(k)+nw(i,j,k) 
      niz(k)=niz(k)+ni(i,j,k) 
      ncn_az(k)=ncn_az(k)+ncn_a(i,j,k)
      nin_az(k)=nin_az(k)+nin_a(i,j,k)
    end do
    end do
  end do

  call hbuf_put('QC',qcz,1.e3*factor_xy)
  call hbuf_put('QI',qiz,1.e3*factor_xy)
  call hbuf_put('QR',qrz,1.e3*factor_xy)
  call hbuf_put('QS',qsz,1.e3*factor_xy)
  call hbuf_put('QG',qgz,1.e3*factor_xy)
  call hbuf_put('NCD',nwz,1.e-6*factor_xy)
  call hbuf_put('NIC',niz,1.e-6*factor_xy)
  call hbuf_put('NCN_a',ncn_az,1.e-6*factor_xy)
  call hbuf_put('NIN_a',nin_az,1.e-6*factor_xy)

!  call hbuf_avg_put('QC',qc,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
!  call hbuf_avg_put('NC',Nc,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
!  call hbuf_avg_put('QI',qc,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
!  call hbuf_avg_put('NI',Nc,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
!  call hbuf_avg_put('QR',qc,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
!  call hbuf_avg_put('NR',Nc,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)

! ...
end subroutine micro_statistics


!-----------------------------------------------------------------------
! Function that computes total water in the domain:
! Don't change this one.

real function total_water()

  use vars, only : nstep,nprint,adz,dz,rho

  integer k,m

  total_water = 0.
  if(mod(nstep,nprint).ne.0) return

  do m=1,nmicro_fields

   if(flag_wmass(m).eq.1) then

    do k=1,nzm
      total_water = total_water + &
       sum(micro_field(1:nx,1:ny,k,m))*adz(k)*dz*rho(k)
    end do

   end if

  end do

end function total_water


double precision function esatw_dble(ttt)
implicit none
double precision :: ttt	! temperature (K)
double precision :: a0,a1,a2,a3,a4,a5,a6,a7,a8 
data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
!	6.105851, 0.4440316, 0.1430341e-1, &
!        0.2641412e-3, 0.2995057e-5, 0.2031998e-7, &
!        0.6936113e-10, 0.2564861e-13,-0.3704404e-15/
     	6.11239921, 0.443987641, 0.142986287e-1, &
       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
double precision :: dt
dt = max(-80.,ttt-273.16)
esatw_dble = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
end function esatw_dble

double precision function esati_dble(ttt)
implicit none
double precision :: ttt	! temperature (K)
double precision :: a0,a1,a2,a3,a4,a5,a6,a7,a8 
data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
	6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/	
double precision :: dt
if(ttt.gt.185.) then
  dt = ttt-273.16
  esati_dble = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
else   ! use some additional interpolation below 184K
  dt = max(-100.,ttt-273.16)
  esati_dble = 0.00763685 + dt*(0.000151069+dt*7.48215e-07)
end if
end function esati_dble


end module microphysics

