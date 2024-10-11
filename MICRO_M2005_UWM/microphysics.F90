module microphysics

! main interface to Morrison microphysics.
! original implementation by Peter Blossey, UW

use params, only: lcond, lsub, fac_cond, fac_sub, ggr

use grid, only: nx,ny,nzm,nz, &  !grid dimensions; nzm = nz-1 # of scalar lvls
     dimx1_s,dimx2_s,dimy1_s,dimy2_s, & ! actual scalar-array dimensions in x,y
     dz, adz, dostatis, masterproc, &
     doSAMconditionals, dosatupdnconditionals

use vars, only: pres, rho, dtn, w, t, tlatqi, condavg_mask, &
#ifdef PNNL_STATS
     ! Heng Xiao: to get total water flux into the unit needed for CLUBB
     ncondavg, condavgname, condavglongname, rhow
#else
     ncondavg, condavgname, condavglongname
#endif /* PNNL_STATS */
use params, only: doprecip, docloud
#ifdef CLUBB
use sgs_params, only: doclubb
#endif

use module_mp_GRAUPEL, only: GRAUPEL_INIT, M2005MICRO_GRAUPEL, &
      doicemicro, &         ! use ice species (snow/cloud ice/graupel)
      dograupel, &          ! use graupel
      dohail, &             ! use graupel
      dosb_warm_rain, &     ! use Seifert & Beheng (2001) warm rain parameterization
      dopredictNc, &        ! prediction of cloud droplet number
      aerosol_mode, &   ! specify two modes of (sulfate) aerosol
      dosubgridw, &         ! input estimate of subgrid w to microphysics
      doarcticicenucl,&     ! use arctic parameter values for ice nucleation
      docloudedgeactivation,&! activate droplets at cloud edges as well as base
      Nc0,            &     ! initial/specified cloud droplet number conc (#/cm3)
      ccnconst, ccnexpnt, & ! parameters for aerosol_mode=1 (powerlaw CCN)
      aer_rm1, aer_rm2, &   ! two modes of aerosol for aerosol_mode=2
      aer_n1, aer_n2, &     ! rm=geometric mean radius (um), n=aerosol conc. (#/cm3)
      aer_sig1, aer_sig2, & ! sig=geom standard deviation of aerosol size distn.
      dofix_pgam, pgam_fixed ! option to specify pgam (exponent of cloud water's gamma distn)

implicit none

logical :: isallocatedMICRO = .false.

integer :: nmicro_fields ! total number of prognostic water vars

#ifdef CLUBB
real, allocatable, target, dimension(:,:,:,:) :: micro_field  ! holds mphys quantities
#else
real, allocatable, dimension(:,:,:,:) :: micro_field  ! holds mphys quantities
#endif
#if defined(LASSO_ENA) && defined(CALC_RADAR_REFL)
real, allocatable, dimension(:,:,:) :: refl  ! holds mphys quantities
#endif

! indices of water quantities in micro_field, e.g. qv = micro_field(:,:,:,iqv)
integer :: iqv, iqcl, iqci, iqr, iqs, iqg, incl, inci, inr, ins, ing
integer :: index_water_vapor ! separate water vapor index used by SAM

#ifdef PNNL_STATS
integer :: indx_qr, indx_qi, indx_qs, indx_qg ! For advection theta-l/qtog stats
#endif /*PNNL_STATS*/

real, allocatable, dimension(:) :: lfac
integer, allocatable, dimension(:) :: flag_wmass, flag_precip, flag_number
integer, allocatable, dimension(:) :: flag_micro3Dout

integer, parameter :: index_cloud_ice = -1 ! historical variable (don't change)

real, allocatable, dimension(:,:,:) :: fluxbmk, fluxtmk !surface/top fluxes
real, allocatable, dimension(:,:,:) :: reffc, reffi
#ifdef CLUBB
real, allocatable, target, dimension(:,:,:) :: cloudliq
#else
real, allocatable, dimension(:,:,:) :: cloudliq
#endif

real, allocatable, dimension(:,:) :: & ! statistical arrays
     mkwle, & ! resolved vertical flux
     mkwsb, & ! SGS vertical flux
     mksed, & ! sedimentation vertical flux
     mkadv, & ! tendency due to vertical advection
     mkdiff, &! tendency due to vertical diffusion
     mklsadv, & ! tendency due to large-scale vertical advection
     mfrac, & ! fraction of domain with microphysical quantity > 1.e-6
     stend, & ! tendency due to sedimentation
     mtend, & ! tendency due to microphysical processes (other than sedimentation)
     mstor, & ! storage terms of microphysical variables 
     trtau    ! optical depths of various species

real, allocatable, dimension(:) :: tmtend

real :: sfcpcp, sfcicepcp

! arrays with names/units for microphysical outputs in statistics.
character*3, allocatable, dimension(:) :: mkname
character*80, allocatable, dimension(:) :: mklongname
character*10, allocatable, dimension(:) :: mkunits
real, allocatable, dimension(:) :: mkoutputscale
logical douse_reffc, douse_reffi

! You can also have some additional, diagnostic, arrays, for example, total
! nonprecipitating cloud water, etc:

!bloss: array which holds temperature tendency due to microphysics
real, allocatable, dimension(:,:,:), SAVE :: tmtend3d
#ifdef CLUBB 
real, dimension(:,:,:), pointer :: conc, qn
#endif

#ifdef UWM_STATS
#include "uwm_stats_declarations.inc"
#endif /*UWM_STATS*/

! Microphysical process rates  +++mhwang
real, allocatable, dimension(:) :: mPRC, mPRA, mPSMLT, mEVPMS, mPRACS, mEVPMG, mPRACG, mPRE, mPGMLT, &
                        mMNUCCC, mPSACWS, mPSACWI, mQMULTS, mQMULTG, mPSACWG, mPGSACW, &
                        mPRD, mPRCI, mPRAI, mQMULTR, mQMULTRG, mMNUCCD, mPRACI, mPRACIS, mEPRD, &
                        mMNUCCR, mPIACR, mPIACRS, mPGRACS, &
                        mPRDS, mEPRDS, mPSACR, &
                        mPRDG, mEPRDG, &
                        mNPRC1, mNRAGG, mNPRACG, mNSUBR, mNSMLTR, &
                        mNGMLTR, mNPRACS, mNNUCCR, mNIACR, mNIACRS, mNGRACS, &
                        mNSMLTS, mNSAGG, mNPRCI, mNSCNG, mNSUBS, mPCC, &
                        mNNUCCC, mNPSACWS, mNPRA, mNPRC, mNPSACWI, &
                        mNPSACWG, mNPRAI, mNMULTS, mNMULTG, mNMULTR, &
                        mNMULTRG, mNNUCCD, mNSUBI, mNGMLTG, mNSUBG, mNACT, &
                        mSIZEFIX_NR, mSIZEFIX_NC, mSIZEFIX_NI, mSIZEFIX_NS, mSIZEFIX_NG, &
                        mNEGFIX_NR, mNEGFIX_NC, mNEGFIX_NI, mNEGFIX_NS, mNEGFIX_NG, &
                        mNIM_MORR_CL, mQC_INST, mQR_INST, mQI_INST, mQS_INST, mQG_INST, &
                        mNC_INST, mNR_INST, mNI_INST, mNS_INST, mNG_INST  

#if defined(UWM_STATS) || defined(LASSO_ENA_3D_MICRO_OUTPUT)
#include "uwm_stats_declarations_3d.inc"
#endif

#ifdef UWM_STATS

! weberjk(UWM), liquid water potential temperature, chi (s_mellor),
! and eta (t_mellor)
real, allocatable, dimension(:,:,:) :: theta_l, chi, eta

#endif /*UWM_STATS*/


CONTAINS

!----------------------------------------------------------------------
!!! Read microphysical options from prm file and allocate variables
!
subroutine micro_setparm()
  use vars
#ifdef CLUBB
  use module_mp_graupel, only: NNUCCD_REDUCE_COEF, NNUCCC_REDUCE_COEF
#endif
  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder
  
   NAMELIST /MICRO_M2005/ &
#ifdef CLUBB
      NNUCCD_REDUCE_COEF, NNUCCC_REDUCE_COEF, &
#endif
#ifdef UWM_STATS
      doensemble_fractions,&  ! Rather than a fixed threshold to compute
                              ! in-'precip' means, variances, etc, use an ensemble of thresholds
      frac_threshold_init, &  ! Initial threshold value
      nfractions         , &  ! Number of fractions to compute

#endif

      doicemicro, &         ! use ice species (snow/cloud ice/graupel)
      dograupel, &          ! use graupel
      dohail, &             ! graupel species has qualities of hail
      dosb_warm_rain, &     ! use Seifert & Beheng (2001) warm rain parameterization in place of KK(2000)
      dopredictNc, &        ! prediction of cloud droplet number
      aerosol_mode, &   ! specify two modes of (sulfate) aerosol
      dosubgridw, &         ! input estimate of subgrid w to microphysics
      doarcticicenucl,&     ! use arctic parameter values for ice nucleation
      docloudedgeactivation,&! activate droplets at cloud edges as well as base
      Nc0,            &     ! initial/specified cloud droplet number conc (#/cm3)
      ccnconst, ccnexpnt, & ! parameters for aerosol_mode=1 (powerlaw CCN)
      aer_rm1, aer_rm2, &   ! two modes of aerosol for aerosol_mode=2
      aer_n1, aer_n2, &     ! rm=geometric mean radius (um), n=aerosol conc. (#/cm3)
      aer_sig1, aer_sig2, & ! sig=geom standard deviation of aerosol size distn.
      dofix_pgam, pgam_fixed, & ! option to specify pgam (exponent of cloud water's gamma distn)
      douse_reffc, &        ! use computed effective radius in radiation computation
      douse_reffi           ! use computed effective ice size in radiation computation

   !bloss: Create dummy namelist, so that we can figure out error code
   !       for a mising namelist.  This lets us differentiate between
   !       missing namelists and those with an error within the namelist.
   NAMELIST /BNCUIODSBJCB/ place_holder
#ifdef SILHS
   integer :: k
#endif

   ! define default values for namelist variables
   doicemicro = .true.        ! use ice
   dograupel = .true.         ! use graupel
   dohail = .false.           ! graupel species has properties of graupel
   dosb_warm_rain = .false.   ! use KK (2000) warm rain scheme by default
   dopredictNc = .false.       ! prognostic cloud droplet number (default = .false. !bloss Apr 09)
#ifdef CLUBB
   dosubgridw = .true.        ! Use clubb's w'^2 for sgs w
   aerosol_mode = 2           ! use lognormal CCN relationship
   doarcticicenucl = .false.  ! use mid-latitude parameters
#else
   aerosol_mode = 1           ! use powerlaw CCN relationship
   dosubgridw=.false.         ! don't bother with estimating w_sgs for now
   doarcticicenucl = .false.  ! use mid-latitude parameters
#endif
   docloudedgeactivation  = .false. ! activate droplets at cloud base, not edges
   douse_reffc = .false.  ! use computed effective radius in rad computations?
   douse_reffi = .false.  ! use computed effective radius in rad computations?

   Nc0 = 100. ! default droplet number concentration
   
   ccnconst = 120.            ! maritime value (/cm3), adapted from Rasmussen 
   ccnexpnt = 0.4             !   et al (2002) by Hugh Morrison et al.  Values
                              !   of 1000. and 0.5 suggested for continental
#ifdef CLUBB
   ! Aerosol for RF02 from Mikhail
   aer_rm1  = 0.011
   aer_sig1 = 1.2
   aer_n1   = 125.
   aer_rm2  = 0.06
   aer_sig2 = 1.7
   aer_n2   = 65.

#else
   aer_rm1 = 0.052           ! two aerosol mode defaults from MPACE (from Hugh)
   aer_sig1 = 2.04
   aer_n1 = 72.2
   aer_rm2 = 1.3
   aer_sig2 = 2.5
   aer_n2 = 1.8
#endif /* CLUBB */
   
   dofix_pgam = .false.
   pgam_fixed = 5. ! middle range value -- corresponds to radius dispersion ~ 0.4

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
  
  !bloss: get error code for missing namelist (by giving the name for
  !       a namelist that doesn't exist in the prm file).
  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  rewind(55) !note that one must rewind before searching for new namelists

  !bloss: read in MICRO_M2005 namelist
  read (55,MICRO_M2005,IOSTAT=ios)

  if (ios.ne.0) then
     !namelist error checking
     if(ios.ne.ios_missing_namelist) then
        write(*,*) '****** ERROR: bad specification in MICRO_M2005 namelist'
        call task_abort()
     elseif(masterproc) then
        write(*,*) '****************************************************'
        write(*,*) '****** No MICRO_M2005 namelist in prm file *********'
        write(*,*) '****************************************************'
     end if
  end if
  close(55)

   if(.not.doicemicro) dograupel=.false.

   if(dohail.and..NOT.dograupel) then
      if(masterproc) write(*,*) 'dograupel must be .true. for dohail to be used.'
      call task_abort()
   end if

#ifndef CLUBB /* Disable this for UWM simulations */
   ! write namelist values out to file for documentation
   if(masterproc) then
      open(unit=55,file='./'//trim(case)//'/'//trim(case)//'_'//trim(caseid)//'.namelists', form='formatted', position='append')
      write (unit=55,nml=MICRO_M2005,IOSTAT=ios)
      write(55,*) ' '
      close(unit=55)
   end if
#endif

   ! scale values of parameters for m2005micro
   aer_rm1 = 1.e-6*aer_rm1 ! convert from um to m
   aer_rm2 = 1.e-6*aer_rm2 
   aer_n1 = 1.e6*aer_n1 ! convert from #/cm3 to #/m3
   aer_n2 = 1.e6*aer_n2
#ifdef MICRO_RESTART /* Save all fields to allow for restarting with micro enabled */
  nmicro_fields = 10
  iqv  = 1   ! total water (vapor + cloud liq) mass mixing ratio [kg H2O / kg dry air]
  incl = 2  ! cloud water number mixing ratio [#/kg dry air]
  iqr  = 3   ! rain mass mixing ratio [kg H2O / kg dry air]
  inr  = 4   ! rain number mixing ratio [#/kg dry air]
  iqci = 5  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
  inci = 6  ! cloud ice number mixing ratio [#/kg dry air]
  iqs  = 7   ! snow mass mixing ratio [kg H2O / kg dry air]
  ins  = 8   ! snow number mixing ratio [#/kg dry air]
  iqg  = 9  ! graupel mass mixing ratio [kg H2O / kg dry air]
  ing  = 10  ! graupel number mixing ratio [#/kg dry air]
#else
  nmicro_fields = 1 ! start with water vapor and cloud water mass mixing ratio
#ifdef CLUBB
  if(docloud.or.doclubb) then
#else
  if(docloud) then
#endif
!bloss/qt     nmicro_fields = nmicro_fields + 1 ! add cloud water mixing ratio
     if(dopredictNc) nmicro_fields = nmicro_fields + 1 ! add cloud water number concentration (if desired)
  end if
  if(doprecip)    nmicro_fields = nmicro_fields + 2 ! add rain mass and number (if desired)
  if(doicemicro)  nmicro_fields = nmicro_fields + 4 ! add snow and cloud ice number and mass (if desired)
  if(dograupel)   nmicro_fields = nmicro_fields + 2 ! add graupel mass and number (if desired)

  ! specify index of various quantities in micro_field array
  !  *** note that not all of these may be used if(.not.doicemicro) ***
  iqv = 1   ! total water (vapor + cloud liq) mass mixing ratio [kg H2O / kg dry air]
!bloss/qt  iqcl = 2  ! cloud water mass mixing ratio [kg H2O / kg dry air]
  
!bloss/qt: cloud liquid water no longer prognosed
  if(dopredictNc) then
     incl = 2  ! cloud water number mixing ratio [#/kg dry air]
     iqr = 3   ! rain mass mixing ratio [kg H2O / kg dry air]
     inr = 4   ! rain number mixing ratio [#/kg dry air]
     iqci = 5  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
     inci = 6  ! cloud ice number mixing ratio [#/kg dry air]
     iqs = 7   ! snow mass mixing ratio [kg H2O / kg dry air]
     ins = 8   ! snow number mixing ratio [#/kg dry air]
     iqg = 9  ! graupel mass mixing ratio [kg H2O / kg dry air]
     ing = 10  ! graupel number mixing ratio [#/kg dry air]
  else
     iqr = 2   ! rain mass mixing ratio [kg H2O / kg dry air]
     inr = 3   ! rain number mixing ratio [#/kg dry air]
     iqci = 4  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
     inci = 5  ! cloud ice number mixing ratio [#/kg dry air]
     iqs = 6   ! snow mass mixing ratio [kg H2O / kg dry air]
     ins = 7   ! snow number mixing ratio [#/kg dry air]
     iqg = 8   ! graupel mass mixing ratio [kg H2O / kg dry air]
     ing = 9  ! graupel number mixing ratio [#/kg dry air]
  end if
#endif /* MICRO_RESTART */
#ifdef PNNL_STATS
if(doprecip) then
  indx_qr = iqr  ! Index for rain water mixing ratio (theta-l/qtog adv. budgets)
endif

if(doicemicro) then
  indx_qi = iqci ! Index for ice mixing ratio (theta-l/qtog adv. budgets)
  indx_qs = iqs  ! Index for snow mixing ratio (theta-l/qtog adv. budgets)
endif

if(dograupel) then
  indx_qg = iqg  ! Index for graupel mixing ratio (theta-l/qtog adv. budgets)
endif
#endif /* PNNL_STATS */
  ! stop if icemicro is specified without precip -- we don't support this right now.
  if((doicemicro).and.(.not.doprecip)) then
     if(masterproc) write(*,*) 'Morrison 2005 Microphysics does not support both doice and .not.doprecip'
     call task_abort()
  end if
  index_water_vapor = iqv ! set SAM water vapor flag

  if(.not.isallocatedMICRO) then
     ! allocate microphysical variables
     allocate(micro_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,nmicro_fields), &
#if defined(LASSO_ENA) && defined(CALC_RADAR_REFL)
          refl(nx,ny,nzm), &
#endif
          fluxbmk(nx,ny,nmicro_fields), fluxtmk(nx,ny,nmicro_fields), &
          reffc(nx,ny,nzm), reffi(nx,ny,nzm), &
          mkwle(nz,nmicro_fields), mkwsb(nz,nmicro_fields), &
          mkadv(nz,nmicro_fields), mkdiff(nz,nmicro_fields), &
          mklsadv(nz,nmicro_fields), &
          stend(nzm,nmicro_fields), mtend(nzm,nmicro_fields), &
          mfrac(nzm,nmicro_fields), trtau(nzm,nmicro_fields), &
          mksed(nzm,nmicro_fields), tmtend(nzm), &
          mstor(nzm,nmicro_fields),  &
          cloudliq(nx,ny,nzm), &
          tmtend3d(nx,ny,nzm), flag_micro3Dout(nmicro_fields), &
          flag_wmass(nmicro_fields), flag_precip(nmicro_fields), &
          flag_number(nmicro_fields), lfac(nmicro_fields), &
          mkname(nmicro_fields), mklongname(nmicro_fields), &
          mkunits(nmicro_fields), mkoutputscale(nmicro_fields), STAT=ierr)

 
! Microphysical process rates +++mhwang 
      allocate(mPRC(nzm), mPRA(nzm), mPSMLT(nzm), mEVPMS(nzm), mPRACS(nzm), mEVPMG(nzm), mPRACG(nzm), mPRE(nzm), mPGMLT(nzm),&
                        mMNUCCC(nzm), mPSACWS(nzm), mPSACWI(nzm), mQMULTS(nzm), mQMULTG(nzm), mPSACWG(nzm), mPGSACW(nzm), & 
                        mPRD(nzm), mPRCI(nzm), mPRAI(nzm), mQMULTR(nzm), mQMULTRG(nzm), mMNUCCD(nzm), mPRACI(nzm), &
                        mPRACIS(nzm), mEPRD(nzm), mMNUCCR(nzm), mPIACR(nzm), mPIACRS(nzm), mPGRACS(nzm), &
                        mPRDS(nzm), mEPRDS(nzm), mPSACR(nzm), mPRDG(nzm), mEPRDG(nzm), mNPRC1(nzm), mNRAGG(nzm), &
                        mNPRACG(nzm), mNSUBR(nzm), mNSMLTR(nzm), mNGMLTR(nzm), mNPRACS(nzm), &
                        mNNUCCR(nzm), mNIACR(nzm), mNIACRS(nzm), mNGRACS(nzm), mNSMLTS(nzm), &
                        mNSAGG(nzm), mNPRCI(nzm), mNSCNG(nzm), mNSUBS(nzm), mPCC(nzm), &
                        mNNUCCC(nzm), mNPSACWS(nzm), mNPRA(nzm), mNPRC(nzm), mNPSACWI(nzm), &
                        mNPSACWG(nzm), mNPRAI(nzm), mNMULTS(nzm), mNMULTG(nzm), mNMULTR(nzm), &
                        mNMULTRG(nzm), mNNUCCD(nzm), mNSUBI(nzm), mNGMLTG(nzm), mNSUBG(nzm), mNACT(nzm), & 
                        mSIZEFIX_NR(nzm), mSIZEFIX_NC(nzm), mSIZEFIX_NI(nzm), mSIZEFIX_NS(nzm), &
                        mSIZEFIX_NG(nzm), mNEGFIX_NR(nzm), mNEGFIX_NC(nzm), mNEGFIX_NI(nzm), &
                        mNEGFIX_NS(nzm), mNEGFIX_NG(nzm), mNIM_MORR_CL(nzm), &
                        mQC_INST(nzm), mQR_INST(nzm), mQI_INST(nzm), mQS_INST(nzm), mQG_INST(nzm), &
                        mNC_INST(nzm), mNR_INST(nzm), mNI_INST(nzm), mNS_INST(nzm), mNG_INST(nzm), &
                        STAT=ierr)

#if defined(UWM_STATS) || defined(LASSO_ENA_3D_MICRO_OUTPUT)
#include "uwm_stats_allocate_3d.inc"
#endif

#ifdef UWM_STATS
! weberjk(UWM), liquid water potential temperature, chi (s_mellor),
! and eta (t_mellor).
      allocate(theta_l(nx,ny,nzm), chi(nx,ny,nzm), eta(nx,ny,nzm), STAT=ierr)

! weberjk(UWM), store fraction names and units to create them within SAM. Create
! 5 elements for the (possibly) 5 prognostic variables (cloud, rain, ice, snow,
! graupel)
      allocate(frac_name(5), frac_longname(5), frac_in_char(5), STAT=ierr)
#endif /*UWM_STATS*/
     if(ierr.ne.0) then
        write(*,*) 'Failed to allocate microphysical arrays on proc ', rank
        call task_abort()
     else
        isallocatedMICRO = .true.
     end if

     ! zero out statistics variables associated with cloud ice sedimentation
     !   in Marat's default SAM microphysics
     tlatqi = 0.

     ! initialize these arrays
     micro_field = 0.
     cloudliq = 0. !bloss/qt: auxially cloud liquid water variable, analogous to qn in MICRO_SAM1MOM
     fluxbmk = 0.
     fluxtmk = 0.
     mkwle = 0.
     mkwsb = 0.
     mkadv = 0.
     mkdiff = 0.
     mklsadv = 0.
     mstor =0.

    ! initialize flag arrays to all mass, no number, no precip
     flag_wmass = 1
     flag_number = 0
     flag_precip = 0
     flag_micro3Dout = 0

  end if

  compute_reffc = douse_reffc
  compute_reffi = douse_reffi
#ifdef CLUBB
  if ( dopredictnc ) then
    conc => micro_field(1:nx,1:ny,1:nzm,incl)
  else
    allocate( conc(nx,ny,nzm) )
  end if
  qn => cloudliq
#endif 
end subroutine micro_setparm

!----------------------------------------------------------------------
!!! Initialize microphysics:
!
! this one is guaranteed to be called by SAM at the 
!   beginning of each run, initial or restart:
subroutine micro_init()

  use vars
  
  implicit none
  
  real, dimension(nzm) :: qc0 !, qi0  !MWSWong:qi0 not used, comment out

! Commented out by dschanen UWM 23 Nov 2009 to avoid a linking error
! real, external :: satadj_water 
  integer :: k
  integer :: i, j, n

#ifndef MICRO_RESTART
  write(6,*) 'microphysics: '
  write(6,*) 'microphysics: dopredictNc', dopredictNc
  write(6,*) 'microphysics: doicemicro', doicemicro
  write(6,*) 'microphysics: doprecip', doprecip
  write(6,*) 'microphysics: dograupel', dograupel
  ! initialize flag arrays
  if(dopredictNc) then
     ! Cloud droplet number concentration is a prognostic variable
     if(doicemicro) then
        if(dograupel) then
          !bloss/qt: qt, Nc, qr, Nr, qi, Ni, qs, Ns, qg, Ng
           flag_wmass  = (/1,0,1,0,1,0,1,0,1,0/)
           flag_precip = (/0,0,1,1,0,0,1,1,1,1/)
           flag_number = (/0,1,0,1,0,1,0,1,0,1/)
        else
          !bloss/qt: qt, Nc, qr, Nr, qi, Ni, qs, Ns
           flag_wmass  = (/1,0,1,0,1,0,1,0/)
           flag_precip = (/0,0,1,1,0,0,1,1/)
           flag_number = (/0,1,0,1,0,1,0,1/)
        end if
     else
        if(doprecip) then
          !bloss/qt: qt, Nc, qr, Nr
           flag_wmass  = (/1,0,1,0/)
           flag_precip = (/0,0,1,1/)
           flag_number = (/0,1,0,1/)
        else
          !bloss/qt: qt, Nc
           flag_wmass  = (/1,0/)
           flag_precip = (/0,0/)
           flag_number = (/0,1/)
        end if
     end if
 else
     ! Cloud droplet number concentration is NOT a prognostic variable
     if(doicemicro) then
        if(dograupel) then
          !bloss/qt: qt, qr, Nr, qi, Ni, qs, Ns, qg, Ng
           flag_wmass  = (/1,1,0,1,0,1,0,1,0/)
           flag_precip = (/0,1,1,0,0,1,1,1,1/)
           flag_number = (/0,0,1,0,1,0,1,0,1/)
        else
          !bloss/qt: qt, qr, Nr, qi, Ni, qs, Ns
           flag_wmass  = (/1,1,0,1,0,1,0/)
           flag_precip = (/0,1,1,0,0,1,1/)
           flag_number = (/0,0,1,0,1,0,1/)
        end if
     else
        if(doprecip) then
          !bloss/qt: qt, qr, Nr
           flag_wmass  = (/1,1,0/)
           flag_precip = (/0,1,1/)
           flag_number = (/0,0,1/)
        else
          !bloss/qt: only total water variable is needed for no-precip, 
          !            fixed droplet number, warm cloud and no cloud simulations.
           flag_wmass  = (/1/)
           flag_precip = (/0/)
           flag_number = (/0/)
        end if
     end if
  end if
#else
flag_wmass  = (/1,0,1,0,1,0,1,0,1,0/)
flag_precip = (/0,0,1,1,0,0,1,1,1,1/)
flag_number = (/0,1,0,1,0,1,0,1,0,1/)
#endif /* MICRO_RESTART */
  ! output all microphysical fields to 3D output files if using more than
  !   just docloud.  Otherwise, rely on basic SAM outputs
#ifdef CLUBB
  if((docloud.OR.doclubb).AND.(doprecip.OR.dopredictNc)) then
#else
  if(docloud.AND.(doprecip.OR.dopredictNc)) then
#endif
     flag_micro3Dout = 1
  end if

  ! initialize factor for latent heat
  lfac(:) = 1. ! use one as default for number species
  lfac(iqv) = lcond
!bloss/qt  if(docloud) lfac(iqcl) = lcond
  if(doprecip) lfac(iqr) = lcond
  if(doicemicro) then
     lfac(iqci) = lsub
     lfac(iqs) = lsub
     if(dograupel) lfac(iqg) = lsub
  end if

  call graupel_init() ! call initialization routine within mphys module

  if(nrestart.eq.0) then

 ! compute initial profiles of liquid water - M.K.
      call satadj_liquid(nzm,tabs0,q0,qc0,pres*100.)

     ! initialize microphysical quantities
     q0 = q0 + qc0
     do k = 1,nzm
        micro_field(:,:,k,iqv) = q0(k)
        cloudliq(:,:,k) = qc0(k)
        tabs(:,:,k) = tabs0(k)
     end do
     if(dopredictNc) then ! initialize concentration somehow...
       do k = 1,nzm
         if(q0(k).gt.0.) then
            micro_field(:,:,k,incl) = 0.5*ccnconst*1.e6
         end if
       end do
     end if
#ifdef CLUBB
     if(docloud.or.doclubb)  call micro_diagnose()   ! leave this line here
#else
     if(docloud) call micro_diagnose()   ! leave this here
#endif


  end if

  ! set up names, units and scales for these microphysical quantities
  mkname(iqv) = 'QTO'
  mklongname(iqv) = 'TOTAL WATER (VAPOR + CLOUD LIQUID)'
  mkunits(iqv) = 'g/kg'
  mkoutputscale(iqv) = 1.e3

#ifdef MICRO_RESTART
  mkname(incl) = 'NC'
  mklongname(incl) = 'CLOUD WATER NUMBER CONCENTRATION'
  mkunits(incl) = '#/cm3'
  mkoutputscale(incl) = 1.e-6

  mkname(iqr) = 'QR'
  mklongname(iqr) = 'RAIN'
  mkunits(iqr) = 'g/kg'
  mkoutputscale(iqr) = 1.e3

  mkname(inr) = 'NR'
  mklongname(inr) = 'RAIN NUMBER CONCENTRATION'
  mkunits(inr) = '#/cm3'
  mkoutputscale(inr) = 1.e-6

  mkname(iqci) = 'QI'
  mklongname(iqci) = 'CLOUD ICE'
  mkunits(iqci) = 'g/kg'
  mkoutputscale(iqci) = 1.e3

  mkname(inci) = 'NI'
  mklongname(inci) = 'CLOUD ICE NUMBER CONCENTRATION'
  mkunits(inci) = '#/cm3'
  mkoutputscale(inci) = 1.e-6

  mkname(iqs) = 'QS'
  mklongname(iqs) = 'SNOW'
  mkunits(iqs) = 'g/kg'
  mkoutputscale(iqs) = 1.e3

  mkname(ins) = 'NS'
  mklongname(ins) = 'SNOW NUMBER CONCENTRATION'
  mkunits(ins) = '#/cm3'
  mkoutputscale(ins) = 1.e-6

  mkname(iqg) = 'QG'
  mklongname(iqg) = 'GRAUPEL'
  mkunits(iqg) = 'g/kg'
  mkoutputscale(iqg) = 1.e3

  mkname(ing) = 'NG'
  mklongname(ing) = 'GRAUPEL NUMBER CONCENTRATION'
  mkunits(ing) = '#/cm3'
  mkoutputscale(ing) = 1.e-6
#else
#ifdef CLUBB
  if(docloud.or.doclubb) then
#else
  if(docloud) then
#endif
!bloss/qt     mkname(iqcl) = 'QC'
!bloss/qt     mklongname(iqcl) = 'CLOUD WATER'
!bloss/qt     mkunits(iqcl) = 'g/kg'
!bloss/qt     mkoutputscale(iqcl) = 1.e3

! We want to output NC even if dopredictNc is false for certainty. 
     if(dopredictNc) then
        mkname(incl) = 'NC'
        mklongname(incl) = 'CLOUD WATER NUMBER CONCENTRATION'
        mkunits(incl) = '#/cm3'
        mkoutputscale(incl) = 1.e-6
     end if
  end if

  if(doprecip) then
     mkname(iqr) = 'QR'
     mklongname(iqr) = 'RAIN'
     mkunits(iqr) = 'g/kg'
     mkoutputscale(iqr) = 1.e3

     mkname(inr) = 'NR'
     mklongname(inr) = 'RAIN NUMBER CONCENTRATION'
     mkunits(inr) = '#/cm3'
     mkoutputscale(inr) = 1.e-6
  end if

  if(doicemicro) then
     mkname(iqci) = 'QI'
     mklongname(iqci) = 'CLOUD ICE'
     mkunits(iqci) = 'g/kg'
     mkoutputscale(iqci) = 1.e3

     mkname(inci) = 'NI'
     mklongname(inci) = 'CLOUD ICE NUMBER CONCENTRATION'
     mkunits(inci) = '#/cm3'
     mkoutputscale(inci) = 1.e-6

     mkname(iqs) = 'QS'
     mklongname(iqs) = 'SNOW'
     mkunits(iqs) = 'g/kg'
     mkoutputscale(iqs) = 1.e3

     mkname(ins) = 'NS'
     mklongname(ins) = 'SNOW NUMBER CONCENTRATION'
     mkunits(ins) = '#/cm3'
     mkoutputscale(ins) = 1.e-6

     if(dograupel) then
        mkname(iqg) = 'QG'
        mklongname(iqg) = 'GRAUPEL'
        mkunits(iqg) = 'g/kg'
        mkoutputscale(iqg) = 1.e3

        mkname(ing) = 'NG'
        mklongname(ing) = 'GRAUPEL NUMBER CONCENTRATION'
        mkunits(ing) = '#/cm3'
        mkoutputscale(ing) = 1.e-6
     end if
end if
#endif /* MICRO_RESTART */

! set mstor to be the inital microphysical mixing ratios    
     do n=1, nmicro_fields
     do k=1, nzm
         !mstor(k, n) = SUM(micro_field(1:nx,1:ny,k,n))
         mstor(k, n) = SUM(dble(micro_field(1:nx,1:ny,k,n))) !MWSWong: cast to double
     end do
     end do

end subroutine micro_init

!----------------------------------------------------------------------
!!! fill-in surface and top boundary fluxes:
!
! Obviously, for liquid/ice water variables those fluxes are zero. They are not zero
! only for water vapor variable and, possibly, for CCN and IN if you have those.

subroutine micro_flux()

use vars, only: fluxbq, fluxtq
#ifdef CLUBB
use sgs_params, only: doclubb, doclubb_sfc_fluxes
#endif

fluxbmk(:,:,:) = 0. ! initialize all fluxes at surface to zero
fluxtmk(:,:,:) = 0. ! initialize all fluxes at top of domain to zero
#ifdef CLUBB
if ( doclubb .and. doclubb_sfc_fluxes ) then
  fluxbmk(:,:,index_water_vapor) = 0.0 ! surface qv (latent heat) flux
else
  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:) ! surface qv (latent heat) flux
end if
#else
fluxbmk(:,:,index_water_vapor) = fluxbq(:,:) ! surface qv (latent heat) flux
#endif
fluxtmk(:,:,index_water_vapor) = fluxtq(:,:) ! top of domain qv flux

end subroutine micro_flux

!----------------------------------------------------------------------
!!! compute local microphysics processes (beyond advection and SGS diffusion):
!
!  This is the place where the condensation/sublimation, accretion, coagulation, freezing,
!  melting, etc., that is  all the microphysics processes except for the spatial transport happen.

! IMPORTANT: You need to use the thermodynamic constants like specific heat, or
! specific heat of condensation, gas constant, etc, the same as in file params.f90
! Also, you should assume that the conservative thermodynamic variable during these
! proceses is the liquid/ice water static energy: t = tabs + gz - Lc (qc+qr) - Ls (qi+qs+qg) 
! It should not be changed during all of your point microphysical processes!

subroutine micro_proc()

use params, only: fac_cond, fac_sub, rgas
use grid, only: z, zi
use vars, only: t,  gamaz, precsfc, precflux, qpfall, tlat, prec_xy, &
     nstep, nstatis, icycle, total_water_prec

#if defined(UWM_STATS) || defined(LASSO_ENA_3D_MICRO_OUTPUT) || (defined(LASSO_ENA) && defined(CALC_RADAR_REFL))
use grid, only: nsave3D, nsave3dstart, nsave3dend
#endif
#ifdef UWM_STATS
!weberjk(UWM), to compute budget statistics on q2, t2, tw, qw, include the
!effect of precipitation (microphysics)    
use vars, only: t2leprec, q2leprec, qwleprec, twleprec, prespot, qsatw
use module_mp_GRAUPEL, only: Nc0
use compute_chi_module, only: compute_chi_eta
use calc_vars_util, only: t2thetal
#endif /*UWM_STATS*/

#ifdef PNNL_STATS
!MWSWong: budget for THL, QTOG, and QTHL
use vars, only: thel2leprec, thelwleprec, thel, thellat, &
                qtog2leprec, qtogwleprec, qtog, &
                qthelleprec
#endif /*PNNL_STATS*/

#ifdef CLUBB
use params, only: docloud, dosmoke
use sgs_params, only: doclubb
use grid, only: nz
use clubb_api_module, only: &
  clubb_at_least_debug_level_api, &
  fill_holes_vertical_api, &
  core_rknd
use clubbvars, only: wp2, cloud_frac, rho_ds_zt, rho_ds_zm ! are used, but not modified here
use vars, only: qcl ! Used here and updated in micro_diagnose
use vars, only: prespot ! exner^-1
use module_mp_GRAUPEL, only: &
  cloud_frac_thresh ! Threshold for using sgs cloud fraction to weight 
                    ! microphysical quantities [%]

#endif
#ifdef SILHS
use clubb_silhs_vars, only: &
  lh_microphys_type, & ! Variables
  lh_microphys_interactive, &
  lh_microphys_non_interactive, &
  lh_microphys_disabled
#endif

real, dimension(nzm) :: &
#if defined(LASSO_ENA) && defined(CALC_RADAR_REFL)
     refl1d, &
#endif 
     tmpqcl, tmpqci, tmpqr, tmpqs, tmpqg, tmpqv, &
     tmpncl, tmpnci, tmpnr, tmpns, tmpng,  &
     tmpw, tmpwsub, tmppres, tmpdz, tmptabs, &
     tmtend1d, &
     mtendqcl, mtendqci, mtendqr, mtendqs, mtendqg, mtendqv, &
     mtendncl, mtendnci, mtendnr, mtendns, mtendng,  &
     stendqcl, stendqci, stendqr, stendqs, stendqg, stendqv, &
     stendncl, stendnci, stendnr, stendns, stendng,  &
     effg1d, effr1d, effs1d, effc1d, effi1d

real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real f0(nzm),df0(nzm)

#ifdef CLUBB
real(kind=core_rknd), dimension(nz) :: &
     qv_clip, qcl_clip

!real, dimension(nzm) :: & These variables were renamed PRC, PRA, and PRE respectively as per ticket #552 and are declared below
!  qr_auto, qr_accr, qr_evap
real, dimension(nzm) :: &
  qr_auto, qr_accr, qr_evap
real, dimension(nzm) :: cloud_frac_in
#endif /*CLUBB*/

#if defined(CLUBB) || defined(UWM_STATS)
real, dimension(nzm) :: rain_vel
#endif


! Microphysical rate +++mwang  
real, dimension(nzm) :: PRC, PRA, PSMLT, EVPMS, PRACS, EVPMG, PRACG, PRE, PGMLT,  & 
                        MNUCCC, PSACWS, PSACWI, QMULTS, QMULTG, PSACWG, PGSACW, & 
                        PRD, PRCI, PRAI, QMULTR, QMULTRG, MNUCCD, PRACI, PRACIS, EPRD, & 
                        MNUCCR, PIACR, PIACRS, PGRACS, & 
                        PRDS, EPRDS, PSACR, & 
                        PRDG, EPRDG, &
                        NPRC1, NRAGG, NPRACG, NSUBR, NSMLTR, NGMLTR, &
                        NPRACS, NNUCCR, NIACR, NIACRS, NGRACS, &
                        NSMLTS, NSAGG, NPRCI, NSCNG, NSUBS, PCC, NNUCCC, &
                        NPSACWS, NPRA, NPRC, NPSACWI, NPSACWG, NPRAI, &
                        NMULTS, NMULTG, NMULTR, NMULTRG, NNUCCD, NSUBI, NGMLTG, NSUBG, NACT, &
                        SIZEFIX_NR, SIZEFIX_NC, SIZEFIX_NI, SIZEFIX_NS, SIZEFIX_NG, &
                        NEGFIX_NR, NEGFIX_NC, NEGFIX_NI, NEGFIX_NS, NEGFIX_NG, &
                        NIM_MORR_CL, QC_INST, QR_INST, QI_INST, QS_INST, QG_INST, &
                        NC_INST, NR_INST, NI_INST, NS_INST, NG_INST 

real, dimension(nzm,nmicro_fields) :: stend1d, mtend1d
real :: tmpc, tmpr, tmpi, tmps, tmpg
integer :: i1, i2, j1, j2, i, j, k, m, n

real(8) :: tmp_total, tmptot

#ifdef UWM_STATS
!weberjk(UWM), declare arrays to compute statistics. Make them less ambiguous
!than 'f' and 'df'
real t_before(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) 
real qt_before(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) 
real t_avg(nzm), t_before_avg(nzm) 
real qt_avg(nzm), qt_avg_before(nzm)
#endif /*UWM_STATS*/

#ifdef PNNL_STATS
!MWSWong:budgets for THL and QTOG
real thel_before(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) 
real thel_avg(nzm), thel_before_avg(nzm)
real qtog_before(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) 
real qtog_avg(nzm), qtog_before_avg(nzm)
#endif /*PNNL_STATS*/

#if defined(LASSO_ENA) && defined(CALC_RADAR_REFL)
logical :: calc_radar_refl
#endif
call t_startf ('micro_proc')

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

#if defined(LASSO_ENA) && defined(CALC_RADAR_REFL)
refl(:,:,:) = 0.0
calc_radar_refl = mod(nstep,nsave3D).eq.0.and.nstep.ge.nsave3Dstart.and.nstep.le.nsave3Dend
#endif

if(dostatis) then ! initialize arrays for statistics
   mfrac(:,:) = 0.
   mtend(:,:) = 0.
   trtau(:,:) = 0.
   qpfall(:)=0.
   tlat(:) = 0.
   tmtend3d(:,:,:) = 0.
#ifdef PNNL_STATS
   thellat(:) = 0.
#endif /*PNNL_STATS*/
! Microphysic process rates +++mhwang 
   mPRC=0.0  
   mPRA=0.0  
   mPSMLT=0.0  
   mEVPMS=0.0  
   mPRACS=0.0  
   mEVPMG=0.0  
   mPRACG=0.0  
   mPRE=0.0  
   mPGMLT=0.0   
   mMNUCCC=0.0  
   mPSACWS=0.0  
   mPSACWI=0.0  
   mQMULTS=0.0  
   mQMULTG=0.0  
   mPSACWG=0.0  
   mPGSACW=0.0  
   mPRD=0.0  
   mPRCI=0.0  
   mPRAI=0.0  
   mQMULTR=0.0  
   mQMULTRG=0.0  
   mMNUCCD=0.0  
   mPRACI=0.0  
   mPRACIS=0.0  
   mEPRD=0.0 
   mMNUCCR=0.0  
   mPIACR=0.0 
   mPIACRS=0.0  
   mPGRACS=0.0 
   mPRDS=0.0  
   mEPRDS=0.0  
   mPSACR=0.0 
   mPRDG=0.0  
   mEPRDG=0.0

   mNPRC1=0.0
   mNRAGG=0.0
   mNPRACG=0.0 
   mNSUBR=0.0 
   mNSMLTR=0.0 
   mNGMLTR=0.0 
   mNPRACS=0.0 
   mNNUCCR=0.0 
   mNIACR=0.0 
   mNIACRS=0.0 
   mNGRACS=0.0
   mNSMLTS=0.0 
   mNSAGG=0.0 
   mNPRCI=0.0 
   mNSCNG=0.0 
   mNSUBS=0.0
   
   mPCC=0.0
   mNNUCCC=0.0 
   mNPSACWS=0.0 
   mNPRA=0.0 
   mNPRC=0.0 
   mNPSACWI=0.0 
   mNPSACWG=0.0 
   mNPRAI=0.0
   mNMULTS=0.0 
   mNMULTG=0.0 
   mNMULTR=0.0 
   mNMULTRG=0.0 
   mNNUCCD=0.0
   mNSUBI=0.0 
   mNGMLTG=0.0 
   mNSUBG=0.0 
   mNACT=0.0 

   mSIZEFIX_NR=0.0 
   mSIZEFIX_NC=0.0
   mSIZEFIX_NI=0.0
   mSIZEFIX_NS=0.0
   mSIZEFIX_NG=0.0
   mNEGFIX_NR=0.0
   mNEGFIX_NC=0.0
   mNEGFIX_NI=0.0
   mNEGFIX_NS=0.0
   mNEGFIX_NG=0.0
   mNIM_MORR_CL = 0.0

   mQC_INST=0.0 
   mQR_INST=0.0 
   mQI_INST=0.0 
   mQS_INST=0.0 
   mQG_INST=0.0
   mNC_INST=0.0 
   mNR_INST=0.0 
   mNI_INST=0.0 
   mNS_INST=0.0 
   mNG_INST=0.0

#if defined(UWM_STATS) || defined(LASSO_ENA_3D_MICRO_OUTPUT)

   mPRE_3D=0.0
   mPRA_3D=0.0
   mPRC_3D=0.0
   mPRACS_3D=0.0
   mMNUCCR_3D=0.0
   mQMULTR_3D=0.0
   mQMULTRG_3D=0.0
   mPIACR_3D=0.0
   mPIACRS_3D=0.0
   mPRACG_3D=0.0
   mPGRACS_3D=0.0
   mPSMLT_3D=0.0
   mPGMLT_3D=0.0
   mQR_INST_3D=0.0

   mSIZEFIX_NR_3D=0.0
   mNEGFIX_NR_3D=0.0
   mNSUBR_3D=0.0
   mNSMLTR_3D=0.0
   mNGMLTR_3D=0.0
   mNPRC1_3D=0.0
   mNPRACS_3D=0.0
   mNNUCCR_3D=0.0
   mNRAGG_3D=0.0
   mNIACR_3D=0.0
   mNIACRS_3D=0.0
   mNPRACG_3D=0.0
   mNGRACS_3D=0.0
   mNR_INST_3D=0.0
   
   rain_vel_3D=0.0
   EFFR_3D=0.0

   mPSACWG_3D=0.0
   mPGSACW_3D=0.0
   mPRDG_3D=0.0
   mEPRDG_3D=0.0
   mPRACI_3D=0.0
   mPSACR_3D=0.0
   mEVPMG_3D=0.0
   mQG_INST_3D=0.0

   mPRAI_3D=0.0
   mPSACWS_3D=0.0
   mPRDS_3D=0.0
   mPRCI_3D=0.0
   mEPRDS_3D=0.0
   mPRACIS_3D=0.0
   mEVPMS_3D=0.0
   mQS_INST_3D=0.0

   mPRD_3D=0.0
   mEPRD_3D=0.0
   mPSACWI_3D=0.0
   mMNUCCC_3D=0.0
   mQMULTS_3D=0.0
   mQMULTG_3D=0.0
   mMNUCCD_3D=0.0
   mQI_INST_3D=0.0

   mSIZEFIX_NG_3D=0.0
   mNEGFIX_NG_3D=0.0
   mNGMLTG_3D=0.0
   mNSCNG_3D=0.0
   mNSUBG_3D=0.0
   mNG_INST_3D=0.0

   mSIZEFIX_NS_3D=0.0
   mNEGFIX_NS_3D=0.0
   mNSMLTS_3D=0.0
   mNSAGG_3D=0.0
   mNPRCI_3D=0.0
   mNSUBS_3D=0.0
   mNS_INST_3D=0.0

   mSIZEFIX_NI_3D=0.0
   mNEGFIX_NI_3D=0.0
   mNNUCCC_3D=0.0
   mNPRAI_3D=0.0
   mNMULTS_3D=0.0
   mNMULTG_3D=0.0
   mNMULTR_3D=0.0
   mNMULTRG_3D=0.0
   mNNUCCD_3D=0.0
   mNSUBI_3D=0.0
   mNIM_MORR_CL_3D=0.0
   mNI_INST_3D=0.0

   mSTENDQR_3D=0.0
   mSTENDQG_3D=0.0
   mSTENDQCI_3D=0.0
   mSTENDQS_3D=0.0
   mSTENDNR_3D=0.0
   mSTENDNCI_3D=0.0
   mSTENDNG_3D=0.0
   mSTENDNS_3D=0.0
#endif

#ifdef UWM_STATS

   !weberjk(UWM) Set 'before' values for budget statistics        
   do k=1,nzm 
     do j=dimy1_s,dimy2_s
       do i=dimx1_s,dimx2_s
         t_before(i,j,k) = t(i,j,k)
         qt_before(i,j,k) = micro_field(i,j,k,iqv)
       end do 
     end do 
   end do 

#endif /*UWM_STATS*/
end if
stend(:,:) = 0.
mksed(:,:) = 0.

!!$if(doprecip) total_water_prec = total_water_prec + total_water()
 
do j = 1,ny
   do i = 1,nx

      ! zero out mixing ratios of microphysical species
      tmpqv(:) = 0.
      tmpqcl(:) = 0.
      tmpncl(:) = 0.
      tmpqr(:) = 0.
      tmpnr(:) = 0.
      tmpqci(:) = 0.
      tmpnci(:) = 0.
      tmpqs(:) = 0.
      tmpns(:) = 0.
      tmpqg(:) = 0.
      tmpng(:) = 0.

#if defined(LASSO_ENA) && defined(CALC_RADAR_REFL)
      refl1d(:) = 0. 
#endif
      ! get microphysical quantities in this grid column
      tmpqv(:) = micro_field(i,j,:,iqv) !bloss/qt: This is total water (qv+qcl)

!bloss/qt: compute below from saturation adjustment.
!bloss/qt      tmpqcl(:) = micro_field(i,j,:,iqcl)
      if(dopredictNc) tmpncl(:) = micro_field(i,j,:,incl)
      if(doprecip) then
         tmpqr(:) = micro_field(i,j,:,iqr)
         tmpnr(:) = micro_field(i,j,:,inr)
      end if

      if(doicemicro) then
         tmpqci(:) = micro_field(i,j,:,iqci)
         tmpnci(:) = micro_field(i,j,:,inci)
         tmpqs(:) = micro_field(i,j,:,iqs)
         tmpns(:) = micro_field(i,j,:,ins)
         if(dograupel) then
            tmpqg(:) = micro_field(i,j,:,iqg)
            tmpng(:) = micro_field(i,j,:,ing)
         end if
      end if

#ifdef PNNL_STATS
      do k=1,nzm 
!MWSWong: add budget for THL and QTOG
         thel_before(i,j,k) = t2thetal( t(i,j,k), gamaz(k), tmpqr(k), &
                                        tmpqci(k), tmpqs(k)+tmpqg(k), &
                                        prespot(k) )
         qtog_before(i,j,k) = micro_field(i,j,k,iqv)
         do n=2,nmicro_fields ! prevents accessing unreferenced memory 
            qtog_before(i,j,k) = qtog_before(i,j,k) + flag_wmass(n)*micro_field(i,j,k,n)
         end do
      end do 
#endif /*PNNL_STATS*/ 
      ! get absolute temperature in this column
      !bloss/qt: before saturation adjustment for liquid,
      !          this is Tcl = T - (L/Cp)*qcl (the cloud liquid water temperature)
      tmptabs(:) = t(i,j,:)  &           ! liquid water-ice static energy over Cp
           - gamaz(:) &                                   ! potential energy
           + fac_cond * (tmpqr(:)) &    ! bloss/qt: liquid latent energy due to rain only
           + fac_sub  * (tmpqci(:) + tmpqs(:) + tmpqg(:)) ! ice latent energy

      tmpdz = adz(:)*dz
!      tmpw = 0.5*(w(i,j,1:nzm) + w(i,j,2:nz))  ! MK: changed for stretched grids 
      tmpw = ((zi(2:nz)-z(1:nzm))*w(i,j,1:nzm)+ &
             (z(1:nzm)-zi(1:nzm))*w(i,j,2:nz))/(zi(2:nz)-zi(1:nzm))
#ifdef CLUBB
      ! Added by dschanen on 4 Nov 2008 to account for w_sgs 
      if ( doclubb .and. dosubgridw ) then
        ! Compute w_sgs.  Formula is consistent with that used with 
        ! TKE from MYJ pbl scheme in WRF (see module_mp_graupel.f90).
        tmpwsub = sqrt( LIN_INT( real( wp2(i,j,2:nz) ), real( wp2(i,j,1:nzm) ), &
                                  zi(2:nz), zi(1:nzm), z(1:nzm) ) )
      else
        tmpwsub = 0.
      end if
#else /* Old code */
      tmpwsub = 0.
#endif
#ifdef CLUBB
      if ( doclubb ) then
        cloud_frac_in(1:nzm) = cloud_frac(i,j,2:nz)
      else
        cloud_frac_in(1:nzm) = 0.0
      end if
#endif

      tmppres(:) = 100.*pres(1:nzm)

      !bloss/qt: saturation adjustment to compute cloud liquid water content.
      !          Note: tmpqv holds qv+qcl on input, qv on output.
      !                tmptabs hold T-(L/Cp)*qcl on input, T on output.
      !                tmpqcl hold qcl on output.
      !                tmppres is unchanged on output, should be in Pa.
#ifdef CLUBB
      ! In the CLUBB case, we want to call the microphysics on sub-saturated grid
      ! boxes and weight by cloud fraction, therefore we use the CLUBB value of 
      ! liquid water. -dschanen 23 Nov 2009
      if ( .not. ( docloud .or. dosmoke ) ) then
        tmpqcl  = cloudliq(i,j,:) ! Liquid updated by CLUBB just prior to this
        tmpqv   = tmpqv - tmpqcl ! Vapor
        tmptabs = tmptabs + fac_cond * tmpqcl ! Update temperature
      else
        call satadj_liquid(nzm,tmptabs,tmpqv,tmpqcl,tmppres)
      end if
#else
      call satadj_liquid(nzm,tmptabs,tmpqv,tmpqcl,tmppres)
#endif

      i1 = 1 ! dummy variables used by WRF convention in subroutine call
      i2 = 1
      j1 = 1
      j2 = 1

      mtendqv = 0.
      mtendqcl = 0.
      mtendqr = 0.
      mtendqci = 0.
      mtendqs = 0.
      mtendqg = 0.
      mtendncl = 0.
      mtendnr = 0.
      mtendnci = 0.
      mtendns = 0.
      mtendng = 0.

      tmtend1d = 0.

      sfcpcp = 0.
      sfcicepcp = 0.

      effc1d(:) = 10. ! default liquid and ice effective radii
      effi1d(:) = 75.

! Zero out microphysical process rates  +++mhwang 
      PRC=0.0 
      PRA=0.0 
      PSMLT=0.0 
      EVPMS=0.0 
      PRACS=0.0 
      EVPMG=0.0 
      PRACG=0.0 
      PRE=0.0 
      PGMLT=0.0 
      MNUCCC=0.0 
      PSACWS=0.0 
      PSACWI=0.0 
      QMULTS=0.0 
      QMULTG=0.0 
      PSACWG=0.0 
      PGSACW=0.0 
      PRD=0.0 
      PRCI=0.0 
      PRAI=0.0 
      QMULTR=0.0 
      QMULTRG=0.0 
      MNUCCD=0.0 
      PRACI=0.0 
      PRACIS=0.0 
      EPRD=0.0 
      MNUCCR=0.0 
      PIACR=0.0 
      PIACRS=0.0 
      PGRACS=0.0 
      PRDS=0.0 
      EPRDS=0.0 
      PSACR=0.0 
      PRDG=0.0 
      EPRDG=0.0 

      NPRC1=0.0 
      NRAGG=0.0 
      NPRACG=0.0 
      NSUBR=0.0 
      NSMLTR=0.0 
      NGMLTR=0.0 
      NPRACS=0.0 
      NNUCCR=0.0 
      NIACR=0.0 
      NIACRS=0.0 
      NGRACS=0.0
      NSMLTS=0.0
      NSAGG=0.0 
      NPRCI=0.0 
      NSCNG=0.0 
      NSUBS=0.0

      PCC=0.0
      NNUCCC=0.0
      NPSACWS=0.0
      NPRA=0.0
      NPRC=0.0
      NPSACWI=0.0
      NPSACWG=0.0
      NPRAI=0.0
      NMULTS=0.0
      NMULTG=0.0
      NMULTR=0.0
      NMULTRG=0.0
      NNUCCD=0.0
      NSUBI=0.0
      NGMLTG=0.0
      NSUBG=0.0
      NACT=0.0

      SIZEFIX_NR=0.0
      SIZEFIX_NC=0.0
      SIZEFIX_NI=0.0
      SIZEFIX_NS=0.0
      SIZEFIX_NG=0.0
      NEGFIX_NR=0.0
      NEGFIX_NC=0.0
      NEGFIX_NI=0.0
      NEGFIX_NS=0.0
      NEGFIX_NG=0.0
      NIM_MORR_CL = 0.0

      QC_INST=0.0
      QR_INST=0.0
      QI_INST=0.0
      QS_INST=0.0
      QG_INST=0.0
      NC_INST=0.0
      NR_INST=0.0
      NI_INST=0.0
      NS_INST=0.0
      NG_INST=0.0


#ifdef CLUBB
      if ( any( tmpqv < 0. ) ) then
        qv_clip(2:nz) = tmpqv(1:nzm)
        qv_clip(1) = 0.0_core_rknd
        if ( clubb_at_least_debug_level_api( 1 ) ) then
          write(0,*) "M2005 has received a negative water vapor"
        end if
        call fill_holes_vertical_api( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qv_clip )
        tmpqv = qv_clip(2:nz)
      end if
      if ( any( tmpqcl < 0. ) ) then
        qcl_clip(2:nz) = tmpqcl(1:nzm)
        qcl_clip(1) = 0.0_core_rknd
        if ( clubb_at_least_debug_level_api( 1 ) ) then
          write(0,*) "M2005 has received a negative liquid water"
        end if
        call fill_holes_vertical_api( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qcl_clip )
        tmpqcl = qcl_clip(2:nz)
      end if

      ! Set autoconversion and accretion rates to 0;  these are diagnostics and
      ! don't feed back into the calculation.
      PRC = 0.
      PRA = 0.
      PRE = 0.
#endif /*CLUBB*/
      
      ! explanation of variable names:
      !   mtend1d: array of 1d profiles of microphysical tendencies (w/ sed.)
      !   stend1d: array of 1d profiles of sedimentation tendencies for q*
      !   tmp**: on input, current value of **.  On output, new value of **.
      !   eff*1d: one-dim. profile of effective raduis for *
      call m2005micro_graupel(&
           mtendqcl,mtendqci,mtendqs,mtendqr, &
           mtendncl,mtendnci,mtendns,mtendnr, &
           tmpqcl,tmpqci,tmpqs,tmpqr, &
           tmpncl,tmpnci,tmpns,tmpnr, &
           tmtend1d,mtendqv, &
           tmptabs,tmpqv,tmppres,rho,tmpdz,tmpw,tmpwsub, &
#if defined(CLUBB) || defined(UWM_STATS)           
           rain_vel,&
#endif /*CLUBB or UWM_STATS*/
           sfcpcp, sfcicepcp, &
           effc1d,effi1d,effs1d,effr1d, &
           dtn, &
           i1,i2, j1,j2, 1,nzm, i1,i2, j1,j2, 1,nzm, &
           mtendqg,mtendng,tmpqg,tmpng,effg1d,stendqg, &
           stendqr,stendqci,stendqs,stendqcl, &
           stendng, stendnr, stendnci, stendns, stendncl, &
#ifdef CLUBB
           cloud_frac_in, & ! cloud_frac added by dschanen UWM
#endif /*CLUBB*/
#if defined(LASSO_ENA) && defined(CALC_RADAR_REFL)
           refl1d, calc_radar_refl, &
#endif
           ! These variables were renamed during the merging for
           ! clubb:ticket:552
           ! qr_auto = PRC, qr_accr = PRA, qr_evap = PRE
           PRC,  PRA,  &
           PSMLT, EVPMS, PRACS, EVPMG, PRACG, PRE, PGMLT,  &
           MNUCCC, PSACWS, PSACWI, QMULTS, QMULTG, PSACWG, PGSACW, &
           PRD, PRCI, PRAI, QMULTR, QMULTRG, MNUCCD, PRACI, PRACIS, EPRD, &
           MNUCCR, PIACR, PIACRS, PGRACS, PRDS, EPRDS, PSACR, &
           PRDG, EPRDG, NPRC1, NRAGG, NPRACG, NSUBR, NSMLTR, NGMLTR, &
           NPRACS, NNUCCR, NIACR, NIACRS, NGRACS, NSMLTS, NSAGG, NPRCI, NSCNG, NSUBS, &
           PCC, NNUCCC, NPSACWS, NPRA, NPRC, NPSACWI, NPSACWG, NPRAI, &
           NMULTS, NMULTG, NMULTR, NMULTRG, NNUCCD, NSUBI, NGMLTG, NSUBG, NACT, &
           SIZEFIX_NR, SIZEFIX_NC, SIZEFIX_NI, SIZEFIX_NS, SIZEFIX_NG, &
           NEGFIX_NR, NEGFIX_NC, NEGFIX_NI, NEGFIX_NS, NEGFIX_NG, &
           NIM_MORR_CL, QC_INST, QR_INST, QI_INST, QS_INST, QG_INST, &
           NC_INST, NR_INST, NI_INST, NS_INST, NG_INST   )
#ifdef CLUBB
      if ( any( tmpqv < 0. ) ) then
        qv_clip(2:nz) = tmpqv(1:nzm)
        qv_clip(1) = 0.0_core_rknd
        if ( clubb_at_least_debug_level_api( 1 ) ) then
          write(0,*) "M2005 has produced a negative water vapor"
        end if
        call fill_holes_vertical_api( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qv_clip )
        tmpqv = qv_clip(2:nz)
      end if
      if ( any( tmpqcl < 0. ) ) then
        qcl_clip(2:nz) = tmpqcl(1:nzm)
        qcl_clip(1) = 0.0_core_rknd
        if ( clubb_at_least_debug_level_api( 1 ) ) then
          write(0,*) "M2005 has produced a negative liquid water"
        end if
        call fill_holes_vertical_api( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qcl_clip )
        tmpqcl = qcl_clip(2:nz)
      end if
#endif /*CLUBB*/

#ifdef UWM_STATS
! We want the tendencies to be explicity: (before - after) / dt
! Therefore, we overwrite the output from Morrison
if(dostatis) then

  mtendqcl(:) = ( tmpqcl(:) - cloudliq(i,j,:) ) / dtn
  mtendqv(:) = ( tmpqv(:) - ( micro_field(i,j,:,iqv) - cloudliq(i,j,:) ) ) / dtn

if(doprecip .EQV. .FALSE.) then
  tmpqcl(:) = tmpqcl(:) + tmpqr(:) ! add rain mass back to cloud water
  tmpncl(:) = tmpncl(:) + tmpnr(:) ! add rain number back to cloud water

  ! zero out rain 
  tmpqr(:) = 0.
  tmpnr(:) = 0.
end if !doprecip is false

  if(dopredictNc) then
    mtendncl(:) = (tmpncl(:) - micro_field(i,j,:,incl) ) / dtn
  end if

  if(doprecip) then
    mtendnr(:) = (tmpnr(:) - micro_field(i,j,:,inr) )/ dtn
    mtendqr(:) = (tmpqr(:) - micro_field(i,j,:,iqr) )/ dtn
  end if

  if(doicemicro) then
    mtendnci(:) = (tmpnci(:) - micro_field(i,j,:,inci) ) / dtn
    mtendqci(:) = (tmpqci(:) - micro_field(i,j,:,iqci) ) / dtn
    
    mtendns(:) = (tmpns(:) - micro_field(i,j,:,ins) ) / dtn
    mtendqs(:) = (tmpqs(:) - micro_field(i,j,:,iqs) ) / dtn
    
    if(dograupel) then
      mtendng(:) = (tmpng(:) - micro_field(i,j,:,ing) ) / dtn
      mtendqg(:) = (tmpqg(:) - micro_field(i,j,:,iqg) ) / dtn
    endif

  endif
endif
#endif /*UWM_STATS*/

     ! update microphysical quantities in this grid column
      if(doprecip) then
         total_water_prec = total_water_prec + sfcpcp

         ! take care of surface precipitation
         precsfc(i,j) = precsfc(i,j) + sfcpcp/dz
         prec_xy(i,j) = prec_xy(i,j) + sfcpcp/dz

         ! update rain
         micro_field(i,j,:,iqr) = tmpqr(:)
         micro_field(i,j,:,inr) = tmpnr(:)
      else
         ! add rain to cloud
         tmpqcl(:) = tmpqcl(:) + tmpqr(:) ! add rain mass back to cloud water
         tmpncl(:) = tmpncl(:) + tmpnr(:) ! add rain number back to cloud water

         ! zero out rain 
         tmpqr(:) = 0.
         tmpnr(:) = 0.

         ! add rain tendencies to cloud
         stendqcl(:) = stendqcl(:) + stendqr(:)
         stendncl(:) = stendncl(:) + stendnr(:)
         mtendqcl(:) = mtendqcl(:) + mtendqr(:)
         mtendncl(:) = mtendncl(:) + mtendnr(:)

         ! zero out rain tendencies
         stendqr(:) = 0.
         stendnr(:) = 0.
         mtendqr(:) = 0.
         mtendnr(:) = 0.
      end if

      !bloss/qt: update total water and cloud liquid.
      !          Note: update of total water moved to after if(doprecip),
      !                  since no precip moves rain --> cloud liq.
      micro_field(i,j,:,iqv) = tmpqv(:) + tmpqcl(:) !bloss/qt: total water
      cloudliq(i,j,:) = tmpqcl(:) !bloss/qt: auxilliary cloud liquid water variable
      if(dopredictNc) micro_field(i,j,:,incl) = tmpncl(:)
      reffc(i,j,:) = effc1d(:)
#if defined(LASSO_ENA) && defined(CALC_RADAR_REFL)
      if (calc_radar_refl) refl(i,j,:) = refl1d(:)
#endif

      if(doicemicro) then
         micro_field(i,j,:,iqci) = tmpqci(:)
         micro_field(i,j,:,inci) = tmpnci(:)
         micro_field(i,j,:,iqs) = tmpqs(:)
         micro_field(i,j,:,ins) = tmpns(:)
         if(dograupel) then
            micro_field(i,j,:,iqg) = tmpqg(:)
            micro_field(i,j,:,ing) = tmpng(:)
         end if
         reffi(i,j,:) = effi1d(:)  
      end if

      !=====================================================
      ! update liquid-ice static energy due to precipitation
      t(i,j,:) = t(i,j,:) &
           - dtn*fac_cond*(stendqcl+stendqr) &
           - dtn*fac_sub*(stendqci+stendqs+stendqg)
      !=====================================================

      if(dostatis) then

#ifdef PNNL_STATS
      !MWSWong: THL budget
      !thel(i,j,1:nzm) = prespot(1:nzm)*( t(i,j,1:nzm) - gamaz(1:nzm)  &
      !                 + fac_cond * micro_field(i,j,:,iqr) &
      !                 + fac_sub  * (micro_field(i,j,:,iqci)+micro_field(i,j,:,iqs)+micro_field(i,j,:,iqg)))
      thel(i,j,1:nzm) = t2thetal( t(i,j,1:nzm), gamaz(1:nzm), tmpqr(1:nzm), &
                                  tmpqci(1:nzm), tmpqs(1:nzm)+tmpqg(1:nzm), &
                                  prespot(1:nzm) )

      thellat(1:nzm) = thellat(1:nzm) + (thel(i,j,1:nzm)-thel_before(i,j,1:nzm))
#endif /*PNNL_STATS*/

!bloss/qt: total water microphysical tendency includes qv and qcl
         mtend(:,iqv) = mtend(:,iqv) + mtendqv + mtendqcl

!bloss/qt         mtend(:,iqcl) = mtend(:,iqcl) + mtendqcl
         if(dopredictNc) mtend(:,incl) = mtend(:,incl) + mtendncl
         if(doprecip) then
            mtend(:,iqr) = mtend(:,iqr) + mtendqr
            mtend(:,inr) = mtend(:,inr) + mtendnr
         end if

         if(doicemicro) then
            mtend(:,iqci) = mtend(:,iqci) + mtendqci
            mtend(:,inci) = mtend(:,inci) + mtendnci
            !bloss            stend(:,inci) = stend(:,inci) + stendnci

            mtend(:,iqs) = mtend(:,iqs) + mtendqs
            mtend(:,ins) = mtend(:,ins) + mtendns
            !bloss            stend(:,ins) = stend(:,ins) + stendns

            if(dograupel) then
               mtend(:,iqg) = mtend(:,iqg) + mtendqg
               mtend(:,ing) = mtend(:,ing) + mtendng
               !bloss            stend(:,ing) = stend(:,ing) + stendng
            end if
         end if

         do n = 1,nmicro_fields
            do k = 1,nzm
               if(micro_field(i,j,k,n).ge.1.e-6) mfrac(k,n) = mfrac(k,n)+1.
            end do
         end do

         ! approximate optical depth = 0.0018*lwp/effrad
         !  integrated up to level at which output
         tmpc = 0.
         tmpr = 0.
         tmpi = 0.
         tmps = 0.
         tmpg = 0.

         do k = 1,nzm
            tmpc = tmpc + 0.0018*rho(k)*dz*adz(k)*tmpqcl(k)/(1.e-20+1.e-6*effc1d(k))
            tmpr = tmpr + 0.0018*rho(k)*dz*adz(k)*tmpqr(k)/(1.e-20+1.e-6*effr1d(k))
            !bloss/qt: put cloud liquid optical depth in trtau(:,iqv)
            trtau(k,iqv) = trtau(k,iqv) + tmpc
            if(doprecip) trtau(k,iqr) = trtau(k,iqr) + tmpr

            if(doicemicro) then
               tmpi = tmpi + 0.0018*rho(k)*dz*adz(k)*tmpqci(k)/(1.e-20+1.e-6*effi1d(k))
               tmps = tmps + 0.0018*rho(k)*dz*adz(k)*tmpqs(k)/(1.e-20+1.e-6*effs1d(k))
               tmpg = tmpg + 0.0018*rho(k)*dz*adz(k)*tmpqg(k)/(1.e-20+1.e-6*effg1d(k))

               trtau(k,iqci) = trtau(k,iqci) + tmpi
               trtau(k,iqs) = trtau(k,iqs) + tmps
#if defined( CLUBB ) || defined( UWM_MISC ) /* Bug fix -dschanen 9 Mar 2012 */
               if ( dograupel ) then
                 trtau(k,iqg) = trtau(k,iqg) + tmpg
               end if
#else
               trtau(k,iqg) = trtau(k,iqg) + tmpg
#endif /* CLUBB */
            end if
         end do

         tlat(1:nzm) = tlat(1:nzm) &
              - dtn*fac_cond*(stendqcl+stendqr) &
              - dtn*fac_sub*(stendqci+stendqs+stendqg)

         qpfall(1:nzm) = qpfall(1:nzm) + dtn*(stendqr+stendqs+stendqg)

         !bloss: temperature tendency (sensible heating) due to phase changes
         tmtend3d(i,j,1:nzm) = tmtend1d(1:nzm)

! Microphysical process rates +++mhwang 
         mPRC=mPRC + PRC 
         mPRA=mPRA + PRA 
         mPSMLT=mPSMLT + PSMLT 
         mEVPMS=mEVPMS + EVPMS 
         mPRACS=mPRACS + PRACS 
         mEVPMG=mEVPMG + EVPMG 
         mPRACG=mPRACG + PRACG 
         mPRE=mPRE + PRE 
         mPGMLT=mPGMLT + PGMLT 
         mMNUCCC=mMNUCCC + MNUCCC 
         mPSACWS=mPSACWS + PSACWS 
         mPSACWI=mPSACWI + PSACWI 
         mQMULTS=mQMULTS + QMULTS 
         mQMULTG=mQMULTG + QMULTG 
         mPSACWG=mPSACWG + PSACWG 
         mPGSACW=mPGSACW + PGSACW 
         mPRD=mPRD + PRD 
         mPRCI=mPRCI + PRCI 
         mPRAI=mPRAI + PRAI 
         mQMULTR=mQMULTR + QMULTR 
         mQMULTRG=mQMULTRG + QMULTRG 
         mMNUCCD=mMNUCCD + MNUCCD 
         mPRACI=mPRACI + PRACI 
         mPRACIS=mPRACIS + PRACIS 
         mEPRD=mEPRD + EPRD  
         mMNUCCR=mMNUCCR + MNUCCR 
         mPIACR=mPIACR + PIACR 
         mPIACRS=mPIACRS + PIACRS 
         mPGRACS=mPGRACS + PGRACS 
         mPRDS=mPRDS + PRDS 
         mEPRDS=mEPRDS + EPRDS 
         mPSACR=mPSACR + PSACR 
         mPRDG=mPRDG + PRDG 
         mEPRDG=mEPRDG + EPRDG 

         mNPRC1=mNPRC1+NPRC1
         mNRAGG=mNRAGG+NRAGG
         mNPRACG=mNPRACG+NPRACG
         mNSUBR=mNSUBR+NSUBR
         mNSMLTR=mNSMLTR+NSMLTR
         mNGMLTR=mNGMLTR+NGMLTR
           
         mNPRACS=mNPRACS+NPRACS
         mNNUCCR=mNNUCCR+NNUCCR
         mNIACR=mNIACR+NIACR
         mNIACRS=mNIACRS+NIACRS
         mNGRACS=mNGRACS+NGRACS
         mNSMLTS=mNSMLTS+NSMLTS
         mNSAGG=mNSAGG+NSAGG
         mNPRCI=mNPRCI+NPRCI
         mNSCNG=mNSCNG+NSCNG
         mNSUBS=mNSUBS+NSUBS
         mPCC=mPCC+PCC
         mNNUCCC=mNNUCCC+NNUCCC
         mNPSACWS=mNPSACWS+NPSACWS
         mNPRA=mNPRA+NPRA
         mNPRC=mNPRC+NPRC
         mNPSACWI=mNPSACWI+NPSACWI
         mNPSACWG=mNPSACWG+NPSACWG
         mNPRAI=mNPRAI+NPRAI
         mNMULTS=mNMULTS+NMULTS
         mNMULTG=mNMULTG+NMULTG
         mNMULTR=mNMULTR+NMULTR
         mNMULTRG=mNMULTRG+NMULTRG
         mNNUCCD=mNNUCCD+NNUCCD
         mNSUBI=mNSUBI+NSUBI
         mNGMLTG=mNGMLTG+NGMLTG
         mNSUBG=mNSUBG+NSUBG
         mNACT=mNACT+NACT
         mSIZEFIX_NR=mSIZEFIX_NR+SIZEFIX_NR
         mSIZEFIX_NC=mSIZEFIX_NC+SIZEFIX_NC
         mSIZEFIX_NI=mSIZEFIX_NI+SIZEFIX_NI
         mSIZEFIX_NS=mSIZEFIX_NS+SIZEFIX_NS
         mSIZEFIX_NG=mSIZEFIX_NG+SIZEFIX_NG
         mNEGFIX_NR=mNEGFIX_NR+NEGFIX_NR
         mNEGFIX_NC=mNEGFIX_NC+NEGFIX_NC
         mNEGFIX_NI=mNEGFIX_NI+NEGFIX_NI
         mNEGFIX_NS=mNEGFIX_NS+NEGFIX_NS
         mNEGFIX_NG=mNEGFIX_NG+NEGFIX_NG
         mNIM_MORR_CL=mNIM_MORR_CL+NIM_MORR_CL    

         mQC_INST=mQC_INST+QC_INST
         mQR_INST=mQR_INST+QR_INST
         mQI_INST=mQI_INST+QI_INST
         mQS_INST=mQS_INST+QS_INST
         mQG_INST=mQG_INST+QG_INST
         mNC_INST=mNC_INST+NC_INST
         mNR_INST=mNR_INST+NR_INST
         mNI_INST=mNI_INST+NI_INST
         mNS_INST=mNS_INST+NS_INST
         mNG_INST=mNG_INST+NG_INST

#if defined(UWM_STATS) || defined(LASSO_ENA_3D_MICRO_OUTPUT)
   !3d Budgets
   mPRE_3D(i,j,:)=mPRE_3D(i,j,:) + PRE
   mPRA_3D(i,j,:)=mPRA_3D(i,j,:) + PRA
   mPRC_3D(i,j,:)=mPRC_3D(i,j,:) + PRC
   mPRACS_3D(i,j,:)=mPRACS_3D(i,j,:) + PRACS
   mMNUCCR_3D(i,j,:)=mMNUCCR_3D(i,j,:) + MNUCCR
   mQMULTR_3D(i,j,:)=mQMULTR_3D(i,j,:) + QMULTR
   mQMULTRG_3D(i,j,:)=mQMULTRG_3D(i,j,:) + QMULTRG
   mPIACR_3D(i,j,:)=mPIACR_3D(i,j,:) + PIACR
   mPIACRS_3D(i,j,:)=mPIACRS_3D(i,j,:) + PIACRS
   mPRACG_3D(i,j,:)=mPRACG_3D(i,j,:) + PRACG
   mPGRACS_3D(i,j,:)=mPGRACS_3D(i,j,:) + PGRACS
   mPSMLT_3D(i,j,:)=mPSMLT_3D(i,j,:) + PSMLT
   mPGMLT_3D(i,j,:)=mPGMLT_3D(i,j,:) + PGMLT
   mQR_INST_3D(i,j,:)= mQR_INST_3D(i,j,:) + QR_INST

   mSIZEFIX_NR_3D(i,j,:)=mSIZEFIX_NR_3D(i,j,:) + SIZEFIX_NR
   mNEGFIX_NR_3D(i,j,:)=mNEGFIX_NR_3D(i,j,:) + NEGFIX_NR
   mNSUBR_3D(i,j,:)=mNSUBR_3D(i,j,:) + NSUBR
   mNSMLTR_3D(i,j,:)=mNSMLTR_3D(i,j,:) + NSMLTR
   mNGMLTR_3D(i,j,:)=mNGMLTR_3D(i,j,:) + NGMLTR
   mNPRC1_3D(i,j,:)=mNPRC1_3D(i,j,:) + NPRC1 
   mNPRACS_3D(i,j,:)=mNPRACS_3D(i,j,:) + NPRACS
   mNNUCCR_3D(i,j,:)=mNNUCCR_3D(i,j,:) + NNUCCR
   mNRAGG_3D(i,j,:)=mNRAGG_3D(i,j,:) + NRAGG
   mNIACR_3D(i,j,:)=mNIACR_3D(i,j,:) + NIACR
   mNIACRS_3D(i,j,:)=mNIACRS_3D(i,j,:) + NIACRS
   mNPRACG_3D(i,j,:)=mNPRACG_3D(i,j,:) + NPRACG
   mNGRACS_3D(i,j,:)=mNGRACS_3D(i,j,:) + NGRACS
   mNR_INST_3D(i,j,:)=mNR_INST_3D(i,j,:) + NR_INST

   ! Additional micro. variables
   rain_vel_3D(i,j,:)=rain_vel_3D(i,j,:) + rain_vel
   EFFR_3D(i,j,:)=EFFR_3D(i,j,:) + effr1d

   mPSACWG_3D(i,j,:)=mPSACWG_3D(i,j,:) + PSACWG
   mPGSACW_3D(i,j,:)=mPGSACW_3D(i,j,:) + PGSACW
   mPRDG_3D(i,j,:)=mPRDG_3D(i,j,:) + PRDG
   mEPRDG_3D(i,j,:)=mEPRDG_3D(i,j,:) + EPRDG
   mPRACI_3D(i,j,:)=mPRACI_3D(i,j,:) + PRACI
   mPSACR_3D(i,j,:)=mPSACR_3D(i,j,:) + PSACR
   mEVPMG_3D(i,j,:)=mEVPMG_3D(i,j,:) + EVPMG
   mQG_INST_3D(i,j,:)=mQG_INST_3D(i,j,:) + QG_INST

   mPRAI_3D(i,j,:)=mPRAI_3D(i,j,:) + PRAI
   mPSACWS_3D(i,j,:)=mPSACWS_3D(i,j,:) + PSACWS
   mPRDS_3D(i,j,:)=mPRDS_3D(i,j,:) + PRDS
   mPRCI_3D(i,j,:)=mPRCI_3D(i,j,:) + PRCI
   mEPRDS_3D(i,j,:)=mEPRDS_3D(i,j,:) + EPRDS
   mPRACIS_3D(i,j,:)=mPRACIS_3D(i,j,:) + PRACIS
   mEVPMS_3D(i,j,:)=mEVPMS_3D(i,j,:) + EVPMS
   mQS_INST_3D(i,j,:)=mQS_INST_3D(i,j,:) + QS_INST

   mPRD_3D(i,j,:)=mPRD_3D(i,j,:) + PRD
   mEPRD_3D(i,j,:)=mEPRD_3D(i,j,:) + EPRD
   mPSACWI_3D(i,j,:)=mPSACWI_3D(i,j,:) + PSACWI
   mMNUCCC_3D(i,j,:)=mMNUCCC_3D(i,j,:) + MNUCCC
   mQMULTS_3D(i,j,:)=mQMULTS_3D(i,j,:) + QMULTS
   mQMULTG_3D(i,j,:)=mQMULTG_3D(i,j,:) + QMULTG
   mMNUCCD_3D(i,j,:)=mMNUCCD_3D(i,j,:) + MNUCCD
   mQI_INST_3D(i,j,:)=mQI_INST_3D(i,j,:) + QI_INST

   mSIZEFIX_NG_3D(i,j,:)=mSIZEFIX_NG_3D(i,j,:) + SIZEFIX_NG
   mNEGFIX_NG_3D(i,j,:)=mNEGFIX_NG_3D(i,j,:) + NEGFIX_NG
   mNGMLTG_3D(i,j,:)=mNGMLTG_3D(i,j,:) + NGMLTG
   mNSCNG_3D(i,j,:)=mNSCNG_3D(i,j,:) + NSCNG
   mNSUBG_3D(i,j,:)=mNSUBG_3D(i,j,:) + NSUBG
   mNG_INST_3D(i,j,:)=mNG_INST_3D(i,j,:) + NG_INST

   mSIZEFIX_NS_3D(i,j,:)=mSIZEFIX_NS_3D(i,j,:) + SIZEFIX_NS
   mNEGFIX_NS_3D(i,j,:)=mNEGFIX_NS_3D(i,j,:) + NEGFIX_NS
   mNSMLTS_3D(i,j,:)=mNSMLTS_3D(i,j,:) + NSMLTS
   mNSAGG_3D(i,j,:)=mNSAGG_3D(i,j,:) + NSAGG
   mNPRCI_3D(i,j,:)=mNPRCI_3D(i,j,:) + NPRCI
   mNSUBS_3D(i,j,:)=mNSUBS_3D(i,j,:) + NSUBS
   mNS_INST_3D(i,j,:)=mNS_INST_3D(i,j,:) + NS_INST

   mSIZEFIX_NI_3D(i,j,:)=mSIZEFIX_NI_3D(i,j,:) + SIZEFIX_NI
   mNEGFIX_NI_3D(i,j,:)=mNEGFIX_NI_3D(i,j,:) + NEGFIX_NI
   mNNUCCC_3D(i,j,:)=mNNUCCC_3D(i,j,:) + NNUCCC
   mNPRAI_3D(i,j,:)=mNPRAI_3D(i,j,:) + NPRAI
   mNMULTS_3D(i,j,:)=mNMULTS_3D(i,j,:) + NMULTS
   mNMULTG_3D(i,j,:)=mNMULTG_3D(i,j,:) + NMULTG
   mNMULTR_3D(i,j,:)=mNMULTR_3D(i,j,:) + NMULTR
   mNMULTRG_3D(i,j,:)=mNMULTRG_3D(i,j,:) + NMULTRG
   mNNUCCD_3D(i,j,:)=mNNUCCD_3D(i,j,:) + NNUCCD
   mNSUBI_3D(i,j,:)=mNSUBI_3D(i,j,:) + NSUBI
   mNIM_MORR_CL_3D(i,j,:)=mNIM_MORR_CL_3D(i,j,:) + NIM_MORR_CL
   mNI_INST_3D(i,j,:)=mNI_INST_3D(i,j,:) + NI_INST


   mSTENDQR_3D(i,j,:)=mSTENDQR_3D(i,j,:) + STENDQR
   mSTENDQCI_3D(i,j,:)=mSTENDQCI_3D(i,j,:) + STENDQCI
   mSTENDQG_3D(i,j,:)=mSTENDQG_3D(i,j,:) + STENDQG
   mSTENDQS_3D(i,j,:)=mSTENDQS_3D(i,j,:) + STENDQS
   mSTENDNR_3D(i,j,:)=mSTENDNR_3D(i,j,:) + STENDNR
   mSTENDNCI_3D(i,j,:)=mSTENDNCI_3D(i,j,:) + STENDNCI
   mSTENDNG_3D(i,j,:)=mSTENDNG_3D(i,j,:) + STENDNG
   mSTENDNS_3D(i,j,:)=mSTENDNS_3D(i,j,:) + STENDNS


#endif /*UWM_STATS*/
    
endif !do statis

      stend(:,iqv) = stend(:,iqv) + stendqcl !bloss/qt: iqcl --> iqv
      if(dopredictNc) stend(:,incl) = stend(:,incl) + stendncl
      if(doprecip) then
         stend(:,iqr) = stend(:,iqr) + stendqr
         stend(:,inr) = stend(:,inr) + stendnr
      end if

      if(doicemicro) then
         stend(:,iqci) = stend(:,iqci) + stendqci
         stend(:,iqs) = stend(:,iqs) + stendqs
         stend(:,inci) = stend(:,inci) + stendnci
         stend(:,ins) = stend(:,ins) + stendns
         if(dograupel) stend(:,iqg) = stend(:,iqg) + stendqg
         if(dograupel) stend(:,ing) = stend(:,ing) + stendng
      end if

   end do ! i = 1,nx
end do ! j = 1,ny

#ifdef UWM_STATS      
do j = 1,ny
   do i = 1,nx
      if (doprecip) then
         tmpqr(:) = micro_field(i,j,:,iqr)
      else
         tmpqr(:) = 0.0
      endif
      ! Compute liquid water potential temperature
      theta_l(i,j,:) = t2thetal( t(i,j,:), gamaz(:), tmpqr(:), &
      ! pass in zero for the ice species ==> tmpqci(:), tmpqs(:) + tmpqg(:), prespot(:) )
                                 0.0,        0.0,                  prespot(:) )
   enddo ! i = 1, nx
enddo ! j = 1, ny

! Compute extended liquid water mixing ratio
call compute_chi_eta( theta_l, micro_field(1:nx,1:ny,1:nzm,iqv), pres, prespot,&
                      chi, eta )

#endif /*UWM_STATS*/

#if defined(ATEX) || defined(DYCOMSRF01) || defined(BOMEX) || defined(HISCALE) 
! do nothing
#elif defined(LASSO_ENA)
#ifdef LASSO_ENA_3D_MICRO_OUTPUT
if(mod(nstep,nsave3D).eq.0.and.nstep.ge.nsave3Dstart.and.nstep.le.nsave3Dend) then
  call write_3d_micro_fields()
  call write_3d_micro_fields_frzmr()
  call write_3d_micro_fields_frznc()
  call write_3d_micro_fields_sed()
endif
#endif
#else
if(mod(nstep,nsave3D).eq.0.and.nstep.ge.nsave3Dstart.and.nstep.le.nsave3Dend) then
  call write_3d_micro_fields()
  call write_3d_micro_fields_frzmr()
  call write_3d_micro_fields_frznc()
  call write_3d_micro_fields_sed()
endif
#endif

#ifdef UWM_STATS      
if(dostatis) then
        !------------------------------------------------------------------------------------------
        !weberjk(UWM), compute the microphysical effects on t2, q2, tw, and qw
        !budgets
        call stat_varscalar(t,t_before,t_avg,t_before_avg,t2leprec)
        call stat_varscalar(micro_field(:,:,:,iqv),qt_before,qt_avg,qt_avg_before,q2leprec)
         
        call setvalue(twleprec,nzm,0.)
        call stat_sw2(t,t_before,twleprec)
         
        call setvalue(qwleprec,nzm,0.)
        call stat_sw2(micro_field(:,:,:,iqv),qt_before,qwleprec)

#ifdef PNNL_STATS
        !MWSWong: add budget for THL
        do k=1,nzm 
          do j=dimy1_s,dimy2_s
            do i=dimx1_s,dimx2_s
              qtog(i,j,k) = micro_field(i,j,k,index_water_vapor)
              do n=2,nmicro_fields ! prevents accessing unreferenced memory 
                qtog(i,j,k) = qtog(i,j,k) + flag_wmass(n)*micro_field(i,j,k,n)
              end do
            end do 
          end do 
        end do 
        call stat_varscalar(thel,thel_before,thel_avg,thel_before_avg,thel2leprec)
        call stat_varscalar(qtog,qtog_before,qtog_avg,qtog_before_avg,qtog2leprec)

        call setvalue(thelwleprec,nzm,0.)
        call stat_sw2(thel,thel_before,thelwleprec)

        call setvalue(qtogwleprec,nzm,0.)
        call stat_sw2(qtog,qtog_before,qtogwleprec)

        !MWSWong: add budget for covariance r'thl' 
        do k=1,nzm
           qthelleprec(k)=0.
           do j=1,ny
              do i=1,nx
                 qthelleprec(k)=qthelleprec(k) + (thel(i,j,k)-thel_avg(k))*(micro_field(i,j,k,iqv)-qt_avg(k)) &
                              - (thel_before(i,j,k)-thel_before_avg(k))*(qt_before(i,j,k)-qt_avg_before(k))
              end do
           end do
           qthelleprec(k)=qthelleprec(k)*(1./(dtn*nx*ny))
        end do

#endif /*PNNL_STATS*/

endif !do statis
#endif /*UWM_STATS*/

! back sedimentation flux out from sedimentation tendencies
tmpc = 0.
do k = 1,nzm
   m = nz-k
   tmpc = tmpc + stend(m,iqv)*rho(m)*dz*adz(m)  !bloss/qt: iqcl --> iqv
   mksed(m,iqv) = tmpc
end do
precflux(1:nzm) = precflux(1:nzm) - mksed(:,iqv)*dtn/dz

if(doprecip) then
   tmpr = 0.
   do k = 1,nzm
      m = nz-k
      tmpr = tmpr + stend(m,iqr)*rho(m)*dz*adz(m)
      mksed(m,iqr) = tmpr
   end do
   precflux(1:nzm) = precflux(1:nzm) - mksed(:,iqr)*dtn/dz
end if

if(doicemicro) then
   tmpi = 0.
   tmps = 0.
   tmpg = 0.
   do k = 1,nzm
      m = nz-k
      tmpi = tmpi + stend(m,iqci)*rho(m)*dz*adz(m)
      tmps = tmps + stend(m,iqs)*rho(m)*dz*adz(m)
#if defined( CLUBB ) || defined( UWM_MISC ) /* Bug fix -dschanen 9 Mar 2012 */
      if ( dograupel ) then
        tmpg = tmpg + stend(m,iqg)*rho(m)*dz*adz(m)
      else
        tmpg = 0.
      end if
#else
      tmpg = tmpg + stend(m,iqg)*rho(m)*dz*adz(m)
#endif
      mksed(m,iqci) = tmpi
      mksed(m,iqs) = tmps
#if defined( CLUBB ) || defined( UWM_MISC ) /* Bug fix -dschanen 9 Mar 2012 */
      if ( dograupel ) then
        mksed(m,iqg) = tmpg
      end if
#else
      mksed(m,iqg) = tmpg
#endif
   end do
#if defined( CLUBB ) || defined( UWM_MISC ) /* Bug fix -dschanen 9 Mar 2012 */
   if ( dograupel ) then
     precflux(1:nzm) = precflux(1:nzm) &
          - (mksed(:,iqci) + mksed(:,iqs) + mksed(:,iqg))*dtn/dz
   else
     precflux(1:nzm) = precflux(1:nzm) &
          - (mksed(:,iqci) + mksed(:,iqs))*dtn/dz
   end if
#else
   precflux(1:nzm) = precflux(1:nzm) &
        - (mksed(:,iqci) + mksed(:,iqs) + mksed(:,iqg))*dtn/dz
#endif
end if

!!$if(doprecip) total_water_prec = total_water_prec - total_water()

#ifdef CLUBB
if (docloud.or.doclubb)  call micro_diagnose()   ! leave this line here
#else
if (docloud)  call micro_diagnose()   ! leave this line here
#endif

call t_stopf ('micro_proc')

end subroutine micro_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and radiation:
!
!  This is the pace where the microphysics field that SAM actually cares about
!  are diagnosed.

subroutine micro_diagnose()

use vars
#ifdef CLUBB
use clubb_api_module, only: &
  clubb_at_least_debug_level_api, & ! Procedure
  fstderr, zero_threshold
implicit none
#endif

real omn, omp
integer i,j,k

! water vapor = total water - cloud liquid
qv(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqv) &
     - cloudliq(1:nx,1:ny,1:nzm)

#ifdef CLUBB
do i = 1, nx
  do j = 1, ny
    do k = 1, nzm
      ! Apply local hole-filling to vapor by converting liquid to vapor. Moist
      ! static energy should be conserved, so updating temperature is not
      ! needed here. -dschanen 31 August 2011
      if ( qv(i,j,k) < zero_threshold ) then
        cloudliq(i,j,k) = cloudliq(i,j,k) + qv(i,j,k)
        qv(i,j,k) = zero_threshold
        if ( cloudliq(i,j,k) < zero_threshold ) then
          if ( clubb_at_least_debug_level_api( 1 ) ) then
            write(fstderr,*) "Total water at", "i =", i, "j =", j, "k =", k, "is negative.", &
              "Applying non-conservative hard clipping."
          end if
          cloudliq(i,j,k) = zero_threshold
        end if ! cloud_liq < 0
      end if ! qv < 0
    end do ! 1.. nzm
  end do ! 1.. ny
end do ! 1.. nx
#endif /* CLUBB */
! cloud liquid water
qcl(1:nx,1:ny,1:nzm) = cloudliq(1:nx,1:ny,1:nzm)

! rain water
if(doprecip) qpl(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqr)

! cloud ice 
if(doicemicro) then
   qci(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqci)

   if(dograupel) then
      qpi(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqs) &
           + micro_field(1:nx,1:ny,1:nzm,iqg)
   else
      qpi(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqs)
   end if
end if

end subroutine micro_diagnose

#ifdef CLUBB
!---------------------------------------------------------------------
subroutine micro_update()

! Description:
! This subroutine essentially does what micro_proc does but does not
! call any microphysics subroutines.  We need to do this for the 
! single-moment bulk microphysics (SAM1MOM) so that CLUBB gets a
! properly updated value of ice fed in. 
!
! -dschanen UWM
!---------------------------------------------------------------------

  ! Update the dynamical core variables (e.g. qv, qcl) with the value in
  ! micro_field.  Diffusion, advection, and other processes are applied to
  ! micro_field but not the variables in vars.f90
  call micro_diagnose()

  return
end subroutine micro_update

!---------------------------------------------------------------------
subroutine micro_adjust( new_qv, new_qc )
! Description:
!   Adjust total water in SAM based on values from CLUBB.
! References:
!   None
!---------------------------------------------------------------------

  use vars, only: qci

  implicit none

  real, dimension(nx,ny,nzm), intent(in) :: &
    new_qv, & ! Water vapor mixing ratio that has been adjusted by CLUBB [kg/kg]
    new_qc    ! Cloud water mixing ratio that has been adjusted by CLUBB [kg/kg]

  ! Total water mixing ratio
  micro_field(1:nx,1:ny,1:nzm,iqv) = new_qv(1:nx,1:ny,1:nzm) &
                                   + new_qc(1:nx,1:ny,1:nzm)

  ! Cloud water mixing ratio
  cloudliq(1:nx,1:ny,1:nzm) = new_qc(1:nx,1:ny,1:nzm) 

  return
end subroutine micro_adjust

#endif /*CLUBB*/
!----------------------------------------------------------------------
!!! functions to compute terminal velocity for precipitating variables:
!
! you need supply functions to compute terminal velocity for all of your 
! precipitating prognostic variables. Note that all functions should
! compute vertical velocity given two microphysics parameters var1, var2, 
! and temperature, and water vapor (single values, not arrays). Var1 and var2 
! are some microphysics variables like water content and concentration.
! Don't change the number of arguments or their meaning!

!!$real function term_vel_qr(qr,nr,tabs,rho)
!!$! .......  
!!$end function term_vel_qr
!!$
!!$real function term_vel_Nr(qr,nr,tabs,rho)
!!$! .......  
!!$end function term_vel_Nr
!!$
!!$real function term_vel_qs(qs,ns,tabs,rho)
!!$! .......  
!!$end function term_vel_qs

! etc.

!----------------------------------------------------------------------
!!! compute sedimentation 
!
!  The perpose of this subroutine is to prepare variables needed to call
! the precip_all() for each of the falling hydrometeor varibles
subroutine micro_precip_fall()

! before calling precip_fall() for each of falling prognostic variables,
! you need to set hydro_type and omega(:,:,:) variables.
! hydro_type can have four values:
! 0 - variable is liquid water mixing ratio
! 1 - hydrometeor is ice mixing ratio
! 2 - hydrometeor is mixture-of-liquid-and-ice mixing ratio. (As in original SAM microphysics).
! 3 - variable is not mixing ratio, but, for example, rain drop concentration
! OMEGA(:,:,:) is used only for hydro_type=2, and is the fraction of liquid phase (0-1).
! for hour hypothetical case, there is no mixed hydrometeor, so omega is not actually used.

integer hydro_type
real omega(nx,ny,nzm) 

integer i,j,k

return ! do not need this routine -- sedimentation done in m2005micro.

!!$! Initialize arrays that accumulate surface precipitation flux
!!$
!!$ if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
!!$   do j=1,ny
!!$    do i=1,nx
!!$     precsfc(i,j)=0.
!!$    end do
!!$   end do
!!$   do k=1,nzm
!!$    precflux(k) = 0.
!!$   end do
!!$ end if
!!$
!!$ do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
!!$    qpfall(k)=0.
!!$    tlat(k) = 0.
!!$ end do
!!$   
!!$! Compute sedimentation of falling variables:
!!$
!!$ hydro_type=0
!!$ call precip_fall(qr, term_vel_qr, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Nr, term_vel_Nr, hydro_type, omega)
!!$ hydro_type=1
!!$ call precip_fall(qs, term_vel_qs, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Ns, term_vel_Ns, hydro_type, omega)
!!$ hydro_type=1
!!$ call precip_fall(qg, term_vel_qg, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Ng, term_vel_Ng, hydro_type, omega)
!!$


end subroutine micro_precip_fall

!----------------------------------------------------------------------
! called when stepout() called

subroutine micro_print()
  implicit none
  integer :: k

  ! print out min/max values of all microphysical variables
  do k=1,nmicro_fields
     call fminmax_print(trim(mkname(k))//':', &
          micro_field(:,:,:,k),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
  end do

end subroutine micro_print

!----------------------------------------------------------------------
!!! Initialize the list of microphysics statistics that will be outputted
!!  to *.stat statistics file

subroutine micro_hbuf_init(namelist,deflist,unitlist,status,average_type,count,microcount)


#ifdef UWM_STATS
use compute_correlation_module, only: &
    corr_avg
integer :: idx_frac_fields ! Index for each hydrometeor fraction
#endif /*UWM_STATS*/
character(*) namelist(*), deflist(*), unitlist(*)
integer status(*),average_type(*),count,microcount, n, ii, jj, ncond

#ifdef UWM_STATS
character*20 name
integer m
character*10 frac_in_char_temp
real curr_frac
#else
character*8 name
#endif /*UWM_STATS*/
character*80 longname
character*10 units

microcount = 0

name = 'QTFLUX'
longname = 'Total (resolved + subgrid) total water (vapor+cloud) flux'
units = 'W/m2'
call add_to_namelist(count,microcount,name,longname,units,0)

do n = 1,nmicro_fields
!bloss/qt   if(n.ne.iqv) then
  ! add mean value of microphysical field to statistics
  !   EXCEPT for water vapor (added in statistics.f90)
  name = trim(mkname(n))
  longname = trim(mklongname(n))
  units = trim(mkunits(n))
  call add_to_namelist(count,microcount,name,longname,units,0)
  if(n.eq.iqv) then
      ! add variance of ONLY total water (vapor + cloud liq) field to statistics
      !   cloud water variance and cloud ice variance
      !   already output in statistics.f90
      name = trim(mkname(n))//'2'
      longname = 'Variance of '//trim(mklongname(n))
      units = '('//trim(mkunits(n))//')^2'
      call add_to_namelist(count,microcount,name,longname,units,0)
   end if

   ! add vertical advective tendency
   name = trim(mkname(n))//'ADV'
   longname = 'Tendency of '//trim(mklongname(n))// &
        ' due to resolved vertical advection'
   units = trim(mkunits(n))//'/day'
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add vertical diffusive tendency
   name = trim(mkname(n))//'DIFF'
   longname = 'Tendency of '//trim(mklongname(n))// &
        ' due to vertical SGS transport'
   units = trim(mkunits(n))//'/day'
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add tendency due to large-scale vertical advection
   name = trim(mkname(n))//'LSADV'
   longname = 'Tendency of '//trim(mklongname(n))// &
        ' due to large-scale vertical advection'
   units = trim(mkunits(n))//'/day'
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add tendency due to microphysical processes
   name = trim(mkname(n))//'MPHY'
   longname = 'Tendency of '//trim(mklongname(n))// &
        ' due to microphysical processes'
   units = trim(mkunits(n))//'/day'
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add vertical diffusive tendency
   name = trim(mkname(n))//'SED'
   longname = 'Tendency of '//trim(mklongname(n))//' due to sedimentation'
   units = trim(mkunits(n))//'/day'
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add storage terms of microphysical variables 
   name = trim(mkname(n))//'STO'
   longname = 'Storage term of '//trim(mklongname(n))
   units = trim(mkunits(n))//'/day'
   call add_to_namelist(count,microcount,name,longname,units,0)

   if(flag_wmass(n).gt.0) then
      ! fluxes output in W/m2 for mass mixing ratios
      units = 'W/m2'
   else
      ! fluxes output in #/m2/s for number concentrations
      units = '#/m2/s'
   end if

   ! add flux of microphysical fields to scalar
   name = trim(mkname(n))//'FLX'
   longname = 'Total flux (Resolved + Subgrid) of '//trim(mklongname(n))
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add subgrid flux of microphysical fields to scalar
   name = trim(mkname(n))//'FLXS'
   longname = 'Subgrid flux of '//trim(mklongname(n))
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add sedimentation flux of microphysical fields to scalar
   name = trim(mkname(n))//'SDFLX'
   longname = 'Sedimentation flux of '//trim(mklongname(n))
   call add_to_namelist(count,microcount,name,longname,units,0)

   if((flag_wmass(n).gt.0).and.(n.ne.iqv)) then
      ! add area fraction of microphysical field to statistics
      name = trim(mkname(n))//'FRAC'
      longname = trim(mklongname(n))//' FRACTION'
      units = '1'
      call add_to_namelist(count,microcount,name,longname,units,0)

      ! add approximate optical depth of hydrometeor fields
      name = 'TAU'//trim(mkname(n))
      longname = 'Approx optical depth of '//trim(mklongname(n))
      units = '1'
      call add_to_namelist(count,microcount,name,longname,units,0)

!bloss (Apr 09): Eliminate this output.  It is unreliable when
!            hydrometeor fractions are variable across processors
!            or in time.  You can still compute this from 
!            TAU* and Q* values in the statistics file.
!bloss      ! add approximate optical depth of hydrometeor fields
!bloss      name = 'EFFR'//trim(mkname(n))
!bloss      longname = 'Effective radius of '//trim(mklongname(n))
!bloss      units = 'microns'
!bloss      call add_to_namelist(count,microcount,name,longname,units,0)

      ! add field which can be used to recover mean effective radius.
      name = trim(mkname(n))//'OEFFR'
      longname = 'Mixing ratio of '//trim(mklongname(n)) &
           //' over effective radius, EFFR = ' &
           //trim(mkname(n))//'/'//trim(mkname(n))//'OEFFR'
      units = 'g/kg/microns'
      call add_to_namelist(count,microcount,name,longname,units,0)
   end if

end do

!bloss/qt: add output for cloud liquid water (not included explicitly in 
!  total water formulation).
call add_to_namelist(count,microcount,'QC', &
     'Cloud liquid water mass mixing ratio', 'g/kg',0)

! add approximate optical depth of cloud liquid water
name = 'TAUQC'
longname = 'Approx optical depth of cloud liquid water'
units = '1'
call add_to_namelist(count,microcount,name,longname,units,0)

! add field which can be used to recover mean effective radius.
name = 'QCOEFFR'
longname = 'Mixing ratio of QC'// &
     ' over effective radius, EFFR = QC/QCOEFFR'
units = 'g/kg/microns'
call add_to_namelist(count,microcount,name,longname,units,0)

!bloss/qt: Can't be computed reliably in total water formulation.
! add temperature tendency (sensible energy) tendency due to mphys
!bloss call add_to_namelist(count,microcount,'QLAT', &
!bloss      'Sensible energy tendency due to phase changes', 'K/day',0)
 
do ncond = 1,ncondavg
   ! add conditional averages of hydrometeor fields
!bloss/qt: Can't be computed reliably in total water formulation.
!bloss   call add_to_namelist(count,microcount,'QLAT' // TRIM(condavgname(ncond)), &
!bloss        'Sensible energy tendency due to phase changes in ' // TRIM(condavglongname(ncond)), &
!bloss        'K/day',ncond)
   !bloss/qt: add conditional averages for water vapor and cloud liquid water
   call add_to_namelist(count,microcount,'QV' // TRIM(condavgname(ncond)), &
        'Water vapor mixing ratio in ' // TRIM(condavglongname(ncond)),'kg/kg',ncond)
   call add_to_namelist(count,microcount,'QC' // TRIM(condavgname(ncond)), &
        'Cloud liquid water mixing ratio in ' // TRIM(condavglongname(ncond)),'kg/kg',ncond)
   do n = 1,nmicro_fields
      call add_to_namelist(count,microcount,trim(mkname(n)) // TRIM(condavgname(ncond)), &
           trim(mklongname(n)) // ' in ' // TRIM(condavglongname(ncond)), &
           trim(mkunits(n)),ncond)
   end do
end do

! add microphysical process rates into model statstics +++mhwang
#include "mhwang_stats_hbuf_init.inc"

#ifndef UWM_STATS
if(masterproc) then
   write(*,*) 'Added ', microcount, ' arrays to statistics for M2005 microphysics'
end if
#endif /*UWM_STATS*/

#ifdef UWM_STATS
#include "uwm_stats_hbuf_init.inc"
#endif /*UWM_STATS*/

end subroutine micro_hbuf_init

!----------------------------------------------------------------------
!!!! Collect microphysics history statistics (vertical profiles)
!! Note that only the fields declared in micro_hbuf_init() are allowed to
! be collected

subroutine micro_statistics()

use vars
#ifdef UWM_STATS
use hbuffer, only: hbuf_put, hbuf_put_level
#else
use hbuffer, only: hbuf_put
#endif /*UWM_STATS*/
use params, only : lcond

#ifdef UWM_STATS
use compute_correlation_module, only: corravg_count, sum_value_xy_domain,& 
                                      mean_xy_domain, variance_xy_domain,&
                                      mean_ip_xy_domain, variance_ip_xy_domain,&
                                      covariance_xy_domain, covariance_ip_xy_domain,&
                                      compute_correlation, undef_corr  

implicit none

!Microphysics correlations and covariances matricies
real, dimension(11,11,nzm) :: micro_correlations ! Magic numbers are the number
                                                 ! of variables in the correlation
!                                                  and covariance matricies
real, dimension(11,11,nzm) :: micro_covarnce 

!Microphysics correlations and covariances matricies in-cloud
real, dimension(11,11,nzm) :: ic_micro_correlations 
real, dimension(11,11,nzm) :: ic_micro_covarnce 

!Microphysics correlations and covariances matricies out of cloud
real, dimension(11,11,nzm) :: oc_micro_correlations 
real, dimension(11,11,nzm) :: oc_micro_covarnce 

! Mean and variance of chi(s_mellor)
real, dimension(nzm) :: domain_mean_chi 
real, dimension(nzm) :: domain_varnce_chi 

! Mean and variance of chi(s_mellor) in-cloud
real, dimension(nzm) :: ic_mean_chi 
real, dimension(nzm) :: ic_varnce_chi 

! Mean and variance of chi(s_mellor) out of cloud
real, dimension(nzm) :: oc_mean_chi 
real, dimension(nzm) :: oc_varnce_chi 

! Vertical velocity, interpolated, domain means, and variances
real, dimension(nx,ny,nzm) :: w_zt !vertical velocity interpolated on the scalar grid
real, dimension(nzm) :: domain_mean_w_zt !domain average of w_zt
real, dimension(nzm) :: domain_varnce_w_zt !domain variance of w_zt

real, dimension(nzm) :: ic_mean_w_zt !domain average of w_zt in-cloud
real, dimension(nzm) :: ic_varnce_w_zt !domain variance of w_zt in-cloud

real, dimension(nzm) :: oc_mean_w_zt !domain average of w_zt out of cloud
real, dimension(nzm) :: oc_varnce_w_zt !domain variance of w_zt out of cloud

!Microphysical domain-wide means and variances 
real, dimension(nmicro_fields,nzm) :: domain_mean_micro !domain averages of micro_fields
real, dimension(nmicro_fields,nzm) :: domain_varnce_micro !domain variance of micro_fields

!Microphysical within-'precip' means and variances 
real, dimension(nmicro_fields,nzm) :: ip_mean_micro ! within-'precip' averages of micro_fields
real, dimension(nmicro_fields,nzm) :: ip_varnce_micro ! within-'precip' of micro_fields

!Microphysical within-cloud means and variances 
real, dimension(nmicro_fields,nzm) :: ic_mean_micro ! within-cloud averages of micro_fields
real, dimension(nmicro_fields,nzm) :: ic_varnce_micro ! within-cloud of micro_fields

!Microphysical out of cloud means and variances 
real, dimension(nmicro_fields,nzm) :: oc_mean_micro ! within-cloud averages of micro_fields
real, dimension(nmicro_fields,nzm) :: oc_varnce_micro ! within-cloud of micro_fields

integer :: idx_s, idx_w, idx_Nc, idx_rr, idx_Nr, idx_ri, idx_Ni,&
           idx_rs, idx_Ns, idx_rg, idx_Ng, micro_indx_start, &
           offset ! Morrison has indices for Nc, rr, etc. This is the offset
                  ! between Morrison's and the correlation/covariance matrix indicies
   
!--------------------------------------------------------------------------------------------------
!Ensemble fractions 
!
!  These are used to see the effect of tolerance values on the in-'precip'
!  variances and fractions. The magic number 5 refers to the number of tolerance
!  values we test. Starting from SAM's defaut 1e-6 [kg/kg] for all micro.
!  species, we modulate the threshold +/- 2 orders of magnitude. 

!  Any variables that are reliant on these thresholds are defined in this
!  section

! Binary mask of in(1) and out(0) of micro. species. Add '1' for 'out of cloud'
! field
real, dimension(nx,ny,nzm,nfrac_fields,nfractions) :: micro_mask
real, dimension(nx,ny,nzm) :: out_cloud_mask

! Number of gridpoints within micro. species. Add '1' for 'out of cloud'
! field
integer, dimension(nzm,nfrac_fields,nfractions) :: micro_sum
integer, dimension(nzm) :: out_cloud_sum

! Micro. fractions. Add '1' for 'out of cloud'
! field
real, dimension(nzm, nfrac_fields, nfractions) :: micro_frac 

! Cloud, Rain, Ice, Snow, Graupel points as defined by the threshhold,
! (thresh_index [kg kg^-1])
integer :: rc_pts, rr_pts, ri_pts, rs_pts, rg_pts, thresh_index

real ::  curr_thresh ! Current threshold value 
integer :: order_of_magnitude ! 10

! Threshold value used to compute in-'precip' means/variance. Set to nfractions
integer :: thresh_out

! Name for each fraction written to disk. 
character*20 name
!--------------------------------------------------------------------------------------------------
#endif /*UWM_STATS*/
real, dimension(nzm) :: tr0, tr2

real tmp(2), factor_xy
integer i,j,k,m, n, ii, jj, nn, ncond

call t_startf ('micro_statistics')

factor_xy = 1./float(nx*ny)

do n = 1,nmicro_fields
   do k = 1,nzm
      tmp(1) = dz
      tmp(2) = dz/dtn
      tr0(k) = SUM(micro_field(1:nx,1:ny,k,n))
      tr2(k) = SUM(micro_field(1:nx,1:ny,k,n)*micro_field(1:nx,1:ny,k,n))
      mkwle(k,n) = mkwle(k,n)*tmp(2)*lfac(n) ! resolved flux
      mkwsb(k,n) = mkwsb(k,n)*tmp(1)*lfac(n) ! subgrid flux
      mksed(k,n) = mksed(k,n)*lfac(n) ! sedimentation flux
#ifdef PNNL_STATS
      ! Heng Xiao: convert vertical flux of total water into g/kg m/s unit needed for CLUBB input
      if (n .eq. iqv) then
        mkwle(k,n) = mkwle(k,n)/rhow(k)/lfac(n)*1.0e3
        mkwsb(k,n) = mkwsb(k,n)/rhow(k)/lfac(n)*1.0e3
      end if
#endif /* PNNL_STATS */

      !mstor(k, n) = SUM(micro_field(1:nx,1:ny,k,n))-mstor(k,n)
      mstor(k, n) = SUM(dble(micro_field(1:nx,1:ny,k,n)))-mstor(k,n) !MWSWong: cast to double
   end do

   if(flag_wmass(n).lt.1) then
      ! remove factor of rho from number concentrations
      tr0(:) = tr0(:)*rho(:)
      tr2(:) = tr2(:)*rho(:)**2
      mkadv(1:nzm,n) = mkadv(1:nzm,n)*rho(:)
      mkdiff(1:nzm,n) = mkdiff(1:nzm,n)*rho(:)
      mtend(1:nzm,n) = mtend(1:nzm,n)*rho(:)
      stend(1:nzm,n) = stend(1:nzm,n)*rho(:)
      mklsadv(1:nzm,n) = mklsadv(1:nzm,n)*rho(:)

      mstor(1:nzm,n) = mstor(1:nzm,n)*rho(:)
   end if

!bloss/qt: output all microphysical fields
!   if(n.ne.iqv) then
   ! mean microphysical field
   call hbuf_put(trim(mkname(n)),tr0,mkoutputscale(n)*factor_xy)
  if(n.eq.iqv) then
      ! variance of microphysical field,  only for QTO (qv+qcl)
      call hbuf_put(trim(mkname(n))//'2',tr2,mkoutputscale(n)**2*factor_xy)
   end if

   ! do not rescale fluxes
   call hbuf_put(trim(mkname(n))//'FLX',mkwle(1,n),factor_xy)
   call hbuf_put(trim(mkname(n))//'FLXS',mkwsb(1,n),factor_xy)
   call hbuf_put(trim(mkname(n))//'SDFLX',mksed(1,n),factor_xy)

   ! tendencies
   call hbuf_put(trim(mkname(n))//'ADV', &
        mkadv(:,n),mkoutputscale(n)*factor_xy*86400./dtn)
   call hbuf_put(trim(mkname(n))//'DIFF', &
        mkdiff(:,n),mkoutputscale(n)*factor_xy*86400./dtn)
   call hbuf_put(trim(mkname(n))//'LSADV', &
        mklsadv(:,n),mkoutputscale(n)*factor_xy*86400.)
   call hbuf_put(trim(mkname(n))//'MPHY', &
#ifdef UWM_STATS
!mtend already includes sedimentation
        (mtend(:,n)-stend(:,n)),mkoutputscale(n)*factor_xy*86400.)
#else
        mtend(:,n),mkoutputscale(n)*factor_xy*86400.)
#endif
   call hbuf_put(trim(mkname(n))//'SED', &
        stend(:,n),mkoutputscale(n)*factor_xy*86400.)

   ! Storage terms
   call hbuf_put(trim(mkname(n))//'STO', &
        mstor(:,n),mkoutputscale(n)*factor_xy*86400./dtn)

   if((flag_wmass(n).gt.0).and.(n.ne.iqv)) then
      ! fractional area of microphysical field > 1.e-6
      call hbuf_put(trim(mkname(n))//'FRAC',mfrac(1,n),factor_xy)

      ! approx optical depth
      call hbuf_put('TAU'//trim(mkname(n)),trtau(:,n),factor_xy)

      !bloss (Apr 09):  This measure of effective radius is unreliable if the 
      !          microphysical fraction is not roughly uniform across
      !          the processors in a MPI run.  As a result, I am
      !          removing it from the outputs.  It is reliable if computed from
      !          the quantities TAU* and Q* in the output file.
!bloss      ! effective radius
!bloss      tr2(:) = 0.
!bloss      if(trtau(1,n).gt.0.) then
!bloss         tr2(1) = 1.e6*0.0018*rho(1)*dz*adz(1)*tr0(1)/trtau(1,n)
!bloss      end if

!bloss      do k = 2,nzm
!bloss         if(trtau(k,n).gt.trtau(k-1,n)) then
!bloss            tr2(k) = 1.e6*0.0018*rho(k)*dz*adz(k)*tr0(k)/(trtau(k,n)-trtau(k-1,n))
!bloss         end if
!bloss      end do
!bloss      call hbuf_put('EFFR'//trim(mkname(n)),tr2,1.)

      !bloss (Apr 09): Make an alternate statistic that can be used
      ! to easily compute the mean effective radius in a consistent
      ! way from optical depth.  This quantity Q*OEFFR is essentially
      ! the layer optical depth scaled in such a way that
      !
      !    EFFR = <Q*> / <Q*OEFFR>
      !
      ! where <.> is a time- and horizontal average.
      tr2(:) = 0.
      tr2(1) = trtau(1,n) / (1.e6*0.0018*rho(1)*dz*adz(1)*1.e-3)
      do k = 2,nzm
            tr2(k) = (trtau(k,n)-trtau(k-1,n)) / (1.e6*0.0018*rho(k)*dz*adz(k)*1.e-3) 
      end do
      call hbuf_put(trim(mkname(n))//'OEFFR',tr2,factor_xy)
      
   end if

   do ncond = 1,ncondavg
      do k = 1,nzm
         tr0(k) = SUM(micro_field(1:nx,1:ny,k,n)*condavg_mask(1:nx,1:ny,k,ncond))
      end do
      if(flag_number(n).eq.1) tr0(:) = tr0(:)*rho(:) ! remove factor of rho from number concentrations
      call hbuf_put(TRIM(mkname(n)) // TRIM(condavgname(ncond)), &
           tr0,mkoutputscale(n))
   end do

end do

!bloss/qt: in total water formulation, fluxes of qv and qcl computed together.
tr0(:) = mkwle(1:nzm,iqv) + mkwsb(1:nzm,iqv) ! qv + qcl tendencies
if(doicemicro) then
   tr0(:) = tr0(:) + mkwle(1:nzm,iqci) + mkwsb(1:nzm,iqci)
end if
call hbuf_put('QTFLUX',tr0,factor_xy)

!bloss/qt: Can't be computed reliably in total water formulation.
!bloss do k = 1,nzm
!bloss    tr0(k) = SUM(tmtend3d(1:nx,1:ny,k))
!bloss end do
!bloss call hbuf_put('QLAT',tr0,factor_xy*86400.)

!bloss do ncond = 1,ncondavg
!bloss    do k = 1,nzm
!bloss       tr0(k) = SUM(tmtend3d(1:nx,1:ny,k)*condavg_mask(1:nx,1:ny,k,ncond))
!bloss    end do
!bloss    call hbuf_put('QLAT' // TRIM(condavgname(ncond)),tr0,86400.)
!bloss end do

!bloss/qt: add separate output for cloud liquid water
!           and approx cloud liquid optical depth.
do k = 1,nzm
  tr0(k) = SUM(cloudliq(1:nx,1:ny,k))
end do
call hbuf_put('QC',tr0,factor_xy*mkoutputscale(iqv))

call hbuf_put('TAUQC',trtau(:,iqv),factor_xy)

!bloss (Apr 09): Make an alternate statistic that can be used
! to easily compute the mean effective radius in a consistent
! way from optical depth.  This quantity Q*OEFFR is essentially
! the layer optical depth scaled in such a way that
!
!    EFFR = <Q*> / <Q*OEFFR>
!
! where <.> is a time- and horizontal average.
tr2(:) = 0.
tr2(1) = trtau(1,iqv) / (1.e6*0.0018*rho(1)*dz*adz(1)*1.e-3)
do k = 2,nzm
  tr2(k) = (trtau(k,iqv)-trtau(k-1,iqv)) / (1.e6*0.0018*rho(k)*dz*adz(k)*1.e-3) 
end do
call hbuf_put('QCOEFFR',tr2,factor_xy)

!bloss/qt: add separate conditional averages for cloud liquid water and vapor.
do ncond = 1,ncondavg
   do k = 1,nzm
      tr0(k) = SUM(cloudliq(1:nx,1:ny,k)*condavg_mask(1:nx,1:ny,k,ncond))
   end do
   call hbuf_put('QC' // TRIM(condavgname(ncond)),tr0,mkoutputscale(iqv))
   do k = 1,nzm
      tr0(k) = SUM((micro_field(1:nx,1:ny,k,iqv)-cloudliq(1:nx,1:ny,k))*condavg_mask(1:nx,1:ny,k,ncond))
   end do
   call hbuf_put('QV' // TRIM(condavgname(ncond)),tr0,mkoutputscale(iqv))
end do

if(dopredictNc) then
  nn = 0.
  tmp(1)=0.
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      if(qcl(i,j,k).gt.0.) then
         tmp(1) = tmp(1) + micro_field(i,j,k,incl)*1.e-6
         nn = nn + 1
       end if
    end do
   end do      
  end do
  if (nn.gt.0) ncmn = ncmn + tmp(1)/dble(nn)
else
  ncmn = Nc0
end if
if(doprecip) then
  nn = 0.
  tmp(1)=0.
  do k=1,nzm
   do j=1,ny
    do i=1,nx 
      if(micro_field(i,j,k,iqr).gt.0.) then 
         tmp(1) = tmp(1) + micro_field(i,j,k,inr)*1.e-6
         nn = nn + 1
       end if
    end do
   end do
  end do
  if (nn.gt.0) then
      nrainy = nrainy + 1
      nrmn = nrmn + tmp(1)/dble(nn)
  end if
else
  nrmn = 0.
end if

#ifdef UWM_STATS
#include "uwm_stats_hbuf_put_prep.inc"
#endif /* UWM_STATS */

! Microphysical process rates +++mhwnag
#include "mhwang_stats_hbuf_put.inc"

#ifdef UWM_STATS
#include "uwm_stats_hbuf_put.inc"
#endif /*UWM_STATS*/

call t_stopf ('micro_statistics')

end subroutine micro_statistics

!-----------------------------------------
subroutine satadj_liquid(nzm,tabs,qt,qc,pres)
  !bloss/qt: Utility routine based on cloud.f90 in 
  !  MICRO_SAM1MOM that was written by Marat Khairoutdinov.
  !  This routine performs a saturation adjustment for
  !  cloud liquid water only using a Newton method.
  !  While 20 iterations are allowed, most often this
  !  routine should exit in five iterations or less.
  !  Only a single calculation of the saturation vapor
  !  pressure is required in subsaturated air.

  use module_mp_GRAUPEL, only: polysvp
  use params, only: cp, lcond, rv, fac_cond
  implicit none

  integer, intent(in) :: nzm
  real, intent(inout), dimension(nzm) :: tabs ! absolute temperature, K
  real, intent(inout), dimension(nzm) :: qt  ! on input: qt; on output: qv
  real, intent(out), dimension(nzm) :: qc ! cloud liquid water, kg/kg
  real, intent(in), dimension(nzm) :: pres ! pressure, Pa

  real tabs1, dtabs, thresh, esat1, qsat1, fff, dfff
  integer k, niter

  integer, parameter :: maxiter = 20

  !bloss/qt: quick saturation adjustment to compute cloud liquid water content.
  do k = 1,nzm
    tabs1 = tabs(k) 
    esat1 = polysvp(tabs1,0)
    qsat1 = 0.622*esat1/ (pres(k) - esat1)
    qc(k) = 0. ! no cloud unless qt > qsat
    
    if (qt(k).gt.qsat1) then

      ! if unsaturated, nothing to do (i.e., qv=qt, T=Tl) --> just exit.
      ! if saturated, do saturation adjustment 
      !    (modeled after Marat's cloud.f90).

      ! generate initial guess based on above calculation of qsat
      dtabs = + fac_cond*MAX(0.,qt(k) - qsat1) &
           / ( 1. + lcond**2*qsat1/(cp*rv*tabs1**2) )
      tabs1 = tabs1 + dtabs
      niter = 1

      ! convergence threshold: min of 0.01K and latent heating due to
      !    condensation of 1% of saturation mixing ratio.
      thresh = MIN(0.01, 0.01*fac_cond*qsat1)

      ! iterate while temperature increment > thresh and niter < maxiter
      do while((ABS(dtabs).GT.thresh) .AND. (niter.lt.maxiter))

        esat1 = polysvp(tabs1,0)
        qsat1 = 0.622*esat1/ (pres(k) - esat1) ! saturation mixing ratio

        fff = tabs(k) - tabs1 + fac_cond*MAX(0.,qt(k) - qsat1)
        dfff = 1. + lcond**2*qsat1/(cp*rv*tabs1**2)
        dtabs = fff/dfff
        tabs1 = tabs1 + dtabs

        niter = niter + 1

      end do

      qc(k) = MAX( 0.,tabs1 - tabs(k) )/fac_cond ! cloud liquid mass mixing ratio
      qt(k) = qt(k) - qc(k) ! This now holds the water vapor mass mixing ratio.
      tabs(k) = tabs1 ! update temperature.
      
      if(niter.gt.maxiter-1) write(*,*) 'Reached iteration limit in satadj_liquid'

    end if ! qt_in > qsat

  end do ! k = 1,nzm

end subroutine satadj_liquid

!-----------------------------------------------------------------------
! Supply function that computes total water in a domain:
!
real(8) function total_water()

  use vars, only : nstep,nprint,adz,dz,rho
  real(8) tmp
  integer i,j,k,m

  total_water = 0.
  do m=1,nmicro_fields
   if(flag_wmass(m).eq.1) then
    do k=1,nzm
      tmp = 0.
      do j=1,ny
        do i=1,nx
          tmp = tmp + micro_field(i,j,k,m)
        end do
      end do
      total_water = total_water + tmp*adz(k)*dz*rho(k)
    end do
   end if
  end do

end function total_water

function Get_reffc() ! liquid water
  real, dimension(nx,ny,nzm) :: Get_reffc
  Get_reffc = reffc
end function Get_reffc

function Get_reffi() ! ice
  real, dimension(nx,ny,nzm) :: Get_reffi
  Get_reffi = reffi
end function Get_reffi

#if defined(CLUBB) || defined(UWM_STATS)

!-------------------------------------------------------------------------------
ELEMENTAL FUNCTION LIN_INT( var_high, var_low, height_high, height_low, height_int )

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
! Modifications: Dave Schanen added the elemental attribute 4 Nov 2008
! References: None

IMPLICIT NONE

! Input Variables
REAL, INTENT(IN):: var_high
REAL, INTENT(IN):: var_low
REAL, INTENT(IN):: height_high
REAL, INTENT(IN):: height_low
REAL, INTENT(IN):: height_int

! Output Variable
REAL:: LIN_INT

LIN_INT = ( var_high - var_low ) / ( height_high - height_low ) &
         * ( height_int - height_low ) + var_low


END FUNCTION LIN_INT
#endif

#if defined(UWM_STATS) || defined(LASSO_ENA_3D_MICRO_OUTPUT)
#include "uwm_stats_3d_out.inc"
#endif /*UWM_STATS*/

end module microphysics



