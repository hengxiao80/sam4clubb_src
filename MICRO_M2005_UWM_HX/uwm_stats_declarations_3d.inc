! weberjk(UWM), 3D QR budget
real, allocatable, dimension(:,:,:) :: mPRE_3D, mPRA_3D, mPRC_3D, mPRACS_3D, mMNUCCR_3D, &
                        mQMULTR_3D, mQMULTRG_3D, mPIACR_3D, mPIACRS_3D, mPRACG_3D, mPGRACS_3D, &
                        mPSMLT_3D, mPGMLT_3D, mQR_INST_3D

! weberjk(UWM), 3D NR budget
real, allocatable, dimension(:,:,:) :: mSIZEFIX_NR_3D, mNEGFIX_NR_3D, mNSUBR_3D, mNSMLTR_3D, &
                        mNGMLTR_3D, mNPRC1_3D, mNPRACS_3D, mNNUCCR_3D, mNRAGG_3D, mNIACR_3D, &
                        mNIACRS_3D, mNPRACG_3D, mNGRACS_3D, mNR_INST_3D


!! MWSWong(PNNL), 3D QG budget (8 additional terms)
real, allocatable, dimension(:,:,:) :: mPSACWG_3D, mPGSACW_3D, mPRDG_3D, mEPRDG_3D, mPRACI_3D, &
                        mPSACR_3D, mEVPMG_3D, mQG_INST_3D

! MWSWong(PNNL), 3D QS budget (8 additional terms)
real, allocatable, dimension(:,:,:) :: mPRAI_3D, mPSACWS_3D, mPRDS_3D, mPRCI_3D, mEPRDS_3D, &
                        mPRACIS_3D, mEVPMS_3D, mQS_INST_3D

! MWSWong(PNNL), 3D QI budget (8 additional terms)
real, allocatable, dimension(:,:,:) :: mPRD_3D, mEPRD_3D, mPSACWI_3D, mMNUCCC_3D,  &
                        mQMULTS_3D, mQMULTG_3D, mMNUCCD_3D, mQI_INST_3D

! MWSWong(PNNL), 3D NG budget (6 additional terms)
real, allocatable, dimension(:,:,:) :: mSIZEFIX_NG_3D, mNEGFIX_NG_3D, mNGMLTG_3D, mNSCNG_3D, &
                        mNSUBG_3D, mNG_INST_3D

! MWSWong(PNNL), 3D NS budget (7 additional terms)
real, allocatable, dimension(:,:,:) :: mSIZEFIX_NS_3D, mNEGFIX_NS_3D, mNSMLTS_3D, mNSAGG_3D, &
                        mNPRCI_3D, mNSUBS_3D, mNS_INST_3D

! MWSWong(PNNL), 3D NI budget (12 additional terms)
real, allocatable, dimension(:,:,:) :: mSIZEFIX_NI_3D, mNEGFIX_NI_3D, mNNUCCC_3D, &
                        mNPRAI_3D, mNMULTS_3D, mNMULTG_3D, mNMULTR_3D, mNMULTRG_3D, mNNUCCD_3D, &
                        mNSUBI_3D, mNIM_MORR_CL_3D, mNI_INST_3D

! MWSWong(PNNL), 3D sedimentation term (8 additional terms)
real, allocatable, dimension(:,:,:) :: mSTENDQR_3D, mSTENDQCI_3D, mSTENDQG_3D, mSTENDQS_3D, &
                        mSTENDNR_3D, mSTENDNCI_3D, mSTENDNG_3D, mSTENDNS_3D

! weberjk(UWM), Ancilliary micro variables
real, allocatable, dimension(:,:,:) :: rain_vel_3D, EFFR_3D