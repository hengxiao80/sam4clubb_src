! weberjk(UWM), 3D QR budget
      allocate(mPRE_3D(nx,ny,nzm), mPRA_3D(nx,ny,nzm), mPRC_3D(nx,ny,nzm), mPRACS_3D(nx,ny,nzm), mMNUCCR_3D(nx,ny,nzm), &
                        mQMULTR_3D(nx,ny,nzm), mQMULTRG_3D(nx,ny,nzm), mPIACR_3D(nx,ny,nzm), mPIACRS_3D(nx,ny,nzm), &
                        mPRACG_3D(nx,ny,nzm), mPGRACS_3D(nx,ny,nzm), mPSMLT_3D(nx,ny,nzm), mPGMLT_3D(nx,ny,nzm), &
                        mQR_INST_3D(nx,ny,nzm), STAT=ierr)

! weberjk(UWM), 3D NR budget
      allocate(mSIZEFIX_NR_3D(nx,ny,nzm), mNEGFIX_NR_3D(nx,ny,nzm), mNSUBR_3D(nx,ny,nzm), mNSMLTR_3D(nx,ny,nzm), &
                        mNGMLTR_3D(nx,ny,nzm), mNPRC1_3D(nx,ny,nzm), mNPRACS_3D(nx,ny,nzm), mNNUCCR_3D(nx,ny,nzm),&
                        mNRAGG_3D(nx,ny,nzm), mNIACR_3D(nx,ny,nzm), mNIACRS_3D(nx,ny,nzm), mNPRACG_3D(nx,ny,nzm),&
                        mNGRACS_3D(nx,ny,nzm), mNR_INST_3D(nx,ny,nzm), STAT=ierr)

! MWSWong(PNNL), 3D QG budget (additional terms)
      allocate(mPSACWG_3D(nx,ny,nzm), mPGSACW_3D(nx,ny,nzm), mPRDG_3D(nx,ny,nzm), mEPRDG_3D(nx,ny,nzm), mPRACI_3D(nx,ny,nzm), &
                        mPSACR_3D(nx,ny,nzm), mEVPMG_3D(nx,ny,nzm), mQG_INST_3D(nx,ny,nzm), STAT=ierr)

! MWSWong(PNNL), 3D QS budget (additional terms)
      allocate(mPRAI_3D(nx,ny,nzm), mPSACWS_3D(nx,ny,nzm), mPRDS_3D(nx,ny,nzm), mPRCI_3D(nx,ny,nzm), mEPRDS_3D(nx,ny,nzm), &
                        mPRACIS_3D(nx,ny,nzm), mEVPMS_3D(nx,ny,nzm), mQS_INST_3D(nx,ny,nzm), STAT=ierr)

! MWSWong(PNNL), 3D QI budget (additional terms)
      allocate(mPRD_3D(nx,ny,nzm), mEPRD_3D(nx,ny,nzm), mPSACWI_3D(nx,ny,nzm), mMNUCCC_3D(nx,ny,nzm),  &
                        mQMULTS_3D(nx,ny,nzm), mQMULTG_3D(nx,ny,nzm), mMNUCCD_3D(nx,ny,nzm), mQI_INST_3D(nx,ny,nzm), STAT=ierr)

! MWSWong(PNNL), 3D NG budget (7 additional terms)
      allocate(mSIZEFIX_NG_3D(nx,ny,nzm), mNEGFIX_NG_3D(nx,ny,nzm), mNGMLTG_3D(nx,ny,nzm), mNSCNG_3D(nx,ny,nzm), &
                        mNSUBG_3D(nx,ny,nzm), mNG_INST_3D(nx,ny,nzm), STAT=ierr)

! MWSWong(PNNL), 3D NS budget (7 additional terms)
      allocate(mSIZEFIX_NS_3D(nx,ny,nzm), mNEGFIX_NS_3D(nx,ny,nzm), mNSMLTS_3D(nx,ny,nzm), mNSAGG_3D(nx,ny,nzm), &
                        mNPRCI_3D(nx,ny,nzm), mNSUBS_3D(nx,ny,nzm), mNS_INST_3D(nx,ny,nzm), STAT=ierr)

! MWSWong(PNNL), 3D NI budget (12 additional terms)
      allocate(mSIZEFIX_NI_3D(nx,ny,nzm), mNEGFIX_NI_3D(nx,ny,nzm), mNNUCCC_3D(nx,ny,nzm), &
                        mNPRAI_3D(nx,ny,nzm), mNMULTS_3D(nx,ny,nzm), mNMULTG_3D(nx,ny,nzm), mNMULTR_3D(nx,ny,nzm), &
                        mNMULTRG_3D(nx,ny,nzm), &
                        mNNUCCD_3D(nx,ny,nzm), mNSUBI_3D(nx,ny,nzm), mNIM_MORR_CL_3D(nx,ny,nzm), mNI_INST_3D(nx,ny,nzm), STAT=ierr)

! MWSWong(PNNL), 3D sedimentation terms (8 additional terms)
      allocate(mSTENDQR_3D(nx,ny,nzm), mSTENDQCI_3D(nx,ny,nzm), mSTENDQG_3D(nx,ny,nzm), mSTENDQS_3D(nx,ny,nzm), &
                        mSTENDNR_3D(nx,ny,nzm), mSTENDNCI_3D(nx,ny,nzm), mSTENDNG_3D(nx,ny,nzm), mSTENDNS_3D(nx,ny,nzm), STAT=ierr)

! weberjk(UWM), 3D NR budget
      allocate(rain_vel_3D(nx,ny,nzm), EFFR_3D(nx,ny,nzm), STAT=ierr)