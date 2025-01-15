integer ::                                                                                         &
  idx_cor_chi_w, idx_cor_chi_Nc, idx_cor_chi_rr, idx_cor_chi_Nr, idx_cor_chi_ri, idx_cor_chi_Ni,   &
  idx_cor_chi_rs, idx_cor_chi_Ns, idx_cor_chi_rg, idx_cor_chi_Ng,                                  &
  idx_cor_w_Nc, idx_cor_w_rr, idx_cor_w_Nr, idx_cor_w_ri, idx_cor_w_Ni, idx_cor_w_rs, idx_cor_w_Ns,&
  idx_cor_w_rg, idx_cor_w_Ng, idx_cor_Nc_rr, idx_cor_Nc_Nr, idx_cor_Nc_ri, idx_cor_Nc_Ni,          &
  idx_cor_Nc_rs, idx_cor_Nc_Ns, idx_cor_Nc_rg, idx_cor_Nc_Ng, idx_cor_rr_Nr, idx_cor_rr_ri,        &
  idx_cor_rr_Ni, idx_cor_rr_rs, idx_cor_rr_Ns, idx_cor_rr_rg, idx_cor_rr_Ng, idx_cor_Nr_ri,        &
  idx_cor_Nr_Ni, idx_cor_Nr_rs, idx_cor_Nr_Ns, idx_cor_Nr_rg, idx_cor_Nr_Ng, idx_cor_ri_Ni,        &
  idx_cor_ri_rs, idx_cor_ri_Ns, idx_cor_ri_rg, idx_cor_ri_Ng, idx_cor_Ni_rs, idx_cor_Ni_Ns,        &
  idx_cor_Ni_rg, idx_cor_Ni_Ng, idx_cor_rs_Ns, idx_cor_rs_rg, idx_cor_rs_Ng, idx_cor_Ns_rg,        &
  idx_cor_Ns_Ng, idx_cor_rg_Ng,                                                                    &
  
  idx_ic_cor_chi_w, idx_ic_cor_chi_Nc, idx_ic_cor_chi_rr, idx_ic_cor_chi_Nr, idx_ic_cor_chi_ri,    &
  idx_ic_cor_chi_Ni, idx_ic_cor_chi_rs, idx_ic_cor_chi_Ns, idx_ic_cor_chi_rg, idx_ic_cor_chi_Ng,   &
  idx_ic_cor_w_Nc, idx_ic_cor_w_rr, idx_ic_cor_w_Nr, idx_ic_cor_w_ri, idx_ic_cor_w_Ni,             &
  idx_ic_cor_w_rs, idx_ic_cor_w_Ns, idx_ic_cor_w_rg, idx_ic_cor_w_Ng, idx_ic_cor_Nc_rr,            &
  idx_ic_cor_Nc_Nr, idx_ic_cor_Nc_ri, idx_ic_cor_Nc_Ni, idx_ic_cor_Nc_rs, idx_ic_cor_Nc_Ns,        &
  idx_ic_cor_Nc_rg, idx_ic_cor_Nc_Ng, idx_ic_cor_rr_Nr, idx_ic_cor_rr_ri, idx_ic_cor_rr_Ni,        &
  idx_ic_cor_rr_rs, idx_ic_cor_rr_Ns, idx_ic_cor_rr_rg, idx_ic_cor_rr_Ng, idx_ic_cor_Nr_ri,        &
  idx_ic_cor_Nr_Ni, idx_ic_cor_Nr_rs, idx_ic_cor_Nr_Ns, idx_ic_cor_Nr_rg, idx_ic_cor_Nr_Ng,        &
  idx_ic_cor_ri_Ni, idx_ic_cor_ri_rs, idx_ic_cor_ri_Ns, idx_ic_cor_ri_rg, idx_ic_cor_ri_Ng,        &
  idx_ic_cor_Ni_rs, idx_ic_cor_Ni_Ns, idx_ic_cor_Ni_rg, idx_ic_cor_Ni_Ng, idx_ic_cor_rs_Ns,        &
  idx_ic_cor_rs_rg, idx_ic_cor_rs_Ng, idx_ic_cor_Ns_rg, idx_ic_cor_Ns_Ng, idx_ic_cor_rg_Ng,        &

  idx_oc_cor_chi_w, idx_oc_cor_chi_Nc, idx_oc_cor_chi_rr, idx_oc_cor_chi_Nr, idx_oc_cor_chi_ri,    &
  idx_oc_cor_chi_Ni, idx_oc_cor_chi_rs, idx_oc_cor_chi_Ns, idx_oc_cor_chi_rg, idx_oc_cor_chi_Ng,   &
  idx_oc_cor_w_Nc, idx_oc_cor_w_rr, idx_oc_cor_w_Nr, idx_oc_cor_w_ri, idx_oc_cor_w_Ni,             &
  idx_oc_cor_w_rs, idx_oc_cor_w_Ns, idx_oc_cor_w_rg, idx_oc_cor_w_Ng, idx_oc_cor_Nc_rr,            &
  idx_oc_cor_Nc_Nr, idx_oc_cor_Nc_ri, idx_oc_cor_Nc_Ni, idx_oc_cor_Nc_rs, idx_oc_cor_Nc_Ns,        &
  idx_oc_cor_Nc_rg, idx_oc_cor_Nc_Ng, idx_oc_cor_rr_Nr, idx_oc_cor_rr_ri, idx_oc_cor_rr_Ni,        &
  idx_oc_cor_rr_rs, idx_oc_cor_rr_Ns, idx_oc_cor_rr_rg, idx_oc_cor_rr_Ng, idx_oc_cor_Nr_ri,        &
  idx_oc_cor_Nr_Ni, idx_oc_cor_Nr_rs, idx_oc_cor_Nr_Ns, idx_oc_cor_Nr_rg, idx_oc_cor_Nr_Ng,        &
  idx_oc_cor_ri_Ni, idx_oc_cor_ri_rs, idx_oc_cor_ri_Ns, idx_oc_cor_ri_rg, idx_oc_cor_ri_Ng,        &
  idx_oc_cor_Ni_rs, idx_oc_cor_Ni_Ns, idx_oc_cor_Ni_rg, idx_oc_cor_Ni_Ng, idx_oc_cor_rs_Ns,        &
  idx_oc_cor_rs_rg, idx_oc_cor_rs_Ng, idx_oc_cor_Ns_rg, idx_oc_cor_Ns_Ng, idx_oc_cor_rg_Ng,        &
  
  !Fractions
  idx_cld_frac, idx_rain_frac, idx_cldic_frac, idx_snw_frac, idx_grpl_frac 


! Default values for the ensemble fractions. 
logical :: doensemble_fractions = .false.  
real :: frac_threshold_init = 1e-4
integer :: nfractions = 1
integer ::  nfrac_fields = 1

character*10, allocatable, dimension(:) :: frac_name
character*80, allocatable, dimension(:) :: frac_longname
character*10, allocatable, dimension(:) :: frac_in_char