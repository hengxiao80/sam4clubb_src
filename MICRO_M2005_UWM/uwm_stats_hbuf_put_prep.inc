! In micro arrays, indicies for hydrometeor
 rc_pts = 1 
 rr_pts = 2
 ri_pts = 3
 rs_pts = 4
 rg_pts = 5
 
 ! Compute in-precip means/variance with the final threshold
 thresh_out = nfractions
 order_of_magnitude = 10

! Current threshold value
curr_thresh = frac_threshold_init

do thresh_index=1,nfractions

  do k=1,nzm
    do j=1,ny
      do i=1,nx
        ! We only want cloud water, not cloud water + rain water.
        if ( cloudliq(i,j,k)  > curr_thresh) then
          micro_mask(i,j,k,rc_pts,thresh_index) = 1
          out_cloud_mask(i,j,k) = 0
        else
          micro_mask(i,j,k,rc_pts,thresh_index) = 0
          out_cloud_mask(i,j,k) = 1
        end if 
        
        if (doprecip) then
          if (micro_field(i,j,k,iqr) > curr_thresh) then
            micro_mask(i,j,k,rr_pts,thresh_index) = 1
          else
            micro_mask(i,j,k,rr_pts,thresh_index) = 0
          end if
          
          if (doicemicro) then
            if (micro_field(i,j,k,iqci) > curr_thresh) then
              micro_mask(i,j,k,ri_pts,thresh_index) = 1
            else
              micro_mask(i,j,k,ri_pts,thresh_index) = 0
            end if

            if (micro_field(i,j,k,iqs) > curr_thresh) then
              micro_mask(i,j,k,rs_pts,thresh_index) = 1
            else
              micro_mask(i,j,k,rs_pts,thresh_index) = 0
            end if
          
            if (dograupel) then
              if (micro_field(i,j,k,iqg) > curr_thresh) then
                micro_mask(i,j,k,rg_pts,thresh_index) = 1
              else
                micro_mask(i,j,k,rg_pts,thresh_index) = 0
              end if
            
            end if!do graupel  
          end if!do icemicro
        end if!doprecip
      
      end do !i
    end do !j
  end do !k
 
! Increase the threshold by a factor of 10
curr_thresh = curr_thresh * order_of_magnitude
end do !thresh_index

! Compute the fractions by summing and finding the mean of the binary mask above.
do n=1,nfrac_fields
  do m=1,nfractions
    micro_sum(:,n,m) = sum_value_xy_domain(nzm,micro_mask(1:nx,1:ny,1:nzm,n,m))
    micro_frac(:,n,m) = mean_xy_domain(nzm,micro_mask(1:nx,1:ny,1:nzm,n,m))
  end do
end do
! Take care of 'out of' cloud points
out_cloud_sum(:) = sum_value_xy_domain( nzm,out_cloud_mask(1:nx,1:ny,1:nzm) )

!Below is used to compute covariances and correlations of microphysical quantities
!=================================================================================

! These are the indicies in the correlation/covariance matrix for each of the
! respective variables.  
idx_s  = 1  ! chi
idx_w  = 2  ! Vertical velocity
if(dopredictNc) then
  idx_Nc = 3  
  idx_rr = 4  
  idx_Nr = 5  
  idx_ri = 6 
  idx_Ni = 7 
  idx_rs = 8  
  idx_Ns = 9  
  idx_rg = 10 
  idx_Ng = 11 
  ! When computing means, variances, covariances, and correlations, we loop
  ! through Morrison's indicies (iqr) and write to our matricies (idx_rr)
  ! These two are related by the 'offset', weberjk(UWM)
  offset = idx_Nc - incl
  micro_indx_start = incl
else
  idx_rr = 3  
  idx_Nr = 4  
  idx_ri = 5 
  idx_Ni = 6 
  idx_rs = 7  
  idx_Ns = 8  
  idx_rg = 9 
  idx_Ng = 10 
  offset = idx_rr - iqr
  micro_indx_start = iqr
endif
  !Interpolate w to the scalar grid
  do k=1,nz-1
    do i=1,nx
      do j = 1,ny
        w_zt(i,j,k) =  LIN_INT( w(i,j,k+1), w(i,j,k), zi(k+1), zi(k), z(k) ) 
      end do
    end do
  end do

  !Find the domain wide means and variances of chi and w_zt
  domain_mean_chi(:) = mean_xy_domain(nzm, chi)
  domain_varnce_chi(:) = variance_xy_domain(nzm, chi, domain_mean_chi)
  
  domain_mean_w_zt(:) = mean_xy_domain(nzm, w_zt)
  domain_varnce_w_zt(:) = variance_xy_domain(nzm, w_zt, domain_mean_w_zt)

  !Find the in-cloud means and variances of chi and w_zt
  ic_mean_chi(:) = mean_ip_xy_domain( nzm, chi, micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                   micro_sum(1:nzm,rc_pts,thresh_out) )
  ic_varnce_chi(:) = variance_ip_xy_domain( nzm, chi, micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                            ic_mean_chi(:), micro_sum(1:nzm,rc_pts,thresh_out) )

  ic_mean_w_zt(:) = mean_ip_xy_domain( nzm, w_zt, micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                   micro_sum(1:nzm,rc_pts,thresh_out) )
  ic_varnce_w_zt(:) = variance_ip_xy_domain( nzm, w_zt, micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                            ic_mean_w_zt(:), micro_sum(1:nzm,rc_pts,thresh_out) )

  !Find the out of cloud means and variances of chi and w_zt
  oc_mean_chi(:) = mean_ip_xy_domain( nzm, chi, out_cloud_mask(1:nx,1:ny,1:nzm),&
                                   out_cloud_sum(1:nzm) )
  oc_varnce_chi(:) = variance_ip_xy_domain( nzm, chi, out_cloud_mask(1:nx,1:ny,1:nzm),&
                                            oc_mean_chi(:), out_cloud_sum(1:nzm) )

  oc_mean_w_zt(:) = mean_ip_xy_domain( nzm, w_zt, out_cloud_mask(1:nx,1:ny,1:nzm),&
                                   out_cloud_sum(1:nzm) )
  oc_varnce_w_zt(:) = variance_ip_xy_domain( nzm, w_zt, out_cloud_mask(1:nx,1:ny,1:nzm),&
                                            oc_mean_w_zt(:), out_cloud_sum(1:nzm) )
  
  do n=micro_indx_start,nmicro_fields
    ! Domain mean and variances stored according to Morrison indicies
    domain_mean_micro(n,:) = mean_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,n))
    domain_varnce_micro(n,:) = variance_xy_domain(nzm, micro_field(1:nx,1:ny,1:nzm,n),&
                                                domain_mean_micro(n,:))

    ! in-cloud means and variances
    ic_mean_micro(n,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,n),&
                                          micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),micro_sum(1:nzm,rc_pts,thresh_out))

    ic_varnce_micro(n,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,n),&
                                          micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),ic_mean_micro(n,:),&
                                          micro_sum(1:nzm,rc_pts,thresh_out))

    ! out of cloud means and variances
    oc_mean_micro(n,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,n),&
                                          out_cloud_mask(1:nx,1:ny,1:nzm),out_cloud_sum(1:nzm) )

    oc_varnce_micro(n,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,n),&
                                          out_cloud_mask(1:nx,1:ny,1:nzm),oc_mean_micro(n,:),&
                                          out_cloud_sum(1:nzm) )
  end do
  !----------------
  ! Within precip means
  !----------------
if(dopredictNc) then
  ip_mean_micro(incl,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,incl),&
                                          micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),micro_sum(1:nzm,rc_pts,thresh_out))
endif

  if(doprecip) then
     
    ip_mean_micro(iqr,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqr),&
                                          micro_mask(1:nx,1:ny,1:nzm,rr_pts,thresh_out),micro_sum(1:nzm,rr_pts,thresh_out))
  
    ip_mean_micro(inr,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,inr),&
                                          micro_mask(1:nx,1:ny,1:nzm,rr_pts,thresh_out),micro_sum(1:nzm,rr_pts,thresh_out))
    if(doicemicro) then  
      
      ip_mean_micro(iqci,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqci),&
                                          micro_mask(1:nx,1:ny,1:nzm,ri_pts,thresh_out),micro_sum(1:nzm,ri_pts,thresh_out))
  
      ip_mean_micro(inci,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,inci),&
                                          micro_mask(1:nx,1:ny,1:nzm,ri_pts,thresh_out),micro_sum(1:nzm,ri_pts,thresh_out))
      
      ip_mean_micro(iqs,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqs),&
                                          micro_mask(1:nx,1:ny,1:nzm,rs_pts,thresh_out),micro_sum(1:nzm,rs_pts,thresh_out))
  
      ip_mean_micro(ins,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,ins),&
                                          micro_mask(1:nx,1:ny,1:nzm,rs_pts,thresh_out),micro_sum(1:nzm,rs_pts,thresh_out))
      if(dograupel) then
        
        ip_mean_micro(iqg,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqg),&
                                          micro_mask(1:nx,1:ny,1:nzm,rg_pts,thresh_out),micro_sum(1:nzm,rg_pts,thresh_out))
  
        ip_mean_micro(ing,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,ing),&
                                          micro_mask(1:nx,1:ny,1:nzm,rg_pts,thresh_out),micro_sum(1:nzm,rg_pts,thresh_out))
      endif !dograupel
    endif !doice
  endif !doprecip
  
  !----------------
  ! Within precip variances
  !----------------
if(dopredictNc) then
  ip_varnce_micro(incl,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,incl),&
                                          micro_mask(1:nx,1:ny,1:nzm,1,thresh_out),ip_mean_micro(incl,:),&
                                          micro_sum(1:nzm,rc_pts,thresh_out))
endif

  if(doprecip) then
    ! Rain 
    ip_varnce_micro(iqr,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqr),&
                                          micro_mask(1:nx,1:ny,1:nzm,2,thresh_out),ip_mean_micro(iqr,:),&
                                          micro_sum(1:nzm,rr_pts,thresh_out))
  
    ip_varnce_micro(inr,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,inr),&
                                          micro_mask(1:nx,1:ny,1:nzm,2,thresh_out),ip_mean_micro(inr,:),&
                                          micro_sum(1:nzm,rr_pts,thresh_out))
    if(doicemicro) then  
      ! Cloud ice
      ip_varnce_micro(iqci,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqci),&
                                          micro_mask(1:nx,1:ny,1:nzm,3,thresh_out),ip_mean_micro(iqci,:),&
                                          micro_sum(1:nzm,ri_pts,thresh_out))
  
      ip_varnce_micro(inci,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,inci),&
                                          micro_mask(1:nx,1:ny,1:nzm,3,thresh_out),ip_mean_micro(inci,:),&
                                          micro_sum(1:nzm,ri_pts,thresh_out))
      ! Snow 
      ip_varnce_micro(iqs,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqs),&
                                          micro_mask(1:nx,1:ny,1:nzm,rs_pts,thresh_out),ip_mean_micro(iqs,:),&
                                          micro_sum(1:nzm,rs_pts,thresh_out))
  
      ip_varnce_micro(ins,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,ins),&
                                          micro_mask(1:nx,1:ny,1:nzm,rs_pts,thresh_out),ip_mean_micro(ins,:),&
                                          micro_sum(1:nzm,rs_pts,thresh_out))
      if(dograupel) then
      ! Graupel 
        ip_varnce_micro(iqg,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqg),&
                                          micro_mask(1:nx,1:ny,1:nzm,rg_pts,thresh_out),ip_mean_micro(iqg,:),&
                                          micro_sum(1:nzm,rg_pts,thresh_out))
  
        ip_varnce_micro(ing,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,ing),&
                                          micro_mask(1:nx,1:ny,1:nzm,rg_pts,thresh_out),ip_mean_micro(ing,:),&
                                          micro_sum(1:nzm,rg_pts,thresh_out))
      endif !dograupel
    endif !doice
  endif !doprecip
 
!-----------------------------------------
! Domain-wide covariances and correlations
!-----------------------------------------

  ! Compute covariances and correlations, first for chi(s_mellor). Manually for
  ! vertical velocity. 
  micro_covarnce(idx_s,idx_w,:) = covariance_xy_domain(nzm,&
                                    chi(1:nx,1:ny,1:nzm), w_zt(1:nx,1:ny,1:nzm),&
                                    domain_mean_chi(:), domain_mean_w_zt(:) )

  micro_correlations(idx_s,idx_w,:) = compute_correlation(nzm, domain_varnce_chi(:),&
                                        domain_varnce_w_zt(:), micro_covarnce(idx_s,idx_w,:) )
  

  ! Compute covariances and correlations, first for chi(s_mellor). Now, loop through
  ! all hydrometeor fields using Morrison indicies, but write in our
  ! covarnce/corr matrix by offset.
  do n = micro_indx_start,nmicro_fields 
      micro_covarnce(idx_s,n+offset,:) = covariance_xy_domain(nzm, chi(1:nx,1:ny,1:nzm), &
                                    micro_field(1:nx,1:ny,1:nzm,n), domain_mean_chi(:), &
                                    domain_mean_micro(n,:) )

      micro_correlations(idx_s,n+offset,:) = compute_correlation(nzm, domain_varnce_chi(:),&
                                  domain_varnce_micro(n,:), micro_covarnce(idx_s,n+offset,:) )
!                                                                                ^ in covarnce matrix
!                                                                                Morrison index is offset
  end do
  
  ! Compute covariances and correlations, now for vertical velocity and all
  ! hydrometeor variables.
  do n = micro_indx_start,nmicro_fields 
      micro_covarnce(idx_w,n+offset,:) = covariance_xy_domain(nzm, w_zt, micro_field(1:nx,1:ny,1:nzm,n),&
                              domain_mean_w_zt, domain_mean_micro(n,:) )

      micro_correlations(idx_w,n+offset,:) = compute_correlation(nzm, domain_varnce_w_zt,&
                                  domain_varnce_micro(n,:), micro_covarnce(idx_w,n+offset,:) )
  end do

  ! Continue to fill the covariance and correlation matricies  
  do m = micro_indx_start, nmicro_fields 
    do n = m, nmicro_fields
  
      micro_covarnce(m+offset,n+offset,:) = covariance_xy_domain(nzm, micro_field(1:nx,1:ny,1:nzm,m), &
          micro_field(1:nx,1:ny,1:nzm,n),domain_mean_micro(m,:), domain_mean_micro(n,:) )

      micro_correlations(m+offset,n+offset,:) = compute_correlation(nzm, domain_varnce_micro(m,:),&
                                  domain_varnce_micro(n,:), micro_covarnce(m+offset,n+offset,:) )
    end do
  end do

!-----------------------------------------
! In-cloud covariances and correlations
!-----------------------------------------

  ic_micro_covarnce(idx_s,idx_w,:) = covariance_ip_xy_domain(nzm,&
                                    chi(1:nx,1:ny,1:nzm), w_zt(1:nx,1:ny,1:nzm),&
                                    micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                    ic_mean_chi(:), ic_mean_w_zt(:),&
                                    micro_sum(:,rc_pts,thresh_out) )

  ic_micro_correlations(idx_s,idx_w,:) = compute_correlation(nzm, ic_varnce_chi(:),&
                                        ic_varnce_w_zt(:), ic_micro_covarnce(idx_s,idx_w,:) )
  

  do n = micro_indx_start,nmicro_fields 
      ic_micro_covarnce(idx_s,n+offset,:) = covariance_ip_xy_domain(nzm,&
                                             chi(1:nx,1:ny,1:nzm), micro_field(1:nx,1:ny,1:nzm,n),&
                                             micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                             ic_mean_chi(:), ic_mean_micro(n,:),&
                                             micro_sum(:,rc_pts,thresh_out) )

      ic_micro_correlations(idx_s,n+offset,:) = compute_correlation(nzm, ic_varnce_chi(:),&
                                        ic_varnce_micro(n,:), ic_micro_covarnce(idx_s,n+offset,:) )
  end do
  
  ! Compute covariances and correlations, now for vertical velocity and all
  ! hydrometeor variables.
  do n = micro_indx_start,nmicro_fields 
      ic_micro_covarnce(idx_w,n+offset,:) = covariance_ip_xy_domain(nzm,&
                                            w_zt, micro_field(1:nx,1:ny,1:nzm,n),&
                                            micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                            ic_mean_w_zt, ic_mean_micro(n,:), &
                                            micro_sum(:,rc_pts,thresh_out) )

      ic_micro_correlations(idx_w,n+offset,:) = compute_correlation(nzm, ic_varnce_w_zt,&
                                  ic_varnce_micro(n,:), ic_micro_covarnce(idx_w,n+offset,:) )
  end do

  ! Continue to fill the covariance and correlation matricies  
  do m = micro_indx_start, nmicro_fields 
    do n = m, nmicro_fields
  
      ic_micro_covarnce(m+offset,n+offset,:) = covariance_ip_xy_domain(nzm,&
                                 micro_field(1:nx,1:ny,1:nzm,m), micro_field(1:nx,1:ny,1:nzm,n),&
                                 micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                 ic_mean_micro(m,:), ic_mean_micro(n,:),&
                                 micro_sum(:,rc_pts,thresh_out) )

      ic_micro_correlations(m+offset,n+offset,:) = compute_correlation(nzm, ic_varnce_micro(m,:),&
                                  ic_varnce_micro(n,:), ic_micro_covarnce(m+offset,n+offset,:) )
    end do
  end do

!-----------------------------------------
! Out of cloud covariances and correlations
!-----------------------------------------

  oc_micro_covarnce(idx_s,idx_w,:) = covariance_ip_xy_domain(nzm,&
                                    chi(1:nx,1:ny,1:nzm), w_zt(1:nx,1:ny,1:nzm),&
                                    out_cloud_mask(1:nx,1:ny,1:nzm),&
                                    oc_mean_chi(:), oc_mean_w_zt(:),&
                                    out_cloud_sum(:) )

  oc_micro_correlations(idx_s,idx_w,:) = compute_correlation(nzm, oc_varnce_chi(:),&
                                        oc_varnce_w_zt(:), oc_micro_covarnce(idx_s,idx_w,:) )
  

  do n = micro_indx_start,nmicro_fields 
      oc_micro_covarnce(idx_s,n+offset,:) = covariance_ip_xy_domain(nzm,&
                                             chi(1:nx,1:ny,1:nzm), micro_field(1:nx,1:ny,1:nzm,n),&
                                             out_cloud_mask(1:nx,1:ny,1:nzm),&
                                             oc_mean_chi(:), oc_mean_micro(n,:),&
                                             out_cloud_sum(:) )

      oc_micro_correlations(idx_s,n+offset,:) = compute_correlation(nzm, oc_varnce_chi(:),&
                                        oc_varnce_micro(n,:), oc_micro_covarnce(idx_s,n+offset,:) )
  end do
  
  ! Compute covariances and correlations, now for vertical velocity and all
  ! hydrometeor variables.
  do n = micro_indx_start,nmicro_fields 
      oc_micro_covarnce(idx_w,n+offset,:) = covariance_ip_xy_domain(nzm,&
                                            w_zt, micro_field(1:nx,1:ny,1:nzm,n),&
                                            out_cloud_mask(1:nx,1:ny,1:nzm),&
                                            oc_mean_w_zt, oc_mean_micro(n,:), &
                                            out_cloud_sum(:) )

      oc_micro_correlations(idx_w,n+offset,:) = compute_correlation(nzm, oc_varnce_w_zt,&
                                  oc_varnce_micro(n,:), oc_micro_covarnce(idx_w,n+offset,:) )
  end do

  ! Continue to fill the covariance and correlation matricies  
  do m = micro_indx_start, nmicro_fields 
    do n = m, nmicro_fields
  
      oc_micro_covarnce(m+offset,n+offset,:) = covariance_ip_xy_domain(nzm,&
                                 micro_field(1:nx,1:ny,1:nzm,m), micro_field(1:nx,1:ny,1:nzm,n),&
                                 out_cloud_mask(1:nx,1:ny,1:nzm),&
                                 oc_mean_micro(m,:), oc_mean_micro(n,:),&
                                 out_cloud_sum(:) )

      oc_micro_correlations(m+offset,n+offset,:) = compute_correlation(nzm, oc_varnce_micro(m,:),&
                                  oc_varnce_micro(n,:), oc_micro_covarnce(m+offset,n+offset,:) )
    end do
  end do