! Write out each hydrometeor fraction
do n = 1,nfrac_fields
  do m = 1,nfractions
    name = trim(frac_name(n))//'_'//trim(adjustl(frac_in_char(m)))

    call hbuf_put(name,micro_frac(:,n,m),1.)

  end do
end do

!----------------------
! Domain-wide mean
!----------------------
  call hbuf_put('chim', domain_mean_chi(:), 1.)
if(dopredictNc) then
  call hbuf_put('Ncm', domain_mean_micro(incl,:), 1.)
endif
!----------------------
! Domain-wide variance
!----------------------
  call hbuf_put('chip2', domain_varnce_chi(:), 1.)

if(dopredictNc) then
  call hbuf_put('Ncp2', domain_varnce_micro(incl,:), 1.)
endif

  if (doprecip) then
    call hbuf_put('Nrp2', domain_varnce_micro(inr,:), 1.)
    call hbuf_put('rrp2', domain_varnce_micro(iqr,:), 1.)
    
    if (doicemicro) then
      call hbuf_put('Nip2', domain_varnce_micro(inci,:), 1.)
      call hbuf_put('rip2', domain_varnce_micro(iqci,:), 1.)
      
      call hbuf_put('Nsp2', domain_varnce_micro(ins,:), 1.)
      call hbuf_put('rsp2', domain_varnce_micro(iqs,:), 1.)
      
      if (dograupel) then  
        call hbuf_put('Ngp2', domain_varnce_micro(ing,:), 1.)
        call hbuf_put('rgp2', domain_varnce_micro(iqg,:), 1.)
      end if
    end if 
  endif

!----------------------
! Within-precip mean
!----------------------
if(dopredictNc) then
  call hbuf_put('Ncm_ip', ip_mean_micro(incl,:), 1.)
endif

  if (doprecip) then
    call hbuf_put('Nrm_ip', ip_mean_micro(inr,:), 1.)
    call hbuf_put('rrm_ip', ip_mean_micro(iqr,:), 1.)
    
    if (doicemicro) then
      call hbuf_put('Nim_ip', ip_mean_micro(inci,:), 1.)
      call hbuf_put('rim_ip', ip_mean_micro(iqci,:), 1.)
      
      call hbuf_put('Nsm_ip', ip_mean_micro(ins,:), 1.)
      call hbuf_put('rsm_ip', ip_mean_micro(iqs,:), 1.)
      
      if (dograupel) then  
        call hbuf_put('Ngm_ip', ip_mean_micro(ing,:), 1.)
        call hbuf_put('rgm_ip', ip_mean_micro(iqg,:), 1.)
      end if
    end if 
  endif


!----------------------
! Within-precip variance
!----------------------
if(dopredictNc) then
  call hbuf_put('Ncp2_ip', ip_varnce_micro(incl,:), 1.)
endif

  if (doprecip) then
    call hbuf_put('Nrp2_ip', ip_varnce_micro(inr,:), 1.)
    call hbuf_put('rrp2_ip', ip_varnce_micro(iqr,:), 1.)
    
    if (doicemicro) then
      call hbuf_put('Nip2_ip', ip_varnce_micro(inci,:), 1.)
      call hbuf_put('rip2_ip', ip_varnce_micro(iqci,:), 1.)
      
      call hbuf_put('Nsp2_ip', ip_varnce_micro(ins,:), 1.)
      call hbuf_put('rsp2_ip', ip_varnce_micro(iqs,:), 1.)
      
      if (dograupel) then  
        call hbuf_put('Ngp2_ip', ip_varnce_micro(ing,:), 1.)
        call hbuf_put('rgp2_ip', ip_varnce_micro(iqg,:), 1.)
      end if
    end if 
  endif

! weberjk(UWM), Microphysical correlations

   do k = 1, nzm, 1
!-------------
! chi
!-------------
      if ( micro_correlations(idx_s,idx_w,k) /= undef_corr ) then
         call hbuf_put_level('corr_chi_w',micro_correlations(idx_s,idx_w,k),1.,k)
         corravg_count(idx_cor_chi_w,k) = corravg_count(idx_cor_chi_w,k) + 1
      else ! micro_correlations(1,incl,k) == undef_corr
         call hbuf_put_level('corr_chi_w',0.0,1.,k)
      endif ! micro_correlations(1,incl,k) /= undef_corr
      
      if(dopredictNc) then
        if ( micro_correlations(idx_s,idx_Nc,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_Nc',micro_correlations(idx_s,idx_Nc,k),1.,k)
           corravg_count(idx_cor_chi_Nc,k) = corravg_count(idx_cor_chi_Nc,k) + 1
        else ! micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('corr_chi_Nc',0.0,1.,k)
        endif ! micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( micro_correlations(idx_s,idx_rr,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_rr',micro_correlations(idx_s,idx_rr,k),1.,k)
           corravg_count(idx_cor_chi_rr,k) = corravg_count(idx_cor_chi_rr,k) + 1
        else ! micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('corr_chi_rr',0.0,1.,k)
        endif ! micro_correlations(1,iqr,k) /= undef_corr
        if ( micro_correlations(idx_s,idx_Nr,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_Nr',micro_correlations(idx_s,idx_Nr,k),1.,k)
           corravg_count(idx_cor_chi_Nr,k) = corravg_count(idx_cor_chi_Nr,k) + 1
        else ! micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('corr_chi_Nr',0.0,1.,k)
        endif ! ! micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( micro_correlations(idx_s,idx_ri,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_ri',micro_correlations(idx_s,idx_ri,k),1.,k)
           corravg_count(idx_cor_chi_ri,k) = corravg_count(idx_cor_chi_ri,k) + 1
        else ! micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('corr_chi_ri',0.0,1.,k)
        endif ! micro_correlations(1,iqci,k) /= undef_corr
        if ( micro_correlations(idx_s,idx_Ni,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_Ni',micro_correlations(idx_s,idx_Ni,k),1.,k)
           corravg_count(idx_cor_chi_Ni,k) = corravg_count(idx_cor_chi_Ni,k) + 1
        else ! micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('corr_chi_Ni',0.0,1.,k)
        endif ! micro_correlations(1,inci,k) /= undef_corr
        if ( micro_correlations(idx_s,idx_rs,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_rs',micro_correlations(idx_s,idx_rs,k),1.,k)
           corravg_count(idx_cor_chi_rs,k) = corravg_count(idx_cor_chi_rs,k) + 1
        else ! micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('corr_chi_rs',0.0,1.,k)
        endif ! micro_correlations(1,iqs,k) /= undef_corr
        if ( micro_correlations(idx_s,idx_Ns,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_Ns',micro_correlations(idx_s,idx_Ns,k),1.,k)
           corravg_count(idx_cor_chi_Ns,k) = corravg_count(idx_cor_chi_Ns,k) + 1
        else ! micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('corr_chi_Ns',0.0,1.,k)
        endif ! micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( micro_correlations(idx_s,idx_rg,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_rg',micro_correlations(idx_s,idx_rg,k),1.,k)
           corravg_count(idx_cor_chi_rg,k) = corravg_count(idx_cor_chi_rg,k) + 1
        else ! micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('corr_chi_rg',0.0,1.,k)
        endif ! micro_correlations(1,iqg,k) /= undef_corr
        if ( micro_correlations(idx_s,idx_Ng,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_Ng',micro_correlations(idx_s,idx_Ng,k),1.,k)
           corravg_count(idx_cor_chi_Ng,k) = corravg_count(idx_cor_chi_Ng,k) + 1
        else ! micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('corr_chi_Ng',0.0,1.,k)
        endif ! micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

!-------------
! Vertical Velocity
!-------------
      if(dopredictNc) then
        if ( micro_correlations(idx_w,idx_Nc,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_Nc',micro_correlations(idx_w,idx_Nc,k),1.,k)
           corravg_count(idx_cor_w_Nc,k) = corravg_count(idx_cor_w_Nc,k) + 1
        else ! micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('corr_w_Nc',0.0,1.,k)
        endif ! micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( micro_correlations(idx_w,idx_rr,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_rr',micro_correlations(idx_w,idx_rr,k),1.,k)
           corravg_count(idx_cor_w_rr,k) = corravg_count(idx_cor_w_rr,k) + 1
        else ! micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('corr_w_rr',0.0,1.,k)
        endif ! micro_correlations(1,iqr,k) /= undef_corr
        if ( micro_correlations(idx_w,idx_Nr,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_Nr',micro_correlations(idx_w,idx_Nr,k),1.,k)
           corravg_count(idx_cor_w_nr,k) = corravg_count(idx_cor_w_nr,k) + 1
        else ! micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('corr_w_Nr',0.0,1.,k)
        endif ! ! micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( micro_correlations(idx_w,idx_ri,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_ri',micro_correlations(idx_w,idx_ri,k),1.,k)
           corravg_count(idx_cor_w_ri,k) = corravg_count(idx_cor_w_ri,k) + 1
        else ! micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('corr_w_ri',0.0,1.,k)
        endif ! micro_correlations(1,iqci,k) /= undef_corr
        if ( micro_correlations(idx_w,idx_Ni,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_Ni',micro_correlations(idx_w,idx_Ni,k),1.,k)
           corravg_count(idx_cor_w_Ni,k) = corravg_count(idx_cor_w_Ni,k) + 1
        else ! micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('corr_w_Ni',0.0,1.,k)
        endif ! micro_correlations(1,inci,k) /= undef_corr
        if ( micro_correlations(idx_w,idx_rs,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_rs',micro_correlations(idx_w,idx_rs,k),1.,k)
           corravg_count(idx_cor_w_rs,k) = corravg_count(idx_cor_w_rs,k) + 1
        else ! micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('corr_w_rs',0.0,1.,k)
        endif ! micro_correlations(1,iqs,k) /= undef_corr
        if ( micro_correlations(idx_w,idx_Ns,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_Ns',micro_correlations(idx_w,idx_Ns,k),1.,k)
           corravg_count(idx_cor_w_Ns,k) = corravg_count(idx_cor_w_Ns,k) + 1
        else ! micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('corr_w_Ns',0.0,1.,k)
        endif ! micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( micro_correlations(idx_w,idx_rg,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_rg',micro_correlations(idx_w,idx_rg,k),1.,k)
           corravg_count(idx_cor_w_rg,k) = corravg_count(idx_cor_w_rg,k) + 1
        else ! micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('corr_w_rg',0.0,1.,k)
        endif ! micro_correlations(1,iqg,k) /= undef_corr
        if ( micro_correlations(idx_w,idx_Ng,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_Ng',micro_correlations(idx_w,idx_Ng,k),1.,k)
           corravg_count(idx_cor_w_Ng,k) = corravg_count(idx_cor_w_Ng,k) + 1
        else ! micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('corr_w_Ng',0.0,1.,k)
        endif ! micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

!-------------
! Cloud droplet number concentration
!-------------
      if (dopredictNc) then
     
        if (doprecip) then
           if ( micro_correlations(idx_Nc,idx_rr,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_rr',micro_correlations(idx_Nc,idx_rr,k),1.,k)
              corravg_count(idx_cor_Nc_rr,k) = corravg_count(idx_cor_Nc_rr,k) + 1
           else ! micro_correlations(incl,iqr,k) == undef_corr
              call hbuf_put_level('corr_Nc_rr',0.0,1.,k)
           endif ! micro_correlations(incl,iqr,k) /= undef_corr
           if ( micro_correlations(idx_Nc,idx_Nr,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_Nr',micro_correlations(idx_Nc,idx_Nr,k),1.,k)
              corravg_count(idx_cor_Nc_nr,k) = corravg_count(idx_cor_Nc_nr,k) + 1
           else ! micro_correlations(incl,inr,k) == undef_corr
              call hbuf_put_level('corr_Nc_Nr',0.0,1.,k)
           endif ! micro_correlations(incl,inr,k) /= undef_corr
        end if
     
        if (doicemicro) then
           if ( micro_correlations(idx_Nc,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_ri',micro_correlations(idx_Nc,idx_ri,k),1.,k) 
              corravg_count(idx_cor_Nc_ri,k) = corravg_count(idx_cor_Nc_ri,k) + 1
           else ! micro_correlations(incl,iqci,k) == undef_corr
              call hbuf_put_level('corr_Nc_ri',0.0,1.,k)
           endif ! micro_correlations(incl,iqci,k) /= undef_corr
           if ( micro_correlations(idx_Nc,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_Ni',micro_correlations(idx_Nc,idx_Ni,k),1.,k)
              corravg_count(idx_cor_Nc_Ni,k) = corravg_count(idx_cor_Nc_Ni,k) + 1
           else ! micro_correlations(incl,inci,k) == undef_corr
              call hbuf_put_level('corr_Nc_Ni',0.0,1.,k)
           endif ! micro_correlations(incl,inci,k) /= undef_corr
           if ( micro_correlations(idx_Nc,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_rs',micro_correlations(idx_Nc,idx_rs,k),1.,k)
              corravg_count(idx_cor_Nc_rs,k) = corravg_count(idx_cor_Nc_rs,k) + 1
           else ! micro_correlations(incl,iqs,k) == undef_corr
              call hbuf_put_level('corr_Nc_rs',0.0,1.,k)
           endif ! micro_correlations(incl,iqs,k) /= undef_corr
           if ( micro_correlations(idx_Nc,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_Ns',micro_correlations(idx_Nc,idx_Ns,k),1.,k)
              corravg_count(idx_cor_Nc_Ns,k) = corravg_count(idx_cor_Nc_Ns,k) + 1
           else ! micro_correlations(incl,ins,k) == undef_corr
              call hbuf_put_level('corr_Nc_Ns',0.0,1.,k)
           endif ! micro_correlations(incl,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( micro_correlations(idx_Nc,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_rg',micro_correlations(idx_Nc,idx_rg,k),1.,k)
              corravg_count(idx_cor_Nc_rg,k) = corravg_count(idx_cor_Nc_rg,k) + 1
           else ! micro_correlations(incl,iqg,k) == undef_corr
              call hbuf_put_level('corr_Nc_rg',0.0,1.,k)
           endif ! micro_correlations(incl,iqg,k) /= undef_corr
           if ( micro_correlations(idx_Nc,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_Ng',micro_correlations(idx_Nc,idx_Ng,k),1.,k)
              corravg_count(idx_cor_Nc_Ng,k) = corravg_count(idx_cor_Nc_Ng,k) + 1
           else ! micro_correlations(incl,ing,k) == undef_corr
              call hbuf_put_level('corr_Nc_Ng',0.0,1.,k)
           endif ! micro_correlations(incl,ing,k) /= undef_corr
        end if   

      end if !dopredictNc
  
!-------------
! Rainwater mixing ratio
!-------------
      if (doprecip) then 
      
         if ( micro_correlations(idx_rr,idx_Nr,k) /= undef_corr ) then
             call hbuf_put_level('corr_rr_Nr',micro_correlations(idx_rr,idx_Nr,k),1.,k)
             corravg_count(idx_cor_rr_nr,k) = corravg_count(idx_cor_rr_nr,k) + 1
          else ! micro_correlations(iqr,inr,k) == undef_corr
             call hbuf_put_level('corr_rr_Nr',0.0,1.,k)
          endif ! micro_correlations(iqr,inr,k) /= undef_corr

        if (doicemicro) then
           if ( micro_correlations(idx_rr,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('corr_rr_ri',micro_correlations(idx_rr,idx_ri,k),1.,k)
              corravg_count(idx_cor_rr_ri,k) = corravg_count(idx_cor_rr_ri,k) + 1
           else ! micro_correlations(iqr,iqci,k) == undef_corr
              call hbuf_put_level('corr_rr_ri',0.0,1.,k)
           endif ! micro_correlations(iqr,iqci,k) /= undef_corr
           if ( micro_correlations(idx_rr,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('corr_rr_Ni',micro_correlations(idx_rr,idx_Ni,k),1.,k)
              corravg_count(idx_cor_rr_Ni,k) = corravg_count(idx_cor_rr_Ni,k) + 1
           else ! micro_correlations(iqr,inci,k) == undef_corr
              call hbuf_put_level('corr_rr_Ni',0.0,1.,k)
           endif ! micro_correlations(iqr,inci,k) /= undef_corr
           if ( micro_correlations(idx_rr,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('corr_rr_rs',micro_correlations(idx_rr,idx_rs,k),1.,k)
              corravg_count(idx_cor_rr_rs,k) = corravg_count(idx_cor_rr_rs,k) + 1
           else ! micro_correlations(iqr,iqs,k) == undef_corr
              call hbuf_put_level('corr_rr_rs',0.0,1.,k)
           endif ! micro_correlations(iqr,iqs,k) /= undef_corr
           if ( micro_correlations(idx_rr,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('corr_rr_Ns',micro_correlations(idx_rr,idx_Ns,k),1.,k)
              corravg_count(idx_cor_rr_Ns,k) = corravg_count(idx_cor_rr_Ns,k) + 1
           else ! micro_correlations(iqr,ins,k) == undef_corr
              call hbuf_put_level('corr_rr_Ns',0.0,1.,k)
           endif ! micro_correlations(iqr,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( micro_correlations(idx_rr,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('corr_rr_rg',micro_correlations(idx_rr,idx_rg,k),1.,k)
              corravg_count(idx_cor_rr_rg,k) = corravg_count(idx_cor_rr_rg,k) + 1
           else ! micro_correlations(iqr,iqg,k) == undef_corr
              call hbuf_put_level('corr_rr_rg',0.0,1.,k)
           endif ! micro_correlations(iqr,iqg,k) /= undef_corr
           if ( micro_correlations(idx_rr,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('corr_rr_Ng',micro_correlations(idx_rr,idx_Ng,k),1.,k)
              corravg_count(idx_cor_rr_Ng,k) = corravg_count(idx_cor_rr_Ng,k) + 1
           else ! micro_correlations(iqr,ing,k) == undef_corr
              call hbuf_put_level('corr_rr_Ng',0.0,1.,k)
           endif ! micro_correlations(iqr,ing,k) /= undef_corr
        end if
!-------------
! Rainwater number concentration
!-------------

        if (doicemicro) then
           if ( micro_correlations(idx_Nr,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nr_ri',micro_correlations(idx_Nr,idx_ri,k),1.,k)
              corravg_count(idx_cor_Nr_ri,k) = corravg_count(idx_cor_Nr_ri,k) + 1
           else ! micro_correlations(inr,iqci,k) == undef_corr
              call hbuf_put_level('corr_Nr_ri',0.0,1.,k)
           endif ! micro_correlations(inr,iqci,k) /= undef_corr
           if ( micro_correlations(idx_Nr,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nr_Ni',micro_correlations(idx_Nr,idx_Ni,k),1.,k)
              corravg_count(idx_cor_Nr_Ni,k) = corravg_count(idx_cor_Nr_Ni,k) + 1
           else ! micro_correlations(inr,inci,k) == undef_corr
              call hbuf_put_level('corr_Nr_Ni',0.0,1.,k)
           endif ! micro_correlations(inr,inci,k) /= undef_corr
           if ( micro_correlations(idx_Nr,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nr_rs',micro_correlations(idx_Nr,idx_rs,k),1.,k)
              corravg_count(idx_cor_Nr_rs,k) = corravg_count(idx_cor_Nr_rs,k) + 1
           else ! micro_correlations(inr,iqs,k) == undef_corr
              call hbuf_put_level('corr_Nr_rs',0.0,1.,k)
           endif ! micro_correlations(inr,iqs,k) /= undef_corr
           if ( micro_correlations(idx_Nr,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nr_Ns',micro_correlations(idx_Nr,idx_Ns,k),1.,k)
              corravg_count(idx_cor_Nr_Ns,k) = corravg_count(idx_cor_Nr_Ns,k) + 1
           else ! micro_correlations(inr,ins,k) == undef_corr
              call hbuf_put_level('corr_Nr_Ns',0.0,1.,k)
           endif ! micro_correlations(inr,ins,k) /= undef_corr
        end if

        if (dograupel) then
        if ( micro_correlations(idx_Nr,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nr_rg',micro_correlations(idx_Nr,idx_rg,k),1.,k)
              corravg_count(idx_cor_Nr_rg,k) = corravg_count(idx_cor_Nr_rg,k) + 1
           else ! micro_correlations(inr,iqg,k) == undef_corr
              call hbuf_put_level('corr_Nr_rg',0.0,1.,k)
           endif ! micro_correlations(inr,iqg,k) /= undef_corr
           if ( micro_correlations(idx_Nr,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nr_Ng',micro_correlations(idx_Nr,idx_Ng,k),1.,k)
              corravg_count(idx_cor_Nr_Ng,k) = corravg_count(idx_cor_Nr_Ng,k) + 1
           else ! micro_correlations(inr,ing,k) == undef_corr
              call hbuf_put_level('corr_Nr_Ng',0.0,1.,k)
           endif ! micro_correlations(inr,ing,k) /= undef_corr
        end if

      end if !doprecip
!-------------
! Cloud-ice mixing ratio
!-------------

      if (doicemicro) then   
         if ( micro_correlations(idx_ri,idx_Ni,k) /= undef_corr ) then
            call hbuf_put_level('corr_ri_Ni',micro_correlations(idx_ri,idx_Ni,k),1.,k)
            corravg_count(idx_cor_ri_Ni,k) = corravg_count(idx_cor_ri_Ni,k) + 1
         else ! micro_correlations(iqci,inci,k) == undef_corr
            call hbuf_put_level('corr_ri_Ni',0.0,1.,k)
         endif ! micro_correlations(iqci,inci,k) /= undef_corr
         if ( micro_correlations(idx_ri,idx_rs,k) /= undef_corr ) then
            call hbuf_put_level('corr_ri_rs',micro_correlations(idx_ri,idx_rs,k),1.,k)
            corravg_count(idx_cor_ri_rs,k) = corravg_count(idx_cor_ri_rs,k) + 1
         else ! micro_correlations(iqci,iqs,k) == undef_corr
            call hbuf_put_level('corr_ri_rs',0.0,1.,k)
         endif ! micro_correlations(iqci,iqs,k) /= undef_corr
         if ( micro_correlations(idx_ri,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('corr_ri_Ns',micro_correlations(idx_ri,idx_Ns,k),1.,k)
            corravg_count(idx_cor_ri_Ns,k) = corravg_count(idx_cor_ri_Ns,k) + 1
         else ! micro_correlations(iqci,ins,k) == undef_corr
            call hbuf_put_level('corr_ri_Ns',0.0,1.,k)
         endif ! micro_correlations(iqci,ins,k) /= undef_corr
      end if   

      if (dograupel) then
         if ( micro_correlations(idx_ri,idx_rg,k) /= undef_corr ) then
            call hbuf_put_level('corr_ri_rg',micro_correlations(idx_ri,idx_rg,k),1.,k)
            corravg_count(idx_cor_ri_rg,k) = corravg_count(idx_cor_ri_rg,k) + 1
         else ! micro_correlations(iqci,iqg,k) == undef_corr
            call hbuf_put_level('corr_ri_rg',0.0,1.,k)
         endif ! micro_correlations(iqci,iqg,k) /= undef_corr
         if ( micro_correlations(idx_ri,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('corr_ri_Ng',micro_correlations(idx_ri,idx_Ng,k),1.,k)
            corravg_count(idx_cor_ri_Ng,k) = corravg_count(idx_cor_ri_Ng,k) + 1
         else ! micro_correlations(iqci,ing,k) == undef_corr
            call hbuf_put_level('corr_ri_Ng',0.0,1.,k)
         endif ! micro_correlations(iqci,ing,k) /= undef_corr
      end if
 
!-------------
! Cloud-ice number concentration
!-------------
      if (doicemicro) then   
         if ( micro_correlations(idx_Ni,idx_rs,k) /= undef_corr ) then
            call hbuf_put_level('corr_Ni_rs',micro_correlations(idx_Ni,idx_rs,k),1.,k)
            corravg_count(idx_cor_Ni_rs,k) = corravg_count(idx_cor_Ni_rs,k) + 1
         else ! micro_correlations(inci,iqs,k) == undef_corr
            call hbuf_put_level('corr_Ni_rs',0.0,1.,k)
         endif ! micro_correlations(inci,iqs,k) /= undef_corr
         if ( micro_correlations(idx_Ni,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('corr_Ni_Ns',micro_correlations(idx_Ni,idx_Ns,k),1.,k)
            corravg_count(idx_cor_Ni_Ns,k) = corravg_count(idx_cor_Ni_Ns,k) + 1
         else ! micro_correlations(inci,ins,k) == undef_corr
            call hbuf_put_level('corr_Ni_Ns',0.0,1.,k)
         endif ! micro_correlations(inci,ins,k) /= undef_corr
      end if

      if (dograupel) then     
         if ( micro_correlations(idx_Ni,idx_rg,k) /= undef_corr ) then
            call hbuf_put_level('corr_Ni_rg',micro_correlations(idx_Ni,idx_rg,k),1.,k)
            corravg_count(idx_cor_Ni_rg,k) = corravg_count(idx_cor_Ni_rg,k) + 1
         else ! micro_correlations(inci,iqg,k) == undef_corr
            call hbuf_put_level('corr_Ni_rg',0.0,1.,k)
         endif ! micro_correlations(inci,iqg,k) /= undef_corr
         if ( micro_correlations(idx_Ni,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('corr_Ni_Ng',micro_correlations(idx_Ni,idx_Ng,k),1.,k)
            corravg_count(idx_cor_Ni_Ng,k) = corravg_count(idx_cor_Ni_Ng,k) + 1
         else ! micro_correlations(inci,ing,k) == undef_corr
            call hbuf_put_level('corr_Ni_Ng',0.0,1.,k)
         endif ! micro_correlations(inci,ing,k) /= undef_corr
      end if 
!-------------
! Snow mixing ratio
!-------------

      if (doicemicro) then
         if ( micro_correlations(idx_rs,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('corr_rs_Ns',micro_correlations(idx_rs,idx_Ns,k),1.,k)
            corravg_count(idx_cor_rs_Ns,k) = corravg_count(idx_cor_rs_Ns,k) + 1
         else ! micro_correlations(iqs,ins,k) == undef_corr
            call hbuf_put_level('corr_rs_Ns',0.0,1.,k)
         endif ! micro_correlations(iqs,ins,k) /= undef_corr
       
          if (dograupel) then   
             if ( micro_correlations(idx_rs,idx_rg,k) /= undef_corr ) then
                call hbuf_put_level('corr_rs_rg',micro_correlations(idx_rs,idx_rg,k),1.,k)
                corravg_count(idx_cor_rs_rg,k) = corravg_count(idx_cor_rs_rg,k) + 1
             else ! micro_correlations(iqs,iqg,k) == undef_corr
                call hbuf_put_level('corr_rs_rg',0.0,1.,k)
             endif ! micro_correlations(iqs,iqg,k) /= undef_corr
             if ( micro_correlations(idx_rs,idx_Ng,k) /= undef_corr ) then
                call hbuf_put_level('corr_rs_Ng',micro_correlations(idx_rs,idx_Ng,k),1.,k)
                corravg_count(idx_cor_rs_Ng,k) = corravg_count(idx_cor_rs_Ng,k) + 1
             else ! micro_correlations(iqs,ing,k) == undef_corr
                call hbuf_put_level('corr_rs_Ng',0.0,1.,k)
             endif ! micro_correlations(iqs,ing,k) /= undef_corr
!-------------
! Snow number concentration
!-------------
             if ( micro_correlations(idx_Ns,idx_rg,k) /= undef_corr ) then
                call hbuf_put_level('corr_Ns_rg',micro_correlations(idx_Ns,idx_rg,k),1.,k)
                corravg_count(idx_cor_Ns_rg,k) = corravg_count(idx_cor_Ns_rg,k) + 1
             else ! micro_correlations(ins,iqg,k) == undef_corr
                call hbuf_put_level('corr_Ns_rg',0.0,1.,k)
             endif ! micro_correlations(ins,iqg,k) /= undef_corr
             if ( micro_correlations(idx_Ns,idx_Ng,k) /= undef_corr ) then
                call hbuf_put_level('corr_Ns_Ng',micro_correlations(idx_Ns,idx_Ng,k),1.,k)
                corravg_count(idx_cor_Ns_Ng,k) = corravg_count(idx_cor_Ns_Ng,k) + 1
             else ! micro_correlations(ins,ing,k) == undef_corr
                call hbuf_put_level('corr_Ns_Ng',0.0,1.,k)
             endif ! micro_correlations(ins,ing,k) /= undef_corr
          end if
      end if

!-------------
! Graupel mixing ratio / number concentration
!-------------
      if (dograupel) then      
        if ( micro_correlations(idx_rg,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('corr_rg_Ng',micro_correlations(idx_rg,idx_Ng,k),1.,k)
            corravg_count(idx_cor_rg_Ng,k) = corravg_count(idx_cor_rg_Ng,k) + 1
         else ! micro_correlations(iqg,ing,k) == undef_corr
            call hbuf_put_level('corr_rg_Ng',0.0,1.,k)
         endif ! micro_correlations(iqg,ing,k) /= undef_corr
      endif

   enddo ! k = 1, nzm, 1

!------------
! In-cloud correlations
!------------
   do k = 1, nzm, 1
!-------------
! chi
!-------------
      if ( ic_micro_correlations(idx_s,idx_w,k) /= undef_corr ) then
         call hbuf_put_level('ic_corr_chi_w',ic_micro_correlations(idx_s,idx_w,k),1.,k)
         corravg_count(idx_ic_cor_chi_w,k) = corravg_count(idx_ic_cor_chi_w,k) + 1
      else ! ic_micro_correlations(1,incl,k) == undef_corr
         call hbuf_put_level('ic_corr_chi_w',0.0,1.,k)
      endif ! ic_micro_correlations(1,incl,k) /= undef_corr
      
      if(dopredictNc) then
        if ( ic_micro_correlations(idx_s,idx_Nc,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_Nc',ic_micro_correlations(idx_s,idx_Nc,k),1.,k)
           corravg_count(idx_ic_cor_chi_Nc,k) = corravg_count(idx_ic_cor_chi_Nc,k) + 1
        else ! ic_micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_Nc',0.0,1.,k)
        endif ! ic_micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( ic_micro_correlations(idx_s,idx_rr,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_rr',ic_micro_correlations(idx_s,idx_rr,k),1.,k)
           corravg_count(idx_ic_cor_chi_rr,k) = corravg_count(idx_ic_cor_chi_rr,k) + 1
        else ! ic_micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_rr',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqr,k) /= undef_corr
        if ( ic_micro_correlations(idx_s,idx_Nr,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_Nr',ic_micro_correlations(idx_s,idx_Nr,k),1.,k)
           corravg_count(idx_ic_cor_chi_Nr,k) = corravg_count(idx_ic_cor_chi_Nr,k) + 1
        else ! ic_micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_Nr',0.0,1.,k)
        endif ! ! ic_micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( ic_micro_correlations(idx_s,idx_ri,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_ri',ic_micro_correlations(idx_s,idx_ri,k),1.,k)
           corravg_count(idx_ic_cor_chi_ri,k) = corravg_count(idx_ic_cor_chi_ri,k) + 1
        else ! ic_micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_ri',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqci,k) /= undef_corr
        if ( ic_micro_correlations(idx_s,idx_Ni,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_Ni',ic_micro_correlations(idx_s,idx_Ni,k),1.,k)
           corravg_count(idx_ic_cor_chi_Ni,k) = corravg_count(idx_ic_cor_chi_Ni,k) + 1
        else ! ic_micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_Ni',0.0,1.,k)
        endif ! ic_micro_correlations(1,inci,k) /= undef_corr
        if ( ic_micro_correlations(idx_s,idx_rs,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_rs',ic_micro_correlations(idx_s,idx_rs,k),1.,k)
           corravg_count(idx_ic_cor_chi_rs,k) = corravg_count(idx_ic_cor_chi_rs,k) + 1
        else ! ic_micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_rs',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqs,k) /= undef_corr
        if ( ic_micro_correlations(idx_s,idx_Ns,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_Ns',ic_micro_correlations(idx_s,idx_Ns,k),1.,k)
           corravg_count(idx_ic_cor_chi_Ns,k) = corravg_count(idx_ic_cor_chi_Ns,k) + 1
        else ! ic_micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_Ns',0.0,1.,k)
        endif ! ic_micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( ic_micro_correlations(idx_s,idx_rg,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_rg',ic_micro_correlations(idx_s,idx_rg,k),1.,k)
           corravg_count(idx_ic_cor_chi_rg,k) = corravg_count(idx_ic_cor_chi_rg,k) + 1
        else ! ic_micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_rg',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqg,k) /= undef_corr
        if ( ic_micro_correlations(idx_s,idx_Ng,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_Ng',ic_micro_correlations(idx_s,idx_Ng,k),1.,k)
           corravg_count(idx_ic_cor_chi_Ng,k) = corravg_count(idx_ic_cor_chi_Ng,k) + 1
        else ! ic_micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_Ng',0.0,1.,k)
        endif ! ic_micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

!-------------
! Vertical Velocity
!-------------
      if(dopredictNc) then
        if ( ic_micro_correlations(idx_w,idx_Nc,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_Nc',ic_micro_correlations(idx_w,idx_Nc,k),1.,k)
           corravg_count(idx_ic_cor_w_Nc,k) = corravg_count(idx_ic_cor_w_Nc,k) + 1
        else ! ic_micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('ic_corr_w_Nc',0.0,1.,k)
        endif ! ic_micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( ic_micro_correlations(idx_w,idx_rr,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_rr',ic_micro_correlations(idx_w,idx_rr,k),1.,k)
           corravg_count(idx_ic_cor_w_rr,k) = corravg_count(idx_ic_cor_w_rr,k) + 1
        else ! ic_micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('ic_corr_w_rr',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqr,k) /= undef_corr
        if ( ic_micro_correlations(idx_w,idx_Nr,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_Nr',ic_micro_correlations(idx_w,idx_Nr,k),1.,k)
           corravg_count(idx_ic_cor_w_nr,k) = corravg_count(idx_ic_cor_w_nr,k) + 1
        else ! ic_micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('ic_corr_w_Nr',0.0,1.,k)
        endif ! ! ic_micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( ic_micro_correlations(idx_w,idx_ri,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_ri',ic_micro_correlations(idx_w,idx_ri,k),1.,k)
           corravg_count(idx_ic_cor_w_ri,k) = corravg_count(idx_ic_cor_w_ri,k) + 1
        else ! ic_micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('ic_corr_w_ri',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqci,k) /= undef_corr
        if ( ic_micro_correlations(idx_w,idx_Ni,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_Ni',ic_micro_correlations(idx_w,idx_Ni,k),1.,k)
           corravg_count(idx_ic_cor_w_Ni,k) = corravg_count(idx_ic_cor_w_Ni,k) + 1
        else ! ic_micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('ic_corr_w_Ni',0.0,1.,k)
        endif ! ic_micro_correlations(1,inci,k) /= undef_corr
        if ( ic_micro_correlations(idx_w,idx_rs,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_rs',ic_micro_correlations(idx_w,idx_rs,k),1.,k)
           corravg_count(idx_ic_cor_w_rs,k) = corravg_count(idx_ic_cor_w_rs,k) + 1
        else ! ic_micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('ic_corr_w_rs',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqs,k) /= undef_corr
        if ( ic_micro_correlations(idx_w,idx_Ns,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_Ns',ic_micro_correlations(idx_w,idx_Ns,k),1.,k)
           corravg_count(idx_ic_cor_w_Ns,k) = corravg_count(idx_ic_cor_w_Ns,k) + 1
        else ! ic_micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('ic_corr_w_Ns',0.0,1.,k)
        endif ! ic_micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( ic_micro_correlations(idx_w,idx_rg,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_rg',ic_micro_correlations(idx_w,idx_rg,k),1.,k)
           corravg_count(idx_ic_cor_w_rg,k) = corravg_count(idx_ic_cor_w_rg,k) + 1
        else ! ic_micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('ic_corr_w_rg',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqg,k) /= undef_corr
        if ( ic_micro_correlations(idx_w,idx_Ng,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_Ng',ic_micro_correlations(idx_w,idx_Ng,k),1.,k)
           corravg_count(idx_ic_cor_w_Ng,k) = corravg_count(idx_ic_cor_w_Ng,k) + 1
        else ! ic_micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('ic_corr_w_Ng',0.0,1.,k)
        endif ! ic_micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

!-------------
! Cloud droplet number concentration
!-------------
      if (dopredictNc) then
     
        if (doprecip) then
           if ( ic_micro_correlations(idx_Nc,idx_rr,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_rr',ic_micro_correlations(idx_Nc,idx_rr,k),1.,k)
              corravg_count(idx_ic_cor_Nc_rr,k) = corravg_count(idx_ic_cor_Nc_rr,k) + 1
           else ! ic_micro_correlations(incl,iqr,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_rr',0.0,1.,k)
           endif ! ic_micro_correlations(incl,iqr,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nc,idx_Nr,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_Nr',ic_micro_correlations(idx_Nc,idx_Nr,k),1.,k)
              corravg_count(idx_ic_cor_Nc_nr,k) = corravg_count(idx_ic_cor_Nc_nr,k) + 1
           else ! ic_micro_correlations(incl,inr,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_Nr',0.0,1.,k)
           endif ! ic_micro_correlations(incl,inr,k) /= undef_corr
        end if
     
        if (doicemicro) then
           if ( ic_micro_correlations(idx_Nc,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_ri',ic_micro_correlations(idx_Nc,idx_ri,k),1.,k) 
              corravg_count(idx_ic_cor_Nc_ri,k) = corravg_count(idx_ic_cor_Nc_ri,k) + 1
           else ! ic_micro_correlations(incl,iqci,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_ri',0.0,1.,k)
           endif ! ic_micro_correlations(incl,iqci,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nc,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_Ni',ic_micro_correlations(idx_Nc,idx_Ni,k),1.,k)
              corravg_count(idx_ic_cor_Nc_Ni,k) = corravg_count(idx_ic_cor_Nc_Ni,k) + 1
           else ! ic_micro_correlations(incl,inci,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_Ni',0.0,1.,k)
           endif ! ic_micro_correlations(incl,inci,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nc,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_rs',ic_micro_correlations(idx_Nc,idx_rs,k),1.,k)
              corravg_count(idx_ic_cor_Nc_rs,k) = corravg_count(idx_ic_cor_Nc_rs,k) + 1
           else ! ic_micro_correlations(incl,iqs,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_rs',0.0,1.,k)
           endif ! ic_micro_correlations(incl,iqs,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nc,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_Ns',ic_micro_correlations(idx_Nc,idx_Ns,k),1.,k)
              corravg_count(idx_ic_cor_Nc_Ns,k) = corravg_count(idx_ic_cor_Nc_Ns,k) + 1
           else ! ic_micro_correlations(incl,ins,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_Ns',0.0,1.,k)
           endif ! ic_micro_correlations(incl,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( ic_micro_correlations(idx_Nc,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_rg',ic_micro_correlations(idx_Nc,idx_rg,k),1.,k)
              corravg_count(idx_ic_cor_Nc_rg,k) = corravg_count(idx_ic_cor_Nc_rg,k) + 1
           else ! ic_micro_correlations(incl,iqg,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_rg',0.0,1.,k)
           endif ! ic_micro_correlations(incl,iqg,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nc,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_Ng',ic_micro_correlations(idx_Nc,idx_Ng,k),1.,k)
              corravg_count(idx_ic_cor_Nc_Ng,k) = corravg_count(idx_ic_cor_Nc_Ng,k) + 1
           else ! ic_micro_correlations(incl,ing,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_Ng',0.0,1.,k)
           endif ! ic_micro_correlations(incl,ing,k) /= undef_corr
        end if   

      end if !dopredictNc
  
!-------------
! Rainwater mixing ratio
!-------------
      if (doprecip) then 
      
         if ( ic_micro_correlations(idx_rr,idx_Nr,k) /= undef_corr ) then
             call hbuf_put_level('ic_corr_rr_Nr',ic_micro_correlations(idx_rr,idx_Nr,k),1.,k)
             corravg_count(idx_ic_cor_rr_nr,k) = corravg_count(idx_ic_cor_rr_nr,k) + 1
          else ! ic_micro_correlations(iqr,inr,k) == undef_corr
             call hbuf_put_level('ic_corr_rr_Nr',0.0,1.,k)
          endif ! ic_micro_correlations(iqr,inr,k) /= undef_corr

        if (doicemicro) then
           if ( ic_micro_correlations(idx_rr,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_rr_ri',ic_micro_correlations(idx_rr,idx_ri,k),1.,k)
              corravg_count(idx_ic_cor_rr_ri,k) = corravg_count(idx_ic_cor_rr_ri,k) + 1
           else ! ic_micro_correlations(iqr,iqci,k) == undef_corr
              call hbuf_put_level('ic_corr_rr_ri',0.0,1.,k)
           endif ! ic_micro_correlations(iqr,iqci,k) /= undef_corr
           if ( ic_micro_correlations(idx_rr,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_rr_Ni',ic_micro_correlations(idx_rr,idx_Ni,k),1.,k)
              corravg_count(idx_ic_cor_rr_Ni,k) = corravg_count(idx_ic_cor_rr_Ni,k) + 1
           else ! ic_micro_correlations(iqr,inci,k) == undef_corr
              call hbuf_put_level('ic_corr_rr_Ni',0.0,1.,k)
           endif ! ic_micro_correlations(iqr,inci,k) /= undef_corr
           if ( ic_micro_correlations(idx_rr,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_rr_rs',ic_micro_correlations(idx_rr,idx_rs,k),1.,k)
              corravg_count(idx_ic_cor_rr_rs,k) = corravg_count(idx_ic_cor_rr_rs,k) + 1
           else ! ic_micro_correlations(iqr,iqs,k) == undef_corr
              call hbuf_put_level('ic_corr_rr_rs',0.0,1.,k)
           endif ! ic_micro_correlations(iqr,iqs,k) /= undef_corr
           if ( ic_micro_correlations(idx_rr,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_rr_Ns',ic_micro_correlations(idx_rr,idx_Ns,k),1.,k)
              corravg_count(idx_ic_cor_rr_Ns,k) = corravg_count(idx_ic_cor_rr_Ns,k) + 1
           else ! ic_micro_correlations(iqr,ins,k) == undef_corr
              call hbuf_put_level('ic_corr_rr_Ns',0.0,1.,k)
           endif ! ic_micro_correlations(iqr,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( ic_micro_correlations(idx_rr,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_rr_rg',ic_micro_correlations(idx_rr,idx_rg,k),1.,k)
              corravg_count(idx_ic_cor_rr_rg,k) = corravg_count(idx_ic_cor_rr_rg,k) + 1
           else ! ic_micro_correlations(iqr,iqg,k) == undef_corr
              call hbuf_put_level('ic_corr_rr_rg',0.0,1.,k)
           endif ! ic_micro_correlations(iqr,iqg,k) /= undef_corr
           if ( ic_micro_correlations(idx_rr,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_rr_Ng',ic_micro_correlations(idx_rr,idx_Ng,k),1.,k)
              corravg_count(idx_ic_cor_rr_Ng,k) = corravg_count(idx_ic_cor_rr_Ng,k) + 1
           else ! ic_micro_correlations(iqr,ing,k) == undef_corr
              call hbuf_put_level('ic_corr_rr_Ng',0.0,1.,k)
           endif ! ic_micro_correlations(iqr,ing,k) /= undef_corr
        end if
!-------------
! Rainwater number concentration
!-------------

        if (doicemicro) then
           if ( ic_micro_correlations(idx_Nr,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nr_ri',ic_micro_correlations(idx_Nr,idx_ri,k),1.,k)
              corravg_count(idx_ic_cor_Nr_ri,k) = corravg_count(idx_ic_cor_Nr_ri,k) + 1
           else ! ic_micro_correlations(inr,iqci,k) == undef_corr
              call hbuf_put_level('ic_corr_Nr_ri',0.0,1.,k)
           endif ! ic_micro_correlations(inr,iqci,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nr,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nr_Ni',ic_micro_correlations(idx_Nr,idx_Ni,k),1.,k)
              corravg_count(idx_ic_cor_Nr_Ni,k) = corravg_count(idx_ic_cor_Nr_Ni,k) + 1
           else ! ic_micro_correlations(inr,inci,k) == undef_corr
              call hbuf_put_level('ic_corr_Nr_Ni',0.0,1.,k)
           endif ! ic_micro_correlations(inr,inci,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nr,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nr_rs',ic_micro_correlations(idx_Nr,idx_rs,k),1.,k)
              corravg_count(idx_ic_cor_Nr_rs,k) = corravg_count(idx_ic_cor_Nr_rs,k) + 1
           else ! ic_micro_correlations(inr,iqs,k) == undef_corr
              call hbuf_put_level('ic_corr_Nr_rs',0.0,1.,k)
           endif ! ic_micro_correlations(inr,iqs,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nr,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nr_Ns',ic_micro_correlations(idx_Nr,idx_Ns,k),1.,k)
              corravg_count(idx_ic_cor_Nr_Ns,k) = corravg_count(idx_ic_cor_Nr_Ns,k) + 1
           else ! ic_micro_correlations(inr,ins,k) == undef_corr
              call hbuf_put_level('ic_corr_Nr_Ns',0.0,1.,k)
           endif ! ic_micro_correlations(inr,ins,k) /= undef_corr
        end if

        if (dograupel) then
        if ( ic_micro_correlations(idx_Nr,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nr_rg',ic_micro_correlations(idx_Nr,idx_rg,k),1.,k)
              corravg_count(idx_ic_cor_Nr_rg,k) = corravg_count(idx_ic_cor_Nr_rg,k) + 1
           else ! ic_micro_correlations(inr,iqg,k) == undef_corr
              call hbuf_put_level('ic_corr_Nr_rg',0.0,1.,k)
           endif ! ic_micro_correlations(inr,iqg,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nr,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nr_Ng',ic_micro_correlations(idx_Nr,idx_Ng,k),1.,k)
              corravg_count(idx_ic_cor_Nr_Ng,k) = corravg_count(idx_ic_cor_Nr_Ng,k) + 1
           else ! ic_micro_correlations(inr,ing,k) == undef_corr
              call hbuf_put_level('ic_corr_Nr_Ng',0.0,1.,k)
           endif ! ic_micro_correlations(inr,ing,k) /= undef_corr
        end if

      end if !doprecip
!-------------
! Cloud-ice mixing ratio
!-------------

      if (doicemicro) then   
         if ( ic_micro_correlations(idx_ri,idx_Ni,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_ri_Ni',ic_micro_correlations(idx_ri,idx_Ni,k),1.,k)
            corravg_count(idx_ic_cor_ri_Ni,k) = corravg_count(idx_ic_cor_ri_Ni,k) + 1
         else ! ic_micro_correlations(iqci,inci,k) == undef_corr
            call hbuf_put_level('ic_corr_ri_Ni',0.0,1.,k)
         endif ! ic_micro_correlations(iqci,inci,k) /= undef_corr
         if ( ic_micro_correlations(idx_ri,idx_rs,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_ri_rs',ic_micro_correlations(idx_ri,idx_rs,k),1.,k)
            corravg_count(idx_ic_cor_ri_rs,k) = corravg_count(idx_ic_cor_ri_rs,k) + 1
         else ! ic_micro_correlations(iqci,iqs,k) == undef_corr
            call hbuf_put_level('ic_corr_ri_rs',0.0,1.,k)
         endif ! ic_micro_correlations(iqci,iqs,k) /= undef_corr
         if ( ic_micro_correlations(idx_ri,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_ri_Ns',ic_micro_correlations(idx_ri,idx_Ns,k),1.,k)
            corravg_count(idx_ic_cor_ri_Ns,k) = corravg_count(idx_ic_cor_ri_Ns,k) + 1
         else ! ic_micro_correlations(iqci,ins,k) == undef_corr
            call hbuf_put_level('ic_corr_ri_Ns',0.0,1.,k)
         endif ! ic_micro_correlations(iqci,ins,k) /= undef_corr
      end if   

      if (dograupel) then
         if ( ic_micro_correlations(idx_ri,idx_rg,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_ri_rg',ic_micro_correlations(idx_ri,idx_rg,k),1.,k)
            corravg_count(idx_ic_cor_ri_rg,k) = corravg_count(idx_ic_cor_ri_rg,k) + 1
         else ! ic_micro_correlations(iqci,iqg,k) == undef_corr
            call hbuf_put_level('ic_corr_ri_rg',0.0,1.,k)
         endif ! ic_micro_correlations(iqci,iqg,k) /= undef_corr
         if ( ic_micro_correlations(idx_ri,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_ri_Ng',ic_micro_correlations(idx_ri,idx_Ng,k),1.,k)
            corravg_count(idx_ic_cor_ri_Ng,k) = corravg_count(idx_ic_cor_ri_Ng,k) + 1
         else ! ic_micro_correlations(iqci,ing,k) == undef_corr
            call hbuf_put_level('ic_corr_ri_Ng',0.0,1.,k)
         endif ! ic_micro_correlations(iqci,ing,k) /= undef_corr
      end if
 
!-------------
! Cloud-ice number concentration
!-------------
      if (doicemicro) then   
         if ( ic_micro_correlations(idx_Ni,idx_rs,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_Ni_rs',ic_micro_correlations(idx_Ni,idx_rs,k),1.,k)
            corravg_count(idx_ic_cor_Ni_rs,k) = corravg_count(idx_ic_cor_Ni_rs,k) + 1
         else ! ic_micro_correlations(inci,iqs,k) == undef_corr
            call hbuf_put_level('ic_corr_Ni_rs',0.0,1.,k)
         endif ! ic_micro_correlations(inci,iqs,k) /= undef_corr
         if ( ic_micro_correlations(idx_Ni,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_Ni_Ns',ic_micro_correlations(idx_Ni,idx_Ns,k),1.,k)
            corravg_count(idx_ic_cor_Ni_Ns,k) = corravg_count(idx_ic_cor_Ni_Ns,k) + 1
         else ! ic_micro_correlations(inci,ins,k) == undef_corr
            call hbuf_put_level('ic_corr_Ni_Ns',0.0,1.,k)
         endif ! ic_micro_correlations(inci,ins,k) /= undef_corr
      end if

      if (dograupel) then     
         if ( ic_micro_correlations(idx_Ni,idx_rg,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_Ni_rg',ic_micro_correlations(idx_Ni,idx_rg,k),1.,k)
            corravg_count(idx_ic_cor_Ni_rg,k) = corravg_count(idx_ic_cor_Ni_rg,k) + 1
         else ! ic_micro_correlations(inci,iqg,k) == undef_corr
            call hbuf_put_level('ic_corr_Ni_rg',0.0,1.,k)
         endif ! ic_micro_correlations(inci,iqg,k) /= undef_corr
         if ( ic_micro_correlations(idx_Ni,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_Ni_Ng',ic_micro_correlations(idx_Ni,idx_Ng,k),1.,k)
            corravg_count(idx_ic_cor_Ni_Ng,k) = corravg_count(idx_ic_cor_Ni_Ng,k) + 1
         else ! ic_micro_correlations(inci,ing,k) == undef_corr
            call hbuf_put_level('ic_corr_Ni_Ng',0.0,1.,k)
         endif ! ic_micro_correlations(inci,ing,k) /= undef_corr
      end if 
!-------------
! Snow mixing ratio
!-------------

      if (doicemicro) then
         if ( ic_micro_correlations(idx_rs,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_rs_Ns',ic_micro_correlations(idx_rs,idx_Ns,k),1.,k)
            corravg_count(idx_ic_cor_rs_Ns,k) = corravg_count(idx_ic_cor_rs_Ns,k) + 1
         else ! ic_micro_correlations(iqs,ins,k) == undef_corr
            call hbuf_put_level('ic_corr_rs_Ns',0.0,1.,k)
         endif ! ic_micro_correlations(iqs,ins,k) /= undef_corr
       
          if (dograupel) then   
             if ( ic_micro_correlations(idx_rs,idx_rg,k) /= undef_corr ) then
                call hbuf_put_level('ic_corr_rs_rg',ic_micro_correlations(idx_rs,idx_rg,k),1.,k)
                corravg_count(idx_ic_cor_rs_rg,k) = corravg_count(idx_ic_cor_rs_rg,k) + 1
             else ! ic_micro_correlations(iqs,iqg,k) == undef_corr
                call hbuf_put_level('ic_corr_rs_rg',0.0,1.,k)
             endif ! ic_micro_correlations(iqs,iqg,k) /= undef_corr
             if ( ic_micro_correlations(idx_rs,idx_Ng,k) /= undef_corr ) then
                call hbuf_put_level('ic_corr_rs_Ng',ic_micro_correlations(idx_rs,idx_Ng,k),1.,k)
                corravg_count(idx_ic_cor_rs_Ng,k) = corravg_count(idx_ic_cor_rs_Ng,k) + 1
             else ! ic_micro_correlations(iqs,ing,k) == undef_corr
                call hbuf_put_level('ic_corr_rs_Ng',0.0,1.,k)
             endif ! ic_micro_correlations(iqs,ing,k) /= undef_corr
!-------------
! Snow number concentration
!-------------
             if ( ic_micro_correlations(idx_Ns,idx_rg,k) /= undef_corr ) then
                call hbuf_put_level('ic_corr_Ns_rg',ic_micro_correlations(idx_Ns,idx_rg,k),1.,k)
                corravg_count(idx_ic_cor_Ns_rg,k) = corravg_count(idx_ic_cor_Ns_rg,k) + 1
             else ! ic_micro_correlations(ins,iqg,k) == undef_corr
                call hbuf_put_level('ic_corr_Ns_rg',0.0,1.,k)
             endif ! ic_micro_correlations(ins,iqg,k) /= undef_corr
             if ( ic_micro_correlations(idx_Ns,idx_Ng,k) /= undef_corr ) then
                call hbuf_put_level('ic_corr_Ns_Ng',ic_micro_correlations(idx_Ns,idx_Ng,k),1.,k)
                corravg_count(idx_ic_cor_Ns_Ng,k) = corravg_count(idx_ic_cor_Ns_Ng,k) + 1
             else ! ic_micro_correlations(ins,ing,k) == undef_corr
                call hbuf_put_level('ic_corr_Ns_Ng',0.0,1.,k)
             endif ! ic_micro_correlations(ins,ing,k) /= undef_corr
          end if
      end if

!-------------
! Graupel mixing ratio / number concentration
!-------------
      if (dograupel) then      
        if ( ic_micro_correlations(idx_rg,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_rg_Ng',ic_micro_correlations(idx_rg,idx_Ng,k),1.,k)
            corravg_count(idx_ic_cor_rg_Ng,k) = corravg_count(idx_ic_cor_rg_Ng,k) + 1
         else ! ic_micro_correlations(iqg,ing,k) == undef_corr
            call hbuf_put_level('ic_corr_rg_Ng',0.0,1.,k)
         endif ! ic_micro_correlations(iqg,ing,k) /= undef_corr
      endif

    end do ! do k = 1, nzm, 1

!------------
! Out of cloud correlations
!------------
   do k = 1, nzm, 1
!-------------
! chi
!-------------
      if ( oc_micro_correlations(idx_s,idx_w,k) /= undef_corr ) then
         call hbuf_put_level('oc_corr_chi_w',oc_micro_correlations(idx_s,idx_w,k),1.,k)
         corravg_count(idx_oc_cor_chi_w,k) = corravg_count(idx_oc_cor_chi_w,k) + 1
      else ! oc_micro_correlations(1,incl,k) == undef_corr
         call hbuf_put_level('oc_corr_chi_w',0.0,1.,k)
      endif ! oc_micro_correlations(1,incl,k) /= undef_corr
      
      if(dopredictNc) then
        if ( oc_micro_correlations(idx_s,idx_Nc,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_Nc',oc_micro_correlations(idx_s,idx_Nc,k),1.,k)
           corravg_count(idx_oc_cor_chi_Nc,k) = corravg_count(idx_oc_cor_chi_Nc,k) + 1
        else ! oc_micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_Nc',0.0,1.,k)
        endif ! oc_micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( oc_micro_correlations(idx_s,idx_rr,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_rr',oc_micro_correlations(idx_s,idx_rr,k),1.,k)
           corravg_count(idx_oc_cor_chi_rr,k) = corravg_count(idx_oc_cor_chi_rr,k) + 1
        else ! oc_micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_rr',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqr,k) /= undef_corr
        if ( oc_micro_correlations(idx_s,idx_Nr,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_Nr',oc_micro_correlations(idx_s,idx_Nr,k),1.,k)
           corravg_count(idx_oc_cor_chi_Nr,k) = corravg_count(idx_oc_cor_chi_Nr,k) + 1
        else ! oc_micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_Nr',0.0,1.,k)
        endif ! ! oc_micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( oc_micro_correlations(idx_s,idx_ri,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_ri',oc_micro_correlations(idx_s,idx_ri,k),1.,k)
           corravg_count(idx_oc_cor_chi_ri,k) = corravg_count(idx_oc_cor_chi_ri,k) + 1
        else ! oc_micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_ri',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqci,k) /= undef_corr
        if ( oc_micro_correlations(idx_s,idx_Ni,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_Ni',oc_micro_correlations(idx_s,idx_Ni,k),1.,k)
           corravg_count(idx_oc_cor_chi_Ni,k) = corravg_count(idx_oc_cor_chi_Ni,k) + 1
        else ! oc_micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_Ni',0.0,1.,k)
        endif ! oc_micro_correlations(1,inci,k) /= undef_corr
        if ( oc_micro_correlations(idx_s,idx_rs,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_rs',oc_micro_correlations(idx_s,idx_rs,k),1.,k)
           corravg_count(idx_oc_cor_chi_rs,k) = corravg_count(idx_oc_cor_chi_rs,k) + 1
        else ! oc_micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_rs',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqs,k) /= undef_corr
        if ( oc_micro_correlations(idx_s,idx_Ns,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_Ns',oc_micro_correlations(idx_s,idx_Ns,k),1.,k)
           corravg_count(idx_oc_cor_chi_Ns,k) = corravg_count(idx_oc_cor_chi_Ns,k) + 1
        else ! oc_micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_Ns',0.0,1.,k)
        endif ! oc_micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( oc_micro_correlations(idx_s,idx_rg,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_rg',oc_micro_correlations(idx_s,idx_rg,k),1.,k)
           corravg_count(idx_oc_cor_chi_rg,k) = corravg_count(idx_oc_cor_chi_rg,k) + 1
        else ! oc_micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_rg',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqg,k) /= undef_corr
        if ( oc_micro_correlations(idx_s,idx_Ng,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_Ng',oc_micro_correlations(idx_s,idx_Ng,k),1.,k)
           corravg_count(idx_oc_cor_chi_Ng,k) = corravg_count(idx_oc_cor_chi_Ng,k) + 1
        else ! oc_micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_Ng',0.0,1.,k)
        endif ! oc_micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

!-------------
! Vertical Velocity
!-------------
      if(dopredictNc) then
        if ( oc_micro_correlations(idx_w,idx_Nc,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_Nc',oc_micro_correlations(idx_w,idx_Nc,k),1.,k)
           corravg_count(idx_oc_cor_w_Nc,k) = corravg_count(idx_oc_cor_w_Nc,k) + 1
        else ! oc_micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('oc_corr_w_Nc',0.0,1.,k)
        endif ! oc_micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( oc_micro_correlations(idx_w,idx_rr,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_rr',oc_micro_correlations(idx_w,idx_rr,k),1.,k)
           corravg_count(idx_oc_cor_w_rr,k) = corravg_count(idx_oc_cor_w_rr,k) + 1
        else ! oc_micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('oc_corr_w_rr',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqr,k) /= undef_corr
        if ( oc_micro_correlations(idx_w,idx_Nr,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_Nr',oc_micro_correlations(idx_w,idx_Nr,k),1.,k)
           corravg_count(idx_oc_cor_w_nr,k) = corravg_count(idx_oc_cor_w_nr,k) + 1
        else ! oc_micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('oc_corr_w_Nr',0.0,1.,k)
        endif ! ! oc_micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( oc_micro_correlations(idx_w,idx_ri,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_ri',oc_micro_correlations(idx_w,idx_ri,k),1.,k)
           corravg_count(idx_oc_cor_w_ri,k) = corravg_count(idx_oc_cor_w_ri,k) + 1
        else ! oc_micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('oc_corr_w_ri',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqci,k) /= undef_corr
        if ( oc_micro_correlations(idx_w,idx_Ni,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_Ni',oc_micro_correlations(idx_w,idx_Ni,k),1.,k)
           corravg_count(idx_oc_cor_w_Ni,k) = corravg_count(idx_oc_cor_w_Ni,k) + 1
        else ! oc_micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('oc_corr_w_Ni',0.0,1.,k)
        endif ! oc_micro_correlations(1,inci,k) /= undef_corr
        if ( oc_micro_correlations(idx_w,idx_rs,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_rs',oc_micro_correlations(idx_w,idx_rs,k),1.,k)
           corravg_count(idx_oc_cor_w_rs,k) = corravg_count(idx_oc_cor_w_rs,k) + 1
        else ! oc_micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('oc_corr_w_rs',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqs,k) /= undef_corr
        if ( oc_micro_correlations(idx_w,idx_Ns,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_Ns',oc_micro_correlations(idx_w,idx_Ns,k),1.,k)
           corravg_count(idx_oc_cor_w_Ns,k) = corravg_count(idx_oc_cor_w_Ns,k) + 1
        else ! oc_micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('oc_corr_w_Ns',0.0,1.,k)
        endif ! oc_micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( oc_micro_correlations(idx_w,idx_rg,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_rg',oc_micro_correlations(idx_w,idx_rg,k),1.,k)
           corravg_count(idx_oc_cor_w_rg,k) = corravg_count(idx_oc_cor_w_rg,k) + 1
        else ! oc_micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('oc_corr_w_rg',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqg,k) /= undef_corr
        if ( oc_micro_correlations(idx_w,idx_Ng,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_Ng',oc_micro_correlations(idx_w,idx_Ng,k),1.,k)
           corravg_count(idx_oc_cor_w_Ng,k) = corravg_count(idx_oc_cor_w_Ng,k) + 1
        else ! oc_micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('oc_corr_w_Ng',0.0,1.,k)
        endif ! oc_micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

!-------------
! Cloud droplet number concentration
!-------------
      if (dopredictNc) then
     
        if (doprecip) then
           if ( oc_micro_correlations(idx_Nc,idx_rr,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_rr',oc_micro_correlations(idx_Nc,idx_rr,k),1.,k)
              corravg_count(idx_oc_cor_Nc_rr,k) = corravg_count(idx_oc_cor_Nc_rr,k) + 1
           else ! oc_micro_correlations(incl,iqr,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_rr',0.0,1.,k)
           endif ! oc_micro_correlations(incl,iqr,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nc,idx_Nr,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_Nr',oc_micro_correlations(idx_Nc,idx_Nr,k),1.,k)
              corravg_count(idx_oc_cor_Nc_nr,k) = corravg_count(idx_oc_cor_Nc_nr,k) + 1
           else ! oc_micro_correlations(incl,inr,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_Nr',0.0,1.,k)
           endif ! oc_micro_correlations(incl,inr,k) /= undef_corr
        end if
     
        if (doicemicro) then
           if ( oc_micro_correlations(idx_Nc,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_ri',oc_micro_correlations(idx_Nc,idx_ri,k),1.,k) 
              corravg_count(idx_oc_cor_Nc_ri,k) = corravg_count(idx_oc_cor_Nc_ri,k) + 1
           else ! oc_micro_correlations(incl,iqci,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_ri',0.0,1.,k)
           endif ! oc_micro_correlations(incl,iqci,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nc,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_Ni',oc_micro_correlations(idx_Nc,idx_Ni,k),1.,k)
              corravg_count(idx_oc_cor_Nc_Ni,k) = corravg_count(idx_oc_cor_Nc_Ni,k) + 1
           else ! oc_micro_correlations(incl,inci,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_Ni',0.0,1.,k)
           endif ! oc_micro_correlations(incl,inci,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nc,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_rs',oc_micro_correlations(idx_Nc,idx_rs,k),1.,k)
              corravg_count(idx_oc_cor_Nc_rs,k) = corravg_count(idx_oc_cor_Nc_rs,k) + 1
           else ! oc_micro_correlations(incl,iqs,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_rs',0.0,1.,k)
           endif ! oc_micro_correlations(incl,iqs,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nc,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_Ns',oc_micro_correlations(idx_Nc,idx_Ns,k),1.,k)
              corravg_count(idx_oc_cor_Nc_Ns,k) = corravg_count(idx_oc_cor_Nc_Ns,k) + 1
           else ! oc_micro_correlations(incl,ins,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_Ns',0.0,1.,k)
           endif ! oc_micro_correlations(incl,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( oc_micro_correlations(idx_Nc,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_rg',oc_micro_correlations(idx_Nc,idx_rg,k),1.,k)
              corravg_count(idx_oc_cor_Nc_rg,k) = corravg_count(idx_oc_cor_Nc_rg,k) + 1
           else ! oc_micro_correlations(incl,iqg,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_rg',0.0,1.,k)
           endif ! oc_micro_correlations(incl,iqg,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nc,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_Ng',oc_micro_correlations(idx_Nc,idx_Ng,k),1.,k)
              corravg_count(idx_oc_cor_Nc_Ng,k) = corravg_count(idx_oc_cor_Nc_Ng,k) + 1
           else ! oc_micro_correlations(incl,ing,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_Ng',0.0,1.,k)
           endif ! oc_micro_correlations(incl,ing,k) /= undef_corr
        end if   

      end if !dopredictNc
  
!-------------
! Rainwater mixing ratio
!-------------
      if (doprecip) then 
      
         if ( oc_micro_correlations(idx_rr,idx_Nr,k) /= undef_corr ) then
             call hbuf_put_level('oc_corr_rr_Nr',oc_micro_correlations(idx_rr,idx_Nr,k),1.,k)
             corravg_count(idx_oc_cor_rr_nr,k) = corravg_count(idx_oc_cor_rr_nr,k) + 1
          else ! oc_micro_correlations(iqr,inr,k) == undef_corr
             call hbuf_put_level('oc_corr_rr_Nr',0.0,1.,k)
          endif ! oc_micro_correlations(iqr,inr,k) /= undef_corr

        if (doicemicro) then
           if ( oc_micro_correlations(idx_rr,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_rr_ri',oc_micro_correlations(idx_rr,idx_ri,k),1.,k)
              corravg_count(idx_oc_cor_rr_ri,k) = corravg_count(idx_oc_cor_rr_ri,k) + 1
           else ! oc_micro_correlations(iqr,iqci,k) == undef_corr
              call hbuf_put_level('oc_corr_rr_ri',0.0,1.,k)
           endif ! oc_micro_correlations(iqr,iqci,k) /= undef_corr
           if ( oc_micro_correlations(idx_rr,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_rr_Ni',oc_micro_correlations(idx_rr,idx_Ni,k),1.,k)
              corravg_count(idx_oc_cor_rr_Ni,k) = corravg_count(idx_oc_cor_rr_Ni,k) + 1
           else ! oc_micro_correlations(iqr,inci,k) == undef_corr
              call hbuf_put_level('oc_corr_rr_Ni',0.0,1.,k)
           endif ! oc_micro_correlations(iqr,inci,k) /= undef_corr
           if ( oc_micro_correlations(idx_rr,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_rr_rs',oc_micro_correlations(idx_rr,idx_rs,k),1.,k)
              corravg_count(idx_oc_cor_rr_rs,k) = corravg_count(idx_oc_cor_rr_rs,k) + 1
           else ! oc_micro_correlations(iqr,iqs,k) == undef_corr
              call hbuf_put_level('oc_corr_rr_rs',0.0,1.,k)
           endif ! oc_micro_correlations(iqr,iqs,k) /= undef_corr
           if ( oc_micro_correlations(idx_rr,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_rr_Ns',oc_micro_correlations(idx_rr,idx_Ns,k),1.,k)
              corravg_count(idx_oc_cor_rr_Ns,k) = corravg_count(idx_oc_cor_rr_Ns,k) + 1
           else ! oc_micro_correlations(iqr,ins,k) == undef_corr
              call hbuf_put_level('oc_corr_rr_Ns',0.0,1.,k)
           endif ! oc_micro_correlations(iqr,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( oc_micro_correlations(idx_rr,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_rr_rg',oc_micro_correlations(idx_rr,idx_rg,k),1.,k)
              corravg_count(idx_oc_cor_rr_rg,k) = corravg_count(idx_oc_cor_rr_rg,k) + 1
           else ! oc_micro_correlations(iqr,iqg,k) == undef_corr
              call hbuf_put_level('oc_corr_rr_rg',0.0,1.,k)
           endif ! oc_micro_correlations(iqr,iqg,k) /= undef_corr
           if ( oc_micro_correlations(idx_rr,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_rr_Ng',oc_micro_correlations(idx_rr,idx_Ng,k),1.,k)
              corravg_count(idx_oc_cor_rr_Ng,k) = corravg_count(idx_oc_cor_rr_Ng,k) + 1
           else ! oc_micro_correlations(iqr,ing,k) == undef_corr
              call hbuf_put_level('oc_corr_rr_Ng',0.0,1.,k)
           endif ! oc_micro_correlations(iqr,ing,k) /= undef_corr
        end if
!-------------
! Rainwater number concentration
!-------------

        if (doicemicro) then
           if ( oc_micro_correlations(idx_Nr,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nr_ri',oc_micro_correlations(idx_Nr,idx_ri,k),1.,k)
              corravg_count(idx_oc_cor_Nr_ri,k) = corravg_count(idx_oc_cor_Nr_ri,k) + 1
           else ! oc_micro_correlations(inr,iqci,k) == undef_corr
              call hbuf_put_level('oc_corr_Nr_ri',0.0,1.,k)
           endif ! oc_micro_correlations(inr,iqci,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nr,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nr_Ni',oc_micro_correlations(idx_Nr,idx_Ni,k),1.,k)
              corravg_count(idx_oc_cor_Nr_Ni,k) = corravg_count(idx_oc_cor_Nr_Ni,k) + 1
           else ! oc_micro_correlations(inr,inci,k) == undef_corr
              call hbuf_put_level('oc_corr_Nr_Ni',0.0,1.,k)
           endif ! oc_micro_correlations(inr,inci,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nr,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nr_rs',oc_micro_correlations(idx_Nr,idx_rs,k),1.,k)
              corravg_count(idx_oc_cor_Nr_rs,k) = corravg_count(idx_oc_cor_Nr_rs,k) + 1
           else ! oc_micro_correlations(inr,iqs,k) == undef_corr
              call hbuf_put_level('oc_corr_Nr_rs',0.0,1.,k)
           endif ! oc_micro_correlations(inr,iqs,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nr,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nr_Ns',oc_micro_correlations(idx_Nr,idx_Ns,k),1.,k)
              corravg_count(idx_oc_cor_Nr_Ns,k) = corravg_count(idx_oc_cor_Nr_Ns,k) + 1
           else ! oc_micro_correlations(inr,ins,k) == undef_corr
              call hbuf_put_level('oc_corr_Nr_Ns',0.0,1.,k)
           endif ! oc_micro_correlations(inr,ins,k) /= undef_corr
        end if

        if (dograupel) then
        if ( oc_micro_correlations(idx_Nr,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nr_rg',oc_micro_correlations(idx_Nr,idx_rg,k),1.,k)
              corravg_count(idx_oc_cor_Nr_rg,k) = corravg_count(idx_oc_cor_Nr_rg,k) + 1
           else ! oc_micro_correlations(inr,iqg,k) == undef_corr
              call hbuf_put_level('oc_corr_Nr_rg',0.0,1.,k)
           endif ! oc_micro_correlations(inr,iqg,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nr,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nr_Ng',oc_micro_correlations(idx_Nr,idx_Ng,k),1.,k)
              corravg_count(idx_oc_cor_Nr_Ng,k) = corravg_count(idx_oc_cor_Nr_Ng,k) + 1
           else ! oc_micro_correlations(inr,ing,k) == undef_corr
              call hbuf_put_level('oc_corr_Nr_Ng',0.0,1.,k)
           endif ! oc_micro_correlations(inr,ing,k) /= undef_corr
        end if

      end if !doprecip
!-------------
! Cloud-ice mixing ratio
!-------------

      if (doicemicro) then   
         if ( oc_micro_correlations(idx_ri,idx_Ni,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_ri_Ni',oc_micro_correlations(idx_ri,idx_Ni,k),1.,k)
            corravg_count(idx_oc_cor_ri_Ni,k) = corravg_count(idx_oc_cor_ri_Ni,k) + 1
         else ! oc_micro_correlations(iqci,inci,k) == undef_corr
            call hbuf_put_level('oc_corr_ri_Ni',0.0,1.,k)
         endif ! oc_micro_correlations(iqci,inci,k) /= undef_corr
         if ( oc_micro_correlations(idx_ri,idx_rs,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_ri_rs',oc_micro_correlations(idx_ri,idx_rs,k),1.,k)
            corravg_count(idx_oc_cor_ri_rs,k) = corravg_count(idx_oc_cor_ri_rs,k) + 1
         else ! oc_micro_correlations(iqci,iqs,k) == undef_corr
            call hbuf_put_level('oc_corr_ri_rs',0.0,1.,k)
         endif ! oc_micro_correlations(iqci,iqs,k) /= undef_corr
         if ( oc_micro_correlations(idx_ri,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_ri_Ns',oc_micro_correlations(idx_ri,idx_Ns,k),1.,k)
            corravg_count(idx_oc_cor_ri_Ns,k) = corravg_count(idx_oc_cor_ri_Ns,k) + 1
         else ! oc_micro_correlations(iqci,ins,k) == undef_corr
            call hbuf_put_level('oc_corr_ri_Ns',0.0,1.,k)
         endif ! oc_micro_correlations(iqci,ins,k) /= undef_corr
      end if   

      if (dograupel) then
         if ( oc_micro_correlations(idx_ri,idx_rg,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_ri_rg',oc_micro_correlations(idx_ri,idx_rg,k),1.,k)
            corravg_count(idx_oc_cor_ri_rg,k) = corravg_count(idx_oc_cor_ri_rg,k) + 1
         else ! oc_micro_correlations(iqci,iqg,k) == undef_corr
            call hbuf_put_level('oc_corr_ri_rg',0.0,1.,k)
         endif ! oc_micro_correlations(iqci,iqg,k) /= undef_corr
         if ( oc_micro_correlations(idx_ri,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_ri_Ng',oc_micro_correlations(idx_ri,idx_Ng,k),1.,k)
            corravg_count(idx_oc_cor_ri_Ng,k) = corravg_count(idx_oc_cor_ri_Ng,k) + 1
         else ! oc_micro_correlations(iqci,ing,k) == undef_corr
            call hbuf_put_level('oc_corr_ri_Ng',0.0,1.,k)
         endif ! oc_micro_correlations(iqci,ing,k) /= undef_corr
      end if
 
!-------------
! Cloud-ice number concentration
!-------------
      if (doicemicro) then   
         if ( oc_micro_correlations(idx_Ni,idx_rs,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_Ni_rs',oc_micro_correlations(idx_Ni,idx_rs,k),1.,k)
            corravg_count(idx_oc_cor_Ni_rs,k) = corravg_count(idx_oc_cor_Ni_rs,k) + 1
         else ! oc_micro_correlations(inci,iqs,k) == undef_corr
            call hbuf_put_level('oc_corr_Ni_rs',0.0,1.,k)
         endif ! oc_micro_correlations(inci,iqs,k) /= undef_corr
         if ( oc_micro_correlations(idx_Ni,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_Ni_Ns',oc_micro_correlations(idx_Ni,idx_Ns,k),1.,k)
            corravg_count(idx_oc_cor_Ni_Ns,k) = corravg_count(idx_oc_cor_Ni_Ns,k) + 1
         else ! oc_micro_correlations(inci,ins,k) == undef_corr
            call hbuf_put_level('oc_corr_Ni_Ns',0.0,1.,k)
         endif ! oc_micro_correlations(inci,ins,k) /= undef_corr
      end if

      if (dograupel) then     
         if ( oc_micro_correlations(idx_Ni,idx_rg,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_Ni_rg',oc_micro_correlations(idx_Ni,idx_rg,k),1.,k)
            corravg_count(idx_oc_cor_Ni_rg,k) = corravg_count(idx_oc_cor_Ni_rg,k) + 1
         else ! oc_micro_correlations(inci,iqg,k) == undef_corr
            call hbuf_put_level('oc_corr_Ni_rg',0.0,1.,k)
         endif ! oc_micro_correlations(inci,iqg,k) /= undef_corr
         if ( oc_micro_correlations(idx_Ni,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_Ni_Ng',oc_micro_correlations(idx_Ni,idx_Ng,k),1.,k)
            corravg_count(idx_oc_cor_Ni_Ng,k) = corravg_count(idx_oc_cor_Ni_Ng,k) + 1
         else ! oc_micro_correlations(inci,ing,k) == undef_corr
            call hbuf_put_level('oc_corr_Ni_Ng',0.0,1.,k)
         endif ! oc_micro_correlations(inci,ing,k) /= undef_corr
      end if 
!-------------
! Snow mixing ratio
!-------------

      if (doicemicro) then
         if ( oc_micro_correlations(idx_rs,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_rs_Ns',oc_micro_correlations(idx_rs,idx_Ns,k),1.,k)
            corravg_count(idx_oc_cor_rs_Ns,k) = corravg_count(idx_oc_cor_rs_Ns,k) + 1
         else ! oc_micro_correlations(iqs,ins,k) == undef_corr
            call hbuf_put_level('oc_corr_rs_Ns',0.0,1.,k)
         endif ! oc_micro_correlations(iqs,ins,k) /= undef_corr
       
          if (dograupel) then   
             if ( oc_micro_correlations(idx_rs,idx_rg,k) /= undef_corr ) then
                call hbuf_put_level('oc_corr_rs_rg',oc_micro_correlations(idx_rs,idx_rg,k),1.,k)
                corravg_count(idx_oc_cor_rs_rg,k) = corravg_count(idx_oc_cor_rs_rg,k) + 1
             else ! oc_micro_correlations(iqs,iqg,k) == undef_corr
                call hbuf_put_level('oc_corr_rs_rg',0.0,1.,k)
             endif ! oc_micro_correlations(iqs,iqg,k) /= undef_corr
             if ( oc_micro_correlations(idx_rs,idx_Ng,k) /= undef_corr ) then
                call hbuf_put_level('oc_corr_rs_Ng',oc_micro_correlations(idx_rs,idx_Ng,k),1.,k)
                corravg_count(idx_oc_cor_rs_Ng,k) = corravg_count(idx_oc_cor_rs_Ng,k) + 1
             else ! oc_micro_correlations(iqs,ing,k) == undef_corr
                call hbuf_put_level('oc_corr_rs_Ng',0.0,1.,k)
             endif ! oc_micro_correlations(iqs,ing,k) /= undef_corr
!-------------
! Snow number concentration
!-------------
             if ( oc_micro_correlations(idx_Ns,idx_rg,k) /= undef_corr ) then
                call hbuf_put_level('oc_corr_Ns_rg',oc_micro_correlations(idx_Ns,idx_rg,k),1.,k)
                corravg_count(idx_oc_cor_Ns_rg,k) = corravg_count(idx_oc_cor_Ns_rg,k) + 1
             else ! oc_micro_correlations(ins,iqg,k) == undef_corr
                call hbuf_put_level('oc_corr_Ns_rg',0.0,1.,k)
             endif ! oc_micro_correlations(ins,iqg,k) /= undef_corr
             if ( oc_micro_correlations(idx_Ns,idx_Ng,k) /= undef_corr ) then
                call hbuf_put_level('oc_corr_Ns_Ng',oc_micro_correlations(idx_Ns,idx_Ng,k),1.,k)
                corravg_count(idx_oc_cor_Ns_Ng,k) = corravg_count(idx_oc_cor_Ns_Ng,k) + 1
             else ! oc_micro_correlations(ins,ing,k) == undef_corr
                call hbuf_put_level('oc_corr_Ns_Ng',0.0,1.,k)
             endif ! oc_micro_correlations(ins,ing,k) /= undef_corr
          end if
      end if

!-------------
! Graupel mixing ratio / number concentration
!-------------
      if (dograupel) then      
        if ( oc_micro_correlations(idx_rg,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_rg_Ng',oc_micro_correlations(idx_rg,idx_Ng,k),1.,k)
            corravg_count(idx_oc_cor_rg_Ng,k) = corravg_count(idx_oc_cor_rg_Ng,k) + 1
         else ! oc_micro_correlations(iqg,ing,k) == undef_corr
            call hbuf_put_level('oc_corr_rg_Ng',0.0,1.,k)
         endif ! oc_micro_correlations(iqg,ing,k) /= undef_corr
      endif

    end do ! do k = 1, nzm, 1

!---------------------------------
! Covariances
!--------------------------------

!---------------------------------
! chi
!--------------------------------

call hbuf_put('covarnce_chi_w',micro_covarnce(idx_s,idx_w,:),1.)

if(dopredictNc) then
  call hbuf_put('covarnce_chi_Nc',micro_covarnce(idx_s,idx_Nc,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('covarnce_chi_rr',micro_covarnce(idx_s,idx_rr,:),1.)
  call hbuf_put('covarnce_chi_Nr',micro_covarnce(idx_s,idx_Nr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('covarnce_chi_ri',micro_covarnce(idx_s,idx_ri,:),1.)
  call hbuf_put('covarnce_chi_Ni',micro_covarnce(idx_s,idx_Ni,:),1.)
  call hbuf_put('covarnce_chi_rs',micro_covarnce(idx_s,idx_rs,:),1.)
  call hbuf_put('covarnce_chi_Ns',micro_covarnce(idx_s,idx_Ns,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('covarnce_chi_rg',micro_covarnce(idx_s,idx_rg,:),1.)
  call hbuf_put('covarnce_chi_Ng',micro_covarnce(idx_s,idx_Ng,:),1.)
end if ! dograupel


!---------------------------------
! Vertical velocity
!--------------------------------
if(dopredictNc) then
  call hbuf_put('covarnce_w_Nc',micro_covarnce(idx_w,idx_Nc,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('covarnce_w_rr',micro_covarnce(idx_w,idx_rr,:),1.)
  call hbuf_put('covarnce_w_Nr',micro_covarnce(idx_w,idx_Nr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('covarnce_w_ri',micro_covarnce(idx_w,idx_ri,:),1.)
  call hbuf_put('covarnce_w_Ni',micro_covarnce(idx_w,idx_Ni,:),1.)
  call hbuf_put('covarnce_w_rs',micro_covarnce(idx_w,idx_rs,:),1.)
  call hbuf_put('covarnce_w_Ns',micro_covarnce(idx_w,idx_Ns,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('covarnce_w_rg',micro_covarnce(idx_w,idx_rg,:),1.)
  call hbuf_put('covarnce_w_Ng',micro_covarnce(idx_w,idx_Ng,:),1.)
end if ! dograupel

!---------------------------------
! Cloud liquid number concentration
!--------------------------------

if(dopredictNc) then  
  if (doprecip) then
    call hbuf_put('covarnce_Nc_rr',micro_covarnce(idx_Nc,idx_rr,:),1.)
    call hbuf_put('covarnce_Nc_Nr',micro_covarnce(idx_Nc,idx_Nr,:),1.)
  end if ! doprecip
     
  if (doicemicro) then
    call hbuf_put('covarnce_Nc_ri',micro_covarnce(idx_Nc,idx_ri,:),1.) 
    call hbuf_put('covarnce_Nc_Ni',micro_covarnce(idx_Nc,idx_Ni,:),1.)
    call hbuf_put('covarnce_Nc_rs',micro_covarnce(idx_Nc,idx_rs,:),1.)
    call hbuf_put('covarnce_Nc_Ns',micro_covarnce(idx_Nc,idx_Ns,:),1.)
  end if ! doicemicro    

  if (dograupel) then
    call hbuf_put('covarnce_Nc_rg',micro_covarnce(idx_Nc,idx_rg,:),1.)
    call hbuf_put('covarnce_Nc_Ng',micro_covarnce(idx_Nc,idx_Ng,:),1.)
  end if ! dograupel

end if !dopredictNc

!---------------------------------
! Rainwater mixing ratio
!--------------------------------
  
if (doprecip) then 
  call hbuf_put('covarnce_rr_Nr',micro_covarnce(idx_rr,idx_Nr,:),1.)

  if (doicemicro) then
    call hbuf_put('covarnce_rr_ri',micro_covarnce(idx_rr,idx_ri,:),1.)
    call hbuf_put('covarnce_rr_Ni',micro_covarnce(idx_rr,idx_Ni,:),1.)
    call hbuf_put('covarnce_rr_rs',micro_covarnce(idx_rr,idx_rs,:),1.)
    call hbuf_put('covarnce_rr_Ns',micro_covarnce(idx_rr,idx_Ns,:),1.)

    if (dograupel) then
      call hbuf_put('covarnce_rr_rg',micro_covarnce(idx_rr,idx_rg,:),1.)
      call hbuf_put('covarnce_rr_Ng',micro_covarnce(idx_rr,idx_Ng,:),1.)
    end if !dograupel

!---------------------------------
! Rainwater number concentration
!--------------------------------

    call hbuf_put('covarnce_Nr_ri',micro_covarnce(idx_Nr,idx_ri,:),1.)
    call hbuf_put('covarnce_Nr_Ni',micro_covarnce(idx_Nr,idx_Ni,:),1.)
    call hbuf_put('covarnce_Nr_rs',micro_covarnce(idx_Nr,idx_rs,:),1.)
    call hbuf_put('covarnce_Nr_Ns',micro_covarnce(idx_Nr,idx_Ns,:),1.)

    if (dograupel) then
      call hbuf_put('covarnce_Nr_rg',micro_covarnce(idx_Nr,idx_rg,:),1.)
      call hbuf_put('covarnce_Nr_Ng',micro_covarnce(idx_Nr,idx_Ng,:),1.)
    end if !dograupel
  
  end if !doicemicro

end if !doprecip

!---------------------------------
! Cloud-ice mixing ratio
!--------------------------------

if (doicemicro) then   
  call hbuf_put('covarnce_ri_Ni',micro_covarnce(idx_ri,idx_Ni,:),1.)
  call hbuf_put('covarnce_ri_rs',micro_covarnce(idx_ri,idx_rs,:),1.)
  call hbuf_put('covarnce_ri_Ns',micro_covarnce(idx_ri,idx_Ns,:),1.)

  if (dograupel) then
    call hbuf_put('covarnce_ri_rg',micro_covarnce(idx_ri,idx_rg,:),1.)
    call hbuf_put('covarnce_ri_Ng',micro_covarnce(idx_ri,idx_Ng,:),1.)
  end if !dograupel

!---------------------------------
! Cloud-ice number concentration
!--------------------------------
 
  call hbuf_put('covarnce_Ni_rs',micro_covarnce(idx_Ni,idx_rs,:),1.)
  call hbuf_put('covarnce_Ni_Ns',micro_covarnce(idx_Ni,idx_Ns,:),1.)
      
  if (dograupel) then     
    call hbuf_put('covarnce_Ni_rg',micro_covarnce(idx_Ni,idx_rg,:),1.)
    call hbuf_put('covarnce_Ni_Ng',micro_covarnce(idx_Ni,idx_Ng,:),1.)
  end if !dograupel

!---------------------------------
! Snow mixing ratio
!--------------------------------

  call hbuf_put('covarnce_rs_Ns',micro_covarnce(idx_rs,idx_Ns,:),1.)
       
  if (dograupel) then   
    call hbuf_put('covarnce_rs_rg',micro_covarnce(idx_rs,idx_rg,:),1.)
    call hbuf_put('covarnce_rs_Ng',micro_covarnce(idx_rs,idx_Ng,:),1.)

!---------------------------------
! Snow number concentration
!--------------------------------
    call hbuf_put('covarnce_Ns_rg',micro_covarnce(idx_Ns,idx_rg,:),1.)
    call hbuf_put('covarnce_Ns_Ng',micro_covarnce(idx_Ns,idx_Ng,:),1.)

!---------------------------------
! Graupel mixing ratio / number concentration
!--------------------------------
    call hbuf_put('covarnce_rg_Ng',micro_covarnce(idx_rg,idx_Ng,:),1.)
  endif !dograupel
endif !doicemicro


!---------------------------------
! In-cloud covariances
!--------------------------------

!---------------------------------
! chi
!--------------------------------

call hbuf_put('ic_covarnce_chi_w',ic_micro_covarnce(idx_s,idx_w,:),1.)

if(dopredictNc) then
  call hbuf_put('ic_covarnce_chi_Nc',ic_micro_covarnce(idx_s,idx_Nc,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('ic_covarnce_chi_rr',ic_micro_covarnce(idx_s,idx_rr,:),1.)
  call hbuf_put('ic_covarnce_chi_Nr',ic_micro_covarnce(idx_s,idx_Nr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('ic_covarnce_chi_ri',ic_micro_covarnce(idx_s,idx_ri,:),1.)
  call hbuf_put('ic_covarnce_chi_Ni',ic_micro_covarnce(idx_s,idx_Ni,:),1.)
  call hbuf_put('ic_covarnce_chi_rs',ic_micro_covarnce(idx_s,idx_rs,:),1.)
  call hbuf_put('ic_covarnce_chi_Ns',ic_micro_covarnce(idx_s,idx_Ns,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('ic_covarnce_chi_rg',ic_micro_covarnce(idx_s,idx_rg,:),1.)
  call hbuf_put('ic_covarnce_chi_Ng',ic_micro_covarnce(idx_s,idx_Ng,:),1.)
end if ! dograupel


!---------------------------------
! Vertical velocity
!--------------------------------
if(dopredictNc) then
  call hbuf_put('ic_covarnce_w_Nc',ic_micro_covarnce(idx_w,idx_Nc,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('ic_covarnce_w_rr',ic_micro_covarnce(idx_w,idx_rr,:),1.)
  call hbuf_put('ic_covarnce_w_Nr',ic_micro_covarnce(idx_w,idx_Nr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('ic_covarnce_w_ri',ic_micro_covarnce(idx_w,idx_ri,:),1.)
  call hbuf_put('ic_covarnce_w_Ni',ic_micro_covarnce(idx_w,idx_Ni,:),1.)
  call hbuf_put('ic_covarnce_w_rs',ic_micro_covarnce(idx_w,idx_rs,:),1.)
  call hbuf_put('ic_covarnce_w_Ns',ic_micro_covarnce(idx_w,idx_Ns,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('ic_covarnce_w_rg',ic_micro_covarnce(idx_w,idx_rg,:),1.)
  call hbuf_put('ic_covarnce_w_Ng',ic_micro_covarnce(idx_w,idx_Ng,:),1.)
end if ! dograupel

!---------------------------------
! Cloud liquid number concentration
!--------------------------------

if(dopredictNc) then  
  if (doprecip) then
    call hbuf_put('ic_covarnce_Nc_rr',ic_micro_covarnce(idx_Nc,idx_rr,:),1.)
    call hbuf_put('ic_covarnce_Nc_Nr',ic_micro_covarnce(idx_Nc,idx_Nr,:),1.)
  end if ! doprecip
     
  if (doicemicro) then
    call hbuf_put('ic_covarnce_Nc_ri',ic_micro_covarnce(idx_Nc,idx_ri,:),1.) 
    call hbuf_put('ic_covarnce_Nc_Ni',ic_micro_covarnce(idx_Nc,idx_Ni,:),1.)
    call hbuf_put('ic_covarnce_Nc_rs',ic_micro_covarnce(idx_Nc,idx_rs,:),1.)
    call hbuf_put('ic_covarnce_Nc_Ns',ic_micro_covarnce(idx_Nc,idx_Ns,:),1.)
  end if ! doicemicro    

  if (dograupel) then
    call hbuf_put('ic_covarnce_Nc_rg',ic_micro_covarnce(idx_Nc,idx_rg,:),1.)
    call hbuf_put('ic_covarnce_Nc_Ng',ic_micro_covarnce(idx_Nc,idx_Ng,:),1.)
  end if ! dograupel

end if !dopredictNc

!---------------------------------
! Rainwater mixing ratio
!--------------------------------
  
if (doprecip) then 
  call hbuf_put('ic_covarnce_rr_Nr',ic_micro_covarnce(idx_rr,idx_Nr,:),1.)

  if (doicemicro) then
    call hbuf_put('ic_covarnce_rr_ri',ic_micro_covarnce(idx_rr,idx_ri,:),1.)
    call hbuf_put('ic_covarnce_rr_Ni',ic_micro_covarnce(idx_rr,idx_Ni,:),1.)
    call hbuf_put('ic_covarnce_rr_rs',ic_micro_covarnce(idx_rr,idx_rs,:),1.)
    call hbuf_put('ic_covarnce_rr_Ns',ic_micro_covarnce(idx_rr,idx_Ns,:),1.)

    if (dograupel) then
      call hbuf_put('ic_covarnce_rr_rg',ic_micro_covarnce(idx_rr,idx_rg,:),1.)
      call hbuf_put('ic_covarnce_rr_Ng',ic_micro_covarnce(idx_rr,idx_Ng,:),1.)
    end if !dograupel

!---------------------------------
! Rainwater number concentration
!--------------------------------

    call hbuf_put('ic_covarnce_Nr_ri',ic_micro_covarnce(idx_Nr,idx_ri,:),1.)
    call hbuf_put('ic_covarnce_Nr_Ni',ic_micro_covarnce(idx_Nr,idx_Ni,:),1.)
    call hbuf_put('ic_covarnce_Nr_rs',ic_micro_covarnce(idx_Nr,idx_rs,:),1.)
    call hbuf_put('ic_covarnce_Nr_Ns',ic_micro_covarnce(idx_Nr,idx_Ns,:),1.)

    if (dograupel) then
      call hbuf_put('ic_covarnce_Nr_rg',ic_micro_covarnce(idx_Nr,idx_rg,:),1.)
      call hbuf_put('ic_covarnce_Nr_Ng',ic_micro_covarnce(idx_Nr,idx_Ng,:),1.)
    end if !dograupel
  
  end if !doicemicro

end if !doprecip

!---------------------------------
! Cloud-ice mixing ratio
!--------------------------------

if (doicemicro) then   
  call hbuf_put('ic_covarnce_ri_Ni',ic_micro_covarnce(idx_ri,idx_Ni,:),1.)
  call hbuf_put('ic_covarnce_ri_rs',ic_micro_covarnce(idx_ri,idx_rs,:),1.)
  call hbuf_put('ic_covarnce_ri_Ns',ic_micro_covarnce(idx_ri,idx_Ns,:),1.)

  if (dograupel) then
    call hbuf_put('ic_covarnce_ri_rg',ic_micro_covarnce(idx_ri,idx_rg,:),1.)
    call hbuf_put('ic_covarnce_ri_Ng',ic_micro_covarnce(idx_ri,idx_Ng,:),1.)
  end if !dograupel

!---------------------------------
! Cloud-ice number concentration
!--------------------------------
 
  call hbuf_put('ic_covarnce_Ni_rs',ic_micro_covarnce(idx_Ni,idx_rs,:),1.)
  call hbuf_put('ic_covarnce_Ni_Ns',ic_micro_covarnce(idx_Ni,idx_Ns,:),1.)
      
  if (dograupel) then     
    call hbuf_put('ic_covarnce_Ni_rg',ic_micro_covarnce(idx_Ni,idx_rg,:),1.)
    call hbuf_put('ic_covarnce_Ni_Ng',ic_micro_covarnce(idx_Ni,idx_Ng,:),1.)
  end if !dograupel

!---------------------------------
! Snow mixing ratio
!--------------------------------

  call hbuf_put('ic_covarnce_rs_Ns',ic_micro_covarnce(idx_rs,idx_Ns,:),1.)
       
  if (dograupel) then   
    call hbuf_put('ic_covarnce_rs_rg',ic_micro_covarnce(idx_rs,idx_rg,:),1.)
    call hbuf_put('ic_covarnce_rs_Ng',ic_micro_covarnce(idx_rs,idx_Ng,:),1.)

!---------------------------------
! Snow number concentration
!--------------------------------
    call hbuf_put('ic_covarnce_Ns_rg',ic_micro_covarnce(idx_Ns,idx_rg,:),1.)
    call hbuf_put('ic_covarnce_Ns_Ng',ic_micro_covarnce(idx_Ns,idx_Ng,:),1.)

!---------------------------------
! Graupel mixing ratio / number concentration
!--------------------------------
    call hbuf_put('ic_covarnce_rg_Ng',ic_micro_covarnce(idx_rg,idx_Ng,:),1.)
  endif !dograupel
endif !doicemicro

!---------------------------------
! Out of cloud covariances
!--------------------------------

!---------------------------------
! chi
!--------------------------------

call hbuf_put('oc_covarnce_chi_w',oc_micro_covarnce(idx_s,idx_w,:),1.)

if(dopredictNc) then
  call hbuf_put('oc_covarnce_chi_Nc',oc_micro_covarnce(idx_s,idx_Nc,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('oc_covarnce_chi_rr',oc_micro_covarnce(idx_s,idx_rr,:),1.)
  call hbuf_put('oc_covarnce_chi_Nr',oc_micro_covarnce(idx_s,idx_Nr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('oc_covarnce_chi_ri',oc_micro_covarnce(idx_s,idx_ri,:),1.)
  call hbuf_put('oc_covarnce_chi_Ni',oc_micro_covarnce(idx_s,idx_Ni,:),1.)
  call hbuf_put('oc_covarnce_chi_rs',oc_micro_covarnce(idx_s,idx_rs,:),1.)
  call hbuf_put('oc_covarnce_chi_Ns',oc_micro_covarnce(idx_s,idx_Ns,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('oc_covarnce_chi_rg',oc_micro_covarnce(idx_s,idx_rg,:),1.)
  call hbuf_put('oc_covarnce_chi_Ng',oc_micro_covarnce(idx_s,idx_Ng,:),1.)
end if ! dograupel


!---------------------------------
! Vertical velocity
!--------------------------------
if(dopredictNc) then
  call hbuf_put('oc_covarnce_w_Nc',oc_micro_covarnce(idx_w,idx_Nc,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('oc_covarnce_w_rr',oc_micro_covarnce(idx_w,idx_rr,:),1.)
  call hbuf_put('oc_covarnce_w_Nr',oc_micro_covarnce(idx_w,idx_Nr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('oc_covarnce_w_ri',oc_micro_covarnce(idx_w,idx_ri,:),1.)
  call hbuf_put('oc_covarnce_w_Ni',oc_micro_covarnce(idx_w,idx_Ni,:),1.)
  call hbuf_put('oc_covarnce_w_rs',oc_micro_covarnce(idx_w,idx_rs,:),1.)
  call hbuf_put('oc_covarnce_w_Ns',oc_micro_covarnce(idx_w,idx_Ns,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('oc_covarnce_w_rg',oc_micro_covarnce(idx_w,idx_rg,:),1.)
  call hbuf_put('oc_covarnce_w_Ng',oc_micro_covarnce(idx_w,idx_Ng,:),1.)
end if ! dograupel

!---------------------------------
! Cloud liquid number concentration
!--------------------------------

if(dopredictNc) then  
  if (doprecip) then
    call hbuf_put('oc_covarnce_Nc_rr',oc_micro_covarnce(idx_Nc,idx_rr,:),1.)
    call hbuf_put('oc_covarnce_Nc_Nr',oc_micro_covarnce(idx_Nc,idx_Nr,:),1.)
  end if ! doprecip
     
  if (doicemicro) then
    call hbuf_put('oc_covarnce_Nc_ri',oc_micro_covarnce(idx_Nc,idx_ri,:),1.) 
    call hbuf_put('oc_covarnce_Nc_Ni',oc_micro_covarnce(idx_Nc,idx_Ni,:),1.)
    call hbuf_put('oc_covarnce_Nc_rs',oc_micro_covarnce(idx_Nc,idx_rs,:),1.)
    call hbuf_put('oc_covarnce_Nc_Ns',oc_micro_covarnce(idx_Nc,idx_Ns,:),1.)
  end if ! doicemicro    

  if (dograupel) then
    call hbuf_put('oc_covarnce_Nc_rg',oc_micro_covarnce(idx_Nc,idx_rg,:),1.)
    call hbuf_put('oc_covarnce_Nc_Ng',oc_micro_covarnce(idx_Nc,idx_Ng,:),1.)
  end if ! dograupel

end if !dopredictNc

!---------------------------------
! Rainwater mixing ratio
!--------------------------------
  
if (doprecip) then 
  call hbuf_put('oc_covarnce_rr_Nr',oc_micro_covarnce(idx_rr,idx_Nr,:),1.)

  if (doicemicro) then
    call hbuf_put('oc_covarnce_rr_ri',oc_micro_covarnce(idx_rr,idx_ri,:),1.)
    call hbuf_put('oc_covarnce_rr_Ni',oc_micro_covarnce(idx_rr,idx_Ni,:),1.)
    call hbuf_put('oc_covarnce_rr_rs',oc_micro_covarnce(idx_rr,idx_rs,:),1.)
    call hbuf_put('oc_covarnce_rr_Ns',oc_micro_covarnce(idx_rr,idx_Ns,:),1.)

    if (dograupel) then
      call hbuf_put('oc_covarnce_rr_rg',oc_micro_covarnce(idx_rr,idx_rg,:),1.)
      call hbuf_put('oc_covarnce_rr_Ng',oc_micro_covarnce(idx_rr,idx_Ng,:),1.)
    end if !dograupel

!---------------------------------
! Rainwater number concentration
!--------------------------------

    call hbuf_put('oc_covarnce_Nr_ri',oc_micro_covarnce(idx_Nr,idx_ri,:),1.)
    call hbuf_put('oc_covarnce_Nr_Ni',oc_micro_covarnce(idx_Nr,idx_Ni,:),1.)
    call hbuf_put('oc_covarnce_Nr_rs',oc_micro_covarnce(idx_Nr,idx_rs,:),1.)
    call hbuf_put('oc_covarnce_Nr_Ns',oc_micro_covarnce(idx_Nr,idx_Ns,:),1.)

    if (dograupel) then
      call hbuf_put('oc_covarnce_Nr_rg',oc_micro_covarnce(idx_Nr,idx_rg,:),1.)
      call hbuf_put('oc_covarnce_Nr_Ng',oc_micro_covarnce(idx_Nr,idx_Ng,:),1.)
    end if !dograupel
  
  end if !doicemicro

end if !doprecip

!---------------------------------
! Cloud-ice mixing ratio
!--------------------------------

if (doicemicro) then   
  call hbuf_put('oc_covarnce_ri_Ni',oc_micro_covarnce(idx_ri,idx_Ni,:),1.)
  call hbuf_put('oc_covarnce_ri_rs',oc_micro_covarnce(idx_ri,idx_rs,:),1.)
  call hbuf_put('oc_covarnce_ri_Ns',oc_micro_covarnce(idx_ri,idx_Ns,:),1.)

  if (dograupel) then
    call hbuf_put('oc_covarnce_ri_rg',oc_micro_covarnce(idx_ri,idx_rg,:),1.)
    call hbuf_put('oc_covarnce_ri_Ng',oc_micro_covarnce(idx_ri,idx_Ng,:),1.)
  end if !dograupel

!---------------------------------
! Cloud-ice number concentration
!--------------------------------
 
  call hbuf_put('oc_covarnce_Ni_rs',oc_micro_covarnce(idx_Ni,idx_rs,:),1.)
  call hbuf_put('oc_covarnce_Ni_Ns',oc_micro_covarnce(idx_Ni,idx_Ns,:),1.)
      
  if (dograupel) then     
    call hbuf_put('oc_covarnce_Ni_rg',oc_micro_covarnce(idx_Ni,idx_rg,:),1.)
    call hbuf_put('oc_covarnce_Ni_Ng',oc_micro_covarnce(idx_Ni,idx_Ng,:),1.)
  end if !dograupel

!---------------------------------
! Snow mixing ratio
!--------------------------------

  call hbuf_put('oc_covarnce_rs_Ns',oc_micro_covarnce(idx_rs,idx_Ns,:),1.)
       
  if (dograupel) then   
    call hbuf_put('oc_covarnce_rs_rg',oc_micro_covarnce(idx_rs,idx_rg,:),1.)
    call hbuf_put('oc_covarnce_rs_Ng',oc_micro_covarnce(idx_rs,idx_Ng,:),1.)

!---------------------------------
! Snow number concentration
!--------------------------------
    call hbuf_put('oc_covarnce_Ns_rg',oc_micro_covarnce(idx_Ns,idx_rg,:),1.)
    call hbuf_put('oc_covarnce_Ns_Ng',oc_micro_covarnce(idx_Ns,idx_Ng,:),1.)

!---------------------------------
! Graupel mixing ratio / number concentration
!--------------------------------
    call hbuf_put('oc_covarnce_rg_Ng',oc_micro_covarnce(idx_rg,idx_Ng,:),1.)
  endif !dograupel
endif !doicemicro
