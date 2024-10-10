   

!----------------------------
! In-cloud  Correlations 
!----------------------------

!----------------------------
! chi(S_mellor) 
!----------------------------
  name = 'corr_chi_w'
  longname = '  Correlation of chi and vertical velocity'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_chi_w = count

if(dopredictNc) then
  name = 'corr_chi_Nc'
  longname = '  Correlation of chi and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_chi_Nc = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'corr_chi_rr'
    longname = '  Correlation of chi and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_chi_rr = count
    
    name = 'corr_chi_Nr'
    longname = '  Correlation of chi and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_chi_Nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'corr_chi_ri'
    longname = '  Correlation of chi and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_chi_ri = count
    
    name = 'corr_chi_Ni'
    longname = '  Correlation of chi and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_chi_Ni = count
    
    name = 'corr_chi_rs'
    longname = '  Correlation of chi and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_chi_rs = count
    
    name = 'corr_chi_Ns'
    longname = '  Correlation of chi and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_chi_Ns = count
    

    if (dograupel) then
    
      name = 'corr_chi_rg'
      longname = '  Correlation of chi and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_chi_rg = count
    
      name = 'corr_chi_Ng'
      longname = '  Correlation of chi and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_chi_Ng = count
    
    end if ! dograupel
  endif ! doicemicro
!----------------------------
! Vertical Velocity 
!----------------------------

if(dopredictNc) then
  name = 'corr_w_Nc'
  longname = '  Correlation of vertical velocity and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_w_Nc = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'corr_w_rr'
    longname = '  Correlation of vertical velocity and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_rr = count
    
    name = 'corr_w_Nr'
    longname = '  Correlation of vertical velocity and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_Nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'corr_w_ri'
    longname = '  Correlation of vertical velocity and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_ri = count
    
    name = 'corr_w_Ni'
    longname = '  Correlation of vertical velocity and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_Ni = count
    
    name = 'corr_w_rs'
    longname = '  Correlation of vertical velocity and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_rs = count
    
    name = 'corr_w_Ns'
    longname = '  Correlation of vertical velocity and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_Ns = count
    

    if (dograupel) then
    
      name = 'corr_w_rg'
      longname = '  Correlation of vertical velocity and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_w_rg = count
    
      name = 'corr_w_Ng'
      longname = '  Correlation of vertical velocity and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_w_Ng = count
    
    end if ! dograupel
  end if ! doicemicro

!----------------------------
! Cloud droplet number concentration 
!----------------------------
if(dopredictNc) then
  if (doprecip) then
        
    name = 'corr_Nc_rr'
    longname = '  Correlation of cloud droplet conc. and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nc_rr = count
        
    name = 'corr_Nc_Nr'
    longname = '  Correlation of cloud droplet conc. and rain drop conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nc_Nr = count
        
  end if ! doprecip
    
  if (doicemicro) then
        
    name = 'corr_Nc_ri'
    longname = '  Correlation of cloud droplet conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nc_ri = count
        
    name = 'corr_Nc_Ni'
    longname = '  Correlation of cloud droplet conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nc_Ni = count
       
    name = 'corr_Nc_rs'
    longname = '  Correlation of cloud droplet conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nc_rs = count
      
    name = 'corr_Nc_Ns'
    longname = '  Correlation of cloud droplet conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nc_Ns = count
        
    if (dograupel) then
        
      name = 'corr_Nc_rg'
      longname = '  Correlation of cloud droplet conc. and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_Nc_rg = count
        
      name = 'corr_Nc_Ng'
      longname = '  Correlation of cloud droplet conc. and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_Nc_Ng = count
        
    end if ! dograupel
  end if ! doicemicro 
end if ! dopredictNc

!----------------------------
! Rainwater mixing ratio 
!----------------------------
if (doprecip) then

    name = 'corr_rr_Nr'
    longname = '  Correlation of rain water mixing ratio and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rr_Nr = count
  
  if(doicemicro) then

    name = 'corr_rr_ri'
    longname = '  Correlation of rain water mixing ratio and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rr_ri = count
        
    name = 'corr_rr_Ni'
    longname = '  Correlation of rain water mixing ratio and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rr_Ni = count
        
    name = 'corr_rr_rs'
    longname = '  Correlation of rain water mixing ratio and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rr_rs = count
        
    name = 'corr_rr_Ns'
    longname = '  Correlation of rain water mixing ratio and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rr_Ns = count
        

    if(dograupel) then

      name = 'corr_rr_rg'
      longname = '  Correlation of rain water mixing ratio and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_rr_rg = count
        
      name = 'corr_rr_Ng'
      longname = '  Correlation of rain water mixing ratio and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_rr_Ng = count
        
    end if !dograupel

!----------------------------
! Rainwater number concentration 
!----------------------------
        
    name = 'corr_Nr_ri'
    longname = '  Correlation of rain drop conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nr_ri = count
        
    name = 'corr_Nr_Ni'
    longname = '  Correlation of rain drop conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nr_Ni = count
        
    name = 'corr_Nr_rs'
    longname = '  Correlation of rain drop conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nr_rs = count
        
    name = 'corr_Nr_Ns'
    longname = '  Correlation of rain drop conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nr_Ns = count
        
      if(dograupel) then

        name = 'corr_Nr_rg'
        longname = '  Correlation of rain drop conc. and graupel mixing ratio'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_cor_Nr_rg = count
        
        name = 'corr_Nr_Ng'
        longname = '  Correlation of rain drop conc. and graupel conc.'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_cor_Nr_Ng = count
        
      end if !dograupel
    end if !doicemicro
end if !doprecip

!----------------------------
! Cloud-ice mixing ratio 
!----------------------------
if (doicemicro) then
    
  name = 'corr_ri_Ni'
  longname = '  Correlation of cloud ice mixing ratio and cloud ice conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_ri_Ni = count
    
  name = 'corr_ri_rs'
  longname = '  Correlation of cloud ice mixing ratio and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_ri_rs = count
    
  name = 'corr_ri_Ns'
  longname = '  Correlation of cloud ice mixing ratio and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_ri_Ns = count
    
  if (dograupel) then

    name = 'corr_ri_rg'
    longname = '  Correlation of cloud ice mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_ri_rg = count
        
    name = 'corr_ri_Ng'
    longname = '  Correlation of cloud ice mixing ratio and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_ri_Ng = count
        
  end if !dograupel

!----------------------------
! Cloud-ice number concentration 
!----------------------------
    
  name = 'corr_Ni_rs'
  longname = '  Correlation of cloud ice concentration and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_Ni_rs = count
    
  name = 'corr_Ni_Ns'
  longname = '  Correlation of cloud ice concentration and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_Ni_Ns = count
    
  if (dograupel) then

    name = 'corr_Ni_rg'
    longname = '  Correlation of cloud ice conc. and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Ni_rg = count
        
    name = 'corr_Ni_Ng'
    longname = '  Correlation of cloud ice conc. and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Ni_Ng = count
        
  end if !dograupel

!----------------------------
! Snow mixing ratio 
!----------------------------
  
  name = 'corr_rs_Ns'
  longname = '  Correlation of snow mixing ratio and snowflake concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_rs_Ns = count
    
  if (dograupel) then

    name = 'corr_rs_rg'
    longname = '  Correlation of snow mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rs_rg = count
        
    name = 'corr_rs_Ng'
    longname = '  Correlation of snow mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rs_Ng = count
        
!----------------------------
! Snow number concentration 
!----------------------------
    name = 'corr_Ns_rg'
    longname = '  Correlation of snowflake concentration and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Ns_rg = count
        
    name = 'corr_Ns_Ng'
    longname = '  Correlation of snowflake concentration and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Ns_Ng = count

!----------------------------
! Graupel mixing ratio / number concentration 
!----------------------------
        
    name = 'corr_rg_Ng'
    longname = '  Correlation of graupel mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rg_Ng = count
        
  end if !dograupel
end if !doicemicro



!---------------------------
! In-cloud ic_correlations 
!----------------------------

!----------------------------
! chi(S_mellor) 
!----------------------------
  name = 'ic_corr_chi_w'
  longname = 'In-Cloud Correlation of chi and vertical velocity'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_chi_w = count

if(dopredictNc) then
  name = 'ic_corr_chi_Nc'
  longname = 'In-Cloud Correlation of chi and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_chi_Nc = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'ic_corr_chi_rr'
    longname = 'In-Cloud Correlation of chi and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_chi_rr = count
    
    name = 'ic_corr_chi_Nr'
    longname = 'In-Cloud Correlation of chi and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_chi_Nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'ic_corr_chi_ri'
    longname = 'In-Cloud Correlation of chi and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_chi_ri = count
    
    name = 'ic_corr_chi_Ni'
    longname = 'In-Cloud Correlation of chi and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_chi_Ni = count
    
    name = 'ic_corr_chi_rs'
    longname = 'In-Cloud Correlation of chi and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_chi_rs = count
    
    name = 'ic_corr_chi_Ns'
    longname = 'In-Cloud Correlation of chi and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_chi_Ns = count
    

    if (dograupel) then
    
      name = 'ic_corr_chi_rg'
      longname = 'In-Cloud Correlation of chi and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_chi_rg = count
    
      name = 'ic_corr_chi_Ng'
      longname = 'In-Cloud Correlation of chi and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_chi_Ng = count
    
    end if ! dograupel
  endif ! doicemicro
!----------------------------
! Vertical Velocity 
!----------------------------

if(dopredictNc) then
  name = 'ic_corr_w_Nc'
  longname = 'In-Cloud Correlation of vertical velocity and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_w_Nc = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'ic_corr_w_rr'
    longname = 'In-Cloud Correlation of vertical velocity and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_w_rr = count
    
    name = 'ic_corr_w_Nr'
    longname = 'In-Cloud Correlation of vertical velocity and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_w_Nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'ic_corr_w_ri'
    longname = 'In-Cloud Correlation of vertical velocity and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_w_ri = count
    
    name = 'ic_corr_w_Ni'
    longname = 'In-Cloud Correlation of vertical velocity and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_w_Ni = count
    
    name = 'ic_corr_w_rs'
    longname = 'In-Cloud Correlation of vertical velocity and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_w_rs = count
    
    name = 'ic_corr_w_Ns'
    longname = 'In-Cloud Correlation of vertical velocity and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_w_Ns = count
    

    if (dograupel) then
    
      name = 'ic_corr_w_rg'
      longname = 'In-Cloud Correlation of vertical velocity and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_w_rg = count
    
      name = 'ic_corr_w_Ng'
      longname = 'In-Cloud Correlation of vertical velocity and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_w_Ng = count
    
    end if ! dograupel
  end if ! doicemicro

!----------------------------
! Cloud droplet number concentration 
!----------------------------
if(dopredictNc) then
  if (doprecip) then
        
    name = 'ic_corr_Nc_rr'
    longname = 'In-Cloud Correlation of cloud droplet conc. and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nc_rr = count
        
    name = 'ic_corr_Nc_Nr'
    longname = 'In-Cloud Correlation of cloud droplet conc. and rain drop conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nc_Nr = count
        
  end if ! doprecip
    
  if (doicemicro) then
        
    name = 'ic_corr_Nc_ri'
    longname = 'In-Cloud Correlation of cloud droplet conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nc_ri = count
        
    name = 'ic_corr_Nc_Ni'
    longname = 'In-Cloud Correlation of cloud droplet conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nc_Ni = count
       
    name = 'ic_corr_Nc_rs'
    longname = 'In-Cloud Correlation of cloud droplet conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nc_rs = count
      
    name = 'ic_corr_Nc_Ns'
    longname = 'In-Cloud Correlation of cloud droplet conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nc_Ns = count
        
    if (dograupel) then
        
      name = 'ic_corr_Nc_rg'
      longname = 'In-Cloud Correlation of cloud droplet conc. and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_Nc_rg = count
        
      name = 'ic_corr_Nc_Ng'
      longname = 'In-Cloud Correlation of cloud droplet conc. and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_Nc_Ng = count
        
    end if ! dograupel
  end if ! doicemicro 
end if ! dopredictNc

!----------------------------
! Rainwater mixing ratio 
!----------------------------
if (doprecip) then

    name = 'ic_corr_rr_Nr'
    longname = 'In-Cloud Correlation of rain water mixing ratio and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rr_Nr = count
  
  if(doicemicro) then

    name = 'ic_corr_rr_ri'
    longname = 'In-Cloud Correlation of rain water mixing ratio and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rr_ri = count
        
    name = 'ic_corr_rr_Ni'
    longname = 'In-Cloud Correlation of rain water mixing ratio and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rr_Ni = count
        
    name = 'ic_corr_rr_rs'
    longname = 'In-Cloud Correlation of rain water mixing ratio and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rr_rs = count
        
    name = 'ic_corr_rr_Ns'
    longname = 'In-Cloud Correlation of rain water mixing ratio and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rr_Ns = count
        

    if(dograupel) then

      name = 'ic_corr_rr_rg'
      longname = 'In-Cloud Correlation of rain water mixing ratio and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_rr_rg = count
        
      name = 'ic_corr_rr_Ng'
      longname = 'In-Cloud Correlation of rain water mixing ratio and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_rr_Ng = count
        
    end if !dograupel

!----------------------------
! Rainwater number concentration 
!----------------------------
        
    name = 'ic_corr_Nr_ri'
    longname = 'In-Cloud Correlation of rain drop conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nr_ri = count
        
    name = 'ic_corr_Nr_Ni'
    longname = 'In-Cloud Correlation of rain drop conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nr_Ni = count
        
    name = 'ic_corr_Nr_rs'
    longname = 'In-Cloud Correlation of rain drop conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nr_rs = count
        
    name = 'ic_corr_Nr_Ns'
    longname = 'In-Cloud Correlation of rain drop conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nr_Ns = count
        
      if(dograupel) then

        name = 'ic_corr_Nr_rg'
        longname = 'In-Cloud Correlation of rain drop conc. and graupel mixing ratio'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_ic_cor_Nr_rg = count
        
        name = 'ic_corr_Nr_Ng'
        longname = 'In-Cloud Correlation of rain drop conc. and graupel conc.'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_ic_cor_Nr_Ng = count
        
      end if !dograupel
    end if !doicemicro
end if !doprecip

!----------------------------
! Cloud-ice mixing ratio 
!----------------------------
if (doicemicro) then
    
  name = 'ic_corr_ri_Ni'
  longname = 'In-Cloud Correlation of cloud ice mixing ratio and cloud ice conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_ri_Ni = count
    
  name = 'ic_corr_ri_rs'
  longname = 'In-Cloud Correlation of cloud ice mixing ratio and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_ri_rs = count
    
  name = 'ic_corr_ri_Ns'
  longname = 'In-Cloud Correlation of cloud ice mixing ratio and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_ri_Ns = count
    
  if (dograupel) then

    name = 'ic_corr_ri_rg'
    longname = 'In-Cloud Correlation of cloud ice mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_ri_rg = count
        
    name = 'ic_corr_ri_Ng'
    longname = 'In-Cloud Correlation of cloud ice mixing ratio and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_ri_Ng = count
        
  end if !dograupel

!----------------------------
! Cloud-ice number concentration 
!----------------------------
    
  name = 'ic_corr_Ni_rs'
  longname = 'In-Cloud Correlation of cloud ice concentration and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_Ni_rs = count
    
  name = 'ic_corr_Ni_Ns'
  longname = 'In-Cloud Correlation of cloud ice concentration and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_Ni_Ns = count
    
  if (dograupel) then

    name = 'ic_corr_Ni_rg'
    longname = 'In-Cloud Correlation of cloud ice conc. and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Ni_rg = count
        
    name = 'ic_corr_Ni_Ng'
    longname = 'In-Cloud Correlation of cloud ice conc. and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Ni_Ng = count
        
  end if !dograupel

!----------------------------
! Snow mixing ratio 
!----------------------------
  
  name = 'ic_corr_rs_Ns'
  longname = 'In-Cloud Correlation of snow mixing ratio and snowflake concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_rs_Ns = count
    
  if (dograupel) then

    name = 'ic_corr_rs_rg'
    longname = 'In-Cloud Correlation of snow mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rs_rg = count
        
    name = 'ic_corr_rs_Ng'
    longname = 'In-Cloud Correlation of snow mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rs_Ng = count
        
!----------------------------
! Snow number concentration 
!----------------------------
    name = 'ic_corr_Ns_rg'
    longname = 'In-Cloud Correlation of snowflake concentration and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Ns_rg = count
        
    name = 'ic_corr_Ns_Ng'
    longname = 'In-Cloud Correlation of snowflake concentration and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Ns_Ng = count

!----------------------------
! Graupel mixing ratio / number concentration 
!----------------------------
        
    name = 'ic_corr_rg_Ng'
    longname = 'In-Cloud Correlation of graupel mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rg_Ng = count
        
  end if !dograupel
end if !doicemicro

!---------------------------
! Out of cloud oc_correlations 
!----------------------------

!----------------------------
! chi(S_mellor) 
!----------------------------
  name = 'oc_corr_chi_w'
  longname = 'In-Cloud Correlation of chi and vertical velocity'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_chi_w = count

if(dopredictNc) then
  name = 'oc_corr_chi_Nc'
  longname = 'In-Cloud Correlation of chi and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_chi_Nc = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'oc_corr_chi_rr'
    longname = 'In-Cloud Correlation of chi and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_chi_rr = count
    
    name = 'oc_corr_chi_Nr'
    longname = 'In-Cloud Correlation of chi and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_chi_Nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'oc_corr_chi_ri'
    longname = 'In-Cloud Correlation of chi and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_chi_ri = count
    
    name = 'oc_corr_chi_Ni'
    longname = 'In-Cloud Correlation of chi and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_chi_Ni = count
    
    name = 'oc_corr_chi_rs'
    longname = 'In-Cloud Correlation of chi and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_chi_rs = count
    
    name = 'oc_corr_chi_Ns'
    longname = 'In-Cloud Correlation of chi and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_chi_Ns = count
    

    if (dograupel) then
    
      name = 'oc_corr_chi_rg'
      longname = 'In-Cloud Correlation of chi and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_chi_rg = count
    
      name = 'oc_corr_chi_Ng'
      longname = 'In-Cloud Correlation of chi and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_chi_Ng = count
    
    end if ! dograupel
  endif ! doicemicro
!----------------------------
! Vertical Velocity 
!----------------------------

if(dopredictNc) then
  name = 'oc_corr_w_Nc'
  longname = 'In-Cloud Correlation of vertical velocity and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_w_Nc = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'oc_corr_w_rr'
    longname = 'In-Cloud Correlation of vertical velocity and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_w_rr = count
    
    name = 'oc_corr_w_Nr'
    longname = 'In-Cloud Correlation of vertical velocity and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_w_Nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'oc_corr_w_ri'
    longname = 'In-Cloud Correlation of vertical velocity and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_w_ri = count
    
    name = 'oc_corr_w_Ni'
    longname = 'In-Cloud Correlation of vertical velocity and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_w_Ni = count
    
    name = 'oc_corr_w_rs'
    longname = 'In-Cloud Correlation of vertical velocity and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_w_rs = count
    
    name = 'oc_corr_w_Ns'
    longname = 'In-Cloud Correlation of vertical velocity and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_w_Ns = count
    

    if (dograupel) then
    
      name = 'oc_corr_w_rg'
      longname = 'In-Cloud Correlation of vertical velocity and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_w_rg = count
    
      name = 'oc_corr_w_Ng'
      longname = 'In-Cloud Correlation of vertical velocity and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_w_Ng = count
    
    end if ! dograupel
  end if ! doicemicro

!----------------------------
! Cloud droplet number concentration 
!----------------------------
if(dopredictNc) then
  if (doprecip) then
        
    name = 'oc_corr_Nc_rr'
    longname = 'In-Cloud Correlation of cloud droplet conc. and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nc_rr = count
        
    name = 'oc_corr_Nc_Nr'
    longname = 'In-Cloud Correlation of cloud droplet conc. and rain drop conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nc_Nr = count
        
  end if ! doprecip
    
  if (doicemicro) then
        
    name = 'oc_corr_Nc_ri'
    longname = 'In-Cloud Correlation of cloud droplet conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nc_ri = count
        
    name = 'oc_corr_Nc_Ni'
    longname = 'In-Cloud Correlation of cloud droplet conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nc_Ni = count
       
    name = 'oc_corr_Nc_rs'
    longname = 'In-Cloud Correlation of cloud droplet conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nc_rs = count
      
    name = 'oc_corr_Nc_Ns'
    longname = 'In-Cloud Correlation of cloud droplet conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nc_Ns = count
        
    if (dograupel) then
        
      name = 'oc_corr_Nc_rg'
      longname = 'In-Cloud Correlation of cloud droplet conc. and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_Nc_rg = count
        
      name = 'oc_corr_Nc_Ng'
      longname = 'In-Cloud Correlation of cloud droplet conc. and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_Nc_Ng = count
        
    end if ! dograupel
  end if ! doicemicro 
end if ! dopredictNc

!----------------------------
! Rainwater mixing ratio 
!----------------------------
if (doprecip) then

    name = 'oc_corr_rr_Nr'
    longname = 'In-Cloud Correlation of rain water mixing ratio and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rr_Nr = count
  
  if(doicemicro) then

    name = 'oc_corr_rr_ri'
    longname = 'In-Cloud Correlation of rain water mixing ratio and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rr_ri = count
        
    name = 'oc_corr_rr_Ni'
    longname = 'In-Cloud Correlation of rain water mixing ratio and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rr_Ni = count
        
    name = 'oc_corr_rr_rs'
    longname = 'In-Cloud Correlation of rain water mixing ratio and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rr_rs = count
        
    name = 'oc_corr_rr_Ns'
    longname = 'In-Cloud Correlation of rain water mixing ratio and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rr_Ns = count
        

    if(dograupel) then

      name = 'oc_corr_rr_rg'
      longname = 'In-Cloud Correlation of rain water mixing ratio and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_rr_rg = count
        
      name = 'oc_corr_rr_Ng'
      longname = 'In-Cloud Correlation of rain water mixing ratio and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_rr_Ng = count
        
    end if !dograupel

!----------------------------
! Rainwater number concentration 
!----------------------------
        
    name = 'oc_corr_Nr_ri'
    longname = 'In-Cloud Correlation of rain drop conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nr_ri = count
        
    name = 'oc_corr_Nr_Ni'
    longname = 'In-Cloud Correlation of rain drop conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nr_Ni = count
        
    name = 'oc_corr_Nr_rs'
    longname = 'In-Cloud Correlation of rain drop conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nr_rs = count
        
    name = 'oc_corr_Nr_Ns'
    longname = 'In-Cloud Correlation of rain drop conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nr_Ns = count
        
      if(dograupel) then

        name = 'oc_corr_Nr_rg'
        longname = 'In-Cloud Correlation of rain drop conc. and graupel mixing ratio'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_oc_cor_Nr_rg = count
        
        name = 'oc_corr_Nr_Ng'
        longname = 'In-Cloud Correlation of rain drop conc. and graupel conc.'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_oc_cor_Nr_Ng = count
        
      end if !dograupel
    end if !doicemicro
end if !doprecip

!----------------------------
! Cloud-ice mixing ratio 
!----------------------------
if (doicemicro) then
    
  name = 'oc_corr_ri_Ni'
  longname = 'In-Cloud Correlation of cloud ice mixing ratio and cloud ice conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_ri_Ni = count
    
  name = 'oc_corr_ri_rs'
  longname = 'In-Cloud Correlation of cloud ice mixing ratio and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_ri_rs = count
    
  name = 'oc_corr_ri_Ns'
  longname = 'In-Cloud Correlation of cloud ice mixing ratio and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_ri_Ns = count
    
  if (dograupel) then

    name = 'oc_corr_ri_rg'
    longname = 'In-Cloud Correlation of cloud ice mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_ri_rg = count
        
    name = 'oc_corr_ri_Ng'
    longname = 'In-Cloud Correlation of cloud ice mixing ratio and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_ri_Ng = count
        
  end if !dograupel

!----------------------------
! Cloud-ice number concentration 
!----------------------------
    
  name = 'oc_corr_Ni_rs'
  longname = 'In-Cloud Correlation of cloud ice concentration and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_Ni_rs = count
    
  name = 'oc_corr_Ni_Ns'
  longname = 'In-Cloud Correlation of cloud ice concentration and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_Ni_Ns = count
    
  if (dograupel) then

    name = 'oc_corr_Ni_rg'
    longname = 'In-Cloud Correlation of cloud ice conc. and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Ni_rg = count
        
    name = 'oc_corr_Ni_Ng'
    longname = 'In-Cloud Correlation of cloud ice conc. and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Ni_Ng = count
        
  end if !dograupel

!----------------------------
! Snow mixing ratio 
!----------------------------
  
  name = 'oc_corr_rs_Ns'
  longname = 'In-Cloud Correlation of snow mixing ratio and snowflake concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_rs_Ns = count
    
  if (dograupel) then

    name = 'oc_corr_rs_rg'
    longname = 'In-Cloud Correlation of snow mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rs_rg = count
        
    name = 'oc_corr_rs_Ng'
    longname = 'In-Cloud Correlation of snow mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rs_Ng = count
        
!----------------------------
! Snow number concentration 
!----------------------------
    name = 'oc_corr_Ns_rg'
    longname = 'In-Cloud Correlation of snowflake concentration and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Ns_rg = count
        
    name = 'oc_corr_Ns_Ng'
    longname = 'In-Cloud Correlation of snowflake concentration and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Ns_Ng = count

!----------------------------
! Graupel mixing ratio / number concentration 
!----------------------------
        
    name = 'oc_corr_rg_Ng'
    longname = 'In-Cloud Correlation of graupel mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rg_Ng = count
        
  end if !dograupel
end if !doicemicro

!----------------------------
! Fractions
!----------------------------


! This section creates the names and the number of hydrometeor fraction fields
! we have
frac_name(1) = 'rc_frac'
frac_longname(1) = 'Cloud liquid water fraction'
nfrac_fields = 1 

if (doprecip) then 

  frac_name(2) = 'rr_frac'
  frac_longname(2) = 'Rain water fraction'
  nfrac_fields = nfrac_fields + 1

  if (doicemicro) then

    frac_name(3) = 'ri_frac'
    frac_longname(3) = 'Cloud ice fraction'
    nfrac_fields = nfrac_fields + 1

    frac_name(4) = 'rs_frac'
    frac_longname(4) = 'Snow fraction'
    nfrac_fields = nfrac_fields + 1

    if (dograupel) then

      frac_name(5) = 'rg_frac'
      frac_longname(5) = 'Graupel fraction'
      nfrac_fields = nfrac_fields + 1

    endif ! dograupel
  endif ! doicemicro
endif ! doprecip

do n = 1,nfrac_fields ! Loop through number of hydrometeor fraction fields
  do m = 1,nfractions ! Loop through the number of fractions the user requests,
                      ! which is set in the namelist
    
    ! The user specifies an initial threshold to determine fractional coverage.
    ! Increment each subsequet fraction an order of magnitude
    curr_frac = (10**(m-1))*frac_threshold_init

    ! Convert the number to a string
    WRITE(frac_in_char_temp, '(ES8.1)' ) curr_frac 

    ! Kludgy, but replace the decimal (.) and minus (-) with 'p' and 'm'
    ! respectively.
    frac_in_char_temp(3:3) = 'p'
    frac_in_char_temp(6:6) = 'm'

    ! Stor in array
    frac_in_char(m) = frac_in_char_temp  

    ! Create hydrometeor fraction name, append the threshold value used. 
    name = trim(frac_name(n))//'_'//trim(adjustl(frac_in_char(m)))
    longname = trim(frac_longname(n))
    write(*,*) name
    units = '[fraction]'
    call add_to_namelist(count,microcount,name,longname,units,0)

  end do ! nfrac_fields
end do ! nfractions

!----------------------------
! Domain-wide means
!----------------------------
name = 'Ncm'
longname = 'Domain-wide mean of Nc'
units = '[(kg/kg)]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'chim'
longname = 'Domain-wide mean of chi'
units = '[(kg/kg)]'
call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Domain-wide variances
!----------------------------

name = 'chip2'
longname = 'Domain-wide variance of cloud liquid number concentration'
units = '[(kg/kg)^2]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'Ncp2'
longname = 'Domain-wide variance of cloud liquid number concentration'
units = '[(#/kg)^2]'
call add_to_namelist(count,microcount,name,longname,units,0)

if (doprecip) then
  
    name = 'Nrp2'
    longname = 'Domain-wide variance of rain number concentration'
    units = '[(#/kg)^2]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'rrp2'
    longname = 'Domain-wide variance of rain mixing ratio'
    units = '[(kg/kg)^2]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    if (doicemicro) then
      
      name = 'Nip2'
      longname = 'Domain-wide variance of cloud ice number concentration'
      units = '[(#/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'rip2'
      longname = 'Domain-wide variance of cloud ice mixing ratio'
      units = '[(kg/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'Nsp2'
      longname = 'Domain-wide variance of snow number concentration'
      units = '[(#/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)

      name = 'rsp2'
      longname = 'Domain-wide variance of snow mixing ratio'
      units = '[(kg/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)

      if (dograupel) then
    
        name = 'Ngp2'
        longname = 'Domain-wide variance of graupel number concentration'
        units = '[(#/kg)^2]'
        call add_to_namelist(count,microcount,name,longname,units,0)

        name = 'rgp2'
        longname = 'Domain-wide variance of graupel mixing ratio'
        units = '[(kg/kg)^2]'
        call add_to_namelist(count,microcount,name,longname,units,0)

      end if !dograupel
    end if !doicemicro
  end if !doprecip


!----------------------------
! In-cloud means and variances
!----------------------------
name = 'chim_ic'
longname = 'In-cloud mean of chi'
units = '[(kg/kg)]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'chip2_ic'
longname = 'In-cloud variance of chi'
units = '[(kg/kg)^2]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'wm_zt_ic'
longname = 'In-cloud mean of w (on zt grid)'
units = '[(m/s)]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'wp2_zt_ic'
longname = 'In-cloud variance of w (on zt grid)'
units = '[(kg/kg)^2]'
call add_to_namelist(count,microcount,name,longname,units,0)



!----------------------------
! In-precip means and variances
!----------------------------

name = 'Ncp2_ip'
longname = 'Within cloud variance of cloud liquid number concentration'
units = '[(#/kg)^2]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'Ncm_ip'
longname = 'Within cloud mean of cloud liquid number concentration'
units = '[(#/kg)]'
call add_to_namelist(count,microcount,name,longname,units,0)

if (doprecip) then
  
    name = 'Nrp2_ip'
    longname = 'Within rain variance of rain number concentration'
    units = '[(#/kg)^2]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'Nrm_ip'
    longname = 'Within rain mean of rain number concentration'
    units = '[(#/kg)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'rrp2_ip'
    longname = 'Within rain variance of rain mixing ratio'
    units = '[(kg/kg)^2]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'rrm_ip'
    longname = 'Within rain mean of rain mixing ratio'
    units = '[(kg/kg)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    if (doicemicro) then
      
      name = 'Nip2_ip'
      longname = 'Within cloudice variance of cloud ice number concentration'
      units = '[(#/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'Nim_ip'
      longname = 'Within cloudice mean of cloud ice number concentration'
      units = '[(#/kg)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'rip2_ip'
      longname = 'Within cloudice variance of cloud ice mixing ratio'
      units = '[(kg/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'rim_ip'
      longname = 'Within cloudice mean of cloud ice mixing ratio'
      units = '[(kg/kg)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'Nsp2_ip'
      longname = 'Within snow variance of snow number concentration'
      units = '[(#/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'Nsm_ip'
      longname = 'Within snow mean of snow number concentration'
      units = '[(#/kg)]'
      call add_to_namelist(count,microcount,name,longname,units,0)

      name = 'rsp2_ip'
      longname = 'Within snow variance of snow mixing ratio'
      units = '[(kg/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'rsm_ip'
      longname = 'Within snow mean of snow mixing ratio'
      units = '[(kg/kg)]'
      call add_to_namelist(count,microcount,name,longname,units,0)

      if (dograupel) then
    
        name = 'Ngp2_ip'
        longname = 'Within graupel variance of graupel number concentration'
        units = '[(#/kg)^2]'
        call add_to_namelist(count,microcount,name,longname,units,0)
        
        name = 'Ngm_ip'
        longname = 'Within graupel mean of graupel number concentration'
        units = '[(#/kg)]'
        call add_to_namelist(count,microcount,name,longname,units,0)

        name = 'rgp2_ip'
        longname = 'Within graupel variance of graupel mixing ratio'
        units = '[(kg/kg)^2]'
        call add_to_namelist(count,microcount,name,longname,units,0)
        
        name = 'rgm_ip'
        longname = 'Within graupel mean of graupel mixing ratio'
        units = '[(kg/kg)]'
        call add_to_namelist(count,microcount,name,longname,units,0)

      end if !dograupel
    end if !doicemicro
  end if !doprecip

!----------------------------
! Covariances
!---------------------------

!----------------------------
! chi(s_mellor)
!---------------------------
name = 'covarnce_chi_w'
longname = 'Covariance of chi and cloud droplet concentration'
units = '[(kg kg^-1)(m s^-1) ]'
call add_to_namelist(count,microcount,name,longname,units,0)

if(dopredictNc) then
  
  name = 'covarnce_chi_Nc'
  longname = 'Covariance of chi and cloud droplet concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'covarnce_chi_rr'
  longname = 'Covariance of chi and rain water mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_chi_Nr'
  longname = 'Covariance of chi and rain drop concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'covarnce_chi_ri'
  longname = 'Covariance of chi and cloud ice mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_chi_Ni'
  longname = 'Covariance of chi and cloud ice concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_chi_rs'
  longname = 'Covariance of chi and snow mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_chi_Ns'
  longname = 'Covariance of chi and snowflake concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'covarnce_chi_rg'
  longname = 'Covariance of chi and graupel mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_chi_Ng'
  longname = 'Covariance of chi and graupel concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

endif !doicemicro

!----------------------------
! Vertical velocity
!---------------------------
if(dopredictNc) then
  
  name = 'covarnce_w_Nc'
  longname = 'Covariance of vertical velocity and cloud droplet concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'covarnce_w_rr'
  longname = 'Covariance of vertical velocity and rain water mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_Nr'
  longname = 'Covariance of vertical velocity and rain drop concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'covarnce_w_ri'
  longname = 'Covariance of vertical velocity and cloud ice mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_Ni'
  longname = 'Covariance of vertical velocity and cloud ice concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_rs'
  longname = 'Covariance of vertical velocity and snow mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_Ns'
  longname = 'Covariance of vertical velocity and snowflake concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'covarnce_w_rg'
  longname = 'Covariance of vertical velocity and graupel mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_Ng'
  longname = 'Covariance of vertical velocity and graupel concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! dograupel

!----------------------------
! Cloud droplet number concentration
!---------------------------

if(dopredictNc) then
    if (doprecip) then
        
      name = 'covarnce_Nc_rr'
      longname = 'Covariance of cloud droplet conc. and rain water mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'covarnce_Nc_Nr'
      longname = 'Covariance of cloud droplet conc. and rain drop conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doprecip
    
    if (doicemicro) then
        
      name = 'covarnce_Nc_ri'
      longname = 'Covariance of cloud droplet conc. and cloud ice mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'covarnce_Nc_Ni'
      longname = 'Covariance of cloud droplet conc. and cloud ice conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'covarnce_Nc_rs'
      longname = 'Covariance of cloud droplet conc. and snow mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'covarnce_Nc_Ns'
      longname = 'Covariance of cloud droplet conc. and snowflake conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doicemicro
    
    if (dograupel) then
        
      name = 'covarnce_Nc_rg'
      longname = 'Covariance of cloud droplet conc. and graupel mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'covarnce_Nc_Ng'
      longname = 'Covariance of cloud droplet conc. and graupel conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)

    end if ! dograupel
end if ! do predictNc

!----------------------------
! Rainwater mixing ratio
!---------------------------
if (doprecip) then

    name = 'covarnce_rr_Nr'
    longname = 'Covariance of rain water mixing ratio and rain drop concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

  if(doicemicro) then

    name = 'covarnce_rr_ri'
    longname = 'Covariance of rain water mixing ratio and cloud ice mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'covarnce_rr_Ni'
    longname = 'Covariance of rain water mixing ratio and cloud ice conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_rr_rs'
    longname = 'Covariance of rain water mixing ratio and snow mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_rr_Ns'
    longname = 'Covariance of rain water mixing ratio and snowflake conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'covarnce_rr_rg'
    longname = 'Covariance of rain water mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'covarnce_rr_Ng'
    longname = 'Covariance of rain water mixing ratio and graupel conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

!----------------------------
! Rainwater number concentration
!---------------------------
  if(doicemicro) then
        
    name = 'covarnce_Nr_ri'
    longname = 'Covariance of rain drop conc. and cloud ice mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'covarnce_Nr_Ni'
    longname = 'Covariance of rain drop conc. and cloud ice conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_Nr_rs'
    longname = 'Covariance of rain drop conc. and snow mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_Nr_Ns'
    longname = 'Covariance of rain drop conc. and snowflake conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'covarnce_Nr_rg'
    longname = 'Covariance of rain drop conc. and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'covarnce_Nr_Ng'
    longname = 'Covariance of rain drop conc. and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doprecip
!----------------------------
! Cloudice mixing ratio
!---------------------------

if (doicemicro) then
    
  name = 'covarnce_ri_Ni'
  longname = 'Covariance of cloud ice mixing ratio and cloud ice conc.'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_ri_rs'
  longname = 'Covariance of cloud ice mixing ratio and snow mixing ratio'
  units = '[(kg kg^-1)(kg kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_ri_Ns'
  longname = 'Covariance of cloud ice mixing ratio and snowflake conc.'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'covarnce_ri_rg'
    longname = 'Covariance of cloud ice mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'covarnce_ri_Ng'
    longname = 'Covariance of cloud ice mixing ratio and graupel conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

!----------------------------
! Cloudice number concentration
!---------------------------
    
  name = 'covarnce_Ni_rs'
  longname = 'Covariance of cloud ice concentration and snow mixing ratio'
  units = '[(# kg^-1)(kg kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_Ni_Ns'
  longname = 'Covariance of cloud ice concentration and snowflake conc.'
  units = '[(# kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'covarnce_Ni_rg'
    longname = 'Covariance of cloud ice conc. and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'covarnce_Ni_Ng'
    longname = 'Covariance of cloud ice conc. and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  end if !dograupel

!----------------------------
! Snow mixing ratio
!---------------------------
    
  name = 'covarnce_rs_Ns'
  longname = 'Covariance of snow mixing ratio and snowflake concentration'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'covarnce_rs_rg'
    longname = 'Covariance of snow mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'covarnce_rs_Ng'
    longname = 'Covariance of snow mixing ratio and graupel concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Snow number concentration
!---------------------------
    
    name = 'covarnce_Ns_rg'
    longname = 'Covariance of snowflake concentration and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'covarnce_Ns_Ng'
    longname = 'Covariance of snowflake concentration and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Graupel mixing ratio / number concentration
!---------------------------
    
    name = 'covarnce_rg_Ng'
    longname = 'Covariance of graupel mixing ratio and graupel concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doicemicro

!----------------------------
! In-cloud covariances
!---------------------------

!----------------------------
! chi(s_mellor)
!---------------------------
name = 'ic_covarnce_chi_w'
longname = 'In-cloud Covariance of chi and cloud droplet concentration'
units = '[(kg kg^-1)(m s^-1) ]'
call add_to_namelist(count,microcount,name,longname,units,0)

if(dopredictNc) then
  
  name = 'ic_covarnce_chi_Nc'
  longname = 'In-cloud Covariance of chi and cloud droplet concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'ic_covarnce_chi_rr'
  longname = 'In-cloud Covariance of chi and rain water mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_chi_Nr'
  longname = 'In-cloud Covariance of chi and rain drop concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'ic_covarnce_chi_ri'
  longname = 'In-cloud Covariance of chi and cloud ice mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_chi_Ni'
  longname = 'In-cloud Covariance of chi and cloud ice concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_chi_rs'
  longname = 'In-cloud Covariance of chi and snow mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_chi_Ns'
  longname = 'In-cloud Covariance of chi and snowflake concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'ic_covarnce_chi_rg'
  longname = 'In-cloud Covariance of chi and graupel mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_chi_Ng'
  longname = 'In-cloud Covariance of chi and graupel concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

endif !doicemicro

!----------------------------
! Vertical velocity
!---------------------------
if(dopredictNc) then
  
  name = 'ic_covarnce_w_Nc'
  longname = 'In-cloud Covariance of vertical velocity and cloud droplet concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'ic_covarnce_w_rr'
  longname = 'In-cloud Covariance of vertical velocity and rain water mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_w_Nr'
  longname = 'In-cloud Covariance of vertical velocity and rain drop concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'ic_covarnce_w_ri'
  longname = 'In-cloud Covariance of vertical velocity and cloud ice mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_w_Ni'
  longname = 'In-cloud Covariance of vertical velocity and cloud ice concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_w_rs'
  longname = 'In-cloud Covariance of vertical velocity and snow mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_w_Ns'
  longname = 'In-cloud Covariance of vertical velocity and snowflake concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'ic_covarnce_w_rg'
  longname = 'In-cloud Covariance of vertical velocity and graupel mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_w_Ng'
  longname = 'In-cloud Covariance of vertical velocity and graupel concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! dograupel

!----------------------------
! Cloud droplet number concentration
!---------------------------

if(dopredictNc) then
    if (doprecip) then
        
      name = 'ic_covarnce_Nc_rr'
      longname = 'In-cloud Covariance of cloud droplet conc. and rain water mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'ic_covarnce_Nc_Nr'
      longname = 'In-cloud Covariance of cloud droplet conc. and rain drop conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doprecip
    
    if (doicemicro) then
        
      name = 'ic_covarnce_Nc_ri'
      longname = 'In-cloud Covariance of cloud droplet conc. and cloud ice mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'ic_covarnce_Nc_Ni'
      longname = 'In-cloud Covariance of cloud droplet conc. and cloud ice conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'ic_covarnce_Nc_rs'
      longname = 'In-cloud Covariance of cloud droplet conc. and snow mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'ic_covarnce_Nc_Ns'
      longname = 'In-cloud Covariance of cloud droplet conc. and snowflake conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doicemicro
    
    if (dograupel) then
        
      name = 'ic_covarnce_Nc_rg'
      longname = 'In-cloud Covariance of cloud droplet conc. and graupel mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'ic_covarnce_Nc_Ng'
      longname = 'In-cloud Covariance of cloud droplet conc. and graupel conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)

    end if ! dograupel
end if ! do predictNc

!----------------------------
! Rainwater mixing ratio
!---------------------------
if (doprecip) then

    name = 'ic_covarnce_rr_Nr'
    longname = 'In-cloud Covariance of rain water mixing ratio and rain drop concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

  if(doicemicro) then

    name = 'ic_covarnce_rr_ri'
    longname = 'In-cloud Covariance of rain water mixing ratio and cloud ice mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'ic_covarnce_rr_Ni'
    longname = 'In-cloud Covariance of rain water mixing ratio and cloud ice conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'ic_covarnce_rr_rs'
    longname = 'In-cloud Covariance of rain water mixing ratio and snow mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'ic_covarnce_rr_Ns'
    longname = 'In-cloud Covariance of rain water mixing ratio and snowflake conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'ic_covarnce_rr_rg'
    longname = 'In-cloud Covariance of rain water mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'ic_covarnce_rr_Ng'
    longname = 'In-cloud Covariance of rain water mixing ratio and graupel conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

!----------------------------
! Rainwater number concentration
!---------------------------
  if(doicemicro) then
        
    name = 'ic_covarnce_Nr_ri'
    longname = 'In-cloud Covariance of rain drop conc. and cloud ice mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'ic_covarnce_Nr_Ni'
    longname = 'In-cloud Covariance of rain drop conc. and cloud ice conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'ic_covarnce_Nr_rs'
    longname = 'In-cloud Covariance of rain drop conc. and snow mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'ic_covarnce_Nr_Ns'
    longname = 'In-cloud Covariance of rain drop conc. and snowflake conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'ic_covarnce_Nr_rg'
    longname = 'In-cloud Covariance of rain drop conc. and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'ic_covarnce_Nr_Ng'
    longname = 'In-cloud Covariance of rain drop conc. and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doprecip
!----------------------------
! Cloudice mixing ratio
!---------------------------

if (doicemicro) then
    
  name = 'ic_covarnce_ri_Ni'
  longname = 'In-cloud Covariance of cloud ice mixing ratio and cloud ice conc.'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_ri_rs'
  longname = 'In-cloud Covariance of cloud ice mixing ratio and snow mixing ratio'
  units = '[(kg kg^-1)(kg kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_ri_Ns'
  longname = 'In-cloud Covariance of cloud ice mixing ratio and snowflake conc.'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'ic_covarnce_ri_rg'
    longname = 'In-cloud Covariance of cloud ice mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'ic_covarnce_ri_Ng'
    longname = 'In-cloud Covariance of cloud ice mixing ratio and graupel conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

!----------------------------
! Cloudice number concentration
!---------------------------
    
  name = 'ic_covarnce_Ni_rs'
  longname = 'In-cloud Covariance of cloud ice concentration and snow mixing ratio'
  units = '[(# kg^-1)(kg kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_Ni_Ns'
  longname = 'In-cloud Covariance of cloud ice concentration and snowflake conc.'
  units = '[(# kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'ic_covarnce_Ni_rg'
    longname = 'In-cloud Covariance of cloud ice conc. and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'ic_covarnce_Ni_Ng'
    longname = 'In-cloud Covariance of cloud ice conc. and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  end if !dograupel

!----------------------------
! Snow mixing ratio
!---------------------------
    
  name = 'ic_covarnce_rs_Ns'
  longname = 'In-cloud Covariance of snow mixing ratio and snowflake concentration'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'ic_covarnce_rs_rg'
    longname = 'In-cloud Covariance of snow mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'ic_covarnce_rs_Ng'
    longname = 'In-cloud Covariance of snow mixing ratio and graupel concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Snow number concentration
!---------------------------
    
    name = 'ic_covarnce_Ns_rg'
    longname = 'In-cloud Covariance of snowflake concentration and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'ic_covarnce_Ns_Ng'
    longname = 'In-cloud Covariance of snowflake concentration and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Graupel mixing ratio / number concentration
!---------------------------
    
    name = 'ic_covarnce_rg_Ng'
    longname = 'In-cloud Covariance of graupel mixing ratio and graupel concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doicemicro

!----------------------------
!  Out of cloud covariances
!---------------------------

!----------------------------
! chi(s_mellor)
!---------------------------
name = 'oc_covarnce_chi_w'
longname = 'In-cloud Covariance of chi and cloud droplet concentration'
units = '[(kg kg^-1)(m s^-1) ]'
call add_to_namelist(count,microcount,name,longname,units,0)

if(dopredictNc) then
  
  name = 'oc_covarnce_chi_Nc'
  longname = 'In-cloud Covariance of chi and cloud droplet concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'oc_covarnce_chi_rr'
  longname = 'In-cloud Covariance of chi and rain water mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_chi_Nr'
  longname = 'In-cloud Covariance of chi and rain drop concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'oc_covarnce_chi_ri'
  longname = 'In-cloud Covariance of chi and cloud ice mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_chi_Ni'
  longname = 'In-cloud Covariance of chi and cloud ice concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_chi_rs'
  longname = 'In-cloud Covariance of chi and snow mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_chi_Ns'
  longname = 'In-cloud Covariance of chi and snowflake concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'oc_covarnce_chi_rg'
  longname = 'In-cloud Covariance of chi and graupel mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_chi_Ng'
  longname = 'In-cloud Covariance of chi and graupel concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

endif !doicemicro

!----------------------------
! Vertical velocity
!---------------------------
if(dopredictNc) then
  
  name = 'oc_covarnce_w_Nc'
  longname = 'In-cloud Covariance of vertical velocity and cloud droplet concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'oc_covarnce_w_rr'
  longname = 'In-cloud Covariance of vertical velocity and rain water mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_w_Nr'
  longname = 'In-cloud Covariance of vertical velocity and rain drop concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'oc_covarnce_w_ri'
  longname = 'In-cloud Covariance of vertical velocity and cloud ice mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_w_Ni'
  longname = 'In-cloud Covariance of vertical velocity and cloud ice concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_w_rs'
  longname = 'In-cloud Covariance of vertical velocity and snow mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_w_Ns'
  longname = 'In-cloud Covariance of vertical velocity and snowflake concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'oc_covarnce_w_rg'
  longname = 'In-cloud Covariance of vertical velocity and graupel mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_w_Ng'
  longname = 'In-cloud Covariance of vertical velocity and graupel concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! dograupel

!----------------------------
! Cloud droplet number concentration
!---------------------------

if(dopredictNc) then
    if (doprecip) then
        
      name = 'oc_covarnce_Nc_rr'
      longname = 'In-cloud Covariance of cloud droplet conc. and rain water mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'oc_covarnce_Nc_Nr'
      longname = 'In-cloud Covariance of cloud droplet conc. and rain drop conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doprecip
    
    if (doicemicro) then
        
      name = 'oc_covarnce_Nc_ri'
      longname = 'In-cloud Covariance of cloud droplet conc. and cloud ice mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'oc_covarnce_Nc_Ni'
      longname = 'In-cloud Covariance of cloud droplet conc. and cloud ice conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'oc_covarnce_Nc_rs'
      longname = 'In-cloud Covariance of cloud droplet conc. and snow mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'oc_covarnce_Nc_Ns'
      longname = 'In-cloud Covariance of cloud droplet conc. and snowflake conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doicemicro
    
    if (dograupel) then
        
      name = 'oc_covarnce_Nc_rg'
      longname = 'In-cloud Covariance of cloud droplet conc. and graupel mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'oc_covarnce_Nc_Ng'
      longname = 'In-cloud Covariance of cloud droplet conc. and graupel conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)

    end if ! dograupel
end if ! do predictNc

!----------------------------
! Rainwater mixing ratio
!---------------------------
if (doprecip) then

    name = 'oc_covarnce_rr_Nr'
    longname = 'In-cloud Covariance of rain water mixing ratio and rain drop concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

  if(doicemicro) then

    name = 'oc_covarnce_rr_ri'
    longname = 'In-cloud Covariance of rain water mixing ratio and cloud ice mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'oc_covarnce_rr_Ni'
    longname = 'In-cloud Covariance of rain water mixing ratio and cloud ice conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'oc_covarnce_rr_rs'
    longname = 'In-cloud Covariance of rain water mixing ratio and snow mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'oc_covarnce_rr_Ns'
    longname = 'In-cloud Covariance of rain water mixing ratio and snowflake conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'oc_covarnce_rr_rg'
    longname = 'In-cloud Covariance of rain water mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'oc_covarnce_rr_Ng'
    longname = 'In-cloud Covariance of rain water mixing ratio and graupel conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

!----------------------------
! Rainwater number concentration
!---------------------------
  if(doicemicro) then
        
    name = 'oc_covarnce_Nr_ri'
    longname = 'In-cloud Covariance of rain drop conc. and cloud ice mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'oc_covarnce_Nr_Ni'
    longname = 'In-cloud Covariance of rain drop conc. and cloud ice conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'oc_covarnce_Nr_rs'
    longname = 'In-cloud Covariance of rain drop conc. and snow mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'oc_covarnce_Nr_Ns'
    longname = 'In-cloud Covariance of rain drop conc. and snowflake conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'oc_covarnce_Nr_rg'
    longname = 'In-cloud Covariance of rain drop conc. and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'oc_covarnce_Nr_Ng'
    longname = 'In-cloud Covariance of rain drop conc. and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doprecip
!----------------------------
! Cloudice mixing ratio
!---------------------------

if (doicemicro) then
    
  name = 'oc_covarnce_ri_Ni'
  longname = 'In-cloud Covariance of cloud ice mixing ratio and cloud ice conc.'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_ri_rs'
  longname = 'In-cloud Covariance of cloud ice mixing ratio and snow mixing ratio'
  units = '[(kg kg^-1)(kg kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_ri_Ns'
  longname = 'In-cloud Covariance of cloud ice mixing ratio and snowflake conc.'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'oc_covarnce_ri_rg'
    longname = 'In-cloud Covariance of cloud ice mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'oc_covarnce_ri_Ng'
    longname = 'In-cloud Covariance of cloud ice mixing ratio and graupel conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

!----------------------------
! Cloudice number concentration
!---------------------------
    
  name = 'oc_covarnce_Ni_rs'
  longname = 'In-cloud Covariance of cloud ice concentration and snow mixing ratio'
  units = '[(# kg^-1)(kg kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_Ni_Ns'
  longname = 'In-cloud Covariance of cloud ice concentration and snowflake conc.'
  units = '[(# kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'oc_covarnce_Ni_rg'
    longname = 'In-cloud Covariance of cloud ice conc. and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'oc_covarnce_Ni_Ng'
    longname = 'In-cloud Covariance of cloud ice conc. and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  end if !dograupel

!----------------------------
! Snow mixing ratio
!---------------------------
    
  name = 'oc_covarnce_rs_Ns'
  longname = 'In-cloud Covariance of snow mixing ratio and snowflake concentration'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'oc_covarnce_rs_rg'
    longname = 'In-cloud Covariance of snow mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'oc_covarnce_rs_Ng'
    longname = 'In-cloud Covariance of snow mixing ratio and graupel concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Snow number concentration
!---------------------------
    
    name = 'oc_covarnce_Ns_rg'
    longname = 'In-cloud Covariance of snowflake concentration and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'oc_covarnce_Ns_Ng'
    longname = 'In-cloud Covariance of snowflake concentration and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Graupel mixing ratio / number concentration
!---------------------------
    
    name = 'oc_covarnce_rg_Ng'
    longname = 'In-cloud Covariance of graupel mixing ratio and graupel concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doicemicro

if(masterproc) then
   write(*,*) 'Added ', microcount, ' arrays to statistics for M2005 microphysics'
end if