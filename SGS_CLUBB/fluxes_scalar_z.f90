subroutine fluxes_scalar_z (f,fluxb,fluxt, &
                          fdiff,flux,f2lediff,f2lediss,fwlediff,doit)

!--------------------------------------------------------------------
! This subroutine is only used to apply the surface fluxes for scalars. 
! This is needed when surface fluxes are applied in the host model in SAM_CLUBB 
! Here tkh is zet to zero so vertical diffusion is not calculated. 
! Minghuai Wang, 2013-02  
!---------------------------------------------------------------------

use grid
use vars, only: rho, rhow
!use sgs, only: tkh
implicit none

! input:	
real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real fluxb(nx,ny)		! bottom flux
real fluxt(nx,ny)		! top flux
real flux(nz)
real fdiff(nz)
real f2lediff(nzm)
real f2lediss(nzm)
real fwlediff(nzm)
logical doit
! Local
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real f0(nzm),df0(nzm),factor_xy
real f2lediss_z(nzm)
real tkh2(0:nxp1, 1-YES3D:nyp1, nzm)     ! eddy conductivity
real r2dx,r2dy,r2dx0,r2dy0,r2dz
integer i,j,k,kb,kc,jb,jc

call t_startf ('fluxes_scalars_z')

tkh2 = 0.0

if(dostatis) then
	
  do k=1,nzm
    do j=dimy1_s,dimy2_s
     do i=dimx1_s,dimx2_s
      df(i,j,k) = f(i,j,k)
     end do
    end do
  end do

endif


if(RUN3D) then
  call diffuse_scalar3D_z (f,fluxb,fluxt,tkh2,rho,rhow,flux)
else  
  call diffuse_scalar2D_z (f,fluxb,fluxt,tkh2,rho,rhow,flux)
endif

if(dostatis) then
	
  do k=1,nzm
    fdiff(k)=0.
    do j=1,ny
     do i=1,nx
      fdiff(k)=fdiff(k)+f(i,j,k)-df(i,j,k)
     end do
    end do
  end do

endif

if(dostatis.and.doit) then
	
  call stat_varscalar(f,df,f0,df0,f2lediff)
  call stat_sw2(f,df,fwlediff)

  factor_xy=1./float(nx*ny)
  r2dx0=1./(2.*dx)
  r2dy0=1./(2.*dy)
  do k=1,nzm
    f2lediss_z(k)=0.
    kc=min(nzm,k+1)
    kb=max(1,k-1)
    r2dz=2./((kc-kb)*(adzw(k+1)+adzw(k))*dz)
    r2dx=r2dx0*sqrt((kc-kb)*dx*r2dz) ! grid anisotropy correction
    r2dy=r2dy0*sqrt((kc-kb)*dx*r2dz)
    f2lediss(k)=0.
    do j=1,ny
     jc=j+YES3D
     jb=j-YES3D
     do i=1,nx
      f2lediss(k)=f2lediss(k)-tkh2(i,j,k)*( &
                       ((f(i,j,kc)-f0(kc)-f(i,j,kb)+f0(kb))*r2dz)**2 )
     end do
    end do
    f2lediss(k)=f2lediss(k)*2.*factor_xy
  end do

endif

call t_stopf ('fluxes_scalars_z')

end subroutine fluxes_scalar_z
