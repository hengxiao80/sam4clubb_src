
subroutine kurant

use vars
use sgs, only: kurant_sgs
use params, only: ncycle_max

implicit none

integer i, j, k
real wm(nz)  ! maximum vertical wind velocity
real uhm(nz) ! maximum horizontal wind velocity
real cfl

call t_startf ('kurant')

ncycle = 1
	
wm(nz)=0.
do k = 1,nzm
 wm(k) = maxval(abs(w(1:nx,1:ny,k)))
 uhm(k) = sqrt(maxval(u(1:nx,1:ny,k)**2+YES3D*v(1:nx,1:ny,k)**2))
end do
w_max=max(w_max,maxval(w(1:nx,1:ny,1:nz)))
u_max=max(u_max,maxval(uhm(1:nzm)))

cfl_adv = 0.
do k=1,nzm
  cfl_adv = max(cfl_adv,uhm(k)*dt*sqrt((1./dx)**2+YES3D*(1./dy)**2), &
                   max(wm(k),wm(k+1))*dt/(dz*adzw(k)) )
end do

call kurant_sgs(cfl_sgs)

if(dompi) then
  cfl = cfl_adv
  call task_max_real(cfl,cfl_adv,1)
  cfl = cfl_sgs
  call task_max_real(cfl,cfl_sgs,1)
end if
cfl = max(cfl_adv,cfl_sgs)
ncycle = max(1,ceiling(cfl/0.7))
if(ncycle.gt.ncycle_max) then
   if(masterproc) print *,'the number of cycles exceeded ', ncycle_max
   call task_abort()
end if

call t_stopf ('kurant')

end subroutine kurant	
