! scalar flux budget: 	

	
subroutine stat_sw1(du,twleproc,qwleproc)

use vars

implicit none
real du(nx,ny,nz,3)
real twleproc(nzm), qwleproc(nzm), swleproc(nzm)
integer i,j,k	
do j=1,ny
 do i=1,nx
  du(i,j,nz,3)=0.
 end do
end do
do k=1,nzm
 do j=1,ny
  do i=1,nx
   twleproc(k)=twleproc(k)+0.5*(t(i,j,k)-t0(k))*(du(i,j,k,3)+du(i,j,k+1,3))
#ifdef PNNL_STATS
   qwleproc(k)=qwleproc(k)+0.5*(qv(i,j,k)+qcl(i,j,k)-(q0(k)-qi0(k)))*(du(i,j,k,3)+du(i,j,k+1,3)) !MWSWong: qv+qc
#else
   qwleproc(k)=qwleproc(k)+0.5*(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)-q0(k))*(du(i,j,k,3)+du(i,j,k+1,3))
#endif
  end do
 end do
end do

end
#ifdef PNNL_STATS
!MWSWong: for thel
subroutine stat_sw0(du,thelwleproc,qtogwleproc)
use vars
implicit none
real du(nx,ny,nz,3)
real thelwleproc(nzm), qtogwleproc(nzm)
integer i,j,k   
do j=1,ny
 do i=1,nx
  du(i,j,nz,3)=0.
 end do
end do
do k=1,nzm
 do j=1,ny
  do i=1,nx
   thelwleproc(k)=thelwleproc(k)+0.5*(thel(i,j,k)-thel0(k))*(du(i,j,k,3)+du(i,j,k+1,3))
   qtogwleproc(k)=qtogwleproc(k)+0.5*(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)+qpl(i,j,k)+qpi(i,j,k) &
                                       -q0(k)-qp0(k))*(du(i,j,k,3)+du(i,j,k+1,3))
  end do
 end do
end do
end
#endif /*PNNL_STATS*/
