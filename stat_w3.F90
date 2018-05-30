#ifdef PNNL_STATS
!Heng Xiao w'w'w' budget
subroutine stat_w3(du,w3le)

use vars
implicit none
real du(nx,ny,nz,3)
real w3le(nz),coef
integer i,j,k
coef = 1./float(nx*ny)
do k=1,nz
 w3le(k)=0.
end do
do j=1,ny
 do i=1,nx
   du(i,j,nz,3) = 0.
 end do
end do
do k=1,nzm
 do j=1,ny
  do i=1,nx
   w3le(k)=w3le(k)+3.0*w(i,j,k)*w(i,j,k)*du(i,j,k,3)
  end do
 end do
 w3le(k)=w3le(k)*coef
end do
do k=1,nzm
  w3le(k)=0.5*(w3le(k)+w3le(k+1))
end do
end
#endif /* PNNL_STATS */



