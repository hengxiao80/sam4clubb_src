subroutine header()
use grid, only: version, version_date
integer          :: values(8)
character        :: date*8, time*10, zone*5
character(len=8) :: cdate          ! System date
character(len=8) :: ctime          ! System time

call date_and_time (date, time, zone, values)
  cdate(1:2) = date(5:6)
  cdate(3:3) = '/'
  cdate(4:5) = date(7:8)
  cdate(6:6) = '/'
  cdate(7:8) = date(3:4)
  ctime(1:2) = time(1:2)
  ctime(3:3) = ':'
  ctime(4:5) = time(3:4)
  ctime(6:6) = ':'
  ctime(7:8) = time(5:6)


write(*,*) '          **************************************************'
write(*,*) '          **************************************************'
write(*,*) '          *                                                *'
write(*,*) '          *      *****          **       ***          **   *'
write(*,*) '          *    ***   ***       ****      ****       ****   *'
write(*,*) '          *   ***     ***     *** ***    *****     *****   *'
write(*,*) '          *    ***           ***   ***   *** **   ** ***   *'
write(*,*) '          *       ***       ***     ***  ***  ** **  ***   *'
write(*,*) '          *          ***    ***     ***  ***   ***   ***   *'
write(*,*) '          *   ***     ***   ***********  ***         ***   *'
write(*,*) '          *    ***   ***    ***     ***  ***         ***   *'
write(*,*) '          *      *****      ***     ***  ***         ***   *'
write(*,*) '          *                                                *'
write(*,*) '          **************************************************'
write(*,*) '          **************************************************'
write(*,*) '          ***       System for Atmospheric Modeling      ***'
write(*,*) '          ***                    SAM                     ***'
write(*,*) '                    Version '//version//' ('//version_date//')  '  
write(*,*) '          **************************************************'
write(*,*) '          ***     (C) Marat Khairoutdinov                ***'
write(*,*) '          *** School of Marine and Atmospheric Sciences  ***'
write(*,*) '          ***       Stony Brook University               ***'
write(*,*) '          **************************************************'
write(*,*) '          ***       The model can be used only with      ***'
write(*,*) '          ***         permission from the author!        ***'
write(*,*) '          ***       No transfer to any third party       ***'
write(*,*) '          ***                is allowed!                 ***'
#ifdef CLUBB
write(*,*) '          **************************************************'
write(*,*) '          ***           DOUBLE-MOMENT VERSION OF         ***'
write(*,*) '          ***           LIN MICROPHYSICS SCHEME          ***'
write(*,*) '          ***   by Vaughan Phillips (GFDL, 2006)         ***'
write(*,*) '          ***  (developed from single-moment version     ***'
write(*,*) '          ***   by Lord et al., 1984, JAS, vol 41, 2836, ***'
write(*,*) '          ***   and Krueger et al, 1995, JAM, vol 34, 281)  '
write(*,*) '          ***  Adapted by Mikhail Ovtchinnikov (PNNL, 2008) '
write(*,*) '          **************************************************'
write(*,*) '          ***       Portions of this code                ***'
write(*,*) '          ***        (C) UW Milwaukee 2008-2014          ***'
write(*,*) '          ***    Developed at the UWM Dept of            ***'
write(*,*) '          ***    Mathematical Sciences under the         ***'
write(*,*) '          ***    direction of Vincent E. Larson          ***'
write(*,*) '          ***             vlarson@uwm.edu                ***'
write(*,*) '          **************************************************'
#endif
write(*,*) '          **************************************************'
write(*,*) '          **************************************************'
write(*,*) '          ***       DATE '//cdate//' TIME '//ctime//'          ***'
write(*,*) '          **************************************************'
write(*,*)

return
end
