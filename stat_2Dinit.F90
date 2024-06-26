! Initialize 2D output

subroutine stat_2Dinit()

use vars
implicit none

     prec_xy(:,:) = 0.
     shf_xy(:,:) = 0.
     lhf_xy(:,:) = 0.
     lwnt_xy(:,:) = 0.
     swnt_xy(:,:) = 0.
     lwntc_xy(:,:) = 0.
     swntc_xy(:,:) = 0.
     pw_xy(:,:) = 0.
     cw_xy(:,:) = 0.
     iw_xy(:,:) = 0.
     cld_xy(:,:) = 0.
     u200_xy(:,:) = 0.
     v200_xy(:,:) = 0.
     usfc_xy(:,:) = 0.
     vsfc_xy(:,:) = 0.
     w500_xy = 0.
     lwns_xy(:,:) = 0.
     swns_xy(:,:) = 0.
     solin_xy(:,:) = 0.
     lwnsc_xy(:,:) = 0.
     swnsc_xy(:,:) = 0.
     qocean_xy(:,:) = 0.
!===================================
! UW ADDITIONS: MOSTLY 2D STATISTICS

    !bloss: store initial profiles for computation of storage terms in budgets
#ifndef UWM_STATS
# /* +++mhwang the following lines cause issues with the storge term in budgets. */
# /* This has been done in statistics.f90. So commented out these four lines. */
     ustor(:) = u0(1:nzm)
     vstor(:) = v0(1:nzm)
     tstor(:) = t0(1:nzm)
     qstor(:) = q0(1:nzm)
#endif

     utendcor(:) = 0.
     vtendcor(:) = 0.

     psfc_xy(:,:) = 0.

     u850_xy(:,:) = 0.
     v850_xy(:,:) = 0.

     swvp_xy(:,:) = 0.

     cloudtopheight(:,:) = 0.
     echotopheight(:,:) = 0.
     cloudtoptemp(:,:) = 0.

! END UW ADDITIONS

end

