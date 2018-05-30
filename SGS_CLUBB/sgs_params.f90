!-------------------------------------------------------------------------------
module sgs_params

! Description:
!   Contains information for running the CLUBB scheme.  Needed to prevent
!   circular dependencies in modules.
!-------------------------------------------------------------------------------
  implicit none

  logical:: dosmagor           ! if true, then use Smagorinsky closure
  logical:: doclubb            ! Enabled the CLUBB parameterization (interactively)
  logical:: doclubb_sfc_fluxes ! Apply the surface fluxes within the CLUBB code rather than SAM
  logical:: doclubbnoninter    ! Enable the CLUBB parameterization (non-interactively)
  integer:: nclubb             ! SAM timesteps per CLUBB timestep

end module sgs_params
