!======================================================================
! lassi_kinds.f90
! Floating-point kind, mathematical/physical constants, and global tags.
! Reference: LASSI_ALGORITHM_SPEC.md sections 1, 13.
!======================================================================
module lassi_kinds
   use, intrinsic :: iso_fortran_env, only: real64, int32, int64
   implicit none
   public

   ! --- floating-point precision ---
   integer, parameter :: rk = real64
   integer, parameter :: ik = int32

   ! --- mathematical / physical constants ---
   real(rk), parameter :: PI    = 3.141592653589793238462643_rk
   real(rk), parameter :: TWOPI = 2.0_rk*PI
   real(rk), parameter :: G_ACC = 9.81_rk     ! gravitational acceleration [m/s^2]

   ! --- numerical safeguards ---
   real(rk), parameter :: EPS_TINY = 1.0e-14_rk
   real(rk), parameter :: EPS_SMALL= 1.0e-9_rk
   real(rk), parameter :: BETA_DRY = 1.0e-12_rk   ! treat |β|<this as zero (Void)
   real(rk), parameter :: BETA_FULL_TOL = 1.0e-12_rk ! treat 1-β<this as full

   ! --- object-kind tags (used in TYPE :: object_t) ---
   integer(ik), parameter :: KIND_INLET   = 1
   integer(ik), parameter :: KIND_OUTLET  = 2
   integer(ik), parameter :: KIND_BUBBLE  = 3
   integer(ik), parameter :: KIND_SECTION = 4
   integer(ik), parameter :: KIND_SLUG    = 5

   ! --- Riemann case identifiers (Tables A.1-A.10) ---
   integer(ik), parameter :: CASE_RR    = 1   ! Rarefaction-Rarefaction
   integer(ik), parameter :: CASE_RS    = 2   ! Rarefaction-Shock
   integer(ik), parameter :: CASE_SR    = 3   ! Shock-Rarefaction
   integer(ik), parameter :: CASE_SS    = 4   ! Shock-Shock
   integer(ik), parameter :: CASE_SSS   = 5   ! Saturated Shock-Shock (β_M=1)
   integer(ik), parameter :: CASE_RV    = 6   ! Rarefaction-Void  (right dry)
   integer(ik), parameter :: CASE_VR    = 7   ! Void-Rarefaction  (left dry)
   integer(ik), parameter :: CASE_RRADB = 8   ! RR with appearing dry bed
   integer(ik), parameter :: CASE_SECSLUG = 9 ! Section-Slug nose
   integer(ik), parameter :: CASE_SLUGSEC =10 ! Slug-Section nose

end module lassi_kinds
