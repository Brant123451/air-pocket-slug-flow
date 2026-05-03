!======================================================================
! lassi_eos.f90
! Gas equation of state (isothermal ideal gas by default).
!======================================================================
module lassi_eos
   use lassi_kinds
   implicit none
   private
   public :: eos_params_t, set_isothermal_air, rho_g, drho_g_dp

   type :: eos_params_t
      real(rk) :: Rs    = 287.05_rk    ! specific gas constant [J/(kg·K)]
      real(rk) :: T     = 293.15_rk    ! temperature [K]
      ! cached:
      real(rk) :: inv_RsT = 1.0_rk/(287.05_rk*293.15_rk)
   end type eos_params_t

contains

   subroutine set_isothermal_air(eos, T)
      type(eos_params_t), intent(inout) :: eos
      real(rk),           intent(in)    :: T
      eos%Rs = 287.05_rk
      eos%T  = T
      eos%inv_RsT = 1.0_rk/(eos%Rs*eos%T)
   end subroutine set_isothermal_air

   pure function rho_g(eos, p) result(r)
      type(eos_params_t), intent(in) :: eos
      real(rk),           intent(in) :: p
      real(rk) :: r
      r = p*eos%inv_RsT
   end function rho_g

   pure function drho_g_dp(eos, p) result(d)
      type(eos_params_t), intent(in) :: eos
      real(rk),           intent(in) :: p
      real(rk) :: d
      d = eos%inv_RsT
   end function drho_g_dp

end module lassi_eos
