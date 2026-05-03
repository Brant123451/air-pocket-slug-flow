!======================================================================
! lassi_geom.f90
! Circular pipe geometric closures and IKH coefficient κ.
! Reference: LASSI_ALGORITHM_SPEC.md sections 1, 8.3, 13.
!
! Convention:
!   β = A_l / A   (liquid holdup, 0..1)
!   γ ∈ [0, π]    (wetted half-angle)  with  β = (γ - sin γ cos γ)/π
!   h_l = R(1 - cos γ)
!   S_l = D γ ;  S_g = D(π-γ) ;  S_i = D sin γ
!   dA_l/dh_l = S_i = D sin γ                                  (LASSI thesis App. D, Table D.1)
!   dα_l/dh_l = (dA_l/dh_l)/A_pipe = 4 sin γ /(π D) = 2 sin γ /(π R)
!   dh_l/dα_l = π D /(4 sin γ) = π R /(2 sin γ)               (inverse, used in κ recovery term)
!======================================================================
module lassi_geom
   use lassi_kinds
   implicit none
   private
   public :: gamma_from_beta, hl_from_gamma, perimeters, dalpha_dhl, ikh_kappa
   public :: pipe_area

contains

   ! ---------------------------------------------------------------- !
   pure function pipe_area(D) result(A)
      real(rk), intent(in) :: D
      real(rk) :: A
      A = 0.25_rk*PI*D*D
   end function pipe_area

   ! ---------------------------------------------------------------- !
   ! Invert  β = (γ - sin γ cos γ)/π  by bisection.
   ! Robust: monotone in γ on [0,π], analytic derivative also available.
   pure function gamma_from_beta(beta) result(gamma)
      real(rk), intent(in) :: beta
      real(rk) :: gamma
      real(rk) :: lo, hi, mid, fmid, b
      integer  :: it
      b = max(0.0_rk, min(1.0_rk, beta))
      if (b <= EPS_TINY) then
         gamma = 0.0_rk; return
      end if
      if (b >= 1.0_rk - EPS_TINY) then
         gamma = PI; return
      end if
      lo = 0.0_rk
      hi = PI
      do it = 1, 80
         mid = 0.5_rk*(lo + hi)
         fmid = (mid - sin(mid)*cos(mid))/PI - b
         if (fmid > 0.0_rk) then
            hi = mid
         else
            lo = mid
         end if
         if (hi - lo < 1.0e-13_rk) exit
      end do
      gamma = 0.5_rk*(lo + hi)
   end function gamma_from_beta

   ! ---------------------------------------------------------------- !
   pure function hl_from_gamma(gamma, D) result(hl)
      real(rk), intent(in) :: gamma, D
      real(rk) :: hl
      hl = 0.5_rk*D*(1.0_rk - cos(gamma))
   end function hl_from_gamma

   ! ---------------------------------------------------------------- !
   ! Wetted perimeters  (per cross section):
   !   S_l = liquid–wall ;  S_g = gas–wall ;  S_i = interfacial width
   pure subroutine perimeters(gamma, D, Sl, Sg, Si)
      real(rk), intent(in)  :: gamma, D
      real(rk), intent(out) :: Sl, Sg, Si
      Sl = D*gamma
      Sg = D*(PI - gamma)
      Si = D*sin(gamma)
   end subroutine perimeters

   ! ---------------------------------------------------------------- !
   ! d(α_l)/d(h_l) = (S_i)/A_pipe = D sin γ /(π D²/4) = 4 sin γ /(π D)
   !              = 2 sin γ /(π R)                                  with R=D/2
   ! Reference: LASSI thesis Appendix D, Table D.1 — dA_l/dh_l = 2R sin θ.
   pure function dalpha_dhl(gamma, D) result(deriv)
      real(rk), intent(in) :: gamma, D
      real(rk) :: deriv
      deriv = 4.0_rk*sin(gamma)/(PI*D)   ! 2 sin γ /(π R) with R=D/2
   end function dalpha_dhl

   ! ---------------------------------------------------------------- !
   ! Inviscid Kelvin–Helmholtz coefficient κ (Eq. 3.9 in spec):
   !   κ = ((ρ_l-ρ_g)/ρ_l) g cos φ · A/(dA_l/dh_l)
   !       − (1/α) (ρ_g/ρ_l)(U_g - U_l)^2
   ! Returns κ; well-posedness ⇔ κ > 0.
   pure function ikh_kappa(beta, D, phi, rho_l, rho_g, Ul, Ug) result(kappa)
      real(rk), intent(in) :: beta, D, phi, rho_l, rho_g, Ul, Ug
      real(rk) :: kappa
      real(rk) :: gamma, alpha, inv_dalpha_dhl
      gamma = gamma_from_beta(beta)
      alpha = max(1.0_rk - beta, EPS_SMALL)
      inv_dalpha_dhl = 1.0_rk/max(dalpha_dhl(gamma, D), EPS_SMALL)
      kappa = ((rho_l - rho_g)/rho_l)*G_ACC*cos(phi)*inv_dalpha_dhl &
              - (1.0_rk/alpha)*(rho_g/rho_l)*(Ug - Ul)**2
   end function ikh_kappa

end module lassi_geom
