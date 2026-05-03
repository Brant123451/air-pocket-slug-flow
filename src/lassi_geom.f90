!======================================================================
! lassi_geom.f90
! Geometric closures for a FLAT-TOP D-SHAPED pipe cross-section and IKH
! coefficient κ.  Reference: LASSI_ALGORITHM_SPEC.md sections 1, 8.3,
! 13, and 14 (the flat-top modification introduced to suppress the
! singular U_F = (Um − β·U_l)/(1 − β) as β→1; see discussion with
! user, May 2026).
!
! Physical geometry
! -----------------
! The nominal pipe diameter is D.  The cross-section is the lower
! portion of a full circle of diameter D, cut by a HORIZONTAL PLATE
! at height h = H_MAX_FRAC · D (default 0.95 D).  So the physical
! pipe is D-shaped: a round bottom welded to a flat roof.  No liquid
! can exist above the roof – β = 1 corresponds to h = H_MAX_FRAC · D
! (not h = D as in the classical circular pipe).
!
! Half-angle parametrisation (same γ as before)
! ---------------------------------------------
! γ ∈ [0, γ_max]           ← γ_max = arccos(1 − 2·H_MAX_FRAC) = 2.690566
! h_l(γ) = R (1 − cos γ)   ← liquid depth from the pipe bottom
! A_l(γ) = R² (γ − sin γ cos γ)   ← cross-section area of the liquid
!                                   (this formula is a geometric
!                                   property of the LOWER CIRCULAR
!                                   SEGMENT and is unaffected by
!                                   whether the pipe is circular or
!                                   flat-top, as long as h ≤ h_max).
!
! Convention
! ----------
!   β = A_l / A_pipe_new      with A_pipe_new = R² · AREA_GAMMA_FACTOR
!   AREA_GAMMA_FACTOR = γ_max − sin γ_max cos γ_max ≈ 3.082867
!   (compare π ≈ 3.141593 for the full circle)
!
! Perimeters (per cross section)
! ------------------------------
!   S_l = D γ                                    ← unchanged (wetted arc)
!   S_g = D (γ_max − γ + sin γ_max)              ← gas arc (γ_max−γ part)
!                                                    + flat roof (sin γ_max)
!   S_i = D sin γ                                ← interface chord
!                                                   (γ < γ_max)
!
! d(A_l)/dh_l = S_i = D sin γ                    ← unchanged identity
! d(α_l)/dh_l = S_i / A_pipe_new
!             = D sin γ / (R² · AREA_GAMMA_FACTOR)
!             = 4 sin γ / (D · AREA_GAMMA_FACTOR)
!======================================================================
module lassi_geom
   use lassi_kinds
   implicit none
   private
   public :: gamma_from_beta, hl_from_gamma, perimeters, dalpha_dhl, ikh_kappa
   public :: pipe_area
   public :: H_MAX_FRAC, GAMMA_MAX, SIN_GAMMA_MAX, AREA_GAMMA_FACTOR

   ! -------------------------------------------------------------------
   ! Flat-top cut parameters.  H_MAX_FRAC = 0.95 puts the roof at
   ! h = 0.95·D, which corresponds to the old circular β ≈ 0.98131
   ! (i.e. exactly the Hooke slug-onset threshold β_init).  Flattening
   ! the roof at this height removes the 0.98..1.00 β band entirely,
   ! so IKH wave crests that would previously have pushed a section's
   ! β to 0.97 (one δt before slug_init) now either (a) hit β ≥ 1 in
   ! the same event and trigger Case A/B immediately, or (b) stay
   ! below 1 with enough headroom (β ≤ 1) that U_F = (Um − β·Ul)/(1−β)
   ! cannot diverge.
   real(rk), parameter :: H_MAX_FRAC        = 0.95_rk
   ! γ_max = arccos(1 − 2·H_MAX_FRAC) = arccos(−0.9)
   real(rk), parameter :: GAMMA_MAX         = 2.6905658417937937_rk
   real(rk), parameter :: SIN_GAMMA_MAX     = 0.4358898943540674_rk
   real(rk), parameter :: COS_GAMMA_MAX     = -0.9_rk
   ! AREA_GAMMA_FACTOR = γ_max − sin γ_max · cos γ_max
   !                  = 2.6905658 + 0.3923009 = 3.0828667
   real(rk), parameter :: AREA_GAMMA_FACTOR = 3.0828667361547884_rk

contains

   ! ---------------------------------------------------------------- !
   ! Pipe cross-section area.  FLAT-TOP D-SHAPE:
   !    A_pipe = R² · AREA_GAMMA_FACTOR
   !           = (D/2)² · AREA_GAMMA_FACTOR
   !           = 0.25 · AREA_GAMMA_FACTOR · D²
   !           ≈ 0.77072 · D²         (vs 0.78540·D² for a full circle)
   pure function pipe_area(D) result(A)
      real(rk), intent(in) :: D
      real(rk) :: A
      A = 0.25_rk*AREA_GAMMA_FACTOR*D*D
   end function pipe_area

   ! ---------------------------------------------------------------- !
   ! Invert  β = (γ − sin γ cos γ) / AREA_GAMMA_FACTOR  by bisection.
   ! Robust: monotone in γ on [0, γ_max].
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
         gamma = GAMMA_MAX; return
      end if
      lo = 0.0_rk
      hi = GAMMA_MAX
      do it = 1, 80
         mid = 0.5_rk*(lo + hi)
         fmid = (mid - sin(mid)*cos(mid))/AREA_GAMMA_FACTOR - b
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
   ! Liquid depth as a function of γ.  SAME formula as a circular
   ! pipe, because the liquid segment is the lower circular segment
   ! regardless of whether the roof is round or flat.
   !    h_l = R · (1 − cos γ) = (D/2)(1 − cos γ)
   pure function hl_from_gamma(gamma, D) result(hl)
      real(rk), intent(in) :: gamma, D
      real(rk) :: hl
      hl = 0.5_rk*D*(1.0_rk - cos(gamma))
   end function hl_from_gamma

   ! ---------------------------------------------------------------- !
   ! Wetted perimeters  (per cross section), FLAT-TOP D-SHAPE:
   !   S_l = liquid–wall arc      = D γ                 (unchanged)
   !   S_g = gas–wall perimeter   = D (γ_max − γ)          arc
   !                              + D sin γ_max            flat roof
   !                              = D (γ_max − γ + sin γ_max)
   !   S_i = interface chord      = D sin γ              (γ < γ_max)
   !
   ! Edge case: when γ → γ_max (β → 1) the interface coincides with the
   ! roof itself; we keep the circular formula D sin γ, which smoothly
   ! approaches D sin γ_max.  Physically this is the limit "liquid
   ! touches the ceiling everywhere" so slug_init will absorb the cell
   ! in the next list-management pass.
   pure subroutine perimeters(gamma, D, Sl, Sg, Si)
      real(rk), intent(in)  :: gamma, D
      real(rk), intent(out) :: Sl, Sg, Si
      Sl = D*gamma
      Sg = D*(GAMMA_MAX - gamma + SIN_GAMMA_MAX)
      Si = D*sin(gamma)
   end subroutine perimeters

   ! ---------------------------------------------------------------- !
   ! d(α_l)/d(h_l) = S_i / A_pipe_new
   !               = (D sin γ) / (R² · AREA_GAMMA_FACTOR)
   !               = (D sin γ) / (D²/4 · AREA_GAMMA_FACTOR)
   !               = 4 sin γ / (D · AREA_GAMMA_FACTOR)
   ! 1.9 % larger than the old circular value 4 sin γ/(π D).
   pure function dalpha_dhl(gamma, D) result(deriv)
      real(rk), intent(in) :: gamma, D
      real(rk) :: deriv
      deriv = 4.0_rk*sin(gamma)/(AREA_GAMMA_FACTOR*D)
   end function dalpha_dhl

   ! ---------------------------------------------------------------- !
   ! Inviscid Kelvin–Helmholtz coefficient κ (Eq. 3.9 in spec):
   !   κ = ((ρ_l-ρ_g)/ρ_l) g cos φ / (dα_l/dh_l)
   !       − (1/α) (ρ_g/ρ_l) (U_g − U_l)^2
   ! Uses the NEW dalpha_dhl (flat-top version) — the formula is
   ! otherwise identical; only the geometric factor in the gravity
   ! stabilisation term shifts by 1.9 %.
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
