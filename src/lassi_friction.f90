!======================================================================
! lassi_friction.f90
! Darcy-Weisbach friction closures, body force F, Bendiksen nose vel.
! Reference: LASSI_ALGORITHM_SPEC.md sections 8.1, 8.2.
!======================================================================
module lassi_friction
   use lassi_kinds
   use lassi_geom
   implicit none
   private
   public :: friction_params_t
   public :: dw_lambda, body_force_F, dF_dUl, bendiksen_nose, wake_effect, u_crit_balance
   public :: solve_equilibrium_Ul

   type :: friction_params_t
      real(rk) :: mu_l    = 1.0e-3_rk       ! liquid viscosity [Pa s]
      real(rk) :: mu_g    = 1.8e-5_rk       ! gas viscosity
      real(rk) :: rough   = 0.0_rk          ! pipe roughness [m]
      real(rk) :: ai_factor= 1.0_rk         ! Andritsos-Hanratty factor multiplier
      logical  :: use_AH = .false.          ! use Andritsos-Hanratty for λ_i  (LASSI thesis default: λ_i = λ_g, page 19)
   end type friction_params_t

contains

   ! ---------------------------------------------------------------- !
   ! Haaland approximation to Colebrook:
   !   1/√λ = -1.8 log10[(ε/D / 3.7)^1.11 + 6.9/Re]
   pure function dw_lambda(Re, eps_over_D) result(lam)
      real(rk), intent(in) :: Re, eps_over_D
      real(rk) :: lam, q
      if (Re < 1.0e-6_rk) then
         lam = 0.0_rk
         return
      end if
      if (Re < 2300.0_rk) then
         lam = 64.0_rk/max(Re, 1.0_rk)        ! laminar
      else
         q = -1.8_rk*log10( (eps_over_D/3.7_rk)**1.11_rk + 6.9_rk/Re )
         lam = 1.0_rk/(q*q)
      end if
   end function dw_lambda

   ! ---------------------------------------------------------------- !
   ! Body force per unit volume of liquid in the LASSI liquid-momentum
   ! equation (M3) — Eq. 3.12 in the thesis (=(3.11) in the SPEC).
   ! F(β, U_l, U_g^S)  is the net non-pressure body+friction force
   ! per unit liquid volume; it is what remains on the RHS of (M3)
   ! after substituting Eq. 3.5 into Eq. 3.3 (see SPEC §2.2-2.3).
   ! Returns F  [N/m^3] so that the right-hand side of (M3) is (β/ρ_l) F.
   pure function body_force_F(beta, Ul, Ug, D, phi, rho_l, rho_g, fp) result(F)
      real(rk),                 intent(in) :: beta, Ul, Ug, D, phi, rho_l, rho_g
      type(friction_params_t),  intent(in) :: fp
      real(rk) :: F

      real(rk) :: gamma, Sl, Sg, Si, A, alpha
      real(rk) :: Dl, Dg, Re_l, Re_g, lam_l, lam_g, lam_i
      real(rk) :: tau_l, tau_g, tau_i

      gamma = gamma_from_beta(beta)
      call perimeters(gamma, D, Sl, Sg, Si)
      A     = pipe_area(D)
      alpha = max(1.0_rk - beta, EPS_SMALL)

      ! hydraulic diameters (no zero-division)
      Dl = 4.0_rk*beta*A / max(Sl, EPS_SMALL)
      Dg = 4.0_rk*alpha*A / max(Sg + Si, EPS_SMALL)

      ! Reynolds and friction factors
      Re_l = rho_l*abs(Ul)*Dl / max(fp%mu_l, EPS_SMALL)
      Re_g = rho_g*abs(Ug - Ul)*Dg / max(fp%mu_g, EPS_SMALL)
      lam_l = dw_lambda(Re_l, fp%rough/D)
      lam_g = dw_lambda(Re_g, fp%rough/D)
      lam_i = lam_g
      if (fp%use_AH) lam_i = lam_g*(1.0_rk + 75.0_rk*beta)*fp%ai_factor

      ! shear stresses (Darcy convention: τ = (λ/8) ρ |U| U, thesis Eq. 2.1–2.2)
      tau_l = 0.125_rk*rho_l*lam_l*abs(Ul)*Ul
      tau_g = 0.125_rk*rho_g*lam_g*abs(Ug)*Ug
      tau_i = 0.125_rk*rho_g*lam_i*abs(Ug - Ul)*(Ug - Ul)

      ! F: liquid wall + interface (gas+liquid sides) - buoyancy-corrected gravity - gas wall feedback
      F = -tau_l*Sl/A &
          + tau_i*Si*(1.0_rk/(alpha*A) + 1.0_rk/A) &
          - ((rho_l - rho_g)/rho_l)*rho_l*beta*G_ACC*sin(phi) &
          - rho_g*Sg*tau_g/(alpha*A)
   end function body_force_F

   ! ---------------------------------------------------------------- !
   ! Analytic Jacobian ∂F/∂U_l of the body-force F above, used by the
   ! linearised liquid-momentum implicit step (Eq. 3.20).
   ! Differentiating term by term (treating λ as Re-independent — the
   ! dλ/dRe contribution is < 1% for fully turbulent flow):
   !   ∂τ_l/∂U_l = (λ_l/4)·ρ_l·|U_l|        (since τ_l = (λ/8)ρ|U|U  → ∂/∂U = (λ/4)ρ|U|)
   !   ∂τ_i/∂U_l = -(λ_i/4)·ρ_g·|U_g-U_l|   (sign flip because ∂(U_g-U_l)/∂U_l = -1)
   !   gravity term: ∂/∂U_l = 0
   !   gas-wall term:  ∂/∂U_l = 0  (gas momentum is decoupled in the M3 reduction)
   ! Result is always ≤ 0, so the implicit denom = 1 - (δt/ρ_l)·∂F/∂U_l > 1
   ! is well-conditioned.
   pure function dF_dUl(beta, Ul, Ug, D, phi, rho_l, rho_g, fp) result(dF)
      real(rk),                 intent(in) :: beta, Ul, Ug, D, phi, rho_l, rho_g
      type(friction_params_t),  intent(in) :: fp
      real(rk) :: dF
      real(rk) :: gamma, Sl, Sg, Si, A, alpha, Dl, Dg, Re_l, Re_g, lam_l, lam_i
      real(rk) :: dtau_l_dU, dtau_i_dU, phi_unused
      phi_unused = phi   ! reserved for future centrifugal terms (no-op)
      gamma = gamma_from_beta(beta)
      call perimeters(gamma, D, Sl, Sg, Si)
      A     = pipe_area(D)
      alpha = max(1.0_rk - beta, EPS_SMALL)
      Dl = 4.0_rk*beta*A / max(Sl, EPS_SMALL)
      Dg = 4.0_rk*alpha*A / max(Sg + Si, EPS_SMALL)
      Re_l = rho_l*abs(Ul)*Dl / max(fp%mu_l, EPS_SMALL)
      Re_g = rho_g*abs(Ug - Ul)*Dg / max(fp%mu_g, EPS_SMALL)
      lam_l = dw_lambda(Re_l, fp%rough/D)
      lam_i = dw_lambda(Re_g, fp%rough/D)
      if (fp%use_AH) lam_i = lam_i*(1.0_rk + 75.0_rk*beta)*fp%ai_factor
      ! ∂τ_l/∂U_l and ∂τ_i/∂U_l (analytic, treating λ as constant in U_l)
      dtau_l_dU =  0.25_rk*lam_l*rho_l*abs(Ul)
      dtau_i_dU = -0.25_rk*lam_i*rho_g*abs(Ug - Ul)
      ! ∂F/∂U_l: contributions from -τ_l·Sl/A and +τ_i·Si·(1/(αA)+1/A)
      dF = -dtau_l_dU*Sl/A + dtau_i_dU*Si*(1.0_rk/(alpha*A) + 1.0_rk/A)
   end function dF_dUl

   ! ---------------------------------------------------------------- !
   ! Thesis page 56 rule: when β < 0.001 the linearised implicit step
   ! Eq. 3.20 becomes unstable (F is extremely sensitive to U_l via the
   ! Sl/A and 1/α wall/interface factors at low holdup).  The thesis
   ! prescribes replacing the time-march by a direct root-find:
   !
   !    "U_l is simply set to its equilibrium value, i.e. to U*_l
   !     ensuring F(β^n_J, U*_l, U^{n+1}_{m,J}) = 0"
   !
   ! We solve F(U_l) = 0 by Newton iteration using the analytic
   ! Jacobian dF/dU_l (which we already have for Eq. 3.20), clipped to
   ! ±U_MAX to keep the iterate bounded during any local non-monotone
   ! excursion (body_force_F is monotonically decreasing in U_l for
   ! normal parameter ranges because both the wall-shear and
   ! interfacial-shear contributions have negative dF/dU_l).
   pure function solve_equilibrium_Ul(beta, Ug, D, phi, rho_l, rho_g, fp, Ul_init) result(Ul_star)
      real(rk),                 intent(in) :: beta, Ug, D, phi, rho_l, rho_g, Ul_init
      type(friction_params_t),  intent(in) :: fp
      real(rk) :: Ul_star
      real(rk) :: F_val, dF_val, Ul_prev
      integer  :: iter
      integer,  parameter :: MAX_ITER = 30
      real(rk), parameter :: TOL_F    = 1.0e-6_rk
      real(rk), parameter :: TOL_U    = 1.0e-6_rk
      real(rk), parameter :: U_MAX    = 20.0_rk
      real(rk), parameter :: DFDU_MIN = 1.0e-12_rk
      Ul_star = max(min(Ul_init, U_MAX), -U_MAX)
      do iter = 1, MAX_ITER
         F_val = body_force_F(beta, Ul_star, Ug, D, phi, rho_l, rho_g, fp)
         if (abs(F_val) < TOL_F) exit
         dF_val = dF_dUl(beta, Ul_star, Ug, D, phi, rho_l, rho_g, fp)
         if (abs(dF_val) < DFDU_MIN) exit           ! degenerate Jacobian
         Ul_prev = Ul_star
         Ul_star = Ul_star - F_val/dF_val
         if (Ul_star >  U_MAX) Ul_star =  U_MAX
         if (Ul_star < -U_MAX) Ul_star = -U_MAX
         if (abs(Ul_star - Ul_prev) < TOL_U) exit
      end do
   end function solve_equilibrium_Ul

   ! ---------------------------------------------------------------- !
   ! Bendiksen (1984) bubble nose velocity in inclined pipes.
   ! LASSI thesis Eq. 2.10–2.11 (paper PNG `48.png`):
   !   low Fr (<3.5):  C_01 = 1.05 + 0.15 sin²φ,  v_01 = (0.35 sinφ + 0.54 cosφ)√(gD)
   !   high Fr (>3.5): C_02 = 1.20 + 0.15 sin²φ,  v_02 = 0.35 sinφ √(gD)
   !   U_b_long      = max(C_01·Um + v_01,  C_02·Um + v_02)
   ! cosφ / sinφ are OUTSIDE the square root (the SPEC v0 had cos/sin INSIDE √, which
   ! agrees only at φ∈{0, π/2}).  abs(sinφ) is used so down-pipe (φ<0) gives the same drift.
   pure function bendiksen_nose(Um, D, phi) result(Un)
      real(rk), intent(in) :: Um, D, phi
      real(rk) :: Un, sqgD, sin2, sphi, cphi, U1, U2, C01, v01, C02, v02
      sqgD = sqrt(G_ACC*D)
      sphi = sin(phi); cphi = cos(phi)
      sin2 = sphi*sphi
      C01 = 1.05_rk + 0.15_rk*sin2
      v01 = (0.35_rk*abs(sphi) + 0.54_rk*max(cphi, 0.0_rk))*sqgD
      C02 = 1.20_rk + 0.15_rk*sin2
      v02 = 0.35_rk*abs(sphi)*sqgD
      U1 = C01*Um + v01
      U2 = C02*Um + v02
      Un = max(U1, U2)
   end function bendiksen_nose

   ! ---------------------------------------------------------------- !
   ! Wake-effect multiplier (Moissis & Griffith 1962), thesis Eq. 2.14:
   !   W_eff = min( 1 + 8 exp(-1.06 L_s/D),  W_eff_cap )
   ! L_s : current slug length, D : pipe diameter, W_eff_cap : user-set ceiling.
   ! Thesis (page 24, paper PNG `51.png`) recommends a cap of 2.0 (the Moissis–Griffith
   ! correlation has not been validated for tiny slugs where the un-capped formula
   ! would give W_eff up to 9).  Setting W_eff_cap = 1.0 disables the wake effect.
   pure function wake_effect(L_s, D, W_eff_cap) result(W_eff)
      real(rk), intent(in) :: L_s, D, W_eff_cap
      real(rk) :: W_eff, W_raw, ratio
      if (D <= 0.0_rk) then
         W_eff = max(1.0_rk, W_eff_cap)
         return
      end if
      ratio = max(L_s, 0.0_rk)/D
      W_raw = 1.0_rk + 8.0_rk*exp(-1.06_rk*ratio)
      W_eff = min(W_raw, max(1.0_rk, W_eff_cap))
      W_eff = max(W_eff, 1.0_rk)
   end function wake_effect

   ! ---------------------------------------------------------------- !
   ! Critical mixture velocity inside a slug (thesis Eq. 2.12, page 21):
   !   ½ λ_l ρ_l U_crit² S_l = ρ_l g sin(φ) A
   ! ⇒ U_crit = sign(sin φ) · √( 2 g |sin φ| A / (λ_l S_l) )
   !
   ! For a fully-liquid slug β=1, S_l = π D, A = π D² / 4, so
   !   U_crit = sign(sinφ) · √( g |sinφ| D / (2 λ_l) ).
   ! Horizontal pipe (sinφ=0): U_crit=0.
   ! Sign convention follows the thesis: down-pipe (φ>0 in thesis,
   ! sinφ>0 here) gives U_crit > 0 (gravity drives the slug forward).
   pure function u_crit_balance(D, phi, lam_l) result(U_crit)
      real(rk), intent(in) :: D, phi, lam_l
      real(rk) :: U_crit, sphi, val
      sphi = sin(phi)
      if (D <= 0.0_rk .or. lam_l <= EPS_SMALL .or. abs(sphi) < EPS_SMALL) then
         U_crit = 0.0_rk
         return
      end if
      val = G_ACC*abs(sphi)*D / (2.0_rk*lam_l)
      U_crit = sign(sqrt(max(val, 0.0_rk)), sphi)
   end function u_crit_balance

end module lassi_friction
