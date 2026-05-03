!======================================================================
! lassi_riemann.f90
! Modified shallow-water Riemann solver covering all 10 LASSI cases
!     RR, RS, SR, SS, SSS, RV, VR, RRADB, Sec-Slug Nose, Slug-Sec Nose
! exactly as Tables A.1-A.10 of the thesis.
!
! All formulas are taken verbatim from LASSI_ALGORITHM_SPEC.md §5.4.4.
!
!  Riemann state stored in TYPE(riemann_state_t):
!     beta_M, U_M             intermediate state
!     U_LL, U_LR              left side of left wave / left side of contact
!     U_RL, U_RR              right side of contact / right side of right wave
!     U_b                     border tracking velocity
!     beta_ML, U_ML           left rarefaction-fan-averaged state
!     beta_MR, U_MR           right rarefaction-fan-averaged state
!     case_id                 CASE_RR, CASE_RS, ...  (lassi_kinds)
!======================================================================
module lassi_riemann
   use lassi_kinds
   use lassi_friction, only: bendiksen_nose, wake_effect
   implicit none
   private
   public :: riemann_state_t
   public :: solve_msw_riemann
   public :: solve_section_slug_nose, solve_slug_section_nose

   type :: riemann_state_t
      real(rk) :: beta_M, U_M
      real(rk) :: U_LL, U_LR, U_RL, U_RR
      real(rk) :: U_b
      real(rk) :: beta_ML, U_ML
      real(rk) :: beta_MR, U_MR
      integer(ik) :: case_id = 0
   end type riemann_state_t

contains

!======================================================================
! Public driver: section-section MSW Riemann solver
!   Inputs:  beta_L, U_L, beta_R, U_R, kappa  (κ averaged at the
!            interface; thesis Eq. 3.12: well-posed iff κ > 0).
!   Output:  full RIEMANN_STATE (Tables A.1-A.8)
!
!   If κ ≤ 0, use a tiny positive κ only to keep this local Riemann
!   calculation hyperbolic; slug initiation itself is handled later by
!   the β > β_init list-management rule.
!======================================================================
   subroutine solve_msw_riemann(beta_L, U_L, beta_R, U_R, kappa, st)
      real(rk),            intent(in)  :: beta_L, U_L, beta_R, U_R, kappa
      type(riemann_state_t), intent(out) :: st

      real(rk) :: kp, sqk
      real(rk) :: beta_M, U_M

      kp  = max(kappa, EPS_SMALL)
      sqk = sqrt(kp)

      ! 1) Dry-bed cases first --------------------------------------------------
      if (beta_R <= BETA_DRY .and. beta_L > BETA_DRY) then
         call fill_RV(beta_L, U_L, kp, sqk, st);  return
      end if
      if (beta_L <= BETA_DRY .and. beta_R > BETA_DRY) then
         call fill_VR(beta_R, U_R, kp, sqk, st);  return
      end if
      if (beta_L <= BETA_DRY .and. beta_R <= BETA_DRY) then
         ! both dry — nothing to do
         st%beta_M = 0.0_rk; st%U_M = 0.5_rk*(U_L + U_R)
         st%U_LL=U_L; st%U_LR=U_L; st%U_RL=U_R; st%U_RR=U_R
         st%U_b = 0.5_rk*(U_L+U_R)
         st%beta_ML=0.0_rk; st%U_ML=U_L
         st%beta_MR=0.0_rk; st%U_MR=U_R
         st%case_id = CASE_RRADB
         return
      end if

      ! 2) Appearing dry bed test (Eq. 3.25, 3.26) -----------------------------
      !    if  U_L + 2 sqrt(kβ_L) < U_R - 2 sqrt(kβ_R)  then β_M = 0
      if (U_L + 2.0_rk*sqk*sqrt(beta_L) < U_R - 2.0_rk*sqk*sqrt(beta_R)) then
         call fill_RRADB(beta_L, U_L, beta_R, U_R, kp, sqk, st);  return
      end if

      ! 3) Find β_M by 1D bisection on the difference of left and right curves -
      !    Curve_L(β_M):  U_L_curve(β_M; β_L, U_L)
      !       if β_M < β_L   ->  rarefaction:  U_M = U_L - 2 sqk (sqrt β_M - sqrt β_L)
      !       if β_M > β_L   ->  shock:        U_M = U_L - (1/√2) sqk (β_M-β_L) sqrt(1/β_M + 1/β_L)
      !    Curve_R(β_M):  same with subscript R, sign flipped (right wave).
      call find_beta_M_bisection(beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M)

      ! 4) SSS test ------------------------------------------------------------
      ! Bisection's upper bound is hi = 1-1e-6, so when β_M is pinned near hi
      ! the true root is past 1 ⇒ Saturated Shock-Shock.
      if (beta_M >= 1.0_rk - 1.0e-4_rk) then
         call fill_SSS(beta_L, U_L, beta_R, U_R, kp, sqk, st);  return
      end if

      ! 5) Dispatch RR / RS / SR / SS -----------------------------------------
      if (beta_M < beta_L .and. beta_M < beta_R) then
         call fill_RR(beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M, st)
      else if (beta_M < beta_L .and. beta_M >= beta_R) then
         call fill_RS(beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M, st)
      else if (beta_M >= beta_L .and. beta_M < beta_R) then
         call fill_SR(beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M, st)
      else
         call fill_SS(beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M, st)
      end if
   end subroutine solve_msw_riemann

!======================================================================
! find_beta_M_bisection
!   Robust 1D bisection on
!      f(β_M) = U_L_curve(β_M) - U_R_curve(β_M)
!   to find the unique intersection (Fig. 3.4).
!======================================================================
   subroutine find_beta_M_bisection(beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M)
      real(rk), intent(in)  :: beta_L, U_L, beta_R, U_R, kp, sqk
      real(rk), intent(out) :: beta_M, U_M

      real(rk) :: lo, hi, mid, fmid, ulc, urc, flo, fhi
      integer  :: it

      lo = max(BETA_DRY, 1.0e-6_rk)
      hi = 1.0_rk - 1.0e-6_rk

      ulc = uL_of_beta(lo,  beta_L, U_L, sqk)
      urc = uR_of_beta(lo,  beta_R, U_R, sqk)
      flo = ulc - urc
      ulc = uL_of_beta(hi,  beta_L, U_L, sqk)
      urc = uR_of_beta(hi,  beta_R, U_R, sqk)
      fhi = ulc - urc

      if (flo*fhi > 0.0_rk) then
         ! No root in (0,1) — caller will detect SSS or RRADB
         if (flo > 0.0_rk) then
            beta_M = hi
         else
            beta_M = lo
         end if
         U_M = 0.5_rk*( uL_of_beta(beta_M,beta_L,U_L,sqk) + uR_of_beta(beta_M,beta_R,U_R,sqk) )
         return
      end if

      do it = 1, 100
         mid  = 0.5_rk*(lo + hi)
         ulc  = uL_of_beta(mid, beta_L, U_L, sqk)
         urc  = uR_of_beta(mid, beta_R, U_R, sqk)
         fmid = ulc - urc
         if (fmid*flo < 0.0_rk) then
            hi = mid; fhi = fmid
         else
            lo = mid; flo = fmid
         end if
         if (hi - lo < 1.0e-12_rk) exit
      end do
      beta_M = 0.5_rk*(lo + hi)
      U_M = 0.5_rk*( uL_of_beta(beta_M,beta_L,U_L,sqk) + uR_of_beta(beta_M,beta_R,U_R,sqk) )
   end subroutine find_beta_M_bisection

   ! Helpers: the L/R curve U vs β
   pure function uL_of_beta(b, beta_L, U_L, sqk) result(u)
      real(rk), intent(in) :: b, beta_L, U_L, sqk
      real(rk) :: u
      if (b <= beta_L) then
         u = U_L - 2.0_rk*sqk*(sqrt(b) - sqrt(beta_L))     ! left rarefaction
      else
         u = U_L - (1.0_rk/sqrt(2.0_rk))*sqk*(b - beta_L)*sqrt(1.0_rk/b + 1.0_rk/beta_L) ! left shock
      end if
   end function uL_of_beta

   pure function uR_of_beta(b, beta_R, U_R, sqk) result(u)
      real(rk), intent(in) :: b, beta_R, U_R, sqk
      real(rk) :: u
      if (b <= beta_R) then
         u = U_R + 2.0_rk*sqk*(sqrt(b) - sqrt(beta_R))     ! right rarefaction
      else
         u = U_R + (1.0_rk/sqrt(2.0_rk))*sqk*(b - beta_R)*sqrt(1.0_rk/b + 1.0_rk/beta_R) ! right shock
      end if
   end function uR_of_beta

!======================================================================
!  Individual case fillers (Tables A.1-A.10)
!======================================================================

   !---- A.1 RR ------------------------------------------------------------
   subroutine fill_RR(beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M, st)
      real(rk), intent(in)  :: beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M
      type(riemann_state_t), intent(out) :: st
      real(rk) :: bL_rar, bR_rar
      st%case_id = CASE_RR
      st%beta_M = beta_M; st%U_M = U_M
      st%U_LL = U_L - sqk*sqrt(beta_L)
      st%U_LR = U_M - sqk*sqrt(beta_M)
      st%U_RL = U_M + sqk*sqrt(beta_M)
      st%U_RR = U_R + sqk*sqrt(beta_R)
      st%U_b  = 0.5_rk*(st%U_LR + st%U_RL)
      bL_rar = (1.0_rk/3.0_rk)*(beta_L + sqrt(beta_L*beta_M) + beta_M)
      bR_rar = (1.0_rk/3.0_rk)*(beta_R + sqrt(beta_R*beta_M) + beta_M)
      st%beta_ML = bL_rar*safe_div(st%U_LR-st%U_LL, st%U_b - st%U_LL) &
                  + beta_M*safe_div(st%U_b - st%U_LR, st%U_b - st%U_LL)
      st%U_ML    = U_L - 2.0_rk*sqk*(sqrt(max(st%beta_ML,0.0_rk)) - sqrt(beta_L))
      st%beta_MR = bR_rar*safe_div(st%U_RR-st%U_RL, st%U_RR - st%U_b) &
                  + beta_M*safe_div(st%U_RL - st%U_b, st%U_RR - st%U_b)
      st%U_MR    = U_R + 2.0_rk*sqk*(sqrt(max(st%beta_MR,0.0_rk)) - sqrt(beta_R))
   end subroutine fill_RR

   !---- A.2 RS ------------------------------------------------------------
   subroutine fill_RS(beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M, st)
      real(rk), intent(in)  :: beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M
      type(riemann_state_t), intent(out) :: st
      real(rk) :: bL_rar, s_shock
      st%case_id = CASE_RS
      st%beta_M = beta_M; st%U_M = U_M
      st%U_LL = U_L - sqk*sqrt(beta_L)
      st%U_LR = U_M - sqk*sqrt(beta_M)
      ! right shock speed s = (β_M U_M - β_R U_R)/(β_M - β_R)
      if (abs(beta_M - beta_R) < EPS_SMALL) then
         s_shock = U_M
      else
         s_shock = (beta_M*U_M - beta_R*U_R)/(beta_M - beta_R)
      end if
      st%U_RL = s_shock; st%U_RR = s_shock
      st%U_b  = s_shock                  ! track right shock
      bL_rar = (1.0_rk/3.0_rk)*(beta_L + sqrt(beta_L*beta_M) + beta_M)
      st%beta_ML = bL_rar*safe_div(st%U_LR-st%U_LL, st%U_b - st%U_LL) &
                  + beta_M*safe_div(st%U_b - st%U_LR, st%U_b - st%U_LL)
      st%U_ML    = U_L - 2.0_rk*sqk*(sqrt(max(st%beta_ML,0.0_rk)) - sqrt(beta_L))
      st%beta_MR = beta_R; st%U_MR = U_R       ! right shock => MR = R
   end subroutine fill_RS

   !---- A.3 SR ------------------------------------------------------------
   subroutine fill_SR(beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M, st)
      real(rk), intent(in)  :: beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M
      type(riemann_state_t), intent(out) :: st
      real(rk) :: bR_rar, s_shock
      st%case_id = CASE_SR
      st%beta_M = beta_M; st%U_M = U_M
      if (abs(beta_M - beta_L) < EPS_SMALL) then
         s_shock = U_M
      else
         s_shock = (beta_M*U_M - beta_L*U_L)/(beta_M - beta_L)
      end if
      st%U_LL = s_shock; st%U_LR = s_shock
      st%U_RL = U_M + sqk*sqrt(beta_M)
      st%U_RR = U_R + sqk*sqrt(beta_R)
      st%U_b  = 0.5_rk*(st%U_LR + st%U_RL)        ! SR has NO fast front (left wave is a SLOW shock,
                                                 ! right wave is a rarefaction); thesis Table A.3:
                                                 ! border tracked through the middle of the intermediate state.
      bR_rar = (1.0_rk/3.0_rk)*(beta_R + sqrt(beta_R*beta_M) + beta_M)
      st%beta_MR = bR_rar*safe_div(st%U_RR-st%U_RL, st%U_RR - st%U_b) &
                  + beta_M*safe_div(st%U_RL - st%U_b, st%U_RR - st%U_b)
      st%U_MR    = U_R + 2.0_rk*sqk*(sqrt(max(st%beta_MR,0.0_rk)) - sqrt(beta_R))
      st%beta_ML = beta_L; st%U_ML = U_L
   end subroutine fill_SR

   !---- A.4 SS ------------------------------------------------------------
   subroutine fill_SS(beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M, st)
      real(rk), intent(in)  :: beta_L, U_L, beta_R, U_R, kp, sqk, beta_M, U_M
      type(riemann_state_t), intent(out) :: st
      real(rk) :: sL, sR
      st%case_id = CASE_SS
      st%beta_M = beta_M; st%U_M = U_M
      ! Degenerate-shock guards: when β_M ≈ β_L (or β_R), the RH speed is
      ! formally 0/0. By L'Hôpital it equals the local fluid velocity.
      if (abs(beta_M - beta_L) < EPS_SMALL) then
         sL = U_M
      else
         sL = (beta_M*U_M - beta_L*U_L)/(beta_M - beta_L)
      end if
      if (abs(beta_M - beta_R) < EPS_SMALL) then
         sR = U_M
      else
         sR = (beta_M*U_M - beta_R*U_R)/(beta_M - beta_R)
      end if
      st%U_LL = sL; st%U_LR = sL
      st%U_RL = sR; st%U_RR = sR
      st%U_b = sR                                ! thesis Table A.4: track the right shock unconditionally
      st%beta_ML = beta_M; st%U_ML = U_M
      st%beta_MR = beta_R; st%U_MR = U_R
   end subroutine fill_SS

   !---- A.5 SSS (β_M = 1) -------------------------------------------------
   subroutine fill_SSS(beta_L, U_L, beta_R, U_R, kp, sqk, st)
      real(rk), intent(in)  :: beta_L, U_L, beta_R, U_R, kp, sqk
      type(riemann_state_t), intent(out) :: st
      real(rk) :: bM, U_M, sL, sR, wL, wR
      bM = 1.0_rk
      wL = sqrt(max(beta_L, BETA_DRY)/max(1.0_rk - beta_L, EPS_SMALL))
      wR = sqrt(max(beta_R, BETA_DRY)/max(1.0_rk - beta_R, EPS_SMALL))
      U_M = (wR*U_R + wL*U_L)/max(wL + wR, EPS_SMALL)
      if (abs(bM - beta_L) < EPS_SMALL) then
         sL = U_M
      else
         sL = (bM*U_M - beta_L*U_L)/(bM - beta_L)
      end if
      if (abs(bM - beta_R) < EPS_SMALL) then
         sR = U_M
      else
         sR = (bM*U_M - beta_R*U_R)/(bM - beta_R)
      end if
      st%case_id = CASE_SSS
      st%beta_M = bM; st%U_M = U_M
      st%U_LL = sL; st%U_LR = sL
      st%U_RL = sR; st%U_RR = sR
      st%U_b  = sR
      st%beta_ML = 1.0_rk; st%U_ML = U_M
      ! thesis Table A.5: with U_b = U_RR the right cell-update strip [U_b, U_RR] has zero length,
      ! so (β_MR, U_MR) is multiplied by zero in (3.28)–(3.29).  Set them to the right outer state
      ! (β_R, U_R) for formal consistency with the thesis Tables.
      st%beta_MR = beta_R; st%U_MR = U_R
   end subroutine fill_SSS

   !---- A.6 RV ------------------------------------------------------------
   subroutine fill_RV(beta_L, U_L, kp, sqk, st)
      real(rk), intent(in)  :: beta_L, U_L, kp, sqk
      type(riemann_state_t), intent(out) :: st
      real(rk) :: bL_rar
      st%case_id = CASE_RV
      st%beta_M = 0.0_rk
      st%U_M    = U_L + 2.0_rk*sqk*sqrt(beta_L)
      st%U_LL = U_L - sqk*sqrt(beta_L)
      st%U_LR = U_L + 2.0_rk*sqk*sqrt(beta_L)
      st%U_RL = st%U_LR
      st%U_RR = st%U_LR
      st%U_b  = st%U_LR
      bL_rar = (1.0_rk/3.0_rk)*beta_L
      st%beta_ML = bL_rar
      st%U_ML    = U_L - 2.0_rk*sqk*(sqrt(max(bL_rar,0.0_rk)) - sqrt(beta_L))
      st%beta_MR = 0.0_rk
      st%U_MR    = st%U_LR
   end subroutine fill_RV

   !---- A.7 VR ------------------------------------------------------------
   subroutine fill_VR(beta_R, U_R, kp, sqk, st)
      real(rk), intent(in)  :: beta_R, U_R, kp, sqk
      type(riemann_state_t), intent(out) :: st
      real(rk) :: bR_rar
      st%case_id = CASE_VR
      st%beta_M = 0.0_rk
      st%U_M    = U_R - 2.0_rk*sqk*sqrt(beta_R)
      st%U_LL = st%U_M; st%U_LR = st%U_M; st%U_RL = st%U_M
      st%U_RR = U_R + sqk*sqrt(beta_R)
      st%U_b  = st%U_M
      bR_rar = (1.0_rk/3.0_rk)*beta_R
      st%beta_MR = bR_rar
      st%U_MR    = U_R + 2.0_rk*sqk*(sqrt(max(bR_rar,0.0_rk)) - sqrt(beta_R))
      st%beta_ML = 0.0_rk
      st%U_ML    = st%U_M
   end subroutine fill_VR

   !---- A.8 RRADB ---------------------------------------------------------
   subroutine fill_RRADB(beta_L, U_L, beta_R, U_R, kp, sqk, st)
      real(rk), intent(in)  :: beta_L, U_L, beta_R, U_R, kp, sqk
      type(riemann_state_t), intent(out) :: st
      real(rk) :: bL_rar, bR_rar
      st%case_id = CASE_RRADB
      st%beta_M = 0.0_rk
      st%U_M    = 0.5_rk*((U_L + 2.0_rk*sqk*sqrt(beta_L)) + (U_R - 2.0_rk*sqk*sqrt(beta_R)))
      st%U_LL = U_L - sqk*sqrt(beta_L)
      st%U_LR = U_L + 2.0_rk*sqk*sqrt(beta_L)
      st%U_RL = U_R - 2.0_rk*sqk*sqrt(beta_R)
      st%U_RR = U_R + sqk*sqrt(beta_R)
      st%U_b  = 0.5_rk*(st%U_LR + st%U_RL)
      bL_rar = (1.0_rk/3.0_rk)*beta_L
      bR_rar = (1.0_rk/3.0_rk)*beta_R
      st%beta_ML = bL_rar*safe_div(st%U_LR-st%U_LL, st%U_b - st%U_LL)
      st%U_ML    = U_L - 2.0_rk*sqk*(sqrt(max(st%beta_ML,0.0_rk)) - sqrt(beta_L))
      st%beta_MR = bR_rar*safe_div(st%U_RR-st%U_RL, st%U_RR - st%U_b)
      st%U_MR    = U_R + 2.0_rk*sqk*(sqrt(max(st%beta_MR,0.0_rk)) - sqrt(beta_R))
   end subroutine fill_RRADB

!======================================================================
! Slug-side closures (Tables A.9, A.10)
!======================================================================

   !---- A.9 Section-Slug Nose: right neighbor is a slug ------------------
   ! L_s      : current slug length (m), for Moissis-Griffith wake effect.
   ! W_eff_cap: ceiling on the wake-effect multiplier (typically 2.0;
   !            set to 1.0 to disable wake correction).
   ! is_front : true ⇒ this bubble is a slug FRONT (rare for the upstream
   !            side; happens only if the slug is being overtaken by the
   !            upstream stratified film, U_ls < U_lb < U_crit). In that
   !            case U_b is set by the mass-conservation steep-front
   !            formula U_F = (U_ls - β_b U_lb)/(1 - β_b)  (thesis §4.5).
   ! is_front=false ⇒ classical slug NOSE: U_b = W_eff·U_bendiksen.
   ! See SPEC §4.6 / Table 2.1 for the turning-point criterion.
   subroutine solve_section_slug_nose(beta_L, U_L, U_mslug, L_s, D, phi, W_eff_cap, is_front, st)
      real(rk), intent(in)  :: beta_L, U_L, U_mslug, L_s, D, phi, W_eff_cap
      logical,  intent(in)  :: is_front
      type(riemann_state_t), intent(out) :: st
      real(rk) :: U_b_use, W_loc, b_use
      if (is_front) then
         b_use   = max(min(beta_L, 1.0_rk - EPS_SMALL), EPS_SMALL)
         U_b_use = safe_div(U_mslug - b_use*U_L, 1.0_rk - b_use)
      else
         W_loc   = wake_effect(L_s, D, W_eff_cap)
         U_b_use = W_loc*bendiksen_nose(U_mslug, D, phi)
      end if
      st%case_id = CASE_SECSLUG
      st%U_LL = U_L
      st%U_LR = U_b_use
      st%U_RL = U_b_use
      st%U_RR = U_b_use
      st%U_b  = U_b_use
      st%beta_M = 1.0_rk; st%U_M = U_mslug
      ! mass balance in fan (β_ML · (U_b - U_LL) = β_slug U_mslug + (... )
      ! Clip fan β to [0,1] – the raw safe_div ratio can go negative or
      ! > 1 when wake-drained upstream film momentarily accelerates past
      ! the bubble speed (U_L > U_b) and then gets remapped with a
      ! negative fan β into the flanking section, collapsing its β to
      ! BETA_DRY and producing the β≈0 "dry patch" seen at t≥11 s.
      st%beta_ML = max(min(safe_div(U_b_use - U_mslug, U_b_use - U_L), 1.0_rk), 0.0_rk)
      st%U_ML    = U_mslug
      st%beta_MR = 1.0_rk; st%U_MR = U_mslug
   end subroutine solve_section_slug_nose

   !---- A.10 Slug-Section Nose: left neighbor is a slug ------------------
   ! Same conventions as A.9; is_front=true picks the steep-front speed.
   ! For a downstream bubble FRONT is the COMMON case (slug eats the
   ! stratified film ahead).  See SPEC §4.5–§4.6 / Table 2.1.
   subroutine solve_slug_section_nose(beta_R, U_R, U_mslug, L_s, D, phi, W_eff_cap, is_front, st)
      real(rk), intent(in)  :: beta_R, U_R, U_mslug, L_s, D, phi, W_eff_cap
      logical,  intent(in)  :: is_front
      type(riemann_state_t), intent(out) :: st
      real(rk) :: U_b_use, W_loc, b_use
      if (is_front) then
         b_use   = max(min(beta_R, 1.0_rk - EPS_SMALL), EPS_SMALL)
         U_b_use = safe_div(U_mslug - b_use*U_R, 1.0_rk - b_use)
      else
         W_loc   = wake_effect(L_s, D, W_eff_cap)
         U_b_use = W_loc*bendiksen_nose(U_mslug, D, phi)
      end if
      st%case_id = CASE_SLUGSEC
      st%U_RR = U_R
      st%U_RL = U_b_use
      st%U_LL = U_b_use
      st%U_LR = U_b_use
      st%U_b  = U_b_use
      st%beta_M = 1.0_rk; st%U_M = U_mslug
      ! Clip fan β to [0,1] (see commentary in solve_section_slug_nose).
      st%beta_MR = max(min(safe_div(U_mslug - U_b_use, U_R - U_b_use), 1.0_rk), 0.0_rk)
      st%U_MR    = U_mslug
      st%beta_ML = 1.0_rk; st%U_ML = U_mslug
   end subroutine solve_slug_section_nose

!======================================================================
   pure function safe_div(a, b) result(r)
      real(rk), intent(in) :: a, b
      real(rk) :: r
      if (abs(b) < EPS_TINY) then
         r = 0.0_rk
      else
         r = a/b
      end if
   end function safe_div

end module lassi_riemann
