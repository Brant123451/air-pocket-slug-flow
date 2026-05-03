!======================================================================
! lassi_voidwave.f90
! LASSI void-wave step: at each border, solve a Riemann problem,
! advect borders Lagrangian-fashion, remap holdup and liquid velocity
! using the constant-state pieces left/right of the new cell boundaries
! (Eqs. 3.27-3.29 of the spec).
!
! Two-pass design:
!   Pass A:  for every interior bubble, solve the Riemann problem and
!            store all wave speeds + intermediate states + new position.
!   Pass B:  for every section, remap (β, U_l) onto the new cell using
!            piecewise constant pieces (β_old, β_ML on left fan, β_MR
!            on right fan).
! Reference: LASSI_ALGORITHM_SPEC.md sections 5.2-5.3.
!======================================================================
module lassi_voidwave
   use lassi_kinds
   use lassi_objects
   use lassi_grid
   use lassi_geom,    only: ikh_kappa
   use lassi_riemann
   use lassi_friction, only: u_crit_balance, wake_effect, bendiksen_nose
   use lassi_listmgmt, only: merge_cfl_breaching_sections, slug_shed, slug_front_eat
   implicit none
   private
   public :: void_wave_step
   public :: border_info_t   ! made public so other modules can introspect

   ! Slug-interior λ_l for the turning-point criterion (thesis page 21,
   ! Eq. 2.12) is now read from grid_t%lam_l_slug (case input
   ! "lam_l_slug", default 0.02 — the value cited on page 21 for fully
   ! turbulent water/air flow).

   type :: border_info_t
      integer(ik) :: bid           = 0
      real(rk)    :: U_b           = 0.0_rk
      real(rk)    :: U_LL          = 0.0_rk
      real(rk)    :: U_LR          = 0.0_rk
      real(rk)    :: U_RL          = 0.0_rk
      real(rk)    :: U_RR          = 0.0_rk
      real(rk)    :: beta_ML       = 0.0_rk
      real(rk)    :: U_ML          = 0.0_rk
      real(rk)    :: beta_MR       = 0.0_rk
      real(rk)    :: U_MR          = 0.0_rk
      real(rk)    :: zR_old        = 0.0_rk
      real(rk)    :: zR_new        = 0.0_rk
      integer(ik) :: case_id       = 0
   end type border_info_t

contains

!======================================================================
   subroutine void_wave_step(g, dt, W_eff)
      type(grid_t), intent(inout) :: g
      real(rk),     intent(in)    :: dt, W_eff

      type(border_info_t), allocatable :: bi(:)
      integer(ik) :: nb
      logical :: merged

      do
         call collect_borders(g, bi, nb, W_eff)
         if (nb == 0) return
         call mark_cfl_breaches(g, bi, dt)
         call merge_cfl_breaching_sections(g, dt, merged)
         deallocate(bi)
         if (.not. merged) exit
      end do

      call collect_borders(g, bi, nb, W_eff)
      if (nb == 0) return
      call slug_shed(g, dt, bi%bid, bi%U_b)
      deallocate(bi)
      call collect_borders(g, bi, nb, W_eff)
      if (nb == 0) return
      call advance_bubble_positions(g, bi, dt)
      call recompute_lengths(g)
      call remap_all_sections(g, bi, dt)
      call sync_owned_bubbles_from_sections(g)
      call slug_front_eat(g, dt)

      deallocate(bi)
   end subroutine void_wave_step

!======================================================================
   subroutine collect_borders(g, bi, nb, W_eff)
      type(grid_t),                       intent(inout) :: g    ! writes Bubble%is_nose
      type(border_info_t), allocatable,    intent(out)   :: bi(:)
      integer(ik),                         intent(out)   :: nb
      real(rk),                            intent(in)    :: W_eff

      type(object_t), pointer :: p, sec1, sec2
      type(riemann_state_t) :: st
      integer(ik) :: ib
      real(rk) :: kappa
      real(rk) :: U_crit, U_ls, U_lb, U_b_nose_test, W_loc
      logical  :: is_front

      ! count
      nb = 0
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_INLET .or. p%kind == KIND_BUBBLE) nb = nb + 1
         p => p%next
      end do
      if (nb == 0) then
         allocate(bi(0)); return
      end if
      allocate(bi(nb))

      ! Slug turning-point critical velocity (thesis Eq. 2.12).
      ! In the K-L 2018 horizontal pipe (sin φ = 0) U_crit = 0, which makes
      ! every slug-section bubble a FRONT whenever U_ls > U_lb (the common
      ! case during slug growth).  For inclined pipes U_crit ≠ 0 selects
      ! between Bendiksen drift and steep-front advection per Table 2.1.
      U_crit = u_crit_balance(g%D, g%phi, g%lam_l_slug)

      ib = 0
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_INLET .or. p%kind == KIND_BUBBLE) then
            ib = ib + 1
            if (p%kind == KIND_INLET) then
               sec1 => p
               sec2 => p%next
            else
               sec1 => p%prev
               sec2 => p%next
            end if

            if (.not. associated(sec1) .or. .not. associated(sec2)) then
               call set_static(st, p)
               p%is_nose = .false.
            else if (sec1%kind == KIND_SECTION .and. sec2%kind == KIND_SECTION) then
               kappa = 0.5_rk*( &
                  ikh_kappa(sec1%beta, g%D, g%phi, g%rho_l, section_rhog(sec1), sec1%Ul, section_ug(sec1)) + &
                  ikh_kappa(sec2%beta, g%D, g%phi, g%rho_l, section_rhog(sec2), sec2%Ul, section_ug(sec2)))
               call solve_msw_riemann(sec1%beta, sec1%Ul, sec2%beta, sec2%Ul, kappa, st)
               p%is_nose = .false.

            else if (sec1%kind == KIND_SECTION .and. sec2%kind == KIND_SLUG) then
               ! W_eff is interpreted as the wake-effect ceiling (W_eff_cap):
               !   W_eff_cap = 1.0  ⇒ disable wake effect (long-slug Bendiksen).
               !   W_eff_cap = 2.0  ⇒ thesis default (Moissis–Griffith capped at 2).
               !
               ! Bubble-slug border (bubble on left, slug on right) — Table 2.1
               ! row.  FRONT only when U_ls < U_lb AND U_ls < U_crit (i.e. the
               ! upstream stratified film overtakes a near-stalled slug).
               ! Otherwise NOSE (Bendiksen drift) — the typical case.
               U_ls = sec2%Um
               U_lb = sec1%Ul
               is_front = (U_ls < U_lb) .and. (U_ls < U_crit)
               ! Self-consistency safeguard for the NOSE branch (SPEC §5.4 /
               ! Table A.9): the mass-balance fan formula
               !    β_ML = (U_nose − U_mslug)/(U_b − U_LL)
               ! requires U_b > U_LL = U_lb, otherwise the raw β_ML goes
               ! negative and `solve_section_slug_nose` clips it to 0,
               ! injecting an unphysical β=0 strip into the upstream
               ! SECTION (the dry-patch artefact at t≥11s in case 4).  If
               ! the wake-drained upstream film is faster than the
               ! Bendiksen-candidate nose speed, the NOSE assumption is
               ! physically broken — the slug is in fact being pushed
               ! along by the upstream film, i.e. behaves as a steep
               ! mass-conservation FRONT.  Force is_front=true so that
               ! solve_section_slug_nose uses U_b = U_F (mass-conserving
               ! front), which makes β_ML = β_L self-consistently and no
               ! clip is needed.
               if (.not. is_front) then
                  W_loc = wake_effect(object_length(sec2), g%D, W_eff)
                  U_b_nose_test = W_loc*bendiksen_nose(U_ls, g%D, g%phi)
                  if (U_lb > U_b_nose_test) is_front = .true.
               end if
               p%is_nose = .not. is_front
               call solve_section_slug_nose(sec1%beta, sec1%Ul, sec2%Um, &
                                            object_length(sec2), g%D, g%phi, &
                                            W_eff, is_front, st)

            else if (sec1%kind == KIND_SLUG .and. sec2%kind == KIND_SECTION) then
               ! Slug-bubble border (slug on left, bubble on right) — Table 2.1
               ! row.  FRONT when U_ls > U_lb AND U_ls > U_crit (the COMMON
               ! case: slug eats the stratified film ahead).
               U_ls = sec1%Um
               U_lb = sec2%Ul
               is_front = (U_ls > U_lb) .and. (U_ls > U_crit)
               ! Self-consistency safeguard for the NOSE branch (SPEC
               ! Table A.10): β_MR = (U_mslug − U_nose)/(U_RR − U_b) needs
               ! U_b < U_RR = U_lb (downstream film slower than nose) for
               ! a non-negative raw β_MR.  If the downstream film is
               ! moving faster than the NOSE candidate (the bubble is
               ! being chased by quick film ahead), the NOSE assumption
               ! is physically broken and the configuration is actually a
               ! steep FRONT.  Force is_front=true so U_b = U_F.  See the
               ! upstream branch above for the analogous rationale.
               if (.not. is_front) then
                  W_loc = wake_effect(object_length(sec1), g%D, W_eff)
                  U_b_nose_test = W_loc*bendiksen_nose(U_ls, g%D, g%phi)
                  if (U_lb < U_b_nose_test) then
                     ! downstream film slower than nose → genuine NOSE, OK
                  else
                     is_front = .true.
                  end if
               end if
               p%is_nose = .not. is_front
               call solve_slug_section_nose(sec2%beta, sec2%Ul, sec1%Um, &
                                            object_length(sec1), g%D, g%phi, &
                                            W_eff, is_front, st)

            else if (sec1%kind == KIND_INLET .and. sec2%kind == KIND_SECTION) then
               st%U_b = 0.0_rk
               st%U_LL = sec1%Ul; st%U_LR = sec1%Ul
               st%U_RL = sec1%Ul; st%U_RR = sec1%Ul
               st%beta_M = sec1%beta; st%U_M = sec1%Ul
               st%beta_ML = sec1%beta; st%U_ML = sec1%Ul
               st%beta_MR = sec1%beta; st%U_MR = sec1%Ul
               st%case_id = 0
               p%is_nose = .false.

            else if (sec1%kind == KIND_SECTION .and. sec2%kind == KIND_OUTLET) then
               ! Transmissive (zero-gradient) outlet: mirror the upstream
               ! section state.  The Lagrangian border then moves with
               ! U_b = sec1%Ul and the list-management step cuts the
               ! protruding piece at z = z_outlet, logging that liquid as
               ! outflow.  This matches a true forward-flow outlet with no
               ! spurious dry-bed rarefaction draining the pipe end.
               st%U_b  = sec1%Ul
               st%U_LL = sec1%Ul; st%U_LR = sec1%Ul
               st%U_RL = sec1%Ul; st%U_RR = sec1%Ul
               st%beta_M = sec1%beta; st%U_M = sec1%Ul
               st%beta_ML = sec1%beta; st%U_ML = sec1%Ul
               st%beta_MR = sec1%beta; st%U_MR = sec1%Ul
               st%case_id = 0
               p%is_nose = .false.

            else if (sec1%kind == KIND_SLUG .and. sec2%kind == KIND_OUTLET) then
               ! Slug reaches the outlet (no intervening section).  The
               ! slug body advances at U_m, so the slug-OUTLET bubble
               ! border must move with U_m too; otherwise the slug is
               ! pinned in the Lagrangian list and never clipped past
               ! Lpipe (silent stalling at outlet, with mass piling up).
               st%U_b  = sec1%Um
               st%U_LL = sec1%Um; st%U_LR = sec1%Um
               st%U_RL = sec1%Um; st%U_RR = sec1%Um
               st%beta_M = 1.0_rk; st%U_M = sec1%Um
               st%beta_ML = 1.0_rk; st%U_ML = sec1%Um
               st%beta_MR = 1.0_rk; st%U_MR = sec1%Um
               st%case_id = 0
               p%is_nose = .false.

            else
               call set_static(st, p)
               p%is_nose = .false.
            end if

            bi(ib)%bid     = p%id
            bi(ib)%U_b     = st%U_b
            bi(ib)%U_LL    = st%U_LL
            bi(ib)%U_LR    = st%U_LR
            bi(ib)%U_RL    = st%U_RL
            bi(ib)%U_RR    = st%U_RR
            bi(ib)%beta_ML = st%beta_ML
            bi(ib)%U_ML    = st%U_ML
            bi(ib)%beta_MR = st%beta_MR
            bi(ib)%U_MR    = st%U_MR
            bi(ib)%zR_old  = p%zR
            bi(ib)%zR_new  = p%zR  ! filled in advance step
            bi(ib)%case_id = st%case_id
         end if
         p => p%next
      end do
   end subroutine collect_borders

!======================================================================
   subroutine advance_bubble_positions(g, bi, dt)
      type(grid_t),        intent(inout) :: g
      type(border_info_t), intent(inout) :: bi(:)
      real(rk),            intent(in)    :: dt
      type(object_t), pointer :: p
      integer(ik) :: ib, n
      logical :: hold_slug_front
      n = size(bi)
      ib = 0
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_INLET .or. p%kind == KIND_BUBBLE) then
            ib = ib + 1
            if (ib > n) exit
            if (p%kind == KIND_INLET) then
               p%zR = 0.0_rk
            else
               hold_slug_front = .false.
               if (.not. p%is_nose .and. associated(p%prev) .and. associated(p%next)) then
                  hold_slug_front = (p%prev%kind == KIND_SLUG .and. p%next%kind == KIND_SECTION)
               end if
               if (.not. hold_slug_front) then
                  p%zR = p%zR + bi(ib)%U_b*dt
                  ! also move the Section that owns this bubble (Section.zR coincides with Bubble.zR)
                  if (associated(p%prev)) then
                     if (p%prev%kind == KIND_SECTION .or. p%prev%kind == KIND_SLUG) then
                        p%prev%zR = p%zR
                     end if
                  end if
               end if
            end if
            bi(ib)%zR_new = p%zR
         end if
         p => p%next
      end do
   end subroutine advance_bubble_positions

   subroutine mark_cfl_breaches(g, bi, dt)
      type(grid_t),        intent(inout) :: g
      type(border_info_t), intent(in)    :: bi(:)
      real(rk),            intent(in)    :: dt
      type(object_t), pointer :: p, left_b, right_b
      integer(ik) :: ib_left, ib_right
      real(rk) :: risk

      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SECTION) p%cfl_breach = .false.
         p => p%next
      end do

      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SECTION) then
            left_b => p%prev
            right_b => p%next
            if (associated(left_b) .and. associated(right_b)) then
               ib_left = lookup_border(bi, left_b%id)
               ib_right = lookup_border(bi, right_b%id)
               if (ib_left > 0 .and. ib_right > 0) then
                  risk = (bi(ib_left)%U_RR - bi(ib_right)%U_LL)*dt
                  p%cfl_breach = (risk > object_length(p))
               end if
            end if
         end if
         p => p%next
      end do
   end subroutine mark_cfl_breaches

!======================================================================
! Remap each Section onto its new domain using a piecewise-constant
! decomposition of the wave fan (β_ML on the left fan, β_MR on the right
! fan, central=old state).  This is a first-order conservative remap
! consistent with (3.28)-(3.29).  In production we should split each fan
! at the actual rarefaction head/tail; the present version uses the
! averaged (Heronian) state computed inside solve_msw_riemann, which
! preserves mass exactly to leading order.
!======================================================================
   subroutine remap_all_sections(g, bi, dt)
      type(grid_t),                  intent(inout) :: g
      type(border_info_t),           intent(in)    :: bi(:)
      real(rk),                      intent(in)    :: dt

      type(object_t), pointer :: p, prev_bub, next_bub
      integer(ik) :: ib_left, ib_right
      real(rk)    :: L_new, dz_left, dz_right, dz_cent
      real(rk)    :: beta_acc, mom_acc

      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SECTION) then
            prev_bub => p%prev   ! left  bubble (or inlet)
            next_bub => p%next   ! right bubble

            L_new = max(p%zR - prev_bub%zR, 0.0_rk)
            if (L_new < EPS_TINY) then
               p%beta = max(p%beta, EPS_TINY)
               p => p%next
               cycle
            end if

            ! Fan widths in THIS section:
            !   left bubble's RIGHT fan (state β_MR) enters from the left,
            !   width = (U_RR - U_b)*dt of the LEFT bubble.
            !   right bubble's LEFT fan (state β_ML) enters from the right,
            !   width = (U_b - U_LL)*dt of the RIGHT bubble.
            ib_left  = lookup_border(bi, prev_bub%id)
            ib_right = lookup_border(bi, next_bub%id)
            dz_left  = 0.0_rk
            dz_right = 0.0_rk
            if (ib_left  > 0) dz_left  = max((bi(ib_left )%U_RR - bi(ib_left )%U_b )*dt, 0.0_rk)
            if (ib_right > 0) dz_right = max((bi(ib_right)%U_b  - bi(ib_right)%U_LL)*dt, 0.0_rk)
            ! If the two fans together would exceed the new cell length
            ! (this can happen if the CFL safety margin is too generous),
            ! shrink them PROPORTIONALLY to keep total = L_new and preserve
            ! mass.  Otherwise the central region preserves the old state.
            if (dz_left + dz_right > L_new) then
               block
                  real(rk) :: scale
                  scale = L_new/max(dz_left + dz_right, EPS_TINY)
                  dz_left  = dz_left*scale
                  dz_right = dz_right*scale
               end block
            end if
            dz_cent  = max(L_new - dz_left - dz_right, 0.0_rk)

            ! conservative mix
            beta_acc = p%beta*dz_cent
            mom_acc  = p%beta*p%Ul*dz_cent
            if (ib_left  > 0) then
               beta_acc = beta_acc + bi(ib_left )%beta_MR*dz_left
               mom_acc  = mom_acc  + bi(ib_left )%beta_MR*bi(ib_left )%U_MR*dz_left
            end if
            if (ib_right > 0) then
               beta_acc = beta_acc + bi(ib_right)%beta_ML*dz_right
               mom_acc  = mom_acc  + bi(ib_right)%beta_ML*bi(ib_right)%U_ML*dz_right
            end if

            p%beta = max(beta_acc/L_new, BETA_DRY)
            if (p%beta > BETA_DRY) then
               p%Ul = max(min(mom_acc/(p%beta*L_new), 20.0_rk), -20.0_rk)
            end if
         end if
         p => p%next
      end do
   end subroutine remap_all_sections

!======================================================================
   pure function lookup_border(bi, bid) result(idx)
      type(border_info_t), intent(in) :: bi(:)
      integer(ik),         intent(in) :: bid
      integer(ik) :: idx, k
      idx = 0
      do k = 1, size(bi)
         if (bi(k)%bid == bid) then
            idx = k; return
         end if
      end do
   end function lookup_border

!======================================================================
   subroutine set_static(st, p)
      type(riemann_state_t), intent(out) :: st
      type(object_t), pointer, intent(in) :: p
      st%U_b = 0.0_rk
      st%U_LL = 0.0_rk; st%U_LR = 0.0_rk; st%U_RL = 0.0_rk; st%U_RR = 0.0_rk
      st%beta_M = p%beta; st%U_M = p%Ul
      st%beta_ML = p%beta; st%U_ML = p%Ul
      st%beta_MR = p%beta; st%U_MR = p%Ul
      st%case_id = 0
   end subroutine set_static

!======================================================================
   function section_ug(p) result(Ug)
      type(object_t), pointer, intent(in) :: p
      real(rk) :: Ug
      real(rk) :: alpha
      type(object_t), pointer :: bub
      alpha = max(1.0_rk - p%beta, EPS_SMALL)
      bub => p%next
      if (associated(bub)) then
         if (bub%kind == KIND_BUBBLE) then
            Ug = bub%Ugs/alpha
            return
         end if
      end if
      Ug = p%Ugs/alpha
   end function section_ug

   function section_rhog(p) result(rhog)
      type(object_t), pointer, intent(in) :: p
      real(rk) :: rhog
      type(object_t), pointer :: bub
      bub => p%next
      if (associated(bub)) then
         if (bub%kind == KIND_BUBBLE) then
            rhog = bub%rhog
            return
         end if
      end if
      rhog = p%rhog
   end function section_rhog

end module lassi_voidwave
