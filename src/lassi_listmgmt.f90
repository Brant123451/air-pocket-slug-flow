!======================================================================
! lassi_listmgmt.f90
! LASSI list management: merge CFL-breaching sections, split too-long
! sections, slug initiation, slug shedding, inlet/outlet treatment.
! Reference: LASSI_ALGORITHM_SPEC.md section 6 and 7.
!======================================================================
module lassi_listmgmt
   use lassi_kinds
   use lassi_objects
   use lassi_grid
   use lassi_friction, only: bendiksen_nose, wake_effect
   implicit none
   private
   public :: list_management, slug_shed, merge_cfl_breaching_sections, slug_front_eat

contains

   ! Local helper: divide with a safe fallback when the denominator is
   ! tiny.  Mirrors lassi_press_mom.safe_div_ufront, duplicated here to
   ! avoid a forward module dependency.
   pure function safe_div_uf(num, den, fallback) result(uf)
      real(rk), intent(in) :: num, den, fallback
      real(rk) :: uf
      if (abs(den) < 1.0e-3_rk) then
         uf = fallback
      else
         uf = num/den
      end if
   end function safe_div_uf

   subroutine list_management(g, dt)
      type(grid_t), intent(inout) :: g
      real(rk),     intent(in)    :: dt
      ! Paper page 47 "Check CFL" step (thesis §3.5, scan `20.png`):
      ! for every section test (U_RR,J − U_LL,J+1)·δt > L_J^n and merge
      ! the section with a neighbour if the criterion is breached.  Run
      ! BEFORE merge_short_sections so the static-length pass can clean
      ! up whatever remains below 0.5·TargetLength.
      call merge_short_sections(g)
      call purge_zero_length_sections(g)
      ! Thesis Fig. 3.2/Fig. 3.9 ordering: existing slugs are tested
      ! for survival and their fronts are advanced before the list
      ! management pass creates new slugs from newly flooded sections.
      ! A newly initiated slug must first pass through the next
      ! Pressure-Momentum/Void-Wave cycle before participating in
      ! front/eating or shedding dynamics.
      if (g%enable_slug_coalescence) call merge_adjacent_slugs(g)   ! optional thesis page 24 coalescence
      call slug_init(g)
      call split_long_sections(g)
      call inlet_treatment(g)
      call outlet_treatment(g)
      call assign_unit_ids(g)
      call recompute_lengths(g)
      call sync_section_gas_from_bubbles(g)
      call sync_owned_bubbles_from_sections(g)
   end subroutine list_management

!======================================================================
! Purge any SECTION whose length has collapsed to zero (e.g. squeezed
! between two SLUGs where the regular merge_short_sections cannot find
! a compatible neighbor SECTION to merge into).  Liquid mass in such a
! section is identically zero, so deletion is conservative.  We delete
! the section together with one of its flanking bubbles to preserve the
! alternating ...-(SEC|SLUG)-BUBBLE-(SEC|SLUG)-... structure.
!======================================================================
   subroutine purge_zero_length_sections(g)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer :: p, p_next, bub_l, bub_r
      real(rk), parameter :: L_EPS = 1.0e-4_rk

      p => g%head
      do while (associated(p))
         p_next => p%next
         if (p%kind == KIND_SECTION) then
            if (object_length(p) < L_EPS) then
               bub_l => p%prev
               bub_r => p%next
               ! prefer to delete the right bubble; fall back to left if
               ! right is OUTLET or absent.
               if (associated(bub_r)) then
                  if (bub_r%kind == KIND_BUBBLE) then
                     p_next => bub_r%next
                     call remove_object(g, bub_r)
                     call remove_object(g, p)
                     g%gas_track_initialized = .false.
                  else if (associated(bub_l)) then
                     if (bub_l%kind == KIND_BUBBLE) then
                        call remove_object(g, bub_l)
                        call remove_object(g, p)
                        g%gas_track_initialized = .false.
                     end if
                  end if
               else if (associated(bub_l)) then
                  if (bub_l%kind == KIND_BUBBLE) then
                     call remove_object(g, bub_l)
                     call remove_object(g, p)
                     g%gas_track_initialized = .false.
                  end if
               end if
            end if
         end if
         p => p_next
      end do
   end subroutine purge_zero_length_sections

!======================================================================
! Paper "Check CFL" step (thesis §3.5 page 47, scan `20.png`, Figure
! 3.9 second box).  The paper criterion is
!
!     (U_RR,J − U_LL,J+1) · δt > L_J^n
!
! where U_RR,J is the rightmost-state velocity of section J and
! U_LL,J+1 is the leftmost-state velocity of section J+1, i.e. the
! post-Riemann far-field velocities on either side of the shared
! bubble border.  In our data layout each section stores only its
! single cell-averaged U_l, so we use a conservative proxy: the
! maximum |velocity| reachable through any neighbouring object is
! taken as an upper bound on |U_RR,J − U_LL,J+1|.  A section whose
! current length is shorter than this "risk distance" is in danger
! of being squeezed to zero before the next void-wave step, so it is
! merged pre-emptively with its closer-β section neighbour (paper
! constraint: a CFL-breaching Section may only be merged with another
! Section, not with a Slug).
!======================================================================
   subroutine merge_cfl_breaching_sections(g, dt, merged)
      type(grid_t), intent(inout) :: g
      real(rk),     intent(in)    :: dt
      logical, optional, intent(out) :: merged
      type(object_t), pointer :: p, p_next, bub_l, bub_r, lneigh, rneigh
      real(rk) :: L_J, L_other
      logical  :: merge_with_left

      if (present(merged)) merged = .false.
      if (dt <= 0.0_rk) return
      p => g%head
      do while (associated(p))
         p_next => p%next
         if (p%kind == KIND_SECTION) then
            L_J = object_length(p)
            bub_l => p%prev
            bub_r => p%next
            lneigh => null(); rneigh => null()
            if (associated(bub_l)) lneigh => bub_l%prev
            if (associated(bub_r)) rneigh => bub_r%next
            if (p%cfl_breach) then
               ! Paper constraint: merge with a Section neighbour only.
               merge_with_left = .false.
               if (associated(lneigh) .and. associated(rneigh)) then
                  if (lneigh%kind == KIND_SECTION .and. rneigh%kind == KIND_SECTION) then
                     if (abs(lneigh%beta - p%beta) <= abs(rneigh%beta - p%beta)) then
                        merge_with_left = .true.
                     end if
                  else if (lneigh%kind == KIND_SECTION) then
                     merge_with_left = .true.
                  else if (rneigh%kind == KIND_SECTION) then
                     merge_with_left = .false.
                  else
                     p => p_next
                     cycle
                  end if
               else if (associated(rneigh)) then
                  if (rneigh%kind /= KIND_SECTION) then
                     p => p_next
                     cycle
                  end if
               else if (associated(lneigh)) then
                  if (lneigh%kind /= KIND_SECTION) then
                     p => p_next
                     cycle
                  end if
                  merge_with_left = .true.
               else
                  p => p_next
                  cycle
               end if
               if (merge_with_left) then
                  L_other = object_length(lneigh)
                  call merge_two(g, lneigh, p, L_other, L_J, bub_l)
                  p_next => lneigh%next
               else
                  L_other = object_length(rneigh)
                  call merge_two(g, p, rneigh, L_J, L_other, bub_r)
                  p_next => p%next
               end if
               if (present(merged)) merged = .true.
               return
            end if
         end if
         p => p_next
      end do
   end subroutine merge_cfl_breaching_sections

!======================================================================
! Merge a section with its closer-β neighbor when it becomes too short.
! Triggers: L_J < 0.5*TargetLength  (thesis §4.7 / Sec. 3.5 page 47).
!======================================================================
   subroutine merge_short_sections(g)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer :: p, p_next, lneigh, rneigh, bub_l, bub_r
      real(rk) :: L_thr, L_J, L_other
      logical  :: merge_with_left

      L_thr = 0.5_rk*g%TargetLength
      p => g%head
      do while (associated(p))
         p_next => p%next
         if (p%kind == KIND_SECTION) then
            L_J = object_length(p)
            if (L_J < L_thr) then
               ! find the two section neighbors
               bub_l => p%prev      ! Bubble to left
               bub_r => p%next      ! Bubble to right
               lneigh => null(); rneigh => null()
               if (associated(bub_l)) lneigh => bub_l%prev
               if (associated(bub_r)) rneigh => bub_r%next

               ! prefer merging with a Section, not a Slug or boundary
               merge_with_left = .false.
               if (associated(lneigh) .and. associated(rneigh)) then
                  if (lneigh%kind == KIND_SECTION .and. rneigh%kind == KIND_SECTION) then
                     ! pick closer-β neighbor
                     if (abs(lneigh%beta - p%beta) <= abs(rneigh%beta - p%beta)) then
                        merge_with_left = .true.
                     end if
                  else if (lneigh%kind == KIND_SECTION) then
                     merge_with_left = .true.
                  else if (rneigh%kind == KIND_SECTION) then
                     merge_with_left = .false.
                  else
                     p => p_next
                     cycle
                  end if
               else if (associated(rneigh)) then
                  if (rneigh%kind == KIND_SECTION) then
                     merge_with_left = .false.
                  else
                     p => p_next
                     cycle
                  end if
               else if (associated(lneigh)) then
                  if (lneigh%kind == KIND_SECTION) then
                     merge_with_left = .true.
                  else
                     p => p_next
                     cycle
                  end if
               else
                  p => p_next
                  cycle
               end if

               if (merge_with_left) then
                  L_other = object_length(lneigh)
                  call merge_two(g, lneigh, p, L_other, L_J, bub_l)
                  ! after merging, lneigh extended over p; p has been deleted
                  p_next => lneigh%next   ! continue past the merged bubble
               else
                  L_other = object_length(rneigh)
                  call merge_two(g, p, rneigh, L_J, L_other, bub_r)
                  ! p has now grown over rneigh; continue past
                  p_next => p%next
               end if
            end if
         end if
         p => p_next
      end do
   end subroutine merge_short_sections

   subroutine merge_two(g, sec_left, sec_right, Ll, Lr, bub_between)
      type(grid_t),               intent(inout) :: g
      type(object_t), pointer,    intent(inout) :: sec_left, sec_right, bub_between
      real(rk),                   intent(in)    :: Ll, Lr
      real(rk) :: bL, bR, uL, uR, total
      bL = sec_left%beta; uL = sec_left%Ul
      bR = sec_right%beta; uR = sec_right%Ul
      total = bL*Ll + bR*Lr
      sec_left%zR  = sec_right%zR
      sec_left%beta = (bL*Ll + bR*Lr)/(Ll + Lr)
      if (total > EPS_SMALL) then
         sec_left%Ul = (bL*Ll*uL + bR*Lr*uR)/total
      end if
      ! remove the in-between bubble and the right section
      call remove_object(g, bub_between)
      call remove_object(g, sec_right)
      g%gas_track_initialized = .false.
   end subroutine merge_two

!======================================================================
! Split a section in two when it grows beyond 1.5 * TargetLength.
! (thesis §4.7 / Sec. 3.5 page 47.)
!======================================================================
   subroutine split_long_sections(g)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer :: p, sec_new, bub_new, p_next
      real(rk) :: L_J, L_thr

      L_thr = 1.5_rk*g%TargetLength
      p => g%head
      do while (associated(p))
         p_next => p%next
         if (p%kind == KIND_SECTION) then
            L_J = object_length(p)
            if (L_J > L_thr) then
               ! create a new bubble + section after p, splitting at zR_left + TargetLength
               allocate(sec_new); allocate(bub_new)
               bub_new%kind = KIND_BUBBLE
               bub_new%zR   = p%prev%zR + g%TargetLength
               bub_new%beta = p%beta
               bub_new%Ul   = p%Ul
               bub_new%pres = p%pres
               bub_new%rhog = p%rhog
               bub_new%Ugs  = p%Ugs
               sec_new%kind = KIND_SECTION
               sec_new%zR   = p%zR    ! original right border
               sec_new%beta = p%beta
               sec_new%Ul   = p%Ul
               sec_new%pres = p%pres
               sec_new%rhog = p%rhog
               sec_new%Ugs  = p%Ugs
               ! shrink p to TargetLength
               p%zR = bub_new%zR
               ! insert: p -> bub_new -> sec_new -> (old p%next)
               call insert_after(g, p, bub_new)
               call insert_after(g, bub_new, sec_new)
               g%gas_track_initialized = .false.
               p_next => sec_new%next
            end if
         end if
         p => p_next
      end do
   end subroutine split_long_sections

!======================================================================
! Slug initiation: replace section with β > β_init by a Slug, conserving
! liquid mass by withdrawing the surplus from neighbours proportionally.
!======================================================================
   subroutine slug_init(g)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer :: p, p_next, slug_left, slug_right, bub_remove, sec_left2, sec_right2
      real(rk) :: L_J, surplus, A, L_slug_old, L_total
      logical  :: trigger_holdup, do_trigger
      integer  :: ntry

      A = 0.25_rk*PI*g%D*g%D
      ntry = 0
      p => g%head
      do while (associated(p))
         p%marked_for_remove = .false.
         p => p%next
      end do
      p => g%head
      do while (associated(p))
         p_next => p%next
         if (p%kind == KIND_SECTION) then
            trigger_holdup = (p%beta > g%beta_init)
            do_trigger = trigger_holdup
         else
            do_trigger = .false.
         end if
         if (do_trigger) then
            L_J = object_length(p)
            surplus = (1.0_rk - p%beta)*L_J*A    ! "missing" liquid to make slug

            ! Look for an adjacent slug across the flanking bubble (paper
            ! page 24): if found, MERGE the flooded section into that slug
            ! rather than creating a new slug.  Mass conservation is
            ! automatic because the absorbed section's β·L·A liquid mass
            ! adds directly to the slug body length.
            slug_left  => null()
            slug_right => null()
            if (associated(p%prev)) then
               if (p%prev%kind == KIND_BUBBLE .and. associated(p%prev%prev)) then
                  if (p%prev%prev%kind == KIND_SLUG) slug_left => p%prev%prev
               end if
            end if
            if (associated(p%next)) then
               if (p%next%kind == KIND_BUBBLE .and. associated(p%next%next)) then
                  if (p%next%next%kind == KIND_SLUG) slug_right => p%next%next
               end if
            end if
            if (associated(slug_left)) then
               ! ---- Case B-left: merge section p into slug_left ----
               ! Layout before:  ... slug_left -- bub_LR -- p (section) -- p%next ...
               ! Layout after :  ... slug_left -- p%next ...
               L_slug_old = object_length(slug_left)
               L_total    = L_slug_old + L_J
               if (L_total > EPS_SMALL) then
                  slug_left%Um = (L_slug_old*slug_left%Um + L_J*p%Ul)/L_total
               end if
               slug_left%zR = p%zR
               slug_left%marked_for_remove = .true.
               ! In Case-B-left p_next (= p%next) is the bubble on
               ! the RIGHT of p, which we KEEP.  Safe to re-use as
               ! the next iteration target.
               bub_remove => p%prev   ! the bubble between slug_left and p
               call remove_object(g, p)
               call remove_object(g, bub_remove)
               ! Invalidate per-Unit gas-mass tracking (see detailed
               ! comment in slug_front_eat): merging a section+bubble
               ! into an existing slug changes the set of bubbles
               ! inside the flanking units without changing n_units
               ! or unit boundary slug IDs, so pm_prepare_gas_track
               ! would otherwise keep a stale tracked mass and drive
               ! the PM pressure correction to its cap.
               g%gas_track_initialized = .false.
               ! Paper §3.7.4 mass conservation: the flooded section's
               ! liquid was β_flooded·L·A, but after merging into the
               ! β=1 slug it contributes L·A.  The surplus
               ! (1-β_flooded)·L·A of "extra" liquid must be withdrawn
               ! from the neighbouring sections (same rule as Case A).
               call distribute_to_neighbors(slug_left, -surplus, A)
               ntry = ntry + 1
               if (ntry > 64) exit
               p => p_next
               cycle
            else if (associated(slug_right)) then
               ! ---- Case B-right: merge section p into slug_right ----
               ! Layout before:  ... p%prev -- p (section) -- bub_RL -- slug_right -- ...
               ! Layout after :  ... p%prev -- slug_right -- ...
               ! Note: p_next (captured at top of loop) == p%next == bub_RL,
               ! which we are about to remove.  After removal, advance to
               ! slug_right (which is now p%prev's next).
               L_slug_old = object_length(slug_right)
               L_total    = L_slug_old + L_J
               if (L_total > EPS_SMALL) then
                  slug_right%Um = (L_slug_old*slug_right%Um + L_J*p%Ul)/L_total
               end if
               ! Case-B-right: the merged slug must absorb the full
               ! length of the flooded section on its LEFT.  Extend
               ! slug_right%zL (which equals slug_right%prev%zR, but
               ! since the intermediate bubble bub_RL will be removed
               ! the list layout automatically gives slug_right%prev
               ! = p%prev after removal; no explicit zL update needed).
               slug_right%marked_for_remove = .true.
               bub_remove => p%next   ! the bubble between p and slug_right
               call remove_object(g, p)
               call remove_object(g, bub_remove)
               ! Invalidate per-Unit gas-mass tracking (see Case-B-left).
               g%gas_track_initialized = .false.
               ! Paper §3.7.4 mass conservation (see Case-B-left).
               call distribute_to_neighbors(slug_right, -surplus, A)
               ! Re-aim p_next at slug_right (still alive) so the loop
               ! continues correctly without dereferencing the freed
               ! bub_remove memory.
               p_next => slug_right
               ntry = ntry + 1
               if (ntry > 64) exit
               p => p_next
               cycle
            end if

            ! ---- Case A: no adjacent slug; create a new slug ----
            if (g%single_slug_only .and. g%birth_patch_done) then
               p => p_next
               cycle
            end if
            p%kind = KIND_SLUG
            p%beta = 1.0_rk
            p%marked_for_remove = .true.
            p%Um = p%Ul
            g%gas_track_initialized = .false.
            if (g%birth_patch_mode > 0 .and. &
                (.not. g%birth_patch_done .or. g%birth_patch_mode >= 10)) then
               call apply_first_birth_patch(g, p)
            end if
            call distribute_to_neighbors(p, -surplus, A)
            ntry = ntry + 1
            if (ntry > 64) exit
         end if
         p => p_next
      end do
      p => g%head
      do while (associated(p))
         p%marked_for_remove = .false.
         p => p%next
      end do
   end subroutine slug_init

   subroutine apply_first_birth_patch(g, slug)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer, intent(inout) :: slug
      type(object_t), pointer :: bL, bR
      real(rk) :: pavg, ravg

      g%birth_patch_done = .true.
      g%birth_patch_slug_id = slug%id
      g%birth_patch_t_birth = g%t_now

      if (g%birth_patch_mode == 1 .or. g%birth_patch_mode == 11) slug%Um = 0.0_rk

      if (g%birth_patch_mode == 2 .or. g%birth_patch_mode == 4 .or. &
          g%birth_patch_mode == 12 .or. g%birth_patch_mode == 14) then
         bL => slug%prev
         bR => slug%next
         if (associated(bL) .and. associated(bR)) then
            if (bL%kind == KIND_BUBBLE .and. bR%kind == KIND_BUBBLE) then
               pavg = 0.5_rk*(bL%pres + bR%pres)
               ravg = 0.5_rk*(bL%rhog + bR%rhog)
               bL%pres = pavg
               bR%pres = pavg
               bL%rhog = ravg
               bR%rhog = ravg
               g%gas_track_initialized = .false.
            end if
         end if
      end if

      g%birth_patch_Um = slug%Um
      g%birth_patch_hold_active = (g%birth_patch_mode == 3 .or. g%birth_patch_mode == 4 .or. &
                                   g%birth_patch_mode == 13 .or. g%birth_patch_mode == 14)
      if (g%birth_patch_mode >= 10 .and. g%birth_patch_n_ids < size(g%birth_patch_ids)) then
         g%birth_patch_n_ids = g%birth_patch_n_ids + 1
         g%birth_patch_ids(g%birth_patch_n_ids) = slug%id
         g%birth_patch_Ums(g%birth_patch_n_ids) = slug%Um
         g%birth_patch_t_births(g%birth_patch_n_ids) = g%t_now
      end if

      write(*,'(A,I0,A,I0,A,F10.6,A,F10.5)') &
         'BIRTH_PATCH mode=', g%birth_patch_mode, ' slug_id=', slug%id, &
         ' t=', g%birth_patch_t_birth, ' Um=', g%birth_patch_Um
   end subroutine apply_first_birth_patch

   subroutine distribute_to_neighbors(slug, dV, A)
      type(object_t), pointer, intent(inout) :: slug
      real(rk),                intent(in)    :: dV, A
      type(object_t), pointer :: lft_bub, rgt_bub, lft_sec, rgt_sec
      real(rk) :: half, Ll, Lr
      lft_bub => slug%prev
      rgt_bub => slug%next
      lft_sec => null(); rgt_sec => null()
      if (associated(lft_bub)) lft_sec => lft_bub%prev
      if (associated(rgt_bub)) rgt_sec => rgt_bub%next
      half = 0.5_rk*dV
      if (associated(lft_sec)) then
         if (lft_sec%kind == KIND_SECTION) then
            Ll = object_length(lft_sec)
            if (Ll > EPS_SMALL) lft_sec%beta = max(lft_sec%beta + half/(Ll*A), BETA_DRY)
         end if
      end if
      if (associated(rgt_sec)) then
         if (rgt_sec%kind == KIND_SECTION) then
            Lr = object_length(rgt_sec)
            if (Lr > EPS_SMALL) rgt_sec%beta = max(rgt_sec%beta + half/(Lr*A), BETA_DRY)
         end if
      end if
   end subroutine distribute_to_neighbors

!======================================================================
! Steep-front absorption (thesis §4.5, paper PNG `47.png`).
!
! After the void-wave step has advanced each Bubble by U_b·δt, a slug-side
! Bubble tagged FRONT (slug eats stratified film ahead) may have overshot
! the next downstream Section (i.e. front_bubble%zR ≥ next_section%zR).
! In LASSI's Lagrangian list this would create a tangled geometry; the
! thesis instead absorbs the overshot sections into the slug body in
! one step ("slug eats several sections in a single δt").
!
! Algorithm:
!   For every SLUG p:
!     bub_R = p%next  (downstream bubble)
!     while  bub_R is FRONT (.not. is_nose)  AND
!            bub_R%zR ≥ next_section%zR - eps:
!         absorb that section's liquid into the slug body (slug%zR <-
!         old next_section%zR), remove (bub_R + section_ahead) from the
!         linked list, and re-aim bub_R at the new downstream bubble.
!
! Mass conservation: the absorbed section contributes β_section·L·A
! of liquid which becomes part of the full-pipe slug body when slug%zR
! is extended.  This implicitly increases the slug's liquid mass by
! exactly the right amount (slug body has β=1, so liquid added per unit
! length is A·ρ_l; the difference (1-β)·L·A of "missing" liquid is
! supplied by surrounding gas-mass conservation handled by the Unit
! gas-mass correction in the next pressure-momentum step).
!
! For NOSE bubbles (the upstream side of the slug, or any FRONT that
! has not actually overshot) this routine is a no-op.
!======================================================================
   subroutine slug_front_eat(g, dt)
      type(grid_t), intent(inout) :: g
      real(rk),     intent(in)    :: dt
      type(object_t), pointer :: p, p_next, bub_R, sec_ahead
      real(rk) :: remaining, z_front, z_cross, d_need, d_take, U_front
      integer :: niter
      real(rk), parameter :: EPS_OVERSHOOT = 1.0e-9_rk
      integer, parameter  :: MAX_EAT_PER_STEP = 100000

      p => g%head
      do while (associated(p))
         p_next => p%next
         if (p%kind == KIND_SLUG) then
            U_front = 0.0_rk
            remaining = 0.0_rk
            niter = 0
            bub_R => p%next
            if (.not. associated(bub_R)) then
               p => p_next
               cycle
            end if
            if (bub_R%kind /= KIND_BUBBLE .or. bub_R%is_nose) then
               p => p_next
               cycle
            end if
            sec_ahead => bub_R%next
            if (.not. associated(sec_ahead)) then
               p => p_next
               cycle
            end if
            if (sec_ahead%kind /= KIND_SECTION) then
               p => p_next
               cycle
            end if
            ! Guard against U_F = (Um − β·Ul)/(1−β) diverging as β → 1:
            ! any section whose holdup is already above β_init will be
            ! converted to a SLUG by the next list_management pass
            ! (slug_init Case A/B).  The hydraulic shock front cannot
            ! consume it in a single δt here — attempting to do so
            ! inflates Um to tens of m/s and spawns the "mega-slug"
            ! artefact observed in case 4.
            if (sec_ahead%beta > g%beta_init - EPS_SMALL) then
               p => p_next
               cycle
            end if
            U_front = max(front_speed(p%Um, sec_ahead%beta, sec_ahead%Ul), 0.0_rk)
            remaining = U_front*max(dt, 0.0_rk)
            if (U_front <= EPS_SMALL) then
               p => p_next
               cycle
            end if
            z_front = bub_R%zR
            do while (remaining > EPS_TINY .and. niter < MAX_EAT_PER_STEP)
               niter = niter + 1
               bub_R => p%next
               if (.not. associated(bub_R)) exit
               if (bub_R%kind /= KIND_BUBBLE) exit
               if (bub_R%is_nose) exit
               sec_ahead => bub_R%next
               if (.not. associated(sec_ahead)) exit
               if (sec_ahead%kind /= KIND_SECTION) exit
               ! Same guard inside the iterative absorption loop.
               if (sec_ahead%beta > g%beta_init - EPS_SMALL) exit
               z_cross = sec_ahead%zR
               d_need = max(z_cross - z_front, 0.0_rk)
               if (d_need <= remaining + EPS_OVERSHOOT) then
                  remaining = remaining - d_need
                  p%zR = z_cross
                  call remove_object(g, sec_ahead)
                  call remove_object(g, bub_R)
                  ! Invalidate per-Unit gas-mass tracking: the unit on
                  ! the slug's downstream side just lost one bubble, so
                  ! the cached tracked mass from the previous PM step no
                  ! longer matches the observed mass.  Without this the
                  ! PM correction loop sees M_obs ≠ M_track every step,
                  ! integrates dp±0.5·p_out per step, and drives the
                  ! flanking bubble pressure to its lower cap (5000 Pa)
                  ! within 20-30 steps — producing the Um≈-13 m/s reverse
                  ! flow and β≈0 dry patches observed in case 4 at t≥9 s.
                  g%gas_track_initialized = .false.
                  if (associated(p%next)) then
                     if (p%next%kind == KIND_BUBBLE) then
                        p%next%zR = z_cross
                        p%next%is_nose = .false.
                     end if
                  end if
                  z_front = z_cross
               else
                  d_take = remaining
                  z_front = z_front + d_take
                  p%zR = z_front
                  bub_R%zR = z_front
                  remaining = 0.0_rk
                  exit
               end if
            end do
            p_next => p%next
         end if
         p => p_next
      end do
   end subroutine slug_front_eat

   pure function front_speed(Um, beta_film, Ul_film) result(UF)
      real(rk), intent(in) :: Um, beta_film, Ul_film
      real(rk) :: UF, b
      b = max(min(beta_film, 1.0_rk - EPS_SMALL), EPS_SMALL)
      UF = (Um - b*Ul_film)/(1.0_rk - b)
   end function front_speed

!======================================================================
! Slug shedding: a slug is suppressed and merged into its left section
! under any of these conditions (thesis page 47 "Slug survives?" +
! Figure 3.9 "For each slug check if the slug survives"):
!   (a) Static collapse: L_slug < 0.05·TargetLength — the body has
!       already vanished (classical collapse branch).
!   (b) Paper's shed>eaten branch: within one δt the tail bubble nose
!       advances faster than the front absorbs film, so the slug body
!       is net-losing liquid and the predicted length L_slug^{n+1} ≤
!       L_min = 0.1·TargetLength.  This catches Um-runaway slugs
!       (e.g. newly-born slugs whose implicit Thomas solve returns a
!       non-physical Um) before they poison downstream sections.
!
!   eaten = β_R·(U_F − U_l,R)          (film ingested at the slug front)
!   shed  = (U_b,nose − U_m)            (liquid shed into the tail bubble)
!   dL/dt = eaten − shed               (slug body length rate of change)
!======================================================================
   subroutine slug_shed(g, dt)
      type(grid_t), intent(inout) :: g
      real(rk),     intent(in)    :: dt
      type(object_t), pointer :: p, p_next, bub_L, bub_R, sec_L, sec_R
      real(rk) :: L, Ls, Ll, beta_new, mom
      real(rk) :: Um, beta_R, Ul_R, U_F, U_nose, W_loc
      real(rk) :: eaten_rate, shed_rate, dLdt
      logical  :: paper_dies
      p => g%head
      do while (associated(p))
         p_next => p%next
         if (p%kind == KIND_SLUG) then
            L = object_length(p)
            ! Paper "shed > eaten" branch
            paper_dies = .false.
            bub_R  => null(); sec_R => null()
            bub_L  => p%prev
            if (associated(p%next)) then
               if (p%next%kind == KIND_BUBBLE) then
                  bub_R => p%next
                  if (associated(bub_R%next)) then
                     if (bub_R%next%kind == KIND_SECTION) sec_R => bub_R%next
                  end if
               end if
            end if
            if (associated(sec_R)) then
               Um     = p%Um
               beta_R = max(min(sec_R%beta, 1.0_rk - EPS_SMALL), EPS_SMALL)
               Ul_R   = sec_R%Ul
               ! Hydraulic-shock front speed (paper §4.5, Eq. 2.12).
               U_F    = safe_div_uf(Um - beta_R*Ul_R, 1.0_rk - beta_R, Um)
               ! Bendiksen-nose speed of the tail bubble (paper §2.2.3).
               W_loc  = wake_effect(L, g%D, 1.0_rk)
               U_nose = W_loc*bendiksen_nose(Um, g%D, g%phi)
               ! Per-unit-length rates (A factors cancel).
               eaten_rate = beta_R*(U_F - Ul_R)
               shed_rate  = U_nose - Um
               dLdt       = eaten_rate - shed_rate
               if (dLdt < 0.0_rk) paper_dies = .true.
            end if
            if (L < 0.05_rk*g%TargetLength .or. paper_dies) then
               bub_L => p%prev
               sec_L => null()
               if (associated(bub_L)) then
                  if (associated(bub_L%prev)) then
                     if (bub_L%prev%kind == KIND_SECTION) sec_L => bub_L%prev
                  end if
               end if
               if (associated(sec_L)) then
                  Ls = max(L, 0.0_rk)
                  Ll = max(object_length(sec_L), 0.0_rk)
                  if (Ll + Ls > EPS_SMALL) then
                     beta_new = (sec_L%beta*Ll + Ls)/(Ll + Ls)
                     mom = sec_L%beta*Ll*sec_L%Ul + Ls*p%Um
                     sec_L%Ul = mom/max(beta_new*(Ll + Ls), EPS_SMALL)
                     sec_L%beta = beta_new
                  end if
                  sec_L%zR = p%zR
                  p_next => p%next
                  call remove_object(g, p)
                  call remove_object(g, bub_L)
                  g%gas_track_initialized = .false.
               else
                  p%kind = KIND_SECTION
                  p%beta = 1.0_rk
                  p%Ul   = p%Um
                  g%gas_track_initialized = .false.
               end if
            end if
         end if
         p => p_next
      end do
   end subroutine slug_shed

!======================================================================
! Merge two adjacent slugs separated by a vanishing bubble into a
! single longer slug.  Triggered when a SLUG-BUBBLE-SLUG pattern has
! the middle bubble shorter than L_M_THR (=0.01·TargetLength): the two
! slugs have effectively coalesced (thesis page 24, "slug coalescence").
!
! Mass-conservation note:
!   When slug_L is extended to slug_R%zR, the new slug body would
!   nominally absorb the FULL central bubble length L_M (β jumps from
!   the bubble's β_M ~ neighbour β to 1).  This injects an extra
!   ρ_l·(1-β_M)·L_M·A of liquid mass into the system.  By restricting
!   triggering to L_M < 0.01·TargetLength = 5 mm at TL=0.5 m, the
!   per-event mass injection is bounded by 0.025 kg (≪ 1 kg/s inflow).
!   The bubble's gas mass is redistributed to the surrounding bubbles
!   in proportion to their lengths.
!
! Velocity:  the merged slug U_m is the length-weighted average:
!             U_m = (L1·Um1 + L2·Um2) / (L1 + L2)
!
! Layout before:  ... bub_LL -- slug_L -- bub_M -- slug_R -- bub_RR ...
! Layout after :  ... bub_LL -- slug_L -- bub_RR ...  (slug_R + bub_M removed,
!                                                      slug_L extended to slug_R%zR)
!======================================================================
   subroutine merge_adjacent_slugs(g)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer :: p, p_next, slug_L, bub_M, slug_R, bub_LL, bub_RR
      real(rk) :: L_M_thr, L1, L2, L_M, A
      real(rk) :: m_g_M, L_LL, L_RR, half
      A = 0.25_rk*PI*g%D*g%D
      L_M_thr = 0.01_rk*g%TargetLength    ! 5 mm at TL=0.5 m
      p => g%head
      do while (associated(p))
         p_next => p%next
         if (p%kind == KIND_SLUG) then
            slug_L => p
            ! Look for slug_L -- bub_M -- slug_R pattern
            if (associated(slug_L%next)) then
               if (slug_L%next%kind == KIND_BUBBLE) then
                  bub_M => slug_L%next
                  if (associated(bub_M%next)) then
                     if (bub_M%next%kind == KIND_SLUG) then
                        slug_R => bub_M%next
                        L_M = object_length(bub_M)
                        if (L_M < L_M_thr) then
                           L1 = max(object_length(slug_L), EPS_SMALL)
                           L2 = max(object_length(slug_R), EPS_SMALL)
                           ! Length-weighted average mixture velocity
                           slug_L%Um = (L1*slug_L%Um + L2*slug_R%Um)/(L1 + L2)
                           ! Extend slug_L to encompass slug_R's range
                           slug_L%zR = slug_R%zR
                           ! Redistribute bub_M's gas mass to flanking bubbles
                           m_g_M = bub_M%rhog*max(1.0_rk - bub_M%beta, EPS_SMALL)*L_M*A
                           bub_LL => null(); bub_RR => null()
                           if (associated(slug_L%prev)) then
                              if (slug_L%prev%kind == KIND_BUBBLE) bub_LL => slug_L%prev
                           end if
                           if (associated(slug_R%next)) then
                              if (slug_R%next%kind == KIND_BUBBLE) bub_RR => slug_R%next
                           end if
                           if (associated(bub_LL) .and. associated(bub_RR)) then
                              L_LL = max(object_length(bub_LL), EPS_SMALL)
                              L_RR = max(object_length(bub_RR), EPS_SMALL)
                              half = 0.5_rk*m_g_M
                              bub_LL%rhog = bub_LL%rhog + half/max(1.0_rk - bub_LL%beta, EPS_SMALL)/(L_LL*A)
                              bub_RR%rhog = bub_RR%rhog + half/max(1.0_rk - bub_RR%beta, EPS_SMALL)/(L_RR*A)
                           else if (associated(bub_LL)) then
                              L_LL = max(object_length(bub_LL), EPS_SMALL)
                              bub_LL%rhog = bub_LL%rhog + m_g_M/max(1.0_rk - bub_LL%beta, EPS_SMALL)/(L_LL*A)
                           else if (associated(bub_RR)) then
                              L_RR = max(object_length(bub_RR), EPS_SMALL)
                              bub_RR%rhog = bub_RR%rhog + m_g_M/max(1.0_rk - bub_RR%beta, EPS_SMALL)/(L_RR*A)
                           end if
                           ! Remove the central bubble and the right slug
                           p_next => slug_R%next   ! advance past slug_R before deleting
                           call remove_object(g, slug_R)
                           call remove_object(g, bub_M)
                           g%gas_track_initialized = .false.
                        end if
                     end if
                  end if
               end if
            end if
         end if
         p => p_next
      end do
   end subroutine merge_adjacent_slugs

!======================================================================
! Inlet treatment:
!   1. Pin INLET at zR = 0.
!   2. Drain any non-{INLET,OUTLET} object whose zR <= 0 (it was pushed
!      backward past the inlet by gas overpressure -- treat it as having
!      exited through the upstream boundary).
!   3. If the first surviving non-INLET object is not a SECTION, insert
!      a fresh inlet-state SECTION (and BUBBLE if needed) to restore the
!      ...SECTION-BUBBLE-SECTION-... alternating invariant.
!   4. If the inlet section grows beyond 2*TargetLength, spawn a fresh
!      inlet-state section + bubble (existing logic).
!======================================================================
   subroutine inlet_treatment(g)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer :: in_obj, sec1, sec_new, p, p_next
      real(rk) :: L
      in_obj => g%head
      if (.not. associated(in_obj)) return
      if (in_obj%kind /= KIND_INLET) return

      ! pin the inlet at z = 0
      in_obj%zR = 0.0_rk

      ! drain any object that drifted past the inlet (zR <= 0); stop at
      ! the first object with zR > 0 (or at the outlet sentinel).
      p => in_obj%next
      do while (associated(p))
         p_next => p%next
         if (p%kind == KIND_OUTLET) exit
         if (p%zR > 0.0_rk) exit
         call remove_object(g, p)
         g%gas_track_initialized = .false.
         p => p_next
      end do

      ! restore alternating SECTION-BUBBLE structure if needed
      p => in_obj%next
      if (associated(p) .and. p%kind == KIND_OUTLET) then
         allocate(sec_new)
         sec_new%kind = KIND_SECTION
         sec_new%zR   = min(g%TargetLength, g%Lpipe)
         sec_new%beta = in_obj%beta
         sec_new%Ul   = in_obj%Ul
         sec_new%pres = in_obj%pres
         sec_new%rhog = in_obj%rhog
         sec_new%Ugs  = in_obj%Ugs
         call insert_after(g, in_obj, sec_new)
         block
            type(object_t), pointer :: bub_new
            allocate(bub_new)
            bub_new%kind = KIND_BUBBLE
            bub_new%zR   = sec_new%zR
            bub_new%beta = sec_new%beta
            bub_new%Ul   = sec_new%Ul
            bub_new%pres = sec_new%pres
            bub_new%rhog = sec_new%rhog
            bub_new%Ugs  = sec_new%Ugs
            call insert_after(g, sec_new, bub_new)
            g%gas_track_initialized = .false.
         end block
      else if (associated(p) .and. p%kind /= KIND_OUTLET) then
         if (p%kind == KIND_SLUG) then
            ! INLET -> SLUG is illegal; insert SECTION + BUBBLE in front
            allocate(sec_new)
            sec_new%kind = KIND_SECTION
            sec_new%zR   = max(p%zR - object_length(p), 0.5_rk*p%zR)
            sec_new%beta = in_obj%beta
            sec_new%Ul   = in_obj%Ul
            sec_new%pres = in_obj%pres
            sec_new%rhog = in_obj%rhog
            sec_new%Ugs  = in_obj%Ugs
            call insert_after(g, in_obj, sec_new)
            block
               type(object_t), pointer :: bub_new
               allocate(bub_new)
               bub_new%kind = KIND_BUBBLE
               bub_new%zR   = sec_new%zR
               bub_new%beta = sec_new%beta
               bub_new%Ul   = sec_new%Ul
               bub_new%pres = sec_new%pres
               bub_new%rhog = sec_new%rhog
               bub_new%Ugs  = sec_new%Ugs
               call insert_after(g, sec_new, bub_new)
               g%gas_track_initialized = .false.
            end block
         else if (p%kind == KIND_BUBBLE) then
            ! INLET -> BUBBLE is illegal; insert SECTION in front
            allocate(sec_new)
            sec_new%kind = KIND_SECTION
            sec_new%zR   = p%zR
            sec_new%beta = in_obj%beta
            sec_new%Ul   = in_obj%Ul
            sec_new%pres = in_obj%pres
            sec_new%rhog = in_obj%rhog
            sec_new%Ugs  = in_obj%Ugs
            call insert_after(g, in_obj, sec_new)
            g%gas_track_initialized = .false.
         end if
      end if

      ! existing spawn logic when the inlet section grows too long
      sec1 => in_obj%next
      if (.not. associated(sec1)) return
      if (sec1%kind /= KIND_SECTION) return
      L = object_length(sec1)
      if (L > 1.5_rk*g%TargetLength) then
         allocate(sec_new)
         sec_new%kind = KIND_SECTION
         sec_new%zR   = in_obj%zR + g%TargetLength
         sec_new%beta = in_obj%beta
         sec_new%Ul   = in_obj%Ul
         sec_new%pres = in_obj%pres
         sec_new%rhog = in_obj%rhog
         sec_new%Ugs  = in_obj%Ugs
         call insert_after(g, in_obj, sec_new)
         block
           type(object_t), pointer :: bub_new
           allocate(bub_new)
           bub_new%kind = KIND_BUBBLE
           bub_new%zR   = sec_new%zR
           bub_new%beta = sec_new%beta
           bub_new%Ul   = sec_new%Ul
           bub_new%pres = sec_new%pres
           bub_new%rhog = sec_new%rhog
           bub_new%Ugs  = sec_new%Ugs
           call insert_after(g, sec_new, bub_new)
           g%gas_track_initialized = .false.
         end block
      end if
   end subroutine inlet_treatment

!======================================================================
! Outlet treatment (minimal, KL 2018 paper §5.4 "p(L,t) = p_atm"):
!   1. Pin the OUTLET object at zR = Lpipe.
!   2. Clip any non-outlet object whose zR drifted past Lpipe back to
!      Lpipe.  The mass between old-zR and Lpipe is the physical
!      outflow through the boundary.
!   3. Remove zero-length non-head/non-tail objects that piled up at
!      the outlet (e.g. a SLUG that fully discharged through Lpipe).
! Nothing else.  The previous version deleted whole slugs as soon as
! the last BUBBLE reached Lpipe (even if the slug itself was still
! several metres upstream), and inserted β=0 buffer sections under
! backward flow; both were over-engineered and unconservative.
!======================================================================
   subroutine outlet_treatment(g)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer :: out_obj, p, p_next
      real(rk) :: zL, L_p
      out_obj => g%tail
      if (.not. associated(out_obj)) return
      if (out_obj%kind /= KIND_OUTLET) return

      ! 1) pin the outlet at Lpipe
      out_obj%zR = g%Lpipe

      ! 2) clip anything that drifted past Lpipe (physical outflow)
      p => g%head
      do while (associated(p))
         p_next => p%next
         if (.not. associated(p, out_obj)) then
            if (p%zR > g%Lpipe) p%zR = g%Lpipe
         end if
         p => p_next
      end do

      ! 3) remove zero-length SECTION/SLUG objects that piled up at outlet.
      !    Skip BUBBLE objects: they are intentional zero-length moving
      !    borders, not removable physical volumes.
      p => g%head
      do while (associated(p))
         p_next => p%next
         if (.not. associated(p, g%head) .and. .not. associated(p, g%tail) .and. &
             p%kind /= KIND_BUBBLE) then
            if (associated(p%prev)) then
               zL = p%prev%zR
            else
               zL = 0.0_rk
            end if
            L_p = p%zR - zL
            if (L_p < EPS_TINY) then
               call remove_object(g, p)
               g%gas_track_initialized = .false.
            end if
         end if
         p => p_next
      end do
   end subroutine outlet_treatment

end module lassi_listmgmt
