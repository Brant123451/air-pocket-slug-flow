!======================================================================
! lassi_main.f90
! LASSI driver:  initial conditions, time loop, snapshot writer.
!======================================================================
program lassi_main
   use lassi_kinds
   use lassi_objects
   use lassi_grid
   use lassi_geom
   use lassi_friction
   use lassi_eos
   use lassi_voidwave
   use lassi_press_mom
   use lassi_listmgmt
   use lassi_diag
   use lassi_io
   implicit none

   type(grid_t)            :: g
   type(case_params_t)     :: cp
   type(friction_params_t) :: fp
   type(eos_params_t)      :: eos

   character(256) :: argv, infile, outdir, fname
   integer        :: nargs, lun_mon, lun_diag, lun_mon_multi, lun_track
   integer        :: step, snap_idx, ios
   logical        :: have_multi_mon
   real(rk)       :: t, dt, t_next_snap

   real(rk) :: Ul0, Ug0, rho_g0
   real(rk) :: M_l_init, M_g_init

   ! ------- parse CLI -----------------------------------------------------
   infile = "INPUT_LASSI.txt"
   outdir = "outputs"
   nargs = command_argument_count()
   if (nargs >= 1) call get_command_argument(1, infile)
   if (nargs >= 2) call get_command_argument(2, outdir)

   ! ------- read input ----------------------------------------------------
   call read_input_file(infile, cp)

   ! create output directory (best effort, Windows)
   call execute_command_line('if not exist "'//trim(outdir)//'" mkdir "'//trim(outdir)//'"', wait=.true.)

   ! ------- physics setup -------------------------------------------------
   fp%mu_l    = cp%mu_l
   fp%mu_g    = cp%mu_g
   fp%rough   = cp%rough
   fp%use_AH  = cp%use_AH
   fp%ai_factor = cp%ai_factor
   call set_isothermal_air(eos, cp%T_gas)
   rho_g0 = rho_g(eos, cp%p_out)
   Ul0 = cp%usl/cp%beta0
   Ug0 = cp%usg/(1.0_rk - cp%beta0)

   ! ------- grid setup ----------------------------------------------------
   g%D     = cp%D
   g%phi   = cp%phi
   g%Lpipe = cp%Lpipe
   g%TargetLength = cp%TargetLength
   g%beta_init    = cp%beta_init
   g%beta_ikh     = cp%beta_ikh
   g%kappa_crit   = cp%kappa_crit
   g%dalpha_dt_w  = cp%dalpha_dt_weight
   g%enable_slug_coalescence = cp%enable_slug_coalescence
   g%picard_max_iter = cp%picard_max_iter
   g%birth_patch_mode = cp%birth_patch_mode
   g%single_slug_only = cp%single_slug_only
   g%rho_l        = cp%rho_l
   g%lam_l_slug   = cp%lam_l_slug
   g%rho_g_ref    = rho_g0
   ! Seed the random generator BEFORE the grid build, so the initial-β
   ! perturbation (and the optional inlet white noise) draw from the
   ! same reproducible sequence keyed to perturb_seed.
   call seed_random(cp%perturb_seed)
   call grid_init_uniform_stratified(g, cp%N_init, cp%beta0, Ul0, Ug0, cp%p_out, rho_g0, &
                                     init_beta_noise=cp%init_beta_noise)

   ! optional dam-break initial condition
   if (trim(cp%ic_mode) == 'dambreak') then
      call apply_dambreak_ic(g, cp)
   else if (trim(cp%ic_mode) == 'single_slug') then
      call apply_single_slug_ic(g, cp)
   end if

   ! tell the inlet what to keep imposing
   g%head%beta = cp%beta0
   g%head%Ul   = Ul0
   g%head%Ugs  = (1.0_rk - cp%beta0)*Ug0
   g%head%pres = cp%p_out

   ! (RNG already seeded above before grid_init_uniform_stratified)

   ! ------- time loop -----------------------------------------------------
   t = 0.0_rk
   step = 0
   snap_idx = 0
   t_next_snap = 0.0_rk

   open(newunit=lun_mon, &
        file=trim(outdir)//'/'//trim(cp%tag)//'_monitor.dat', &
        status='replace', action='write', iostat=ios)
   if (ios == 0) call write_monitor_header(lun_mon)

   ! Multi-point monitor (D.2): only opened if at least one x_monitors[k] >= 0
   have_multi_mon = .false.
   if (any(cp%x_monitors >= 0.0_rk)) have_multi_mon = .true.
   lun_mon_multi = 0
   if (have_multi_mon) then
      open(newunit=lun_mon_multi, &
           file=trim(outdir)//'/'//trim(cp%tag)//'_monitor_multi.dat', &
           status='replace', action='write', iostat=ios)
      if (ios == 0) then
         call write_monitor_multi_header(lun_mon_multi, cp%x_monitors)
      else
         lun_mon_multi = 0
      end if
   end if

   open(newunit=lun_diag, &
        file=trim(outdir)//'/'//trim(cp%tag)//'_diag.dat', &
        status='replace', action='write', iostat=ios)
   if (ios /= 0) lun_diag = 0
   if (lun_diag /= 0) call write_diag_header(lun_diag)

   ! Per-slug-track CSV (D.3)
   lun_track = 0
   if (cp%enable_slug_track) then
      open(newunit=lun_track, &
           file=trim(outdir)//'/'//trim(cp%tag)//'_slug_track.csv', &
           status='replace', action='write', iostat=ios)
      if (ios == 0) then
         call write_slug_track_header(lun_track)
      else
         lun_track = 0
      end if
   end if

   ! Slug-onset PM diagnostics (D.4): logs A_slug / D_slug0 break-down
   ! per slug per PM step inside the [pm_diag_t_lo, pm_diag_t_hi] window
   ! so the offending Eq. 3.18 term can be identified offline.  Always
   ! created; activation is via cp%pm_diag_t_lo/hi (default -1 = off).
   block
      integer :: lun_diag_pm
      open(newunit=lun_diag_pm, &
           file=trim(outdir)//'/'//trim(cp%tag)//'_pm_diag.csv', &
           status='replace', action='write', iostat=ios)
      if (ios == 0) then
         write(lun_diag_pm,'(A)') &
            '# t  picard_iter  slug_id  Lslug  Um_old  beta_R  U_l_R  U_F  U_b_back  '// &
            'rhoL_dt  conv_F  conv_B  fric  hydrostat  gravity  intf  '// &
            'A_slug  D_slug0  Um_new  pL  pR'
         g%pm_diag_lun = lun_diag_pm
      end if
      g%pm_diag_t_lo = cp%pm_diag_t_lo
      g%pm_diag_t_hi = cp%pm_diag_t_hi
   end block

   M_l_init = total_liquid_mass(g)
   M_g_init = total_gas_mass(g)

   write(*,'(A)') '====== LASSI start ======'
   write(*,'(A,A)') 'Tag    : ', trim(cp%tag)
   write(*,'(A,F8.4,A)') 'Lpipe  : ', cp%Lpipe, ' m'
   write(*,'(A,F8.5,A)') 'D      : ', cp%D, ' m'
   write(*,'(A,F8.4,A)') 't_end  : ', cp%t_end, ' s'
   write(*,'(A,F8.4)')   'beta0  : ', cp%beta0
   write(*,'(A,F8.4,A,F8.4,A)') 'Inlet U_l, U_g : ', Ul0, ', ', Ug0, ' m/s'
   write(*,'(A,I0)')     'N_init : ', cp%N_init
   write(*,'(A)') '------------------------'

   do while (t < cp%t_end)
      ! 1) compute timestep
      call compute_dt(g, cp%CFL, dt)
      if (t + dt > cp%t_end) dt = cp%t_end - t

      ! Push current time into grid_t so PM diagnostics can window-test.
      g%t_now = t

      ! 2) pressure-momentum step
      call pressure_momentum_step(g, dt, fp, eos, cp%p_out, cp%W_eff)

      ! 3) re-impose inlet boundary state, optionally perturbed.
      !    The INLET object is treated as the thesis inlet bubble by
      !    lassi_voidwave; the first Section is updated only through the
      !    smooth inlet-border remap, not overwritten directly.
      if (cp%usl > 0.0_rk) then
         block
            real(rk) :: rb, ru, beta_p, ul_p, ug_p
            rb = 0.0_rk; ru = 0.0_rk
            if (trim(cp%perturb_mode) == 'sine') then
               rb = sin(2.0_rk*PI*t/max(cp%perturb_period, EPS_SMALL))
               ru = rb
            else
               if (cp%perturb_beta /= 0.0_rk) call random_number(rb)
               if (cp%perturb_Ul   /= 0.0_rk) call random_number(ru)
               rb = 2.0_rk*rb - 1.0_rk
               ru = 2.0_rk*ru - 1.0_rk
            end if
            beta_p = max(0.05_rk, min(0.95_rk, cp%beta0 + cp%perturb_beta*rb))
            ul_p   = cp%usl/max(beta_p, EPS_SMALL) + cp%perturb_Ul*ru
            ug_p   = cp%usg/max(1.0_rk - beta_p, EPS_SMALL)
            g%head%beta = beta_p
            g%head%Ul   = ul_p
            g%head%Ugs  = (1.0_rk - beta_p)*ug_p
         end block
      end if

      ! 4) void wave step (Riemann-driven)
      call void_wave_step(g, dt, cp%W_eff)

      ! 5) topology cleanup
      call list_management(g, dt)

      ! 6) snapshot / monitor / diagnostics
      t = t + dt
      step = step + 1
      if (lun_mon /= 0) call write_monitor(lun_mon, t, g, cp%x_monitor)
      if (lun_mon_multi /= 0) call write_monitor_multi(lun_mon_multi, t, g, cp%x_monitors)
      if (t >= t_next_snap) then
         snap_idx = snap_idx + 1
         write(fname, '(A,A,A,I5.5,A)') trim(outdir)//'/', trim(cp%tag), '_snap_', snap_idx, '.dat'
         call write_snapshot(g, t, trim(fname))
         if (lun_diag  /= 0) call write_diag_step(lun_diag, t, g)
         if (lun_track /= 0) call write_slug_track(lun_track, t, g)
         t_next_snap = t_next_snap + cp%dt_out
         write(*,'(A,F10.5,A,I8,A,F8.5,A,I0)') &
            't=', t, '  step=', step, '  dt=', dt, '  Nobj=', count_objects(g%head)
      end if
   end do

   ! ------- final mass-balance summary -----------------------------------
   block
      real(rk) :: M_l_end, M_g_end, dM_l, dM_g, fill_in_l, fill_in_g, A_pipe
      A_pipe = pipe_area(cp%D)
      M_l_end = total_liquid_mass(g)
      M_g_end = total_gas_mass(g)
      dM_l = M_l_end - M_l_init
      dM_g = M_g_end - M_g_init
      ! Expected inlet inflow over the whole run (constant volumetric rate)
      fill_in_l = cp%rho_l * cp%usl * A_pipe * cp%t_end
      fill_in_g = rho_g0   * cp%usg * A_pipe * cp%t_end
      write(*,'(A)') '------ mass balance ------'
      write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') '  Δm_liq  = ', dM_l, ' kg  (expected inflow =', fill_in_l, ' kg)'
      write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') '  Δm_gas  = ', dM_g, ' kg  (expected inflow =', fill_in_g, ' kg)'
   end block

   if (lun_mon  /= 0) close(lun_mon)
   if (lun_mon_multi /= 0) close(lun_mon_multi)
   if (lun_diag /= 0) close(lun_diag)
   if (lun_track /= 0) close(lun_track)
   call grid_destroy(g)
   write(*,'(A)') '====== LASSI done ======'

contains

   subroutine seed_random(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: a(:)
      call random_seed(size=n)
      allocate(a(n))
      do i = 1, n
         a(i) = seed + i
      end do
      call random_seed(put=a)
      deallocate(a)
   end subroutine seed_random

   subroutine apply_dambreak_ic(g, cp)
      type(grid_t),        intent(inout) :: g
      type(case_params_t), intent(in)    :: cp
      type(object_t), pointer :: p
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SECTION .or. p%kind == KIND_BUBBLE) then
            if (p%zR <= cp%x_dam) then
               p%beta = cp%beta_left
            else
               p%beta = cp%beta_right
            end if
            p%Ul = 0.0_rk
            p%Ugs = 0.0_rk
         end if
         p => p%next
      end do
   end subroutine apply_dambreak_ic

   ! Single-slug IC for isolated slug-front transport tests:
   !   sections with z_R <= x_slug_L  → β=beta_left   (full-supply upstream)
   !   sections in (x_slug_L, x_slug_R] → kind=KIND_SLUG, β=1, Um=Um_slug_init
   !   sections with z_L >= x_slug_R  → β=beta_right (stratified film ahead)
   ! Owned bubbles inherit the section's β/Ul (LASSI grid invariant).
   subroutine apply_single_slug_ic(g, cp)
      type(grid_t),        intent(inout) :: g
      type(case_params_t), intent(in)    :: cp
      type(object_t), pointer :: p
      real(rk) :: zL, zR_p
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SECTION) then
            zL = p%zR - p%L
            zR_p = p%zR
            if (zR_p > cp%x_slug_L .and. zL < cp%x_slug_R) then
               ! section overlaps the prescribed slug interval → make it a slug
               p%kind = KIND_SLUG
               p%beta = 1.0_rk
               p%Um   = cp%Um_slug_init
               p%Ul   = cp%Um_slug_init
            else if (zR_p <= cp%x_slug_L) then
               p%beta = cp%beta_left
            else
               p%beta = cp%beta_right
            end if
         end if
         p => p%next
      end do
      ! Refresh each bubble's β / Ul from the section it owns (immediately to the left).
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_BUBBLE .and. associated(p%prev)) then
            select case (p%prev%kind)
            case (KIND_SECTION)
               p%beta = p%prev%beta
               p%Ul   = p%prev%Ul
            case (KIND_SLUG)
               p%beta = 1.0_rk
               p%Ul   = p%prev%Um
            end select
         end if
         p => p%next
      end do
      ! Reassign Unit IDs after slug placement (units span between slug edges).
      call assign_unit_ids(g)
   end subroutine apply_single_slug_ic

   subroutine compute_dt(g, CFL, dt)
      type(grid_t),  intent(in)  :: g
      real(rk),      intent(in)  :: CFL
      real(rk),      intent(out) :: dt
      type(object_t), pointer :: p
      real(rk) :: c, vmax, Lloc, kappa
      real(rk) :: Ug
      ! Thesis page 47, Figure 3.9, "Update Front" step (verbatim):
      !
      !   "For each slug update the front position and the liquid mass
      !    inside the slug.  The slug front can eat several sections
      !    in a single time step and is NOT limited by the CFL
      !    criterion."
      !
      ! Therefore compute_dt must NOT include a slug-front CFL based on
      ! U_F = (Um − β·Ul)/(1 − β).  The only CFL constraints in the
      ! thesis are:
      !   (1) section-section border CFL   (U_RR − U_LL)·δt < L_J
      !       implemented as a MERGE trigger in list_management, not
      !       as a dt limiter (thesis §3.7.3 / page 47 "Check CFL").
      !   (2) slug body transport CFL       |U_m|·δt < L_slug
      !       (dispersive wave inside the slug).
      !
      ! So here we only apply:
      !   - MSW wave-speed CFL in every SECTION
      !   - slug body CFL in every SLUG
      ! and let slug_front_eat iterate over several sections per δt as
      ! prescribed by the thesis.  Adding a U_F-based CFL here causes
      ! dt to collapse to O(10 μs) whenever an IKH wave crest legitimately
      ! pushes a slug-facing section β toward 1 (see 2026-05-03
      ! case-4 N=1000 run where removing this limiter was necessary to
      ! recover a stable multi-slug regime).
      vmax = 1.0e-3_rk
      dt   = 1.0_rk
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SECTION) then
            Ug = section_ugs_for_dt(p)/max(1.0_rk - p%beta, EPS_SMALL)
            kappa = ikh_kappa(p%beta, g%D, g%phi, g%rho_l, section_rhog_for_dt(p), p%Ul, Ug)
            ! Thesis Eq. 3.12: wave speeds c = ±√(κβ) are real only when
            ! κ > 0 (well-posed).  In the IKH ill-posed regime κ ≤ 0 the
            ! characteristic speeds become imaginary; the void-wave
            ! Riemann already routes such borders to the SSS branch, so
            ! here we simply drop them from the CFL estimate rather than
            ! clamp κ to a non-physical positive value.
            if (kappa > 0.0_rk) then
               c = sqrt(kappa*p%beta)
            else
               c = 0.0_rk
            end if
            vmax = max(abs(p%Ul) + c, vmax)
            Lloc = max(object_length(p), 1.0e-4_rk)
            dt = min(dt, CFL*Lloc/(abs(p%Ul) + c + 1.0e-3_rk))
         else if (p%kind == KIND_SLUG) then
            Lloc = max(object_length(p), 1.0e-4_rk)
            dt = min(dt, CFL*Lloc/(abs(p%Um) + 1.0e-3_rk))
         end if
         p => p%next
      end do
      ! safety bounds
      dt = max(min(dt, 0.05_rk), 1.0e-6_rk)
   end subroutine compute_dt

   function section_rhog_for_dt(sec) result(rhog)
      type(object_t), pointer, intent(in) :: sec
      real(rk) :: rhog
      type(object_t), pointer :: bub
      rhog = sec%rhog
      bub => sec%next
      if (associated(bub)) then
         if (bub%kind == KIND_BUBBLE) rhog = bub%rhog
      end if
   end function section_rhog_for_dt

   function section_ugs_for_dt(sec) result(ugs)
      type(object_t), pointer, intent(in) :: sec
      real(rk) :: ugs
      type(object_t), pointer :: bub
      ugs = sec%Ugs
      bub => sec%next
      if (associated(bub)) then
         if (bub%kind == KIND_BUBBLE) ugs = bub%Ugs
      end if
   end function section_ugs_for_dt

end program lassi_main
