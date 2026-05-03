!======================================================================
! lassi_grid.f90
! Linked-list construction and editing primitives for LASSI.
!======================================================================
module lassi_grid
   use lassi_kinds
   use lassi_objects
   implicit none
   private
   public :: grid_t
   public :: grid_init_uniform_stratified
   public :: grid_destroy
   public :: insert_after, insert_before, remove_object
   public :: walk_print
   public :: assign_unit_ids
   public :: total_liquid_volume, total_gas_volume
   public :: object_length, recompute_lengths
   public :: sync_section_gas_from_bubbles, sync_owned_bubbles_from_sections

   type :: grid_t
      type(object_t), pointer :: head => null()
      type(object_t), pointer :: tail => null()
      integer(ik) :: next_id = 0
      ! pipe physics
      real(rk) :: D    = 0.078_rk
      real(rk) :: phi  = 0.0_rk
      real(rk) :: Lpipe= 36.0_rk
      ! design parameters
      real(rk) :: TargetLength = 0.5_rk      ! [m] design cell length
      real(rk) :: beta_init    = 0.98_rk     ! slug-initiation threshold
      real(rk) :: dalpha_dt_w  = 0.0_rk      ! weight of dα/dt term in pressure-CV g_J
      logical  :: enable_slug_coalescence = .true.    ! merge_adjacent_slugs gate
      integer  :: picard_max_iter = 3                  ! A.1 Picard outer iterations
                                                       ! (1 = legacy sequential, 2+ = active)
      real(rk) :: rho_l        = 1000.0_rk
      real(rk) :: lam_l_slug   = 0.02_rk     ! Darcy friction factor inside the slug body
                                              ! (thesis page 21 — fully turbulent water/air;
                                              ! used by the void-wave turning-point criterion
                                              ! Eq. 2.12 and by pm_update_slug_momentum).
      real(rk) :: rho_g_ref    = 1.2_rk      ! reference gas density for IKH κ evaluation
                                              ! (set from EOS at startup; used only as fallback
                                              ! when section gas density cannot be inferred)
      integer(ik) :: gas_track_n_units = 0
      logical :: gas_track_initialized = .false.
      real(rk), allocatable :: gas_track_mass(:)
      integer(ik), allocatable :: gas_track_left_id(:), gas_track_right_id(:)
      ! Diagnostic logging window for the slug-onset acoustic-transient
      ! analysis (zero-cost when not active).  When pm_diag_lun > 0 and
      ! pm_diag_t_lo ≤ current_time ≤ pm_diag_t_hi, slug_momentum_coeffs
      ! and pm_block_solve dump the per-slug A_slug / D_slug0 break-down
      ! plus the new Um after solve, so the offending term in Eq. 3.18
      ! at slug-birth can be identified offline.
      real(rk) :: t_now       = 0.0_rk
      real(rk) :: pm_diag_t_lo = -1.0_rk
      real(rk) :: pm_diag_t_hi = -1.0_rk
      integer  :: pm_diag_lun  = 0
      integer  :: birth_patch_mode = 0
      logical  :: single_slug_only = .false.
      logical  :: birth_patch_done = .false.
      logical  :: birth_patch_hold_active = .false.
      integer(ik) :: birth_patch_slug_id = 0
      integer(ik) :: birth_patch_n_ids = 0
      integer(ik) :: birth_patch_ids(64) = 0
      real(rk) :: birth_patch_Ums(64) = 0.0_rk
      real(rk) :: birth_patch_t_births(64) = -1.0_rk
      real(rk) :: birth_patch_Um = 0.0_rk
      real(rk) :: birth_patch_t_birth = -1.0_rk
   end type grid_t

contains

   ! ------------------------------------------------------------------ !
   ! Build a grid:  Inlet, [N Sections separated by Bubbles], Outlet.
   ! All sections start with constant β0, U_l0; bubbles with p0, rho_g0,
   ! and Ugs0 = (alpha) * U_g0.
   ! ------------------------------------------------------------------ !
   subroutine grid_init_uniform_stratified(g, N, beta0, Ul0, Ug0, p0, rho_g0, &
                                            init_beta_noise)
      type(grid_t),  intent(inout) :: g
      integer,       intent(in)    :: N
      real(rk),      intent(in)    :: beta0, Ul0, Ug0, p0, rho_g0
      real(rk),      intent(in), optional :: init_beta_noise   ! amplitude (e.g. 0.01)

      type(object_t), pointer :: in_obj, out_obj, prev_obj, sec, bub
      real(rk) :: dz, z
      real(rk) :: noise_amp, rb, beta_loc
      integer  :: i

      noise_amp = 0.0_rk
      if (present(init_beta_noise)) noise_amp = init_beta_noise

      ! reset list
      call grid_destroy(g)

      dz = g%Lpipe/real(N, rk)
      z  = 0.0_rk

      ! inlet
      allocate(in_obj)
      g%next_id = g%next_id + 1
      in_obj%id    = g%next_id
      in_obj%kind  = KIND_INLET
      in_obj%zR    = 0.0_rk
      in_obj%L     = 0.0_rk
      in_obj%beta  = beta0
      in_obj%Ul    = Ul0
      in_obj%pres  = p0
      in_obj%rhog  = rho_g0
      in_obj%Ugs   = (1.0_rk - beta0)*Ug0
      g%head => in_obj
      prev_obj => in_obj

      do i = 1, N
         ! Each section gets a one-shot random β perturbation with
         ! amplitude noise_amp (uniform in [-noise_amp, +noise_amp]).
         ! This is the seed required by the IKH instability mechanism
         ! (Renault & Nydal 2006 §3.5/§4.1; KL 2018 §5.4 alludes to
         ! "any perturbations in inlet flow rates / physical properties"
         ! without quantifying it).  The U_l/U_g superficial velocities
         ! at the inlet are kept fixed so that mass flux remains imposed
         ! exactly; only the holdup β is perturbed at t=0.
         if (noise_amp > 0.0_rk) then
            call random_number(rb)
            rb = 2.0_rk*rb - 1.0_rk    ! [-1, +1]
            beta_loc = max(0.05_rk, min(0.95_rk, beta0 + noise_amp*rb))
         else
            beta_loc = beta0
         end if

         ! section
         allocate(sec)
         g%next_id = g%next_id + 1
         sec%id    = g%next_id
         sec%kind  = KIND_SECTION
         sec%zR    = z + dz
         sec%L     = dz
         sec%beta  = beta_loc
         sec%Ul    = Ul0
         sec%pres  = p0
         sec%rhog  = rho_g0
         sec%Ugs   = (1.0_rk - beta_loc)*Ug0
         call link_left_right(prev_obj, sec)
         prev_obj => sec

         ! bubble (right border of section)
         allocate(bub)
         g%next_id = g%next_id + 1
         bub%id    = g%next_id
         bub%kind  = KIND_BUBBLE
         bub%zR    = z + dz
         bub%L     = 0.0_rk
         bub%beta  = beta_loc
         bub%Ul    = Ul0
         bub%pres  = p0
         bub%rhog  = rho_g0
         bub%Ugs   = (1.0_rk - beta_loc)*Ug0
         call link_left_right(prev_obj, bub)
         prev_obj => bub

         z = z + dz
      end do

      ! outlet
      allocate(out_obj)
      g%next_id = g%next_id + 1
      out_obj%id    = g%next_id
      out_obj%kind  = KIND_OUTLET
      out_obj%zR    = g%Lpipe
      out_obj%L     = 0.0_rk
      out_obj%beta  = 0.0_rk          ! separator at constant pressure
      out_obj%Ul    = 0.0_rk
      out_obj%pres  = p0
      out_obj%rhog  = rho_g0
      out_obj%Ugs   = 0.0_rk
      call link_left_right(prev_obj, out_obj)
      g%tail => out_obj

      call assign_unit_ids(g)
   end subroutine grid_init_uniform_stratified

   ! ------------------------------------------------------------------ !
   subroutine grid_destroy(g)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer :: p, q
      p => g%head
      do while (associated(p))
         q => p%next
         deallocate(p)
         p => q
      end do
      g%head => null()
      g%tail => null()
      g%next_id = 0
      if (allocated(g%gas_track_mass)) deallocate(g%gas_track_mass)
      if (allocated(g%gas_track_left_id)) deallocate(g%gas_track_left_id)
      if (allocated(g%gas_track_right_id)) deallocate(g%gas_track_right_id)
      g%gas_track_n_units = 0
      g%gas_track_initialized = .false.
   end subroutine grid_destroy

   ! ------------------------------------------------------------------ !
   subroutine insert_after(g, anchor, newobj)
      type(grid_t),               intent(inout) :: g
      type(object_t), pointer,    intent(inout) :: anchor
      type(object_t), pointer,    intent(inout) :: newobj
      type(object_t), pointer :: rgt
      g%next_id = g%next_id + 1
      newobj%id = g%next_id
      rgt => anchor%next
      newobj%prev => anchor
      newobj%next => rgt
      anchor%next => newobj
      if (associated(rgt)) then
         rgt%prev => newobj
      else
         g%tail => newobj
      end if
   end subroutine insert_after

   ! ------------------------------------------------------------------ !
   subroutine insert_before(g, anchor, newobj)
      type(grid_t),               intent(inout) :: g
      type(object_t), pointer,    intent(inout) :: anchor
      type(object_t), pointer,    intent(inout) :: newobj
      type(object_t), pointer :: lft
      g%next_id = g%next_id + 1
      newobj%id = g%next_id
      lft => anchor%prev
      newobj%next => anchor
      newobj%prev => lft
      anchor%prev => newobj
      if (associated(lft)) then
         lft%next => newobj
      else
         g%head => newobj
      end if
   end subroutine insert_before

   ! ------------------------------------------------------------------ !
   subroutine remove_object(g, obj)
      type(grid_t), intent(inout)        :: g
      type(object_t), pointer, intent(inout) :: obj
      type(object_t), pointer :: lft, rgt
      lft => obj%prev
      rgt => obj%next
      if (associated(lft)) lft%next => rgt
      if (associated(rgt)) rgt%prev => lft
      if (associated(obj, g%head)) g%head => rgt
      if (associated(obj, g%tail)) g%tail => lft
      deallocate(obj)
      obj => null()
   end subroutine remove_object

   ! ------------------------------------------------------------------ !
   ! A Unit = pipe portion bounded by two consecutive Slugs (or by inlet
   ! / outlet at the ends). Used in gas-mass correction.
   subroutine assign_unit_ids(g)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer :: p
      integer(ik) :: uid
      uid = 1
      p => g%head
      do while (associated(p))
         p%unit_id = uid
         if (p%kind == KIND_SLUG) uid = uid + 1
         p => p%next
      end do
   end subroutine assign_unit_ids

   ! ------------------------------------------------------------------ !
   function object_length(p) result(L)
      type(object_t), pointer, intent(in) :: p
      real(rk) :: L
      type(object_t), pointer :: lft
      L = 0.0_rk
      if (.not. associated(p)) return
      if (p%kind == KIND_BUBBLE .or. p%kind == KIND_INLET .or. p%kind == KIND_OUTLET) then
         L = 0.0_rk
         return
      end if
      lft => p%prev
      if (associated(lft)) then
         L = max(p%zR - lft%zR, 0.0_rk)
      else
         L = p%zR
      end if
   end function object_length

   subroutine recompute_lengths(g)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer :: p
      p => g%head
      do while (associated(p))
         p%L = object_length(p)
         p => p%next
      end do
   end subroutine recompute_lengths

   subroutine sync_section_gas_from_bubbles(g)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer :: p, bub
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SECTION) then
            bub => p%next
            if (associated(bub)) then
               if (bub%kind == KIND_BUBBLE) then
                  p%pres = bub%pres
                  p%rhog = bub%rhog
                  p%Ugs  = bub%Ugs
               end if
            end if
         end if
         p => p%next
      end do
   end subroutine sync_section_gas_from_bubbles

   subroutine sync_owned_bubbles_from_sections(g)
      type(grid_t), intent(inout) :: g
      type(object_t), pointer :: p, bub
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SECTION) then
            bub => p%next
            if (associated(bub)) then
               if (bub%kind == KIND_BUBBLE) then
                  bub%beta = p%beta
                  bub%Ul   = p%Ul
               end if
            end if
         end if
         p => p%next
      end do
   end subroutine sync_owned_bubbles_from_sections

   ! ------------------------------------------------------------------ !
   subroutine walk_print(g, lun)
      type(grid_t), intent(in) :: g
      integer,      intent(in) :: lun
      type(object_t), pointer  :: p
      character(8) :: kindstr
      p => g%head
      do while (associated(p))
         select case (p%kind)
         case (KIND_INLET);   kindstr = 'INLET'
         case (KIND_OUTLET);  kindstr = 'OUTLET'
         case (KIND_BUBBLE);  kindstr = 'BUBBLE'
         case (KIND_SECTION); kindstr = 'SECTION'
         case (KIND_SLUG);    kindstr = 'SLUG'
         case default;        kindstr = '?'
         end select
         write(lun, '(I6,1X,A8,1X,F12.6,1X,F12.6,1X,F8.5,1X,F10.5,1X,F10.2,1X,F8.5)') &
            p%id, kindstr, p%zR, object_length(p), p%beta, p%Ul, p%pres, p%rhog
         p => p%next
      end do
   end subroutine walk_print

   ! ------------------------------------------------------------------ !
   function total_liquid_volume(g) result(V)
      type(grid_t), intent(in) :: g
      real(rk) :: V
      type(object_t), pointer :: p
      real(rk) :: A, L
      A = 0.25_rk*PI*g%D*g%D
      V = 0.0_rk
      p => g%head
      do while (associated(p))
         L = object_length(p)
         if (p%kind == KIND_SECTION) V = V + p%beta*L*A
         if (p%kind == KIND_SLUG)    V = V + L*A
         p => p%next
      end do
   end function total_liquid_volume

   function total_gas_volume(g) result(V)
      type(grid_t), intent(in) :: g
      real(rk) :: V
      type(object_t), pointer :: p
      real(rk) :: A, L
      A = 0.25_rk*PI*g%D*g%D
      V = 0.0_rk
      p => g%head
      do while (associated(p))
         L = object_length(p)
         if (p%kind == KIND_SECTION) V = V + (1.0_rk - p%beta)*L*A
         p => p%next
      end do
   end function total_gas_volume

end module lassi_grid
