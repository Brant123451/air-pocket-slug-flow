!======================================================================
! lassi_io.f90
! Simple input parsing and snapshot writer for LASSI.
!======================================================================
module lassi_io
   use lassi_kinds
   use lassi_objects
   use lassi_grid
   implicit none
   private
   public :: case_params_t, read_input_file, write_snapshot
   public :: write_monitor_header, write_monitor
   public :: write_monitor_multi_header, write_monitor_multi
   public :: write_slug_track_header, write_slug_track

   type :: case_params_t
      character(256) :: tag      = "lassi_default"
      real(rk) :: D        = 0.078_rk
      real(rk) :: Lpipe    = 36.0_rk
      real(rk) :: phi      = 0.0_rk          ! [rad]
      real(rk) :: rho_l    = 1000.0_rk
      real(rk) :: mu_l     = 1.14e-3_rk
      real(rk) :: mu_g     = 1.79e-5_rk
      real(rk) :: rough    = 0.0_rk
      real(rk) :: T_gas    = 293.15_rk
      ! Andritsos-Hanratty interfacial-friction enhancement (thesis Eq. 2.5,
      ! page 19): λ_i = λ_g·(1 + 75·β)·ai_factor.  Default OFF (use_AH=F)
      ! makes λ_i = λ_g (smooth-interface limit).
      logical  :: use_AH    = .false.
      real(rk) :: ai_factor = 1.0_rk

      real(rk) :: beta0    = 0.5_rk
      real(rk) :: usl      = 1.0_rk          ! superficial liquid velocity at inlet
      real(rk) :: usg      = 2.0_rk          ! superficial gas velocity at inlet
      real(rk) :: p_out    = 1.0e5_rk

      ! optional dam-break initial condition (ic_mode=='dambreak')
      character(32) :: ic_mode  = "uniform"
      real(rk)      :: x_dam    = 18.0_rk
      real(rk)      :: beta_left  = 0.8_rk
      real(rk)      :: beta_right = 0.2_rk
      ! optional single-slug initial condition (ic_mode=='single_slug')
      ! Promotes every section in [x_slug_L, x_slug_R] to a SLUG with
      ! Um = Um_slug_init.  Sections with z_R <= x_slug_L receive
      ! β = beta_left (e.g. 0.99 for "full upstream supply"); sections
      ! with z_L >= x_slug_R receive β = beta_right (the stratified
      ! film the slug front will eat).  Useful for isolated slug-front
      ! transport tests independent of slug-init dynamics.
      real(rk)      :: x_slug_L     = 2.0_rk
      real(rk)      :: x_slug_R     = 2.5_rk
      real(rk)      :: Um_slug_init = 1.0_rk
      real(rk)      :: wave_x0      = 9.0_rk
      real(rk)      :: wave_sigma   = 0.25_rk
      real(rk)      :: wave_amp     = 0.0_rk

      integer  :: N_init   = 200
      real(rk) :: TargetLength = 0.5_rk
      real(rk) :: beta_init    = 0.98_rk
      ! dα/dt time-derivative weight in pressure-CV g_J (thesis Eq. 3.13
      ! "geometric corrections", SPEC line 158).  Mathematically this term
      ! belongs in g_J at weight 1.0, but because β^{n+1} is unavailable
      ! during the staggered (p-only) tridiagonal solve we must extrapolate
      ! dα/dt from the past two β samples — and that explicit extrapolation
      ! has been observed to corrupt liquid mass conservation in K-L case4
      ! (Δm_liq drops from +50 kg to −10 kg at full weight).  The default
      ! 0.0 disables the term; positive values up to ~0.25 are physically
      ! meaningful and should be revisited once the full block-tridiagonal
      ! [p, m, U_m] coupling (A.1) is in place.
      real(rk) :: dalpha_dt_weight = 0.0_rk
      ! Slug-coalescence (thesis page 24): merge two slugs separated by
      ! a vanishingly short bubble.
      logical  :: enable_slug_coalescence = .true.
      ! Per-slug-track CSV (D.3): writes one row per active slug per
      ! snapshot interval, listing (t, slug_id, zL, zR, L, U_m).  Used
      ! offline to build slug-length and slug-frequency PDFs.
      logical  :: enable_slug_track = .false.
      ! Slug-onset PM diagnostic window (D.4).  When pm_diag_t_lo ≥ 0
      ! and t ∈ [pm_diag_t_lo, pm_diag_t_hi] the per-slug A_slug/D_slug0
      ! breakdown plus solved Um are appended to <tag>_pm_diag.csv at
      ! every PM Picard iteration.  Default −1 / −1 ⇒ off.
      real(rk) :: pm_diag_t_lo = -1.0_rk
      real(rk) :: pm_diag_t_hi = -1.0_rk
      ! A.1 Picard outer iteration on (p, m, U_m) self-consistency.
      ! Default 3 = active Picard iteration, which converges p^{n+1}
      ! and U_m^{n+1} together through repeated pressure solve + scalar
      ! slug-momentum update.  Any value 2+ gives a closer approximation
      ! to the Fig. 3.3 block solve; set 1 only for legacy staggered
      ! comparison studies.
      integer  :: picard_max_iter = 3
      integer  :: birth_patch_mode = 0
      logical  :: single_slug_only = .false.
      real(rk) :: W_eff        = 1.0_rk
      real(rk) :: lam_l_slug   = 0.02_rk     ! Darcy λ_l inside slug body (thesis page 21)

      real(rk) :: t_end    = 20.0_rk
      real(rk) :: CFL      = 0.4_rk
      real(rk) :: dt_out   = 0.5_rk
      real(rk) :: x_monitor= 30.0_rk
      ! Optional extra monitor stations (D.2).  Sentinel -1.0 = unused.
      ! Parsed from CSV input "x_monitors = 9.0,18.0,27.0,30.0".
      ! The primary x_monitor stays in the legacy single-point monitor
      ! file; extras land in a separate "<tag>_monitor_multi.dat".
      real(rk) :: x_monitors(8) = -1.0_rk

      ! optional inlet perturbation (white noise on β_in & U_l_in)
      real(rk) :: perturb_beta = 0.0_rk     ! amplitude (e.g. 0.01)
      real(rk) :: perturb_Ul   = 0.0_rk     ! amplitude (e.g. 0.05)
      integer  :: perturb_seed = 12345
      character(16) :: perturb_mode = "noise"
      real(rk) :: perturb_period = 2.0_rk

      ! one-shot initial β perturbation amplitude for IKH instability
      ! seeding (thesis §3.5/§4.1; Renault & Nydal 2006 §3.6.2 alludes
      ! to a "small perturbation").  Default 0.0 ⇒ exact uniform initial
      ! state (which a conservative scheme will leave untouched).
      ! Recommended: 0.01 (±1 % of β0) at t=0, with subsequent inlet
      ! perturb_beta=0 / perturb_Ul=0.
      real(rk) :: init_beta_noise = 0.0_rk
   end type case_params_t

contains

   subroutine read_input_file(fname, cp)
      character(*),         intent(in)    :: fname
      type(case_params_t),  intent(inout) :: cp
      integer :: ios, lun
      character(256) :: line, key, val
      integer :: ieq

      open(newunit=lun, file=trim(fname), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,*) 'WARNING: cannot open ', trim(fname), ' — using defaults.'
         return
      end if
      do
         read(lun, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle
         ieq = index(line, '=')
         if (ieq < 2) cycle
         key = adjustl(line(1:ieq-1))
         val = adjustl(line(ieq+1:))
         call assign_kv(cp, key, val)
      end do
      close(lun)
   end subroutine read_input_file

   subroutine assign_kv(cp, k, v)
      type(case_params_t), intent(inout) :: cp
      character(*),        intent(in)    :: k, v
      character(64) :: kk
      kk = trim(adjustl(k))
      select case (trim(kk))
      case('tag');           cp%tag = trim(v)
      case('D');             read(v,*) cp%D
      case('Lpipe');         read(v,*) cp%Lpipe
      case('phi_deg');       read(v,*) cp%phi; cp%phi = cp%phi*PI/180.0_rk
      case('rho_l');         read(v,*) cp%rho_l
      case('mu_l');          read(v,*) cp%mu_l
      case('mu_g');          read(v,*) cp%mu_g
      case('rough');         read(v,*) cp%rough
      case('T_gas');         read(v,*) cp%T_gas
      case('use_AH');        cp%use_AH = parse_logical(v)
      case('ai_factor');     read(v,*) cp%ai_factor
      case('beta0');         read(v,*) cp%beta0
      case('usl');           read(v,*) cp%usl
      case('usg');           read(v,*) cp%usg
      case('p_out');         read(v,*) cp%p_out
      case('N_init');        read(v,*) cp%N_init
      case('TargetLength');  read(v,*) cp%TargetLength
      case('beta_init');     read(v,*) cp%beta_init
      case('dalpha_dt_weight'); read(v,*) cp%dalpha_dt_weight
      case('enable_slug_coalescence'); cp%enable_slug_coalescence = parse_logical(v)
      case('enable_slug_track'); cp%enable_slug_track = parse_logical(v)
      case('pm_diag_t_lo');  read(v,*) cp%pm_diag_t_lo
      case('pm_diag_t_hi');  read(v,*) cp%pm_diag_t_hi
      case('picard_max_iter'); read(v,*) cp%picard_max_iter
      case('birth_patch_mode'); read(v,*) cp%birth_patch_mode
      case('single_slug_only'); cp%single_slug_only = parse_logical(v)
      case('W_eff');         read(v,*) cp%W_eff
      case('lam_l_slug');    read(v,*) cp%lam_l_slug
      case('t_end');         read(v,*) cp%t_end
      case('CFL');           read(v,*) cp%CFL
      case('dt_out');        read(v,*) cp%dt_out
      case('x_monitor');     read(v,*) cp%x_monitor
      case('x_monitors');    call parse_real_csv(v, cp%x_monitors)
      case('ic_mode');       cp%ic_mode = trim(v)
      case('x_dam');         read(v,*) cp%x_dam
      case('beta_left');     read(v,*) cp%beta_left
      case('beta_right');    read(v,*) cp%beta_right
      case('x_slug_L');      read(v,*) cp%x_slug_L
      case('x_slug_R');      read(v,*) cp%x_slug_R
      case('Um_slug_init');  read(v,*) cp%Um_slug_init
      case('wave_x0');       read(v,*) cp%wave_x0
      case('wave_sigma');    read(v,*) cp%wave_sigma
      case('wave_amp');      read(v,*) cp%wave_amp
      case('perturb_beta');  read(v,*) cp%perturb_beta
      case('perturb_Ul');    read(v,*) cp%perturb_Ul
      case('perturb_seed');  read(v,*) cp%perturb_seed
      case('perturb_mode');  cp%perturb_mode = trim(v)
      case('perturb_period'); read(v,*) cp%perturb_period
      case('init_beta_noise'); read(v,*) cp%init_beta_noise
      case default
         write(*,*) 'INFO: unrecognised key "',trim(kk),'"'
      end select
   end subroutine assign_kv

   ! ------------------------------------------------------------------ !
   ! Parse a YAML/INI-style boolean: T/F, True/False, .true./.false.,
   ! 1/0, yes/no.  Falls back to .false. on unrecognised input.
   pure function parse_logical(s) result(b)
      character(*), intent(in) :: s
      logical :: b
      character(:), allocatable :: t
      integer :: i, n
      ! lowercase trim
      n = len_trim(s)
      allocate(character(len=n) :: t)
      do i = 1, n
         if (s(i:i) >= 'A' .and. s(i:i) <= 'Z') then
            t(i:i) = char(ichar(s(i:i)) + 32)
         else
            t(i:i) = s(i:i)
         end if
      end do
      select case (trim(t))
      case ('t','true','.true.','1','yes','y','on')
         b = .true.
      case default
         b = .false.
      end select
   end function parse_logical

   ! ------------------------------------------------------------------ !
   ! Parse a comma-separated list of reals into the front of array arr.
   ! Trailing arr entries that do not receive a value are left at their
   ! pre-call value (sentinel −1 in our usage).  Used for x_monitors.
   subroutine parse_real_csv(s, arr)
      character(*), intent(in)    :: s
      real(rk),     intent(inout) :: arr(:)
      integer :: i0, i1, n, k, ios
      real(rk) :: val
      character(:), allocatable :: tok
      n  = len_trim(s)
      i0 = 1
      k  = 0
      do while (i0 <= n .and. k < size(arr))
         i1 = index(s(i0:n), ',')
         if (i1 == 0) then
            tok = adjustl(s(i0:n))
            if (len_trim(tok) > 0) then
               read(tok,*,iostat=ios) val
               if (ios == 0) then
                  k = k + 1
                  arr(k) = val
               end if
            end if
            exit
         else
            tok = adjustl(s(i0:i0+i1-2))
            if (len_trim(tok) > 0) then
               read(tok,*,iostat=ios) val
               if (ios == 0) then
                  k = k + 1
                  arr(k) = val
               end if
            end if
            i0 = i0 + i1
         end if
      end do
   end subroutine parse_real_csv

   ! ------------------------------------------------------------------ !
   subroutine write_snapshot(g, t, fname)
      type(grid_t), intent(in) :: g
      real(rk),     intent(in) :: t
      character(*), intent(in) :: fname
      integer :: lun, ios
      type(object_t), pointer :: p
      character(8) :: kindstr
      open(newunit=lun, file=trim(fname), status='replace', action='write', iostat=ios)
      if (ios /= 0) return
      write(lun,'(A)') '# t=, id, kind, zR, L, beta, Ul, Um, p, rho_g, Ugs'
      write(lun,'(A,ES16.8)') '# t = ', t
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
         write(lun, '(I8,1X,A8,1X,F12.6,1X,F12.6,1X,F8.5,1X,F10.5,1X,F10.5,1X,F10.2,1X,F8.5,1X,F10.5)') &
            p%id, kindstr, p%zR, object_length(p), p%beta, p%Ul, p%Um, p%pres, p%rhog, p%Ugs
         p => p%next
      end do
      close(lun)
   end subroutine write_snapshot

   ! ------------------------------------------------------------------ !
   subroutine write_monitor_header(lun)
      integer, intent(in) :: lun
      write(lun,'(A)') '# t   beta_at_xm   Ul_at_xm   slug_count_in_pipe   max_beta'
   end subroutine write_monitor_header

   subroutine write_monitor(lun, t, g, x_m)
      integer,     intent(in) :: lun
      real(rk),    intent(in) :: t, x_m
      type(grid_t),intent(in) :: g
      real(rk) :: beta_m, ul_m, max_beta
      integer  :: nslug
      type(object_t), pointer :: p
      real(rk) :: zL
      ! find object containing x_m
      beta_m = 0.0_rk; ul_m = 0.0_rk; max_beta = 0.0_rk; nslug = 0
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SECTION .or. p%kind == KIND_SLUG) then
            if (associated(p%prev)) then
               zL = p%prev%zR
            else
               zL = 0.0_rk
            end if
            if (x_m >= zL .and. x_m <= p%zR) then
               if (p%kind == KIND_SECTION) then
                  beta_m = p%beta; ul_m = p%Ul
               else
                  beta_m = 1.0_rk; ul_m = p%Um
               end if
            end if
            if (p%beta > max_beta) max_beta = p%beta
         end if
         if (p%kind == KIND_SLUG) nslug = nslug + 1
         p => p%next
      end do
      write(lun,'(F12.5,1X,F8.5,1X,F10.5,1X,I6,1X,F8.5)') t, beta_m, ul_m, nslug, max_beta
   end subroutine write_monitor

   ! ------------------------------------------------------------------ !
   ! Multi-station monitor (D.2).  Writes "t  β_1 Ul_1  β_2 Ul_2  ..."
   ! for every x_m in xm_arr that is ≥ 0 (negative values are sentinels
   ! and skipped).  Header is written by write_monitor_multi_header.
   subroutine write_monitor_multi_header(lun, xm_arr)
      integer,  intent(in) :: lun
      real(rk), intent(in) :: xm_arr(:)
      integer :: k, n_active
      character(2048) :: hdr
      n_active = 0
      hdr = '# t'
      do k = 1, size(xm_arr)
         if (xm_arr(k) >= 0.0_rk) then
            n_active = n_active + 1
            write(hdr, '(A,A,F0.3,A,F0.3)') trim(hdr), '   beta@', xm_arr(k), '  Ul@', xm_arr(k)
         end if
      end do
      if (n_active > 0) write(lun,'(A)') trim(hdr)
   end subroutine write_monitor_multi_header

   ! ------------------------------------------------------------------ !
   ! Per-slug tracking (D.3).  Each call writes one row per active SLUG,
   ! recording: t, slug-id, zL (= prev%zR), zR, L, U_m, β_R (= 1).
   ! Combined with the unique slug%id, post-processors can reconstruct
   ! each slug's birth-time, length history, and front/back velocities,
   ! then derive slug length/frequency PDFs.
   subroutine write_slug_track_header(lun)
      integer, intent(in) :: lun
      write(lun,'(A)') '# t   slug_id   zL[m]   zR[m]   L[m]   Um[m/s]'
   end subroutine write_slug_track_header

   subroutine write_slug_track(lun, t, g)
      integer,      intent(in) :: lun
      real(rk),     intent(in) :: t
      type(grid_t), intent(in) :: g
      type(object_t), pointer :: p
      real(rk) :: zL, L
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SLUG) then
            if (associated(p%prev)) then
               zL = p%prev%zR
            else
               zL = 0.0_rk
            end if
            L = max(p%zR - zL, 0.0_rk)
            write(lun,'(F12.5,1X,I8,1X,F12.5,1X,F12.5,1X,F12.5,1X,F10.5)') &
               t, p%id, zL, p%zR, L, p%Um
         end if
         p => p%next
      end do
   end subroutine write_slug_track

   subroutine write_monitor_multi(lun, t, g, xm_arr)
      integer,      intent(in) :: lun
      real(rk),     intent(in) :: t, xm_arr(:)
      type(grid_t), intent(in) :: g
      type(object_t), pointer :: p
      real(rk) :: beta_m, ul_m, zL, x_m
      integer  :: k, n_active
      character(2048) :: line
      character(64)   :: tok
      n_active = 0
      do k = 1, size(xm_arr)
         if (xm_arr(k) >= 0.0_rk) n_active = n_active + 1
      end do
      if (n_active == 0) return
      write(line,'(F12.5)') t
      do k = 1, size(xm_arr)
         if (xm_arr(k) < 0.0_rk) cycle
         x_m = xm_arr(k)
         beta_m = 0.0_rk; ul_m = 0.0_rk
         p => g%head
         do while (associated(p))
            if (p%kind == KIND_SECTION .or. p%kind == KIND_SLUG) then
               if (associated(p%prev)) then
                  zL = p%prev%zR
               else
                  zL = 0.0_rk
               end if
               if (x_m >= zL .and. x_m <= p%zR) then
                  if (p%kind == KIND_SECTION) then
                     beta_m = p%beta; ul_m = p%Ul
                  else
                     beta_m = 1.0_rk; ul_m = p%Um
                  end if
                  exit
               end if
            end if
            p => p%next
         end do
         write(tok,'(1X,F8.5,1X,F10.5)') beta_m, ul_m
         line = trim(line)//trim(tok)
      end do
      write(lun,'(A)') trim(line)
   end subroutine write_monitor_multi

end module lassi_io
