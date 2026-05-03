!======================================================================
! lassi_press_mom.f90
! Implicit Pressure-Momentum step (FULL gas-momentum form).
!
! Reference: thesis Ch.3 §3.5, Eqs. (3.13)-(3.20); spec §4.
!
! Variables on the staggered grid: at every interior bubble J=1..nb,
!   p_J^{n+1}  (pressure at bubble J)
!   m_J^{n+1} := (rho_g U_g^S)_J^{n+1}  (gas mass flux through bubble J)
!
! Two coupled discrete equations per bubble (Eqs. 3.14 + 3.16-3.17):
!
!   Pressure CV (Eq. 3.14):
!       p_J^{n+1} = z_J ( m_{J-1}^{n+1} - m_J^{n+1} ) + g_J
!     with  z_J = dt * A_pipe / ( V_{g,J} * d(rho_g)/dp )
!           g_J = p_J^n   (extra dV_g/dt corrections lumped into next step)
!
!   Gas momentum (Eq. 3.16, after substituting Eq. 3.14):
!       a_J m_J^{n+1} - b_J m_{J+1}^{n+1} - c_J m_{J-1}^{n+1}
!         + (alpha_J/L_J)( p_{J+1}^{n+1} - p_J^{n+1} )  =  d_J
!     with  a_J = 1/dt + b_J + c_J + (1/8) S_g/A_g lambda_g |U_g|
!           b_J = -min(U_g_{J+1}-U_b_{J+1},0)/L_J + alpha_{J+1}/L_J
!           c_J =  max(U_g_{J  }-U_b_{J  },0)/L_J + alpha_J    /L_J
!           d_J = m_J^n/dt + b_J m_{J+1}^n + c_J m_{J-1}^n
!                 + (1/8) S_i/A_g lambda_i |U_g - U_l|
!                 - rho_g_J g sin(phi)
!
! Substituting p_J^{n+1} (Eq. 3.14) into the (alpha_J/L_J)(p_{J+1}-p_J)
! term yields a clean TRIDIAGONAL in m alone:
!
!       A_J m_J - B_J m_{J+1} - C_J m_{J-1} = D_J
!     with  A_J = a_J + (alpha_J/L_J)*(z_J + z_{J+1})
!           B_J = b_J + (alpha_J/L_J)* z_{J+1}
!           C_J = c_J + (alpha_J/L_J)* z_J
!           D_J = d_J - (alpha_J/L_J)*( g_{J+1} - g_J )
!
! Boundary conditions:
!   * Inlet  J=1: m_0 = m_in given, drop C_1 m_0 to RHS, set C_1=0.
!   * Outlet J=N: p_{N+1} = p_outlet (Dirichlet), so g_{N+1}=p_outlet,
!                 z_{N+1}=0; gas mass flux exits zero-gradient (m_{N+1}=m_N),
!                 absorbed by replacing A_N <- A_N - B_N, B_N=0.
!
! After Thomas, m^{n+1} is known; recover p^{n+1} via Eq. 3.14.
!
! Per-Unit gas-mass correction (Eq. 3.19):
!   For each Unit (a pipe portion bounded by two consecutive slugs, or
!   by the inlet/outlet at the ends), enforce discrete gas-mass balance
!   by adding a uniform pressure correction
!     dp = ( M_obs - M_track ) / ( V_unit * d(rho_g)/dp )
!   to every bubble inside that Unit, then refresh rho_g.
!
! Slug momentum (Eq. 3.18) and the post-PM slug Front/Nose decision are
! NOT yet integrated into the simultaneous Thomas system; for the K-L
! 2018 stratified-inflow case there are no pre-existing slugs so the
! current treatment is sufficient to drive IKH-mediated slug formation.
! The naive slug U_m relaxation is kept for already-formed slugs and is
! flagged TODO for later upgrade.
!
! Liquid velocity update (Eq. 3.20) follows the implicit linearisation
! (1 - dt/rho_l * dF/dU_l) U_l^{n+1} = U_l^n + dt/rho_l * F^n
!======================================================================
module lassi_press_mom
   use lassi_kinds
   use lassi_objects
   use lassi_grid
   use lassi_friction
   use lassi_eos
   use lassi_geom,  only: pipe_area, gamma_from_beta, perimeters, hl_from_gamma
   use lassi_thomas
   implicit none
   private
   public :: pressure_momentum_step

   ! Bendiksen weighting coefficient (W_eff in spec); 1.0 is the
   ! commonly-used default cited by Bendiksen (1984).
   real(rk), parameter :: bendiksen_W = 1.0_rk

   ! D.4 diagnostic break-down of the slug momentum coefficients written
   ! by slug_momentum_coeffs and consumed by pm_block_solve so the
   ! per-slug A_slug / D_slug0 contributions can be logged offline at
   ! slug-onset events.  Single-threaded ⇒ a module-level scratch
   ! variable is safe.  Reset to zero by slug_momentum_coeffs every call.
   type :: slug_diag_t
      real(rk) :: Lslug      = 0.0_rk
      real(rk) :: Um_old     = 0.0_rk
      real(rk) :: beta_R     = 0.0_rk
      real(rk) :: U_l_R      = 0.0_rk
      real(rk) :: U_F        = 0.0_rk
      real(rk) :: U_b_back   = 0.0_rk
      real(rk) :: rhoL_dt    = 0.0_rk      ! ρ L_total / dt   (in A_slug)
      real(rk) :: conv_F     = 0.0_rk      ! ρ (U_F − Um) U_l_R   (in D_slug)
      real(rk) :: conv_B     = 0.0_rk      ! ρ (U_b − Um)         (in A_slug, ×Um implicit)
      real(rk) :: fric       = 0.0_rk      ! 0.5 λ/D L_s ρ |Um|   (in A_slug, ×Um implicit)
      real(rk) :: hydrostat  = 0.0_rk      ! −g ρ (h_R − h_L)     (in D_slug)
      real(rk) :: gravity    = 0.0_rk      ! −ρ g L_s sin φ        (in D_slug)
      real(rk) :: intf       = 0.0_rk      ! interfacial-front     (in D_slug)
   end type slug_diag_t

   type(slug_diag_t) :: pm_last_diag

   ! Per-bubble workspace (allocated once per call to pressure_momentum_step)
   type :: pm_state_t
      ! pointers
      type(ptr_wrap_t), allocatable :: bub_obj(:)
      type(ptr_wrap_t), allocatable :: sec_left(:), sec_right(:)
      ! state at time n (note: pres and m_old are overwritten with n+1
      ! values during the PM step; A.1 Picard outer iteration uses the
      ! _save backups below to restore n-state at the start of each iter)
      real(rk), allocatable :: pres(:)        ! p_J^n  (overwritten -> p_J^{n+1})
      real(rk), allocatable :: m_old(:)       ! (rho_g U_g^S)_J^n (overwritten -> m^{n+1})
      real(rk), allocatable :: rhog_old(:)    ! rho_g_J^n
      real(rk), allocatable :: pres_save(:)   ! immutable p_J^n backup for Picard
      real(rk), allocatable :: m_save(:)      ! immutable m_J^n backup for Picard
      ! geometry
      real(rk), allocatable :: alpha_b(:)     ! alpha at bubble (½(L+R) average)
      real(rk), allocatable :: alpha_L(:)     ! alpha of left  neighbour
      real(rk), allocatable :: alpha_R(:)     ! alpha of right neighbour
      real(rk), allocatable :: dalpha_dt_L(:) ! ∂α_L/∂t (from β_prev history)
      real(rk), allocatable :: dalpha_dt_R(:) ! ∂α_R/∂t (from β_prev history)
      real(rk), allocatable :: half_L(:), half_R(:)  ! half-section lengths around bubble
      real(rk), allocatable :: V_g(:)         ! gas volume in pressure CV
      real(rk), allocatable :: L_b(:)         ! distance between adjacent bubble centres
      ! kinematics
      real(rk), allocatable :: U_g(:)         ! gas velocity at bubble
      real(rk), allocatable :: U_l(:)         ! liquid velocity at bubble (½(L+R))
      real(rk), allocatable :: U_b(:)         ! border velocity (≈ U_l for stratified)
      ! discrete coefficients
      real(rk), allocatable :: zJ(:), gJ(:)   ! Eq. 3.14
      real(rk), allocatable :: aJ(:), bJ(:), cJ(:), dJ(:)   ! Eq. 3.16-3.17
      ! tridiagonal A_J m_J - B_J m_{J+1} - C_J m_{J-1} = D_J
      real(rk), allocatable :: tri_a(:), tri_b(:), tri_c(:), tri_d(:)
   end type pm_state_t

contains

   subroutine pressure_momentum_step(g, dt, fp, eos, p_outlet, W_eff)
      type(grid_t),            intent(inout) :: g
      real(rk),                intent(in)    :: dt
      type(friction_params_t), intent(in)    :: fp
      type(eos_params_t),      intent(in)    :: eos
      real(rk),                intent(in)    :: p_outlet, W_eff

      type(pm_state_t) :: pm
      integer(ik)      :: nb, ib
      real(rk)         :: A_pipe, gamma_eos, m_in

      ! A.1 Picard fixed-point iteration: alternate (p, m) tridiagonal
      ! solve and slug U_m^{n+1} update until self-consistent.  This is
      ! the iterative counterpart of the thesis Fig. 3.3 banded
      ! block-tridiagonal solve.  At convergence (max|ΔU_m| < tol),
      ! p^{n+1} and U_m^{n+1} are consistent.  Empirically the
      ! Picard residual on K-L case4 oscillates around 5–8 mm/s at
      ! steady state, so the tolerance is set to 1 cm/s.
      ! Controlled by g%picard_max_iter (case input picard_max_iter):
      !   1  ⇒ legacy sequential staggered solve (back-compatible)
      !   2+ ⇒ active Picard outer iteration
      real(rk), parameter :: PICARD_TOL_UM  = 0.01_rk    ! [m/s]
      real(rk), parameter :: PICARD_OMEGA   = 0.7_rk     ! under-relaxation
      ! Divergence safety net: at iter k≥2, if any slug's |Um^{(k)}| has
      ! grown by more than PICARD_DIVERGE_FACTOR × |Um^{(k-1)}| OR by
      ! more than PICARD_DIVERGE_ABS [m/s], the Picard map is presumed
      ! non-contractive (slug-onset acoustic regime).  In that case the
      ! solver restores every slug Um to its sane iter-1 value (saved
      ! into Um_safe before under-relaxation) and exits the loop.
      ! Without this, Picard drives U_m to ±10 m/s through pressure
      ! excursions of ±MPa scale, producing the spurious slug Um
      ! oscillation observed in K-L case-4 N=200.
      real(rk), parameter :: PICARD_DIVERGE_FACTOR = 3.0_rk   ! [-]
      real(rk), parameter :: PICARD_DIVERGE_ABS    = 5.0_rk   ! [m/s]
      integer :: picard_iter, n_slug, max_picard_iter
      real(rk) :: max_dUm
      real(rk), allocatable :: Um_prev(:)
      real(rk), allocatable :: Um_safe(:)
      real(rk), allocatable :: gas_track_old(:)
      logical :: gas_track_ready, picard_diverged
      type(object_t), pointer :: qslug

      A_pipe    = pipe_area(g%D)
      gamma_eos = drho_g_dp(eos, p_outlet)         ! isothermal -> constant

      !-----------------------------------------------------------------
      ! 1) Count interior bubbles
      !-----------------------------------------------------------------
      nb = count_kind(g%head, KIND_BUBBLE)
      if (nb == 0) return

      call pm_alloc(pm, nb)

      !-----------------------------------------------------------------
      ! 2) Gather n-state at every bubble (only n-state; U_m updates
      !    happen inside the Picard loop and are read directly from g).
      !-----------------------------------------------------------------
      call pm_gather_state(g, pm, nb)

      ! Snapshot n-state into _save backups so that each Picard iteration
      ! can restore (m_old, pres) before re-assembling the tridiagonal.
      pm%pres_save(:) = pm%pres(:)
      pm%m_save(:)    = pm%m_old(:)

      ! Snapshot slug U_m^n for residual computation
      n_slug = count_kind(g%head, KIND_SLUG)
      if (n_slug > 0) then
         allocate(Um_prev(n_slug))
         allocate(Um_safe(n_slug))
         Um_safe(:) = 0.0_rk
      end if
      picard_diverged = .false.

      ! Freeze every slug's U_m^n into slug%Um_n (paper Eq. 3.18 needs the
      ! explicit U^n in U_F, U_b, |U_m| and the RHS  ρ_l L/dt · U_m^n  to
      ! stay constant during the linear PM solve).  Recomputing these
      ! quantities from the running slug%Um (the iterate of the outer
      ! Picard loop) makes D_slug0 nonlinear in U_m^{(k)} and was the
      ! actual root cause of the Picard divergence at slug onset.
      !
      ! Same logic applies to the interfacial-shear term `intf` in
      ! Eq. 3.18: it depends on the bubble's ρ_g and U_g^S at the
      ! slug front, both at time n.  pm_block_solve writes the just-
      ! solved n+1 values back into every bubble's rhog/Ugs at the
      ! end of each Picard inner step, so without separate _n
      ! snapshots the next inner step's slug_momentum_coeffs would
      ! see a moving target — observed in case-4 as `intf` swinging
      ! by 5×10^7 Pa between Picard iters at slug-cluster onset
      ! (which then drives U_m^{n+1} to ±400 m/s).
      qslug => g%head
      do while (associated(qslug))
         if (qslug%kind == KIND_SLUG) qslug%Um_n = qslug%Um
         if (qslug%kind == KIND_BUBBLE) then
            qslug%rhog_n = qslug%rhog
            qslug%Ugs_n  = qslug%Ugs
         end if
         qslug => qslug%next
      end do

      !=================================================================
      ! Picard outer loop — converges (p, m, U_m) self-consistently
      !=================================================================
      max_dUm = huge(1.0_rk)
      picard_iter = 0
      max_picard_iter = 1
      gas_track_ready = .false.
      picard_loop: do picard_iter = 1, max_picard_iter

      ! Restore n-state at the start of each Picard iteration (since the
      ! previous iter overwrote pm%pres and pm%m_old with n+1 values).
      pm%pres(:)  = pm%pres_save(:)
      pm%m_old(:) = pm%m_save(:)

      ! Save current slug U_m for residual check at end of iteration
      if (n_slug > 0) then
         qslug => g%head; ib = 0
         do while (associated(qslug))
            if (qslug%kind == KIND_SLUG) then
               ib = ib + 1
               Um_prev(ib) = qslug%Um
            end if
            qslug => qslug%next
         end do
      end if

      !-----------------------------------------------------------------
      ! 3) Geometry / kinematics at every bubble
      !-----------------------------------------------------------------
      call pm_compute_geometry(pm, nb, A_pipe, dt)
      if (.not. gas_track_ready) then
         call pm_prepare_gas_track(g, pm, nb, A_pipe, gas_track_old)
         gas_track_ready = .true.
      end if
      call pm_compute_kinematics(g, pm, nb, W_eff)

      !-----------------------------------------------------------------
      ! 4) Pressure-CV coefficients (Eq. 3.14):  z_J, g_J
      !
      !   Gas-mass conservation in the J-th pressure CV
      !     γ V_g (dp/dt) = A (m_{J-1} - m_J) - ρ_g (dV_g/dt)
      !   where  V_g = ½(α_{L} L_{L} + α_{R} L_{R})  (SPEC line 157),
      !          dV_g/dt = ½ A · [(α_L − α_R)·U_b      ← border-motion term
      !                          + L_L · dα_L/dt        ← α(t) on the left
      !                          + L_R · dα_R/dt]       ← α(t) on the right
      !          (SPEC line 158 "geometric corrections").
      !
      !   Discretising in time → Eq. 3.14:
      !     p_J^{n+1} = z_J (m_{J-1}^{n+1} - m_J^{n+1}) + g_J
      !   with
      !     z_J = δt·A / (γ V_g),
      !     g_J = p_J^n
      !          − ½·z_J·ρ_g·(α_L − α_R)·U_b,J          (border-motion)
      !          − (z_J·ρ_g)/(2·V_g/A)·(L_L·dα_L/dt
      !                                + L_R·dα_R/dt)·V_g/A
      !       ≡  p_J^n
      !          − ½·z_J·ρ_g·[(α_L − α_R)·U_b
      !                       + (L_L·dα_L/dt + L_R·dα_R/dt)]
      !-----------------------------------------------------------------
      do ib = 1, nb
         pm%zJ(ib) = dt * A_pipe / max(pm%V_g(ib)*gamma_eos, EPS_SMALL)
         ! dα/dt extrapolation is necessarily explicit (no β^{n+1} available
         ! during the p-only tridiagonal solve), so we apply an under-
         ! relaxation factor DALPHA_REL = 0.25 to the time-derivative
         ! contribution.  This balances physical fidelity (the term IS
         ! present in thesis Eq. 3.13) against the explicit-extrapolation
         ! artefacts that would otherwise corrupt mass conservation.
         pm%gJ(ib) = pm%pres(ib) &
                     - 0.5_rk*pm%zJ(ib)*pm%rhog_old(ib)*( &
                          (pm%alpha_L(ib) - pm%alpha_R(ib))*pm%U_b(ib) &
                        + g%dalpha_dt_w*( pm%half_L(ib)*pm%dalpha_dt_L(ib) &
                                        + pm%half_R(ib)*pm%dalpha_dt_R(ib) ) )
      end do

      !-----------------------------------------------------------------
      ! 5) Gas-momentum coefficients (Eq. 3.16-3.17):  a_J, b_J, c_J, d_J
      !-----------------------------------------------------------------
      call pm_compute_momentum_coeffs(g, pm, nb, dt, fp)

      m_in = pm_inlet_mass_flux(g)
      call pm_block_solve(g, pm, nb, dt, fp, eos, p_outlet, W_eff, m_in, picard_iter)
      call pm_apply_unit_pressure_correction(g, pm, nb, dt, A_pipe, gamma_eos, eos, m_in, gas_track_old)

      !-----------------------------------------------------------------
      ! Picard convergence test: max |U_m^{(k+1)} - U_m^{(k)}| < tol
      ! Under-relaxation slug%Um ← ω·U_m_new + (1-ω)·U_m_prev is only
      ! applied when Picard is active AND the loop will iterate again
      ! (i.e. residual not yet converged AND picard_iter < max_picard_iter).
      ! When max_picard_iter = 1 (legacy sequential), this block reduces
      ! to a residual computation only and leaves slug%Um untouched.
      !-----------------------------------------------------------------
      if (n_slug > 0) then
         max_dUm = 0.0_rk
         qslug => g%head; ib = 0
         do while (associated(qslug))
            if (qslug%kind == KIND_SLUG) then
               ib = ib + 1
               max_dUm = max(max_dUm, abs(qslug%Um - Um_prev(ib)))
               ! Save iter-1 raw block-solve Um as the trusted fallback
               ! used when Picard later diverges (see safety net below).
               if (picard_iter == 1) Um_safe(ib) = qslug%Um
               ! Divergence test: from iter 2 onwards, if Um has grown
               ! by more than PICARD_DIVERGE_FACTOR × |Um_prev| AND
               ! by more than PICARD_DIVERGE_ABS [m/s] relative to the
               ! iteration-start value, the fixed-point map is treated
               ! as non-contractive and we fall back to Um_safe.
               if (picard_iter >= 2) then
                  if (abs(qslug%Um - Um_prev(ib)) > PICARD_DIVERGE_ABS .and. &
                      abs(qslug%Um - Um_prev(ib)) > &
                          PICARD_DIVERGE_FACTOR*abs(Um_prev(ib))) then
                     picard_diverged = .true.
                  end if
               end if
            end if
            qslug => qslug%next
         end do
      else
         max_dUm = 0.0_rk    ! no slugs ⇒ trivially converged
      end if

      ! Picard divergence rescue: restore every slug Um to the safe
      ! iter-1 value and stop iterating.  The (p, m) field will remain
      ! consistent with that Um after one extra Picard iter; for
      ! simplicity we accept the latest (p, m) since they were solved
      ! against the iter-1 Um_safe (within ω·dUm).
      if (picard_diverged .and. n_slug > 0) then
         qslug => g%head; ib = 0
         do while (associated(qslug))
            if (qslug%kind == KIND_SLUG) then
               ib = ib + 1
               qslug%Um = Um_safe(ib)
            end if
            qslug => qslug%next
         end do
         exit picard_loop
      end if

      if (max_dUm < PICARD_TOL_UM) exit picard_loop
      if (picard_iter >= max_picard_iter) exit picard_loop

      ! Apply under-relaxation only when about to take another iteration
      if (n_slug > 0 .and. max_picard_iter > 1) then
         qslug => g%head; ib = 0
         do while (associated(qslug))
            if (qslug%kind == KIND_SLUG) then
               ib = ib + 1
               qslug%Um = PICARD_OMEGA*qslug%Um + (1.0_rk - PICARD_OMEGA)*Um_prev(ib)
            end if
            qslug => qslug%next
         end do
      end if
      end do picard_loop

      if (n_slug > 0) deallocate(Um_prev)
      if (allocated(Um_safe)) deallocate(Um_safe)
      if (gas_track_ready) call pm_commit_gas_track(g, pm, nb, dt, A_pipe, m_in, gas_track_old)
      if (allocated(gas_track_old)) deallocate(gas_track_old)

      !-----------------------------------------------------------------
      ! 12) Liquid velocity update for SECTIONS only (Eq. 3.20).
      !-----------------------------------------------------------------
      call update_liquid_velocity(g, dt, fp)

      !-----------------------------------------------------------------
      ! 13) Update SECTION β_prev history for the next-step dα/dt term
      !     in the pressure-CV g_J coefficient (B.1 / SPEC line 158).
      !-----------------------------------------------------------------
      block
         type(object_t), pointer :: q
         q => g%head
         do while (associated(q))
            if (q%kind == KIND_SECTION) q%beta_prev = q%beta
            q => q%next
         end do
      end block

      call pm_dealloc(pm)
      call sync_section_gas_from_bubbles(g)
      call sync_owned_bubbles_from_sections(g)
   end subroutine pressure_momentum_step

!======================================================================
! Workspace allocation / deallocation
!======================================================================
   subroutine pm_alloc(pm, nb)
      type(pm_state_t), intent(inout) :: pm
      integer(ik),      intent(in)    :: nb
      allocate(pm%bub_obj(nb), pm%sec_left(nb), pm%sec_right(nb))
      allocate(pm%pres(nb), pm%m_old(nb), pm%rhog_old(nb))
      allocate(pm%pres_save(nb), pm%m_save(nb))
      allocate(pm%alpha_b(nb), pm%alpha_L(nb), pm%alpha_R(nb))
      allocate(pm%dalpha_dt_L(nb), pm%dalpha_dt_R(nb))
      allocate(pm%half_L(nb), pm%half_R(nb))
      allocate(pm%V_g(nb), pm%L_b(nb))
      allocate(pm%U_g(nb), pm%U_l(nb), pm%U_b(nb))
      allocate(pm%zJ(nb), pm%gJ(nb))
      allocate(pm%aJ(nb), pm%bJ(nb), pm%cJ(nb), pm%dJ(nb))
      allocate(pm%tri_a(nb), pm%tri_b(nb), pm%tri_c(nb), pm%tri_d(nb))
   end subroutine pm_alloc

   subroutine pm_dealloc(pm)
      type(pm_state_t), intent(inout) :: pm
      if (allocated(pm%bub_obj))    deallocate(pm%bub_obj)
      if (allocated(pm%sec_left))   deallocate(pm%sec_left)
      if (allocated(pm%sec_right))  deallocate(pm%sec_right)
      if (allocated(pm%pres))       deallocate(pm%pres)
      if (allocated(pm%m_old))      deallocate(pm%m_old)
      if (allocated(pm%rhog_old))   deallocate(pm%rhog_old)
      if (allocated(pm%pres_save))  deallocate(pm%pres_save)
      if (allocated(pm%m_save))     deallocate(pm%m_save)
      if (allocated(pm%alpha_b))    deallocate(pm%alpha_b)
      if (allocated(pm%alpha_L))    deallocate(pm%alpha_L)
      if (allocated(pm%alpha_R))    deallocate(pm%alpha_R)
      if (allocated(pm%dalpha_dt_L)) deallocate(pm%dalpha_dt_L)
      if (allocated(pm%dalpha_dt_R)) deallocate(pm%dalpha_dt_R)
      if (allocated(pm%half_L))     deallocate(pm%half_L)
      if (allocated(pm%half_R))     deallocate(pm%half_R)
      if (allocated(pm%V_g))        deallocate(pm%V_g)
      if (allocated(pm%L_b))        deallocate(pm%L_b)
      if (allocated(pm%U_g))        deallocate(pm%U_g)
      if (allocated(pm%U_l))        deallocate(pm%U_l)
      if (allocated(pm%U_b))        deallocate(pm%U_b)
      if (allocated(pm%zJ))         deallocate(pm%zJ)
      if (allocated(pm%gJ))         deallocate(pm%gJ)
      if (allocated(pm%aJ))         deallocate(pm%aJ)
      if (allocated(pm%bJ))         deallocate(pm%bJ)
      if (allocated(pm%cJ))         deallocate(pm%cJ)
      if (allocated(pm%dJ))         deallocate(pm%dJ)
      if (allocated(pm%tri_a))      deallocate(pm%tri_a)
      if (allocated(pm%tri_b))      deallocate(pm%tri_b)
      if (allocated(pm%tri_c))      deallocate(pm%tri_c)
      if (allocated(pm%tri_d))      deallocate(pm%tri_d)
   end subroutine pm_dealloc

!======================================================================
! Step 2: gather n-state
!======================================================================
   subroutine pm_gather_state(g, pm, nb)
      type(grid_t),     intent(in)    :: g
      type(pm_state_t), intent(inout) :: pm
      integer(ik),      intent(in)    :: nb
      type(object_t), pointer :: p
      integer(ik) :: ib
      ib = 0
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_BUBBLE) then
            ib = ib + 1
            pm%bub_obj(ib)%p   => p
            pm%sec_left(ib)%p  => p%prev
            pm%sec_right(ib)%p => p%next
            pm%pres(ib)     = p%pres
            pm%rhog_old(ib) = p%rhog
            pm%m_old(ib)    = p%rhog * p%Ugs
         end if
         p => p%next
      end do
   end subroutine pm_gather_state

!======================================================================
! Step 3a: geometry per pressure CV (alpha, V_g, L_b, dα/dt)
!
! dα/dt = -(β^n - β^{n-1})/dt   for SECTION neighbours with valid β_prev.
! For SLUG neighbours dα/dt is set to 0 (α = EPS, no temporal change).
! On the very first PM step β_prev = -1 (sentinel) ⇒ dα/dt = 0.
!
! The dα/dt estimate is capped to ±DBETA_DT_CAP [1/s] to suppress the
! O(100 [1/s]) transients that occur when a section's β jumps from ~0.5
! to >β_init in a single PM step (slug onset event).  Without a cap the
! corresponding (L · dα/dt) contribution to dV_g/dt overwhelms the
! border-motion (α_L−α_R)·U_b term and breaks liquid mass conservation.
! The dα/dt contribution is applied ONLY when |Δβ| < DBETA_SMOOTH (i.e.
! during smooth wave evolution); during slug-onset events Δβ exceeds
! this threshold and the term is set to zero, on the grounds that:
!   (i) the leading-order (α_L−α_R)·U_b term already captures the
!       border-motion-driven volume change, which is the dominant
!       physics during smooth slug front advection;
!   (ii) the algebraic β jump from a slug-init event is non-physical
!       (it represents a topology-change artefact, not a fluid dα/dt),
!       so including it as a fluid-dynamic source would corrupt the
!       gas-mass balance.
! With DBETA_SMOOTH = 0.05 ordinary KH growth (Δβ < 0.05 per step) is
! retained while slug-onset spikes are filtered out.
!======================================================================
   subroutine pm_compute_geometry(pm, nb, A_pipe, dt)
      type(pm_state_t), intent(inout) :: pm
      integer(ik),      intent(in)    :: nb
      real(rk),         intent(in)    :: A_pipe, dt
      type(object_t), pointer :: sec_L, sec_R
      integer(ik) :: ib
      real(rk)    :: aL, aR, hL, hR, daLdt, daRdt, dbeta
      real(rk), parameter :: DBETA_SMOOTH = 0.05_rk   ! threshold for "smooth" β change
      do ib = 1, nb
         sec_L => pm%sec_left(ib)%p
         sec_R => pm%sec_right(ib)%p
         aL = 0.5_rk; aR = 0.5_rk
         hL = 0.0_rk; hR = 0.0_rk
         daLdt = 0.0_rk; daRdt = 0.0_rk
         if (associated(sec_L)) then
            select case (sec_L%kind)
            case (KIND_SECTION)
               aL = max(1.0_rk - sec_L%beta, EPS_SMALL)
               hL = 0.5_rk*max(object_length(sec_L), EPS_SMALL)
               if (sec_L%beta_prev >= 0.0_rk .and. dt > EPS_SMALL) then
                  dbeta = sec_L%beta - sec_L%beta_prev
                  if (abs(dbeta) < DBETA_SMOOTH) daLdt = -dbeta/dt
               end if
            case (KIND_SLUG)
               aL = EPS_SMALL
               hL = 0.5_rk*max(object_length(sec_L), EPS_SMALL)
            end select
         end if
         if (associated(sec_R)) then
            select case (sec_R%kind)
            case (KIND_SECTION)
               aR = max(1.0_rk - sec_R%beta, EPS_SMALL)
               hR = 0.5_rk*max(object_length(sec_R), EPS_SMALL)
               if (sec_R%beta_prev >= 0.0_rk .and. dt > EPS_SMALL) then
                  dbeta = sec_R%beta - sec_R%beta_prev
                  if (abs(dbeta) < DBETA_SMOOTH) daRdt = -dbeta/dt
               end if
            case (KIND_SLUG)
               aR = EPS_SMALL
               hR = 0.5_rk*max(object_length(sec_R), EPS_SMALL)
            end select
         end if
         pm%alpha_L(ib) = aL
         pm%alpha_R(ib) = aR
         pm%dalpha_dt_L(ib) = daLdt
         pm%dalpha_dt_R(ib) = daRdt
         pm%half_L(ib)  = hL
         pm%half_R(ib)  = hR
         pm%alpha_b(ib) = 0.5_rk*(aL + aR)
         pm%V_g(ib)     = (aL*hL + aR*hR)*A_pipe
         pm%L_b(ib)     = max(hL + hR, EPS_SMALL)
      end do
   end subroutine pm_compute_geometry

!======================================================================
! Step 3b: kinematics at each bubble (U_g, U_l, U_b)
!======================================================================
   subroutine pm_compute_kinematics(g, pm, nb, W_eff_cap)
      type(grid_t),     intent(in)    :: g
      type(pm_state_t), intent(inout) :: pm
      integer(ik),      intent(in)    :: nb
      real(rk),         intent(in)    :: W_eff_cap
      type(object_t), pointer :: sec_L, sec_R
      integer(ik) :: ib
      real(rk)    :: uL, uR, U_crit, U_ls, U_lb, beta_use, W_loc
      logical     :: is_front
      U_crit = u_crit_balance(g%D, g%phi, g%lam_l_slug)
      do ib = 1, nb
         sec_L => pm%sec_left(ib)%p
         sec_R => pm%sec_right(ib)%p
         uL = 0.0_rk; uR = 0.0_rk
         if (associated(sec_L)) then
            select case (sec_L%kind)
            case (KIND_SECTION); uL = sec_L%Ul
            case (KIND_SLUG);    uL = sec_L%Um
            case default;        uL = sec_L%Ul
            end select
         end if
         if (associated(sec_R)) then
            select case (sec_R%kind)
            case (KIND_SECTION); uR = sec_R%Ul
            case (KIND_SLUG);    uR = sec_R%Um
            case default;        uR = sec_R%Ul
            end select
         end if
         pm%U_l(ib) = 0.5_rk*(uL + uR)
         pm%U_b(ib) = pm%U_l(ib)
         if (associated(sec_L) .and. associated(sec_R)) then
            if (sec_L%kind == KIND_SECTION .and. sec_R%kind == KIND_SLUG) then
               U_ls = sec_R%Um
               U_lb = sec_L%Ul
               is_front = (U_ls < U_lb) .and. (U_ls < U_crit)
               if (is_front) then
                  beta_use = max(min(sec_L%beta, 1.0_rk - EPS_SMALL), EPS_SMALL)
                  pm%U_b(ib) = safe_div_ufront(U_ls - beta_use*U_lb, 1.0_rk - beta_use, U_ls)
               else
                  W_loc = wake_effect(object_length(sec_R), g%D, W_eff_cap)
                  pm%U_b(ib) = W_loc*bendiksen_nose(U_ls, g%D, g%phi)
               end if
            else if (sec_L%kind == KIND_SLUG .and. sec_R%kind == KIND_SECTION) then
               U_ls = sec_L%Um
               U_lb = sec_R%Ul
               is_front = (U_ls > U_lb) .and. (U_ls > U_crit)
               if (is_front) then
                  beta_use = max(min(sec_R%beta, 1.0_rk - EPS_SMALL), EPS_SMALL)
                  pm%U_b(ib) = safe_div_ufront(U_ls - beta_use*U_lb, 1.0_rk - beta_use, U_ls)
               else
                  W_loc = wake_effect(object_length(sec_L), g%D, W_eff_cap)
                  pm%U_b(ib) = W_loc*bendiksen_nose(U_ls, g%D, g%phi)
               end if
            end if
         end if
         pm%U_g(ib) = pm%bub_obj(ib)%p%Ugs / max(pm%alpha_b(ib), EPS_SMALL)
      end do
   end subroutine pm_compute_kinematics

!======================================================================
! Step 5: gas-momentum coefficients (Eq. 3.17)
!     a_J = 1/dt + b_J + c_J + (1/8) S_g/A_g lambda_g |U_g|
!     b_J = -min(U_g_{J+1}-U_b_{J+1},0)/L_J + alpha_{J+1}/L_J
!     c_J =  max(U_g_J     -U_b_J     ,0)/L_J + alpha_J    /L_J
!     d_J = m_J^n/dt + (1/8) S_i/A_g lambda_i |U_g - U_l| - rho_g_J g sin(phi)
!
!   NOTE: thesis Eq. 3.17 d_J does NOT contain b_J m^n_{J+1} + c_J m^n_{J-1};
!   those terms enter the LHS of the implicit gas momentum equation
!   Eq. 3.16 (a_J m^{n+1}_J = b_J m^{n+1}_{J+1} + c_J m^{n+1}_{J-1} + d_J)
!   and are advanced explicitly in m_star via the predictor step.
!======================================================================
   subroutine pm_compute_momentum_coeffs(g, pm, nb, dt, fp)
      type(grid_t),            intent(in)    :: g
      type(pm_state_t),        intent(inout) :: pm
      integer(ik),             intent(in)    :: nb
      real(rk),                intent(in)    :: dt
      type(friction_params_t), intent(in)    :: fp
      integer(ik) :: ib
      real(rk) :: A_pipe, alpha, beta_avg, gamma_geom, Sl, Sg, Si
      real(rk) :: A_g, Dg, Re_g, lam_g, lam_i
      real(rk) :: alpha_left, alpha_right
      real(rk) :: U_g_right, U_b_right
      A_pipe = pipe_area(g%D)
      do ib = 1, nb
         beta_avg = max(1.0_rk - pm%alpha_b(ib), EPS_SMALL)
         alpha    = pm%alpha_b(ib)
         A_g      = alpha*A_pipe
         gamma_geom = gamma_from_beta(beta_avg)
         call perimeters(gamma_geom, g%D, Sl, Sg, Si)
         Dg    = 4.0_rk*A_g / max(Sg + Si, EPS_SMALL)
         Re_g  = pm%rhog_old(ib)*abs(pm%U_g(ib) - pm%U_l(ib))*Dg / max(fp%mu_g, EPS_SMALL)
         lam_g = dw_lambda(Re_g, fp%rough/g%D)
         if (fp%use_AH) then
            lam_i = lam_g*(1.0_rk + 75.0_rk*beta_avg)*fp%ai_factor   ! Andritsos-Hanratty (optional)
         else
            lam_i = lam_g                                            ! LASSI thesis default (page 19)
         end if
         ! upwind neighbour data
         if (ib < nb) then
            U_g_right   = pm%U_g(ib+1)
            U_b_right   = pm%U_b(ib+1)
            alpha_right = pm%alpha_b(ib+1)
         else
            U_g_right   = pm%U_g(ib)
            U_b_right   = pm%U_b(ib)
            alpha_right = pm%alpha_b(ib)
         end if
         alpha_left = pm%alpha_b(ib)
         pm%bJ(ib) = (-min(U_g_right - U_b_right, 0.0_rk) + alpha_right)/pm%L_b(ib)
         pm%cJ(ib) = ( max(pm%U_g(ib) - pm%U_b(ib), 0.0_rk) + alpha_left )/pm%L_b(ib)
         pm%aJ(ib) = 1.0_rk/dt + pm%bJ(ib) + pm%cJ(ib) &
                     + 0.125_rk*Sg/max(A_g, EPS_SMALL)*lam_g*abs(pm%U_g(ib))
         pm%dJ(ib) = pm%m_old(ib)/dt &
                     + 0.125_rk*Si/max(A_g, EPS_SMALL)*lam_i*abs(pm%U_g(ib) - pm%U_l(ib)) &
                     - pm%rhog_old(ib)*G_ACC*sin(g%phi)
      end do
   end subroutine pm_compute_momentum_coeffs

   subroutine pm_block_solve(g, pm, nb, dt, fp, eos, p_outlet, W_eff_cap, m_in, iter_label)
      type(grid_t),            intent(inout) :: g
      type(pm_state_t),        intent(inout) :: pm
      integer(ik),             intent(in)    :: nb
      real(rk),                intent(in)    :: dt, p_outlet, W_eff_cap, m_in
      type(friction_params_t), intent(in)    :: fp
      type(eos_params_t),      intent(in)    :: eos
      integer,                 intent(in)    :: iter_label    ! D.4 diag label
      type(ptr_wrap_t), allocatable :: slug_obj(:)
      integer(ik), allocatable :: p_col(:), m_col(:), slug_col(:)
      integer(ik), allocatable :: diag_pL_col(:), diag_pR_col(:)
      logical, allocatable :: slug_row_added(:)
      real(rk), allocatable :: A(:,:), rhs(:)
      integer(ik), allocatable :: block_size(:)
      integer(ik) :: n_slug, nvar, col, row, i, k, idx, idxL, idxR
      integer(ik) :: nblock, sum_blocks
      type(object_t), pointer :: p, sec_L, sec_R
      real(rk) :: kJ, A_slug, D_slug0
      type(slug_diag_t), allocatable :: diag_arr(:)
      logical :: do_log

      n_slug = count_kind(g%head, KIND_SLUG)
      do_log = (g%pm_diag_lun /= 0 .and. n_slug > 0 .and. &
                g%pm_diag_t_lo >= 0.0_rk .and. &
                g%t_now >= g%pm_diag_t_lo .and. g%t_now <= g%pm_diag_t_hi)
      allocate(p_col(nb), m_col(nb))
      if (n_slug > 0) then
         allocate(slug_obj(n_slug), slug_col(n_slug), slug_row_added(n_slug))
         allocate(diag_pL_col(n_slug), diag_pR_col(n_slug))
         allocate(diag_arr(n_slug))
         slug_col(:) = 0
         slug_row_added(:) = .false.
         diag_pL_col(:) = 0
         diag_pR_col(:) = 0
         p => g%head
         k = 0
         do while (associated(p))
            if (p%kind == KIND_SLUG) then
               k = k + 1
               slug_obj(k)%p => p
            end if
            p => p%next
         end do
      else
         allocate(slug_obj(0), slug_col(0), slug_row_added(0))
         allocate(diag_pL_col(0), diag_pR_col(0))
         allocate(diag_arr(0))
      end if

      col = 0
      do i = 1, nb
         col = col + 1; p_col(i) = col
         col = col + 1; m_col(i) = col
         sec_R => pm%sec_right(i)%p
         if (associated(sec_R)) then
            if (sec_R%kind == KIND_SLUG) then
               idx = find_slug_idx(slug_obj, n_slug, sec_R)
               if (idx > 0 .and. slug_col(idx) == 0) then
                  col = col + 1
                  slug_col(idx) = col
               end if
            end if
         end if
      end do
      do idx = 1, n_slug
         if (slug_col(idx) == 0) then
            col = col + 1
            slug_col(idx) = col
         end if
      end do
      nvar = col
      allocate(A(nvar,nvar), rhs(nvar))
      A(:,:) = 0.0_rk
      rhs(:) = 0.0_rk

      row = 0
      do i = 1, nb
         sec_L => pm%sec_left(i)%p
         sec_R => pm%sec_right(i)%p

         row = row + 1
         if (i == nb) then
            call add_coeff(A, row, p_col(i), 1.0_rk)
            rhs(row) = p_outlet
         else if (associated(sec_R) .and. sec_R%kind == KIND_SLUG) then
            idx = find_slug_idx(slug_obj, n_slug, sec_R)
            call add_coeff(A, row, p_col(i), 1.0_rk)
            if (i > 1) then
               call add_coeff(A, row, m_col(i-1), -pm%zJ(i))
               rhs(row) = pm%gJ(i)
            else
               rhs(row) = pm%gJ(i) + pm%zJ(i)*m_in
            end if
            if (idx > 0) call add_coeff(A, row, slug_col(idx), pm%zJ(i)*pm%rhog_old(i))
         else if (associated(sec_L) .and. sec_L%kind == KIND_SLUG) then
            idx = find_slug_idx(slug_obj, n_slug, sec_L)
            call add_coeff(A, row, p_col(i), 1.0_rk)
            call add_coeff(A, row, m_col(i), pm%zJ(i))
            if (idx > 0) call add_coeff(A, row, slug_col(idx), -pm%zJ(i)*pm%rhog_old(i))
            rhs(row) = pm%gJ(i)
         else
            call add_coeff(A, row, p_col(i), 1.0_rk)
            call add_coeff(A, row, m_col(i), pm%zJ(i))
            if (i > 1) then
               call add_coeff(A, row, m_col(i-1), -pm%zJ(i))
               rhs(row) = pm%gJ(i)
            else
               rhs(row) = pm%gJ(i) + pm%zJ(i)*m_in
            end if
         end if

         row = row + 1
         kJ = pm%alpha_b(i)/max(pm%L_b(i), EPS_SMALL)
         call add_coeff(A, row, m_col(i), pm%aJ(i))
         if (i > 1) then
            call add_coeff(A, row, m_col(i-1), -pm%cJ(i))
            rhs(row) = pm%dJ(i)
         else
            rhs(row) = pm%dJ(i) + pm%cJ(i)*m_in
         end if
         if (i < nb) then
            call add_coeff(A, row, m_col(i+1), -pm%bJ(i))
            call add_coeff(A, row, p_col(i+1), kJ)
         else
            call add_coeff(A, row, m_col(i), -pm%bJ(i))
            rhs(row) = rhs(row) - kJ*p_outlet
         end if
         call add_coeff(A, row, p_col(i), -kJ)

         if (associated(sec_R)) then
            if (sec_R%kind == KIND_SLUG) then
               idx = find_slug_idx(slug_obj, n_slug, sec_R)
               if (idx > 0) then
                  row = row + 1
                  call slug_momentum_coeffs(g, sec_R, dt, fp, W_eff_cap, A_slug, D_slug0)
                  if (do_log) diag_arr(idx) = pm_last_diag
                  call add_coeff(A, row, slug_col(idx), A_slug)
                  idxL = find_bubble_idx(pm, nb, sec_R%prev)
                  idxR = find_bubble_idx(pm, nb, sec_R%next)
                  if (idxL > 0) then
                     call add_coeff(A, row, p_col(idxL), -1.0_rk)
                     if (do_log) diag_pL_col(idx) = p_col(idxL)
                  end if
                  if (idxR > 0) then
                     call add_coeff(A, row, p_col(idxR), 1.0_rk)
                     if (do_log) diag_pR_col(idx) = p_col(idxR)
                  else
                     rhs(row) = rhs(row) - p_outlet
                  end if
                  rhs(row) = rhs(row) + D_slug0
                  slug_row_added(idx) = .true.
               end if
            end if
         end if
      end do

      do idx = 1, n_slug
         if (.not. slug_row_added(idx)) then
            row = row + 1
            call add_coeff(A, row, slug_col(idx), 1.0_rk)
            rhs(row) = slug_obj(idx)%p%Um
         end if
      end do

      if (row /= nvar) then
         do i = row + 1, nvar
            call add_coeff(A, i, i, 1.0_rk)
            rhs(i) = 0.0_rk
         end do
      end if

      allocate(block_size(nb+1))
      sum_blocks = 0
      do i = 1, nb
         block_size(i) = 2
         sec_R => pm%sec_right(i)%p
         if (associated(sec_R)) then
            if (sec_R%kind == KIND_SLUG) then
               idx = find_slug_idx(slug_obj, n_slug, sec_R)
               if (idx > 0 .and. slug_col(idx) > 0) block_size(i) = 3
            end if
         end if
         sum_blocks = sum_blocks + block_size(i)
      end do
      nblock = nb
      if (sum_blocks < nvar) then
         nblock = nb + 1
         block_size(nblock) = nvar - sum_blocks
         sum_blocks = nvar
      end if
      if (sum_blocks == nvar) then
         call solve_block_thomas(A, rhs, block_size(1:nblock), nblock)
      else
         call solve_banded_dense(A, rhs, nvar, min(block_bandwidth(A, nvar) + 8, nvar-1))
      end if

      do i = 1, nb
         pm%pres(i) = rhs(p_col(i))
         pm%m_old(i) = rhs(m_col(i))
      end do
      pm%pres(nb) = p_outlet
      do idx = 1, n_slug
         slug_obj(idx)%p%Um = rhs(slug_col(idx))
      end do
      do i = 1, nb
         pm%bub_obj(i)%p%pres = pm%pres(i)
         pm%bub_obj(i)%p%rhog = rho_g(eos, pm%bub_obj(i)%p%pres)
         pm%bub_obj(i)%p%Ugs  = pm%m_old(i)/max(pm%bub_obj(i)%p%rhog, EPS_SMALL)
      end do

      ! D.4 emit per-slug PM diagnostic (one row per slug per pm_block_solve call).
      if (do_log) then
         block
            real(rk) :: pL_val, pR_val, A_slug_val, D_slug_val
            real(rk) :: VgL_val, VgR_val, zJL_val, zJR_val, gJL_val, gJR_val
            real(rk) :: mL_val, mR_val
            real(rk) :: Um_new
            do idx = 1, n_slug
               if (.not. slug_row_added(idx)) cycle
               pL_val = 0.0_rk; pR_val = p_outlet
               VgL_val = 0.0_rk; VgR_val = 0.0_rk
               zJL_val = 0.0_rk; zJR_val = 0.0_rk
               gJL_val = 0.0_rk; gJR_val = 0.0_rk
               mL_val = 0.0_rk; mR_val = 0.0_rk
               if (diag_pL_col(idx) > 0) pL_val = rhs(diag_pL_col(idx))
               if (diag_pR_col(idx) > 0) pR_val = rhs(diag_pR_col(idx))
               if (associated(slug_obj(idx)%p%prev)) then
                  i = find_bubble_idx(pm, nb, slug_obj(idx)%p%prev)
                  if (i > 0) then
                     VgL_val = pm%V_g(i); zJL_val = pm%zJ(i); gJL_val = pm%gJ(i)
                     mL_val = rhs(m_col(i))
                  end if
               end if
               if (associated(slug_obj(idx)%p%next)) then
                  i = find_bubble_idx(pm, nb, slug_obj(idx)%p%next)
                  if (i > 0) then
                     VgR_val = pm%V_g(i); zJR_val = pm%zJ(i); gJR_val = pm%gJ(i)
                     mR_val = rhs(m_col(i))
                  end if
               end if
               Um_new = rhs(slug_col(idx))
               A_slug_val = diag_arr(idx)%rhoL_dt &
                          + diag_arr(idx)%conv_B &
                          + diag_arr(idx)%fric
               D_slug_val = diag_arr(idx)%rhoL_dt * diag_arr(idx)%Um_old &
                          + diag_arr(idx)%conv_F &
                          + diag_arr(idx)%hydrostat &
                          + diag_arr(idx)%gravity &
                          + diag_arr(idx)%intf
               write(g%pm_diag_lun, '(F11.6,1X,I3,1X,I8,18(1X,ES13.5E2))') &
                  g%t_now, iter_label, slug_obj(idx)%p%id, &
                  diag_arr(idx)%Lslug, diag_arr(idx)%Um_old, &
                  diag_arr(idx)%beta_R, diag_arr(idx)%U_l_R, &
                  diag_arr(idx)%U_F, diag_arr(idx)%U_b_back, &
                  diag_arr(idx)%rhoL_dt, diag_arr(idx)%conv_F, &
                  diag_arr(idx)%conv_B, diag_arr(idx)%fric, &
                  diag_arr(idx)%hydrostat, diag_arr(idx)%gravity, &
                  diag_arr(idx)%intf, A_slug_val, D_slug_val, Um_new, &
                  pL_val, pR_val
            end do
            flush(g%pm_diag_lun)
         end block
      end if

      deallocate(A, rhs, p_col, m_col, slug_obj, slug_col, slug_row_added, block_size)
      if (allocated(diag_pL_col)) deallocate(diag_pL_col)
      if (allocated(diag_pR_col)) deallocate(diag_pR_col)
      if (allocated(diag_arr))    deallocate(diag_arr)
   end subroutine pm_block_solve

   integer(ik) function pm_count_units(pm, nb) result(n_units)
      type(pm_state_t), intent(in) :: pm
      integer(ik),      intent(in) :: nb
      integer(ik) :: ib_last, ib_scan

      n_units = 0
      ib_scan = 1
      do while (ib_scan <= nb)
         n_units = n_units + 1
         ib_last = ib_scan
         do while (ib_last < nb)
            if (associated(pm%sec_right(ib_last)%p)) then
               if (pm%sec_right(ib_last)%p%kind == KIND_SLUG) exit
            end if
            ib_last = ib_last + 1
         end do
         ib_scan = ib_last + 1
      end do
   end function pm_count_units

   subroutine pm_collect_unit_state(pm, nb, A_pipe, n_units, M_obs, V_unit, left_id, right_id)
      type(pm_state_t), intent(in) :: pm
      integer(ik),      intent(in) :: nb
      real(rk),         intent(in) :: A_pipe
      integer(ik),      intent(out) :: n_units
      real(rk), allocatable, intent(out) :: M_obs(:), V_unit(:)
      integer(ik), allocatable, intent(out) :: left_id(:), right_id(:)
      integer(ik) :: ib, ib_first, ib_last, ib_scan, iu

      n_units = pm_count_units(pm, nb)
      allocate(M_obs(n_units), V_unit(n_units), left_id(n_units), right_id(n_units))
      M_obs = 0.0_rk
      V_unit = 0.0_rk
      left_id = -1
      right_id = -2

      ib_scan = 1
      iu = 0
      do while (ib_scan <= nb)
         iu = iu + 1
         ib_first = ib_scan
         if (associated(pm%sec_left(ib_first)%p)) then
            if (pm%sec_left(ib_first)%p%kind == KIND_SLUG) left_id(iu) = pm%sec_left(ib_first)%p%id
         end if

         ib_last = ib_first
         do while (ib_last < nb)
            if (associated(pm%sec_right(ib_last)%p)) then
               if (pm%sec_right(ib_last)%p%kind == KIND_SLUG) exit
            end if
            ib_last = ib_last + 1
         end do

         if (associated(pm%sec_right(ib_last)%p)) then
            if (pm%sec_right(ib_last)%p%kind == KIND_SLUG) right_id(iu) = pm%sec_right(ib_last)%p%id
         end if

         do ib = ib_first, ib_last
            M_obs(iu)  = M_obs(iu)  + pm%bub_obj(ib)%p%rhog * pm%alpha_b(ib) * pm%L_b(ib) * A_pipe
            V_unit(iu) = V_unit(iu) + pm%alpha_b(ib) * pm%L_b(ib) * A_pipe
         end do
         ib_scan = ib_last + 1
      end do
   end subroutine pm_collect_unit_state

   subroutine pm_prepare_gas_track(g, pm, nb, A_pipe, track_old)
      type(grid_t),     intent(inout) :: g
      type(pm_state_t), intent(in)    :: pm
      integer(ik),      intent(in)    :: nb
      real(rk),         intent(in)    :: A_pipe
      real(rk), allocatable, intent(out) :: track_old(:)
      real(rk), allocatable :: M_obs(:), V_unit(:)
      integer(ik), allocatable :: left_id(:), right_id(:)
      integer(ik) :: n_units, iu
      logical :: reset_track

      call pm_collect_unit_state(pm, nb, A_pipe, n_units, M_obs, V_unit, left_id, right_id)
      reset_track = .not. g%gas_track_initialized
      if (.not. allocated(g%gas_track_mass)) reset_track = .true.
      if (.not. allocated(g%gas_track_left_id)) reset_track = .true.
      if (.not. allocated(g%gas_track_right_id)) reset_track = .true.
      if (g%gas_track_n_units /= n_units) reset_track = .true.
      if (.not. reset_track) then
         do iu = 1, n_units
            if (g%gas_track_left_id(iu) /= left_id(iu)) reset_track = .true.
            if (g%gas_track_right_id(iu) /= right_id(iu)) reset_track = .true.
         end do
      end if

      if (reset_track) then
         if (allocated(g%gas_track_mass)) deallocate(g%gas_track_mass)
         if (allocated(g%gas_track_left_id)) deallocate(g%gas_track_left_id)
         if (allocated(g%gas_track_right_id)) deallocate(g%gas_track_right_id)
         allocate(g%gas_track_mass(n_units), g%gas_track_left_id(n_units), g%gas_track_right_id(n_units))
         g%gas_track_mass(:) = M_obs(:)
         g%gas_track_left_id(:) = left_id(:)
         g%gas_track_right_id(:) = right_id(:)
         g%gas_track_n_units = n_units
         g%gas_track_initialized = .true.
      end if

      allocate(track_old(n_units))
      track_old(:) = g%gas_track_mass(:)
      deallocate(M_obs, V_unit, left_id, right_id)
   end subroutine pm_prepare_gas_track

   subroutine pm_compute_unit_targets(pm, nb, dt, A_pipe, m_in, track_old, target_mass)
      type(pm_state_t), intent(in) :: pm
      integer(ik),      intent(in) :: nb
      real(rk),         intent(in) :: dt, A_pipe, m_in
      real(rk),         intent(in) :: track_old(:)
      real(rk), allocatable, intent(out) :: target_mass(:)
      integer(ik) :: ib_first, ib_last, ib_scan, iu, n_units
      real(rk) :: m_unit_in, m_unit_out
      logical :: left_at_pipe_inlet, right_at_pipe_outlet
      logical :: left_blocked_by_slug, right_blocked_by_slug

      n_units = pm_count_units(pm, nb)
      allocate(target_mass(n_units))
      target_mass = 0.0_rk

      ib_scan = 1
      iu = 0
      do while (ib_scan <= nb)
         iu = iu + 1
         ib_first = ib_scan
         left_blocked_by_slug = .false.
         if (associated(pm%sec_left(ib_first)%p)) then
            if (pm%sec_left(ib_first)%p%kind == KIND_SLUG) left_blocked_by_slug = .true.
         end if
         left_at_pipe_inlet = (ib_first == 1) .and. (.not. left_blocked_by_slug)

         ib_last = ib_first
         do while (ib_last < nb)
            if (associated(pm%sec_right(ib_last)%p)) then
               if (pm%sec_right(ib_last)%p%kind == KIND_SLUG) exit
            end if
            ib_last = ib_last + 1
         end do

         right_blocked_by_slug = .false.
         if (associated(pm%sec_right(ib_last)%p)) then
            if (pm%sec_right(ib_last)%p%kind == KIND_SLUG) right_blocked_by_slug = .true.
         end if
         right_at_pipe_outlet = (ib_last == nb) .and. (.not. right_blocked_by_slug)

         if (left_at_pipe_inlet) then
            m_unit_in = m_in
         else
            m_unit_in = 0.0_rk
         end if
         if (right_at_pipe_outlet) then
            m_unit_out = pm%m_save(ib_last)
         else
            m_unit_out = 0.0_rk
         end if

         if (iu <= size(track_old)) target_mass(iu) = max(0.0_rk, &
            track_old(iu) + dt*(m_unit_in - m_unit_out)*A_pipe)
         ib_scan = ib_last + 1
      end do
   end subroutine pm_compute_unit_targets

   subroutine pm_commit_gas_track(g, pm, nb, dt, A_pipe, m_in, track_old)
      type(grid_t),     intent(inout) :: g
      type(pm_state_t), intent(in)    :: pm
      integer(ik),      intent(in)    :: nb
      real(rk),         intent(in)    :: dt, A_pipe, m_in
      real(rk),         intent(in)    :: track_old(:)
      real(rk), allocatable :: M_obs(:), V_unit(:), target_mass(:)
      integer(ik), allocatable :: left_id(:), right_id(:)
      integer(ik) :: n_units

      call pm_collect_unit_state(pm, nb, A_pipe, n_units, M_obs, V_unit, left_id, right_id)
      if (size(track_old) == n_units) then
         call pm_compute_unit_targets(pm, nb, dt, A_pipe, m_in, track_old, target_mass)
      else
         allocate(target_mass(n_units))
         target_mass(:) = M_obs(:)
      end if

      if (allocated(g%gas_track_mass)) deallocate(g%gas_track_mass)
      if (allocated(g%gas_track_left_id)) deallocate(g%gas_track_left_id)
      if (allocated(g%gas_track_right_id)) deallocate(g%gas_track_right_id)
      allocate(g%gas_track_mass(n_units), g%gas_track_left_id(n_units), g%gas_track_right_id(n_units))
      g%gas_track_mass(:) = target_mass(:)
      g%gas_track_left_id(:) = left_id(:)
      g%gas_track_right_id(:) = right_id(:)
      g%gas_track_n_units = n_units
      g%gas_track_initialized = .true.

      deallocate(M_obs, V_unit, left_id, right_id, target_mass)
   end subroutine pm_commit_gas_track

   subroutine pm_add_unit_correction_to_gJ(pm, nb, dt, A_pipe, gamma_eos, m_in, p_outlet, track_old)
      type(pm_state_t), intent(inout) :: pm
      integer(ik),      intent(in)    :: nb
      real(rk),         intent(in)    :: dt, A_pipe, gamma_eos, m_in, p_outlet
      real(rk),         intent(in)    :: track_old(:)
      integer(ik) :: ib, ib_first, ib_last, ib_scan, iu, n_units
      real(rk) :: dp
      real(rk), allocatable :: M_obs(:), V_unit(:), target_mass(:)
      integer(ik), allocatable :: left_id(:), right_id(:)

      if (nb < 1) return
      call pm_collect_unit_state(pm, nb, A_pipe, n_units, M_obs, V_unit, left_id, right_id)
      if (size(track_old) == n_units) then
         call pm_compute_unit_targets(pm, nb, dt, A_pipe, m_in, track_old, target_mass)
      else
         allocate(target_mass(n_units))
         target_mass(:) = M_obs(:)
      end if

      ib_scan = 1
      iu = 0
      do while (ib_scan <= nb)
         iu = iu + 1
         ib_first = ib_scan

         ib_last = ib_first
         do while (ib_last < nb)
            if (associated(pm%sec_right(ib_last)%p)) then
               if (pm%sec_right(ib_last)%p%kind == KIND_SLUG) exit
            end if
            ib_last = ib_last + 1
         end do
         if (V_unit(iu) < EPS_SMALL) then
            ib_scan = ib_last + 1
            cycle
         end if
         dp = (target_mass(iu) - M_obs(iu))/(V_unit(iu)*gamma_eos)
         if (abs(dp) > 0.5_rk*p_outlet) dp = sign(0.5_rk*p_outlet, dp)

         do ib = ib_first, ib_last
            pm%gJ(ib) = pm%gJ(ib) + dp
         end do
         ib_scan = ib_last + 1
      end do
      deallocate(M_obs, V_unit, target_mass, left_id, right_id)
   end subroutine pm_add_unit_correction_to_gJ

   subroutine pm_apply_unit_pressure_correction(g, pm, nb, dt, A_pipe, gamma_eos, eos, m_in, track_old)
      type(grid_t),       intent(inout) :: g
      type(pm_state_t),   intent(inout) :: pm
      integer(ik),        intent(in)    :: nb
      real(rk),           intent(in)    :: dt, A_pipe, gamma_eos, m_in
      type(eos_params_t), intent(in)    :: eos
      real(rk),           intent(in)    :: track_old(:)
      integer(ik) :: ib, ib_first, ib_last, ib_scan, iu, n_units
      real(rk) :: dp
      real(rk), allocatable :: M_obs(:), V_unit(:), target_mass(:)
      integer(ik), allocatable :: left_id(:), right_id(:)

      if (nb < 1) return
      call pm_collect_unit_state(pm, nb, A_pipe, n_units, M_obs, V_unit, left_id, right_id)
      if (size(track_old) == n_units) then
         call pm_compute_unit_targets(pm, nb, dt, A_pipe, m_in, track_old, target_mass)
      else
         allocate(target_mass(n_units))
         target_mass(:) = M_obs(:)
      end if

      ib_scan = 1
      iu = 0
      do while (ib_scan <= nb)
         iu = iu + 1
         ib_first = ib_scan
         ib_last = ib_first
         do while (ib_last < nb)
            if (associated(pm%sec_right(ib_last)%p)) then
               if (pm%sec_right(ib_last)%p%kind == KIND_SLUG) exit
            end if
            ib_last = ib_last + 1
         end do

         if (V_unit(iu) >= EPS_SMALL) then
            dp = (target_mass(iu) - M_obs(iu))/(V_unit(iu)*gamma_eos)
            do ib = ib_first, ib_last
               pm%pres(ib) = pm%pres(ib) + dp
               pm%bub_obj(ib)%p%pres = pm%pres(ib)
               pm%bub_obj(ib)%p%rhog = rho_g(eos, pm%bub_obj(ib)%p%pres)
               pm%bub_obj(ib)%p%Ugs  = pm%m_old(ib)/max(pm%bub_obj(ib)%p%rhog, EPS_SMALL)
            end do
         end if
         ib_scan = ib_last + 1
      end do

      deallocate(M_obs, V_unit, target_mass, left_id, right_id)
   end subroutine pm_apply_unit_pressure_correction

   subroutine slug_momentum_coeffs(g, slug, dt, fp, W_eff_cap, A_slug, D_slug0)
      type(grid_t),            intent(in)  :: g
      type(object_t), pointer, intent(in)  :: slug
      real(rk),                intent(in)  :: dt, W_eff_cap
      type(friction_params_t), intent(in)  :: fp
      real(rk),                intent(out) :: A_slug, D_slug0
      type(object_t), pointer :: bub_L, bub_R, sec_feed_L, sec_feed_R
      real(rk) :: U_m_n, Lslug, L_total, L_bL, L_bR, lam_l
      real(rk) :: U_b_back, U_F, W_loc, beta_R, U_l_R, h_L, h_R
      real(rk) :: A_pipe, alpha_R, gamma_geom, Sl, Sg, Si, A_g, Dg, Re_g, lam_i
      real(rk) :: rhog_R, Ug_R, interfacial_front

      ! Paper Eq. 3.18 (thesis page 36, scan 9.png) freezes every U_m on
      ! the RHS at time n.  In particular U_F^n, U_b^n, |U_m|^n and the
      ! ρ_l L/dt · U_m^n  inertia term must NOT recompute from the
      ! running outer-iterate slug%Um (which would make D_slug0 nonlinear
      ! in the Picard variable and was the actual source of the
      ! divergence we observed).  The PM step's entry point now snapshots
      ! slug%Um into slug%Um_n exactly once per timestep; we read it back
      ! here so the LASSI Eq. 3.18 row is genuinely linear in
      ! (p^{n+1}, U_m^{n+1}) and the simultaneous Thomas solve converges
      ! in a single pass (paper Figure 3.3).
      U_m_n = slug%Um_n
      Lslug = max(object_length(slug), EPS_SMALL)
      bub_L => slug%prev
      bub_R => slug%next
      sec_feed_L => null(); sec_feed_R => null()
      if (associated(bub_L)) then
         if (associated(bub_L%prev)) then
            if (bub_L%prev%kind == KIND_SECTION) sec_feed_L => bub_L%prev
         end if
      end if
      if (associated(bub_R)) then
         if (associated(bub_R%next)) then
            if (bub_R%next%kind == KIND_SECTION) sec_feed_R => bub_R%next
         end if
      end if
      L_bL = 0.0_rk; L_bR = 0.0_rk
      if (associated(sec_feed_L)) L_bL = 0.5_rk*max(object_length(sec_feed_L), EPS_SMALL)
      if (associated(sec_feed_R)) L_bR = 0.5_rk*max(object_length(sec_feed_R), EPS_SMALL)
      L_total = Lslug + L_bL + L_bR
      W_loc = wake_effect(Lslug, g%D, W_eff_cap)
      U_b_back = W_loc*bendiksen_nose(U_m_n, g%D, g%phi)
      beta_R = 1.0_rk
      U_l_R = U_m_n
      h_R = 0.0_rk
      rhog_R = g%rho_g_ref
      Ug_R = U_l_R
      interfacial_front = 0.0_rk
      if (associated(sec_feed_R)) then
         beta_R = max(min(sec_feed_R%beta, 1.0_rk - EPS_SMALL), EPS_SMALL)
         U_l_R = sec_feed_R%Ul
         h_R = h_l_of_beta(sec_feed_R%beta, g%D)
         alpha_R = max(1.0_rk - beta_R, EPS_SMALL)
         ! Paper Eq. 3.18 explicit terms ⇒ frozen at timestep n.  Read
         ! the bubble's ρ_g^n / U_g^{S,n} from the _n snapshots taken
         ! by pressure_momentum_step at PM entry, NOT from the live
         ! bubble%rhog / bubble%Ugs (which pm_block_solve overwrites
         ! with the n+1 iterate at the end of each Picard inner step).
         rhog_R = section_rhog_n(sec_feed_R)
         Ug_R = section_ugs_n(sec_feed_R)/alpha_R
         A_pipe = pipe_area(g%D)
         A_g = alpha_R*A_pipe
         gamma_geom = gamma_from_beta(beta_R)
         call perimeters(gamma_geom, g%D, Sl, Sg, Si)
         Dg = 4.0_rk*A_g/max(Sg + Si, EPS_SMALL)
         Re_g = rhog_R*abs(Ug_R - U_l_R)*Dg/max(fp%mu_g, EPS_SMALL)
         lam_i = dw_lambda(Re_g, fp%rough/g%D)
         if (fp%use_AH) lam_i = lam_i*(1.0_rk + 75.0_rk*beta_R)*fp%ai_factor
         interfacial_front = 0.125_rk*Si/max(A_g, EPS_SMALL)*lam_i*rhog_R &
                              *(Ug_R - U_l_R)*abs(Ug_R - U_l_R)*U_l_R
      end if
      U_F = safe_div_ufront(U_m_n - beta_R*U_l_R, 1.0_rk - beta_R, U_m_n)
      h_L = 0.0_rk
      if (associated(sec_feed_L)) h_L = h_l_of_beta(sec_feed_L%beta, g%D)
      lam_l = g%lam_l_slug
      A_slug = g%rho_l*L_total/dt &
               + g%rho_l*(U_b_back - U_m_n) &
               + 0.5_rk*lam_l/g%D * Lslug * g%rho_l * abs(U_m_n)
      ! Paper Eq. 3.18 (thesis page 36, scan 9.png) — ALL U_m on the RHS
      ! are at time n:
      !   ρ_l L/dt · U_m^{n+1} = (p_L − p_R) − g ρ_l (h_R − h_L)
      !                          + ρ_l (U_F^n − U_m^n) · U_{l,R}^n     ← explicit
      !                          − ρ_l (U_b^n − U_m^n) · U_m^{n+1}     ← implicit on Um^{n+1}
      !                          − λ_l/(2D) L_s |U_m|^n · U_m^{n+1}    ← implicit on Um^{n+1}
      !                          + ρ_l L/dt · U_m^n                    ← explicit RHS
      !                          + d_p
      ! With U_m^n frozen at the timestep-entry value, A_slug and
      ! D_slug0 are CONSTANTS w.r.t. the Picard iterate and the row is
      ! a true linear equation in (p_L, p_R, U_m^{n+1}).
      D_slug0 = g%rho_l*L_total/dt * U_m_n &
                + g%rho_l*(U_F - U_m_n)*U_l_R &
                - g%rho_l*G_ACC*(h_R - h_L) &
                - g%rho_l*G_ACC*sin(g%phi)*Lslug &
                + interfacial_front

      ! D.4 break-down for the per-slug PM diagnostic CSV.  All terms in
      ! [Pa] (i.e. directly comparable to A_slug × Um or D_slug0).  Sign
      ! convention follows Eq. 3.18 with positive contributions to the
      ! rhs of  A_slug · Um^{n+1} − pL + pR = D_slug.
      pm_last_diag%Lslug     = Lslug
      pm_last_diag%Um_old    = U_m_n
      pm_last_diag%beta_R    = beta_R
      pm_last_diag%U_l_R     = U_l_R
      pm_last_diag%U_F       = U_F
      pm_last_diag%U_b_back  = U_b_back
      pm_last_diag%rhoL_dt   = g%rho_l*L_total/dt
      pm_last_diag%conv_F    = g%rho_l*(U_F - U_m_n)*U_l_R
      pm_last_diag%conv_B    = g%rho_l*(U_b_back - U_m_n)
      pm_last_diag%fric      = 0.5_rk*lam_l/g%D * Lslug * g%rho_l * abs(U_m_n)
      pm_last_diag%hydrostat = -g%rho_l*G_ACC*(h_R - h_L)
      pm_last_diag%gravity   = -g%rho_l*G_ACC*sin(g%phi)*Lslug
      pm_last_diag%intf      = interfacial_front
   end subroutine slug_momentum_coeffs

   subroutine add_coeff(A, row, col, val)
      real(rk), intent(inout) :: A(:,:)
      integer(ik), intent(in) :: row, col
      real(rk), intent(in) :: val
      if (row >= 1 .and. row <= size(A,1) .and. col >= 1 .and. col <= size(A,2)) then
         A(row,col) = A(row,col) + val
      end if
   end subroutine add_coeff

   integer(ik) function find_slug_idx(slug_obj, n_slug, obj) result(idx)
      type(ptr_wrap_t), intent(in) :: slug_obj(:)
      integer(ik), intent(in) :: n_slug
      type(object_t), pointer, intent(in) :: obj
      integer(ik) :: k
      idx = 0
      if (.not. associated(obj)) return
      do k = 1, n_slug
         if (associated(slug_obj(k)%p, obj)) then
            idx = k
            return
         end if
      end do
   end function find_slug_idx

   integer(ik) function find_bubble_idx(pm, nb, obj) result(idx)
      type(pm_state_t), intent(in) :: pm
      integer(ik), intent(in) :: nb
      type(object_t), pointer, intent(in) :: obj
      integer(ik) :: k
      idx = 0
      if (.not. associated(obj)) return
      do k = 1, nb
         if (associated(pm%bub_obj(k)%p, obj)) then
            idx = k
            return
         end if
      end do
   end function find_bubble_idx

   integer(ik) function block_bandwidth(A, n) result(bw)
      real(rk), intent(in) :: A(:,:)
      integer(ik), intent(in) :: n
      integer(ik) :: i, j
      bw = 0
      do i = 1, n
         do j = 1, n
            if (abs(A(i,j)) > 0.0_rk) bw = max(bw, abs(i-j))
         end do
      end do
   end function block_bandwidth

   subroutine solve_block_thomas(A, rhs, block_size, nblock)
      real(rk),    intent(inout) :: A(:,:), rhs(:)
      integer(ik), intent(in)    :: block_size(:)
      integer(ik), intent(in)    :: nblock
      integer(ik), allocatable :: bstart(:)
      integer(ik) :: iblk, k, r, j, l, nK, nN, sK, sN, nvar
      real(rk), allocatable :: D(:,:), coeff(:), sol(:)
      real(rk) :: accum

      nvar = size(rhs)
      allocate(bstart(nblock))
      bstart(1) = 1
      do iblk = 2, nblock
         bstart(iblk) = bstart(iblk-1) + block_size(iblk-1)
      end do

      do k = 1, nblock - 1
         nK = block_size(k)
         nN = block_size(k+1)
         sK = bstart(k)
         sN = bstart(k+1)
         allocate(D(nK,nK), coeff(nK))
         D(:,:) = A(sK:sK+nK-1, sK:sK+nK-1)
         do r = 0, nN - 1
            coeff(:) = A(sN+r, sK:sK+nK-1)
            if (all(abs(coeff) <= EPS_TINY)) cycle
            call solve_dense_small(transpose(D), coeff, nK)
            do j = sK, nvar
               accum = 0.0_rk
               do l = 1, nK
                  accum = accum + coeff(l)*A(sK+l-1,j)
               end do
               A(sN+r,j) = A(sN+r,j) - accum
            end do
            rhs(sN+r) = rhs(sN+r) - dot_product(coeff, rhs(sK:sK+nK-1))
            A(sN+r, sK:sK+nK-1) = 0.0_rk
         end do
         deallocate(D, coeff)
      end do

      do k = nblock, 1, -1
         nK = block_size(k)
         sK = bstart(k)
         allocate(D(nK,nK), sol(nK))
         D(:,:) = A(sK:sK+nK-1, sK:sK+nK-1)
         sol(:) = rhs(sK:sK+nK-1)
         if (k < nblock) then
            nN = block_size(k+1)
            sN = bstart(k+1)
            do r = 1, nK
               sol(r) = sol(r) - dot_product(A(sK+r-1, sN:sN+nN-1), rhs(sN:sN+nN-1))
            end do
         end if
         call solve_dense_small(D, sol, nK)
         rhs(sK:sK+nK-1) = sol(:)
         deallocate(D, sol)
      end do

      deallocate(bstart)
   end subroutine solve_block_thomas

   subroutine solve_dense_small(Ain, x, n)
      real(rk),    intent(in)    :: Ain(:,:)
      real(rk),    intent(inout) :: x(:)
      integer(ik), intent(in)    :: n
      real(rk), allocatable :: M(:,:), b(:)
      integer(ik) :: i, j, k
      real(rk) :: pivot, factor, accum

      allocate(M(n,n), b(n))
      M(:,:) = Ain(1:n,1:n)
      b(:) = x(1:n)

      do k = 1, n - 1
         pivot = M(k,k)
         if (abs(pivot) < EPS_TINY) then
            if (pivot < 0.0_rk) then
               pivot = -EPS_TINY
            else
               pivot = EPS_TINY
            end if
            M(k,k) = pivot
         end if
         do i = k + 1, n
            if (abs(M(i,k)) <= EPS_TINY) cycle
            factor = M(i,k)/pivot
            M(i,k) = 0.0_rk
            do j = k + 1, n
               M(i,j) = M(i,j) - factor*M(k,j)
            end do
            b(i) = b(i) - factor*b(k)
         end do
      end do

      do i = n, 1, -1
         accum = b(i)
         do j = i + 1, n
            accum = accum - M(i,j)*x(j)
         end do
         pivot = M(i,i)
         if (abs(pivot) < EPS_TINY) then
            if (pivot < 0.0_rk) then
               pivot = -EPS_TINY
            else
               pivot = EPS_TINY
            end if
         end if
         x(i) = accum/pivot
      end do

      deallocate(M, b)
   end subroutine solve_dense_small

   subroutine solve_banded_dense(A, rhs, n, bw)
      real(rk), intent(inout) :: A(:,:), rhs(:)
      integer(ik), intent(in) :: n, bw
      integer(ik) :: k, i, j, jmax
      real(rk) :: pivot, factor
      do k = 1, n-1
         pivot = A(k,k)
         if (abs(pivot) < EPS_TINY) then
            if (pivot < 0.0_rk) then
               A(k,k) = -EPS_TINY
            else
               A(k,k) = EPS_TINY
            end if
            pivot = A(k,k)
         end if
         do i = k+1, min(n, k+bw)
            if (abs(A(i,k)) <= EPS_TINY) cycle
            factor = A(i,k)/pivot
            A(i,k) = 0.0_rk
            jmax = min(n, k+bw)
            do j = k+1, jmax
               A(i,j) = A(i,j) - factor*A(k,j)
            end do
            rhs(i) = rhs(i) - factor*rhs(k)
         end do
      end do
      do i = n, 1, -1
         pivot = A(i,i)
         if (abs(pivot) < EPS_TINY) then
            if (pivot < 0.0_rk) then
               pivot = -EPS_TINY
            else
               pivot = EPS_TINY
            end if
         end if
         jmax = min(n, i+bw)
         do j = i+1, jmax
            rhs(i) = rhs(i) - A(i,j)*rhs(j)
         end do
         rhs(i) = rhs(i)/pivot
      end do
   end subroutine solve_banded_dense

!======================================================================
! Step 6: build explicit predictor m_J^* by FREEZING m at neighbours.
!
!   Treating b_J m_{J+1} and c_J m_{J-1} explicitly at time n on the
!   LHS together with d_J from Eq. 3.17 (which itself contains b_J m^n_{J+1}
!   + c_J m^n_{J-1}), the discrete momentum equation simplifies to:
!
!     a_J m_J^* = d_J + (b_J m_{J+1}^n + c_J m_{J-1}^n)
!
!   so m_J^* is the predictor BEFORE the implicit pressure-gradient
!   coupling is applied; the pressure-gradient term is then attached
!   in the next step via
!
!     m_J^{n+1} = m_J^* - kappa_J ( p_{J+1}^{n+1} - p_J^{n+1} )
!
!   with kappa_J = (1/a_J) * alpha_J / L_J  (note: kappa_J is now a
!   quotient because a_J carries 1/dt + advection + friction).
!======================================================================
   subroutine pm_compute_m_star(pm, nb, dt, m_in)
      type(pm_state_t), intent(inout) :: pm
      integer(ik),      intent(in)    :: nb
      real(rk),         intent(in)    :: dt, m_in
      integer(ik) :: ib
      real(rk) :: m_neigh_R, m_neigh_L
      real(rk), allocatable :: m_star(:)
      allocate(m_star(nb))
      do ib = 1, nb
         if (ib < nb) then
            m_neigh_R = pm%m_old(ib+1)
         else
            m_neigh_R = pm%m_old(ib)
         end if
         if (ib > 1) then
            m_neigh_L = pm%m_old(ib-1)
         else
            m_neigh_L = m_in
         end if
         ! m_J^* solves a_J m_J^* = d_J + b_J m_{J+1}^n + c_J m_{J-1}^n
         ! (the b_J m_{J+1}^n + c_J m_{J-1}^n in d_J doubles the explicit
         !  inflow to mimic Crank-Nicolson advection: a_J m^{n+1} =
         !  b_J m^{n+1}_{J+1} + c_J m^{n+1}_{J-1} + d_J  with both n and n+1
         !  contributions, freeze the n+1 contributions at the n predictor.)
         m_star(ib) = (pm%dJ(ib) + pm%bJ(ib)*m_neigh_R + pm%cJ(ib)*m_neigh_L) &
                       / max(pm%aJ(ib), EPS_SMALL)
      end do
      ! reuse pm%m_old to carry m_star into recovery step
      do ib = 1, nb
         pm%m_old(ib) = m_star(ib)
      end do
      deallocate(m_star)
   end subroutine pm_compute_m_star

!======================================================================
! Direct thesis-style gas-momentum solve after substituting Eq. 3.14:
!     p_J = z_J (m_{J-1} - m_J) + g_J
! into
!     a_J m_J - b_J m_{J+1} - c_J m_{J-1}
!       + (alpha_J/L_J)(p_{J+1} - p_J) = d_J.
! This yields a tridiagonal system in m only.
!======================================================================
   subroutine pm_assemble_mflux_tridiag(pm, nb, p_outlet, m_in)
      type(pm_state_t), intent(inout) :: pm
      integer(ik),      intent(in)    :: nb
      real(rk),         intent(in)    :: p_outlet, m_in
      integer(ik) :: ib
      real(rk)    :: kJ, rhs

      do ib = 1, nb
         kJ = pm%alpha_b(ib)/max(pm%L_b(ib), EPS_SMALL)
         if (ib < nb) then
            pm%tri_a(ib) = -(pm%cJ(ib) + kJ*pm%zJ(ib))
            pm%tri_b(ib) =   pm%aJ(ib) + kJ*(pm%zJ(ib) + pm%zJ(ib+1))
            pm%tri_c(ib) = -(pm%bJ(ib) + kJ*pm%zJ(ib+1))
            rhs = pm%dJ(ib) - kJ*(pm%gJ(ib+1) - pm%gJ(ib))
         else
            pm%tri_a(ib) = -(pm%cJ(ib) + kJ*pm%zJ(ib))
            pm%tri_b(ib) =   pm%aJ(ib) - pm%bJ(ib) + kJ*pm%zJ(ib)
            pm%tri_c(ib) = 0.0_rk
            rhs = pm%dJ(ib) - kJ*(p_outlet - pm%gJ(ib))
         end if
         if (ib == 1) then
            rhs = rhs - pm%tri_a(ib)*m_in
            pm%tri_a(ib) = 0.0_rk
         end if
         pm%tri_d(ib) = rhs
      end do
   end subroutine pm_assemble_mflux_tridiag

   subroutine pm_solve_tridiag_mflux(pm, nb)
      type(pm_state_t), intent(inout) :: pm
      integer(ik),      intent(in)    :: nb
      real(rk), allocatable :: a_T(:), b_T(:), c_T(:), d_T(:)
      integer(ik) :: ib
      allocate(a_T(nb), b_T(nb), c_T(nb), d_T(nb))
      a_T(:) = pm%tri_a(:)
      b_T(:) = pm%tri_b(:)
      c_T(:) = pm%tri_c(:)
      d_T(:) = pm%tri_d(:)
      a_T(1) = 0.0_rk
      c_T(nb) = 0.0_rk
      call thomas_solve(a_T, b_T, c_T, d_T, nb)
      do ib = 1, nb
         pm%m_old(ib) = d_T(ib)
         if (pm%m_old(ib) /= pm%m_old(ib)) pm%m_old(ib) = 0.0_rk
         if (abs(pm%m_old(ib)) > 1.0e6_rk) pm%m_old(ib) = sign(1.0e6_rk, pm%m_old(ib))
      end do
      deallocate(a_T, b_T, c_T, d_T)
   end subroutine pm_solve_tridiag_mflux

   subroutine pm_recover_pressure_from_mflux(pm, nb, m_in, p_outlet, eos)
      type(pm_state_t),   intent(inout) :: pm
      integer(ik),        intent(in)    :: nb
      real(rk),           intent(in)    :: m_in, p_outlet
      type(eos_params_t), intent(in)    :: eos
      integer(ik) :: ib
      real(rk) :: m_left
      do ib = 1, nb
         if (ib == 1) then
            m_left = m_in
         else
            m_left = pm%m_old(ib-1)
         end if
         pm%pres(ib) = bounded_pressure(pm%zJ(ib)*(m_left - pm%m_old(ib)) + pm%gJ(ib), p_outlet)
      end do
      pm%pres(nb) = p_outlet
      do ib = 1, nb
         pm%bub_obj(ib)%p%pres = pm%pres(ib)
         pm%bub_obj(ib)%p%rhog = rho_g(eos, pm%bub_obj(ib)%p%pres)
         pm%bub_obj(ib)%p%Ugs  = pm%m_old(ib)/max(pm%bub_obj(ib)%p%rhog, EPS_SMALL)
      end do
   end subroutine pm_recover_pressure_from_mflux

!======================================================================
! Step 7: tridiagonal in p^{n+1}
!
!   Substituting m_J^{n+1} = m_J^* - kappa_J ( p_{J+1} - p_J ) into
!   Eq. 3.14:
!
!       p_J^{n+1} = z_J [ (m_{J-1}^* - kappa_{J-1} (p_J - p_{J-1}))
!                       - (m_J^*     - kappa_J     (p_{J+1} - p_J)) ] + g_J
!
!   yields the symmetric tridiagonal
!
!     -z_J kappa_{J-1} p_{J-1} + (1 + z_J(kappa_{J-1} + kappa_J)) p_J
!         - z_J kappa_J     p_{J+1}  =  z_J (m_{J-1}^* - m_J^*) + g_J
!
!   Notation:
!     z_J     = dt * A_pipe / (V_g_J * d(rho_g)/dp)        Eq. 3.14
!     kappa_J = alpha_J / (L_J * a_J)                     (gas-mom inversion)
!======================================================================
   subroutine pm_assemble_tridiag(pm, nb, p_outlet, m_in, dt)
      type(pm_state_t), intent(inout) :: pm
      integer(ik),      intent(in)    :: nb
      real(rk),         intent(in)    :: p_outlet, m_in, dt
      integer(ik) :: ib
      real(rk)    :: kJ_im1, kJ
      real(rk)    :: m_star_left
      do ib = 1, nb
         kJ = pm%alpha_b(ib)/(pm%L_b(ib) * max(pm%aJ(ib), EPS_SMALL))
         ! kappa_{J-1} for the left coupling
         if (ib > 1) then
            kJ_im1 = pm%alpha_b(ib-1)/(pm%L_b(ib-1) * max(pm%aJ(ib-1), EPS_SMALL))
         else
            kJ_im1 = 0.0_rk          ! inlet treats m_0 as Neumann
         end if
         ! left m_star for the RHS (z_J*(m_left - m_J))
         if (ib > 1) then
            m_star_left = pm%m_old(ib-1)
         else
            m_star_left = m_in
         end if
         pm%tri_a(ib) = -pm%zJ(ib)*kJ_im1
         pm%tri_b(ib) =  1.0_rk + pm%zJ(ib)*(kJ_im1 + kJ)
         pm%tri_c(ib) = -pm%zJ(ib)*kJ
         pm%tri_d(ib) =  pm%zJ(ib)*(m_star_left - pm%m_old(ib)) + pm%gJ(ib)
      end do
      ! Inlet (J=1): subdiagonal already 0 because kJ_im1=0 above.
      pm%tri_a(1) = 0.0_rk
      ! Outlet (J=N): p_{N+1} = p_outlet (Dirichlet) -> move to RHS.
      pm%tri_d(nb) = pm%tri_d(nb) - pm%tri_c(nb)*p_outlet
      pm%tri_c(nb) = 0.0_rk
   end subroutine pm_assemble_tridiag

!======================================================================
! Step 8: Thomas solve in pressure
!======================================================================
   subroutine pm_solve_tridiag_pressure(pm, nb)
      type(pm_state_t), intent(inout) :: pm
      integer(ik),      intent(in)    :: nb
      real(rk), allocatable :: a_T(:), b_T(:), c_T(:), d_T(:)
      integer(ik) :: ib
      allocate(a_T(nb), b_T(nb), c_T(nb), d_T(nb))
      do ib = 1, nb
         a_T(ib) = pm%tri_a(ib)
         b_T(ib) = pm%tri_b(ib)
         c_T(ib) = pm%tri_c(ib)
         d_T(ib) = pm%tri_d(ib)
      end do
      a_T(1)  = 0.0_rk
      c_T(nb) = 0.0_rk
      call thomas_solve(a_T, b_T, c_T, d_T, nb)
      do ib = 1, nb
         pm%pres(ib) = d_T(ib)
      end do
      deallocate(a_T, b_T, c_T, d_T)
   end subroutine pm_solve_tridiag_pressure

!======================================================================
! Step 9: recover m^{n+1} from m_star and pressure gradient, write back
!======================================================================
   subroutine pm_recover_mass_flux(pm, nb, m_in, p_outlet, dt, eos)
      type(pm_state_t),   intent(inout) :: pm
      integer(ik),        intent(in)    :: nb
      real(rk),           intent(in)    :: m_in, p_outlet, dt
      type(eos_params_t), intent(in)    :: eos
      integer(ik) :: ib
      real(rk)    :: kJ, p_R
      do ib = 1, nb
         pm%pres(ib) = bounded_pressure(pm%pres(ib), p_outlet)
      end do
      do ib = 1, nb
         kJ = pm%alpha_b(ib)/(pm%L_b(ib) * max(pm%aJ(ib), EPS_SMALL))
         if (ib < nb) then
            p_R = pm%pres(ib+1)
         else
            p_R = p_outlet
         end if
         ! m^{n+1}_J = m_J^* - kappa_J ( p_{J+1} - p_J )
         pm%m_old(ib) = pm%m_old(ib) - kJ*(p_R - pm%pres(ib))
         if (pm%m_old(ib) /= pm%m_old(ib)) pm%m_old(ib) = 0.0_rk
         if (abs(pm%m_old(ib)) > 1.0e6_rk) pm%m_old(ib) = sign(1.0e6_rk, pm%m_old(ib))
      end do
      ! enforce Dirichlet p at outlet bubble
      pm%pres(nb) = p_outlet
      ! write back
      do ib = 1, nb
         pm%pres(ib) = bounded_pressure(pm%pres(ib), p_outlet)
         pm%bub_obj(ib)%p%pres = pm%pres(ib)
         pm%bub_obj(ib)%p%rhog = rho_g(eos, pm%bub_obj(ib)%p%pres)
         pm%bub_obj(ib)%p%Ugs  = pm%m_old(ib)/max(pm%bub_obj(ib)%p%rhog, EPS_SMALL)
      end do
   end subroutine pm_recover_mass_flux

!======================================================================
! Step 9: per-Unit gas-mass correction (Eq. 3.19)
!
! For each Unit (a contiguous span of sections+bubbles between two
! consecutive slugs, or between inlet/outlet at the ends), enforce
! discrete gas-mass balance:
!     M_obs   = sum_{J in unit} rho_g_J * alpha_J * L_J * A
!     M_track = M_obs(t=0) + dt * (m_in - m_out) * A    (cumulative)
!     dp      = (M_obs - M_track) / (V_unit * d(rho_g)/dp)
!   then add dp to every bubble in the Unit and refresh rho_g / Ugs.
!
! For the typical no-slug initial state the whole pipe is ONE Unit
! and the boundary fluxes are m_in (inlet) and m_N (outlet bubble).
! The cumulative tracker is held on g%head's auxiliary integer in a
! poor-man's persistence: we compute the ONE-STEP residual instead of
! a long-run cumulative tracker (i.e. enforce dM/dt = inflow - outflow)
! which is exactly the rebalanced statement of Eq. 3.19 for the
! incremental form.
!======================================================================
   subroutine pm_unit_gas_mass_correction(g, pm, nb, dt, A_pipe, gamma_eos, eos, m_in, p_outlet)
      type(grid_t),       intent(inout) :: g
      type(pm_state_t),   intent(inout) :: pm
      integer(ik),        intent(in)    :: nb
      real(rk),           intent(in)    :: dt, A_pipe, gamma_eos, m_in, p_outlet
      type(eos_params_t), intent(in)    :: eos
      integer(ik) :: ib, ib_first, ib_last, ib_scan
      real(rk)    :: M_obs, M_pred, V_unit, dp
      real(rk)    :: m_unit_in, m_unit_out
      logical     :: left_at_pipe_inlet, right_at_pipe_outlet
      logical     :: left_blocked_by_slug, right_blocked_by_slug

      if (nb < 1) return

      ib_scan = 1
      do while (ib_scan <= nb)
         !----------------------------------------------------------------
         ! Find the boundaries [ib_first .. ib_last] of the next Unit.
         ! A Unit is a maximal run of bubbles whose enclosing sections do
         ! NOT include a slug; slugs split the pipe into independent gas
         ! regions because no gas mass crosses a liquid slug.
         !----------------------------------------------------------------
         ib_first = ib_scan
         left_blocked_by_slug = .false.
         if (associated(pm%sec_left(ib_first)%p)) then
            if (pm%sec_left(ib_first)%p%kind == KIND_SLUG) left_blocked_by_slug = .true.
         end if
         left_at_pipe_inlet = (ib_first == 1) .and. (.not. left_blocked_by_slug)

         ! advance ib_last forward until we hit a slug or the pipe end
         ib_last = ib_first
         do while (ib_last < nb)
            if (associated(pm%sec_right(ib_last)%p)) then
               if (pm%sec_right(ib_last)%p%kind == KIND_SLUG) exit
            end if
            ib_last = ib_last + 1
         end do

         right_blocked_by_slug = .false.
         if (associated(pm%sec_right(ib_last)%p)) then
            if (pm%sec_right(ib_last)%p%kind == KIND_SLUG) right_blocked_by_slug = .true.
         end if
         right_at_pipe_outlet = (ib_last == nb) .and. (.not. right_blocked_by_slug)

         !----------------------------------------------------------------
         ! Volume / observed mass / predicted mass for this Unit.
         !----------------------------------------------------------------
         M_obs  = 0.0_rk
         M_pred = 0.0_rk
         V_unit = 0.0_rk
         do ib = ib_first, ib_last
            M_obs  = M_obs  + pm%bub_obj(ib)%p%rhog * pm%alpha_b(ib) * pm%L_b(ib) * A_pipe
            M_pred = M_pred + pm%rhog_old(ib)        * pm%alpha_b(ib) * pm%L_b(ib) * A_pipe
            V_unit = V_unit + pm%alpha_b(ib) * pm%L_b(ib) * A_pipe
         end do

         if (V_unit < EPS_SMALL) then
            ib_scan = ib_last + 1
            cycle
         end if

         ! Boundary fluxes:
         !   inlet side : m_in (global inlet) if Unit reaches the pipe
         !                inlet, else 0 (slug blocks gas)
         !   outlet side: pm%m_old(ib_last) if Unit reaches the pipe
         !                outlet, else 0 (slug blocks gas)
         if (left_at_pipe_inlet) then
            m_unit_in = m_in
         else
            m_unit_in = 0.0_rk
         end if
         if (right_at_pipe_outlet) then
            m_unit_out = pm%m_old(ib_last)
         else
            m_unit_out = 0.0_rk
         end if

         M_pred = M_pred + dt*(m_unit_in - m_unit_out)*A_pipe

         dp = (M_pred - M_obs) / (V_unit * gamma_eos)
         if (abs(dp) > 0.5_rk*p_outlet) dp = sign(0.5_rk*p_outlet, dp)

         do ib = ib_first, ib_last
            pm%bub_obj(ib)%p%pres = bounded_pressure(pm%bub_obj(ib)%p%pres + dp, p_outlet)
            pm%bub_obj(ib)%p%rhog = rho_g(eos, pm%bub_obj(ib)%p%pres)
            pm%bub_obj(ib)%p%Ugs  = pm%m_old(ib)/max(pm%bub_obj(ib)%p%rhog, EPS_SMALL)
         end do

         ib_scan = ib_last + 1
      end do
   end subroutine pm_unit_gas_mass_correction

!======================================================================
! Helpers for inlet mass flux
!======================================================================
!======================================================================
! Slug momentum update — semi-implicit Eq. 3.18 (thesis §3.5, page 47).
!
! The full block-Thomas form of Eq. 3.18 (with all convective and
! pressure-thrust terms solved simultaneously with the gas pressure
! tridiagonal, paper Fig. 3.3) couples [p, m, U_m] in one banded matrix.
! Here we exploit that the gas pressures p^{n+1} have ALREADY been
! solved by the preceding tridiagonal Thomas pass, so we can decouple
! and update U_m^{n+1} per-slug as a SCALAR semi-implicit equation:
!
!   A_slug · U_m^{n+1}  =  D_slug                       (SPEC §6.3)
!
!   A_slug = ρ_l L_total / δt
!          + ρ_l (U_b^n − U_m^n)               ← back-side convective inflow
!          + (λ_l/(2D)) L_slug ρ_l |U_m^n|     ← wall friction (linearised)
!
!   D_slug = ρ_l L_total / δt · U_m^n
!          − ρ_l (U_F^n − U_m^n) U_l,R^n        ← front-side convective inflow
!          − (p_R^{n+1} − p_L^{n+1})            ← pressure thrust (using just-solved p^{n+1})
!          − g ρ_l (h_R − h_L)                  ← hydrostatic level change across slug (SPEC line 209)
!          − ρ_l g sin(φ) L_total               ← gravitational body force
!
! where  U_b^n  = W_eff·U_bendiksen(U_m^n, D, φ)         (slug back nose, Eq. 2.10–2.14)
!        U_F^n  = (U_m^n − β_R·U_l,R^n)/(1 − β_R)        (steep-front mass cons., §4.5)
!        U_l,R  = sec_feed_R%Ul                          (downstream stratified film velocity)
!        L_total = L_slug + ½(L_feed_L + L_feed_R)       (effective inertial length)
!        h_R, h_L : liquid level in feed sections (= h_l(β)).
!
! This scalar Picard form keeps the Eq. 3.18 terms literal; the remaining
! difference from the thesis is algebraic organization, not added
! pressure-thrust caps or velocity floors.
!======================================================================
   subroutine pm_update_slug_momentum(g, dt, fp, W_eff_cap)
      type(grid_t),            intent(inout) :: g
      real(rk),                intent(in)    :: dt
      type(friction_params_t), intent(in)    :: fp
      real(rk),                intent(in)    :: W_eff_cap
      type(object_t), pointer :: p, bub_L, bub_R, sec_feed_L, sec_feed_R
      real(rk) :: U_m_n, U_m_new
      real(rk) :: Lslug, L_total, L_bL, L_bR
      real(rk) :: pL, pR
      real(rk) :: lam_l
      real(rk) :: U_b_back, U_F, W_loc
      real(rk) :: beta_R, U_l_R, h_L, h_R
      real(rk) :: A_slug, D_slug
      real(rk) :: A_pipe, alpha_R, gamma_geom, Sl, Sg, Si, A_g, Dg, Re_g, lam_i
      real(rk) :: rhog_R, Ug_R, interfacial_front

      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SLUG) then
            ! Paper Eq. 3.18 freezes U_m^n: read from the per-slug snapshot
            ! (set in pressure_momentum_step entry).  See slug_momentum_coeffs.
            U_m_n   = p%Um_n
            Lslug   = max(object_length(p), EPS_SMALL)

            bub_L => p%prev
            bub_R => p%next

            sec_feed_L => null(); sec_feed_R => null()
            if (associated(bub_L)) then
               if (associated(bub_L%prev)) then
                  if (bub_L%prev%kind == KIND_SECTION) sec_feed_L => bub_L%prev
               end if
            end if
            if (associated(bub_R)) then
               if (associated(bub_R%next)) then
                  if (bub_R%next%kind == KIND_SECTION) sec_feed_R => bub_R%next
               end if
            end if

            ! Effective inertial length L_total = L_slug + ½(L_L + L_R)
            L_bL = 0.0_rk; L_bR = 0.0_rk
            if (associated(sec_feed_L)) L_bL = 0.5_rk*max(object_length(sec_feed_L), EPS_SMALL)
            if (associated(sec_feed_R)) L_bR = 0.5_rk*max(object_length(sec_feed_R), EPS_SMALL)
            L_total = Lslug + L_bL + L_bR

            ! Slug back-nose Bendiksen velocity (Eq. 2.10–2.14)
            W_loc    = wake_effect(Lslug, g%D, W_eff_cap)
            U_b_back = W_loc*bendiksen_nose(U_m_n, g%D, g%phi)

            ! Steep-front velocity (mass conservation, §4.5)
            ! and downstream stratified-film state.
            beta_R = 1.0_rk
            U_l_R  = U_m_n
            h_R    = 0.0_rk
            rhog_R = g%rho_g_ref
            Ug_R   = U_l_R
            interfacial_front = 0.0_rk
            if (associated(sec_feed_R)) then
               beta_R = max(min(sec_feed_R%beta, 1.0_rk - EPS_SMALL), EPS_SMALL)
               U_l_R  = sec_feed_R%Ul
               h_R    = h_l_of_beta(sec_feed_R%beta, g%D)
               alpha_R = max(1.0_rk - beta_R, EPS_SMALL)
               rhog_R = section_rhog(sec_feed_R)
               Ug_R = section_ugs(sec_feed_R)/alpha_R
               A_pipe = pipe_area(g%D)
               A_g = alpha_R*A_pipe
               gamma_geom = gamma_from_beta(beta_R)
               call perimeters(gamma_geom, g%D, Sl, Sg, Si)
               Dg = 4.0_rk*A_g/max(Sg + Si, EPS_SMALL)
               Re_g = rhog_R*abs(Ug_R - U_l_R)*Dg/max(fp%mu_g, EPS_SMALL)
               lam_i = dw_lambda(Re_g, fp%rough/g%D)
               if (fp%use_AH) lam_i = lam_i*(1.0_rk + 75.0_rk*beta_R)*fp%ai_factor
               interfacial_front = 0.125_rk*Si/max(A_g, EPS_SMALL)*lam_i*rhog_R &
                                    *(Ug_R - U_l_R)*abs(Ug_R - U_l_R)*U_l_R
            end if
            ! U_F = (U_m^n − β_R U_l,R^n)/(1 − β_R) — steep front absorbing the film
            U_F = safe_div_ufront(U_m_n - beta_R*U_l_R, 1.0_rk - beta_R, U_m_n)

            ! Upstream liquid level
            h_L = 0.0_rk
            if (associated(sec_feed_L)) h_L = h_l_of_beta(sec_feed_L%beta, g%D)

            ! Pressure thrust (using just-solved p^{n+1} at flanking bubbles).
            pL = 0.0_rk; pR = 0.0_rk
            if (associated(bub_L)) pL = bub_L%pres
            if (associated(bub_R)) pR = bub_R%pres
            ! Assemble Eq. 3.18 scalar  A_slug · U_m^{n+1} = D_slug
            ! Paper (Renault & Nydal Part 1) Eq. 3.18, page 70.png:
            !   ρ_l L/dt · U_m^{n+1} = (p_L − p_R) − g ρ_l (h_R − h_L)
            !                          + ρ_l (U_F − U_m^n) · U_{l,R}^n     ← +
            !                          + ρ_l (U_b − U_m^n) · U_m^{n+1}
            !                          − λ_l/(2D) L_s U_m^{n+1} |U_m^n|
            ! The (U_F − U_m) · U_{l,R} front-eating term is POSITIVE
            ! in the paper; a stale minus sign here used to override
            ! the block-solve answer with a wrong sign and triggered
            ! the runaway Um observed during multi-slug formation
            ! clusters.  slug_momentum_coeffs() at line ~1089 already
            ! uses the correct sign — keep this scalar wrapper in lock
            ! step with that block coefficient assembly.
            lam_l = g%lam_l_slug
            A_slug = g%rho_l*L_total/dt &
                   + g%rho_l*(U_b_back - U_m_n) &
                   + 0.5_rk*lam_l/g%D * Lslug * g%rho_l * abs(U_m_n)
            D_slug = g%rho_l*L_total/dt * U_m_n &
                   + g%rho_l*(U_F - U_m_n)*U_l_R &
                   - (pR - pL) &
                   - g%rho_l*G_ACC*(h_R - h_L) &
                   - g%rho_l*G_ACC*sin(g%phi)*Lslug &
                   + interfacial_front

            if (abs(A_slug) < EPS_SMALL) then
               U_m_new = U_m_n
            else
               U_m_new = D_slug / A_slug
            end if

            p%Um = U_m_new
         end if
         p => p%next
      end do
   end subroutine pm_update_slug_momentum

   pure function pm_inlet_mass_flux(g) result(m_in)
      type(grid_t), intent(in) :: g
      real(rk) :: m_in
      m_in = 0.0_rk
      if (associated(g%head)) then
         if (g%head%kind == KIND_INLET) then
            m_in = g%head%rhog * g%head%Ugs
         end if
      end if
   end function pm_inlet_mass_flux

   pure function pm_inlet_mass_flux_helper(g) result(m_in)
      type(grid_t), intent(in) :: g
      real(rk) :: m_in
      m_in = pm_inlet_mass_flux(g)
   end function pm_inlet_mass_flux_helper


!======================================================================
! Update liquid velocity in every Section / Slug using Eq. 3.20 (linearised)
!======================================================================
   subroutine update_liquid_velocity(g, dt, fp)
      type(grid_t),            intent(inout) :: g
      real(rk),                intent(in)    :: dt
      type(friction_params_t), intent(in)    :: fp
      type(object_t), pointer :: p
      real(rk) :: F0, dFdU, U_old, U_new, alpha, U_g, denom, Ul_per
      ! Thesis page 56 prescribes a 3-tier liquid-velocity update
      ! depending on the local holdup β:
      !
      !   β > BETA_LOW_THR    : standard Eq. 3.20 linearised implicit
      !                         Euler on F.
      !   β < BETA_LOW_THR    : F is too sensitive to U_l for the
      !     and β >= BETA_DIFF  linearisation to be stable; instead
      !                         set U_l to the equilibrium value U*_l
      !                         such that F(β, U*_l, U_g) = 0
      !                         (Newton iteration in
      !                         solve_equilibrium_Ul).
      !   β < BETA_DIFF       : "for practical reasons, U_l is simply
      !                         set to zero when the holdup is
      !                         differential" (page 56, last
      !                         paragraph).
      !
      ! Without this 3-tier rule the Eq. 3.20 denominator
      !   1 − (δt/ρ_l)·∂F/∂U_l
      ! becomes very large at low β (∂F/∂U_l blows up via Sl/A and
      ! 1/α factors), yet the truncation-introduced cap to ±20 m/s
      ! still leaves U_l pinned to ±20 m/s in those cells.  That cap
      ! then propagates next step into the void wave and PM solves,
      ! producing the global Ul-blow-up at slug birth observed in the
      ! 2026-05-03 N=1000 run.
      real(rk), parameter :: BETA_LOW_THR = 1.0e-3_rk
      real(rk), parameter :: BETA_DIFF    = 1.0e-6_rk

      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SECTION) then
            alpha = max(1.0_rk - p%beta, EPS_SMALL)
            U_g   = section_ugs(p)/alpha
            U_old = p%Ul
            if (p%beta < BETA_DIFF) then
               ! "differential" holdup — U_l set to zero (page 56)
               p%Ul = 0.0_rk
            else if (p%beta < BETA_LOW_THR) then
               ! Low-β branch: solve F(β, U_l, U_g) = 0 directly
               ! instead of the linearised implicit step.
               p%Ul = solve_equilibrium_Ul(p%beta, U_g, g%D, g%phi, &
                                           g%rho_l, section_rhog(p), &
                                           fp, U_old)
            else
               ! Liquid velocity update inside a section uses the
               ! LASSI reduced liquid-momentum equation (M3), Eq. 3.12 in
               ! the thesis.  In (M3) the explicit gas-pressure-gradient
               ! term has been absorbed analytically into the IKH
               ! coefficient kappa via the substitution Eq. 3.5 (see
               ! SPEC §2.2-2.3); the dynamic kappa-driven wave is
               ! handled by the void-wave Riemann step in lassi_voidwave.
               ! Here we therefore advance U_l only by the body force F
               ! (friction + reduced-gravity) using the linearised
               ! implicit step Eq. 3.20:
               !     U_l^{n+1} = (U_l^n + (δt/ρ_l)(F^n - U_l^n·dF/dU_l))
               !                 / (1 - (δt/ρ_l) dF/dU_l)
               ! with dF/dU_l from the analytic Jacobian dF_dUl.
               F0    = body_force_F(p%beta, U_old, U_g, g%D, g%phi, g%rho_l, section_rhog(p), fp)
               dFdU  = dF_dUl     (p%beta, U_old, U_g, g%D, g%phi, g%rho_l, section_rhog(p), fp)
               denom = 1.0_rk - (dt/g%rho_l)*dFdU
               if (abs(denom) < EPS_SMALL) denom = sign(EPS_SMALL, denom)
               U_new = (U_old + (dt/g%rho_l)*F0 - U_old*(dt/g%rho_l)*dFdU)/denom
               ! cap the velocity to avoid blow-up (defensive — should
               ! be unreachable for β > BETA_LOW_THR, but kept as a
               ! belt-and-braces guard).
               Ul_per = 20.0_rk
               U_new = max(min(U_new, Ul_per), -Ul_per)
               p%Ul = U_new
            end if
         end if
         ! Slug U_m is updated by pm_update_slug_momentum using p^{n+1}.
         p => p%next
      end do
   end subroutine update_liquid_velocity

!======================================================================
! Helpers for the slug-momentum scalar update
!======================================================================

   pure function h_l_of_beta(beta, D) result(h)
      real(rk), intent(in) :: beta, D
      real(rk) :: h, gam
      gam = gamma_from_beta(beta)
      h   = hl_from_gamma(gam, D)
   end function h_l_of_beta

   function section_rhog_n(sec) result(rhog)
      type(object_t), pointer, intent(in) :: sec
      real(rk) :: rhog
      type(object_t), pointer :: bub
      rhog = 0.0_rk
      bub => sec%next
      if (associated(bub)) then
         if (bub%kind == KIND_BUBBLE) rhog = bub%rhog_n
      end if
   end function section_rhog_n

   function section_ugs_n(sec) result(ugs)
      type(object_t), pointer, intent(in) :: sec
      real(rk) :: ugs
      type(object_t), pointer :: bub
      ugs = 0.0_rk
      bub => sec%next
      if (associated(bub)) then
         if (bub%kind == KIND_BUBBLE) ugs = bub%Ugs_n
      end if
   end function section_ugs_n

   function section_rhog(sec) result(rhog)
      type(object_t), pointer, intent(in) :: sec
      real(rk) :: rhog
      type(object_t), pointer :: bub
      rhog = sec%rhog
      bub => sec%next
      if (associated(bub)) then
         if (bub%kind == KIND_BUBBLE) rhog = bub%rhog
      end if
   end function section_rhog

   function section_ugs(sec) result(ugs)
      type(object_t), pointer, intent(in) :: sec
      real(rk) :: ugs
      type(object_t), pointer :: bub
      ugs = sec%Ugs
      bub => sec%next
      if (associated(bub)) then
         if (bub%kind == KIND_BUBBLE) ugs = bub%Ugs
      end if
   end function section_ugs

   pure function safe_div_ufront(num, den, fallback) result(uf)
      real(rk), intent(in) :: num, den, fallback
      real(rk) :: uf
      if (abs(den) < 1.0e-3_rk) then
         uf = fallback
      else
         uf = num/den
      end if
   end function safe_div_ufront

   pure function bounded_pressure(p, p_ref) result(p_ok)
      real(rk), intent(in) :: p, p_ref
      real(rk) :: p_ok
      p_ok = p
      if (p_ok /= p_ok) p_ok = p_ref
      if (p_ok < 0.05_rk*p_ref) p_ok = 0.05_rk*p_ref
      if (p_ok > 20.0_rk*p_ref) p_ok = 20.0_rk*p_ref
   end function bounded_pressure

end module lassi_press_mom
