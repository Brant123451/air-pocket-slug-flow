!======================================================================
! lassi_diag.f90
! Per-Unit gas-mass and total liquid-mass diagnostics for LASSI.
!
! A "Unit" is a pipe portion bounded by two consecutive Slugs (or by
! the inlet/outlet at the ends).  Inside a Unit, the total gas mass is
! exactly conserved in the continuous LASSI model (no slug allows gas
! to cross), so the discrete sum is a sensitive consistency check.
!======================================================================
module lassi_diag
   use lassi_kinds
   use lassi_objects
   use lassi_grid
   use lassi_geom, only: pipe_area
   implicit none
   private
   public :: write_diag_header, write_diag_step
   public :: total_liquid_mass, total_gas_mass, max_unit_id
   public :: gas_mass_per_unit

contains

   ! ------------------------------------------------------------------ !
   function total_liquid_mass(g) result(M_l)
      type(grid_t), intent(in) :: g
      real(rk) :: M_l
      type(object_t), pointer :: p
      real(rk) :: A, L
      A = pipe_area(g%D)
      M_l = 0.0_rk
      p => g%head
      do while (associated(p))
         L = object_length(p)
         if (p%kind == KIND_SECTION) M_l = M_l + g%rho_l*p%beta*L*A
         if (p%kind == KIND_SLUG)    M_l = M_l + g%rho_l*L*A
         p => p%next
      end do
   end function total_liquid_mass

   ! ------------------------------------------------------------------ !
   function total_gas_mass(g) result(M_g)
      type(grid_t), intent(in) :: g
      real(rk) :: M_g
      type(object_t), pointer :: p
      real(rk) :: A, L
      A = pipe_area(g%D)
      M_g = 0.0_rk
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SECTION) then
            L = object_length(p)
            M_g = M_g + gas_density_for_section(p)*(1.0_rk - p%beta)*L*A
         end if
         p => p%next
      end do
   end function total_gas_mass

   ! ------------------------------------------------------------------ !
   function max_unit_id(g) result(maxu)
      type(grid_t), intent(in) :: g
      integer(ik) :: maxu
      type(object_t), pointer :: p
      maxu = 0
      p => g%head
      do while (associated(p))
         if (p%unit_id > maxu) maxu = p%unit_id
         p => p%next
      end do
   end function max_unit_id

   ! ------------------------------------------------------------------ !
   subroutine gas_mass_per_unit(g, M_per_unit, n_units)
      type(grid_t),         intent(in)    :: g
      real(rk), allocatable, intent(out)  :: M_per_unit(:)
      integer(ik),           intent(out)  :: n_units
      type(object_t), pointer :: p
      real(rk) :: A, L
      integer(ik) :: u
      A = pipe_area(g%D)
      n_units = max_unit_id(g)
      if (n_units < 1) n_units = 1
      allocate(M_per_unit(n_units))
      M_per_unit = 0.0_rk
      p => g%head
      do while (associated(p))
         if (p%kind == KIND_SECTION) then
            L = object_length(p)
            u = p%unit_id
            if (u >= 1 .and. u <= n_units) then
               M_per_unit(u) = M_per_unit(u) + gas_density_for_section(p)*(1.0_rk - p%beta)*L*A
            end if
         end if
         p => p%next
      end do
   end subroutine gas_mass_per_unit

   ! ------------------------------------------------------------------ !
   subroutine write_diag_header(lun)
      integer, intent(in) :: lun
      write(lun, '(A)') '# t  total_liquid_mass[kg]  total_gas_mass[kg]  n_units  max_unit_gas_mass[kg]'
   end subroutine write_diag_header

   subroutine write_diag_step(lun, t, g)
      integer,      intent(in) :: lun
      real(rk),     intent(in) :: t
      type(grid_t), intent(in) :: g
      real(rk), allocatable :: M_per_unit(:)
      integer(ik) :: nu
      real(rk) :: M_l, M_g, M_max_u
      M_l = total_liquid_mass(g)
      M_g = total_gas_mass(g)
      call gas_mass_per_unit(g, M_per_unit, nu)
      M_max_u = 0.0_rk
      if (allocated(M_per_unit)) then
         if (size(M_per_unit) > 0) M_max_u = maxval(M_per_unit)
         deallocate(M_per_unit)
      end if
      write(lun, '(F12.5,1X,ES16.8,1X,ES16.8,1X,I6,1X,ES16.8)') t, M_l, M_g, nu, M_max_u
   end subroutine write_diag_step

   function gas_density_for_section(sec) result(rhog)
      type(object_t), pointer, intent(in) :: sec
      real(rk) :: rhog
      type(object_t), pointer :: bub
      rhog = sec%rhog
      bub => sec%next
      if (associated(bub)) then
         if (bub%kind == KIND_BUBBLE) rhog = bub%rhog
      end if
   end function gas_density_for_section

end module lassi_diag
