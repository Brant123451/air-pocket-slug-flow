!======================================================================
! lassi_objects.f90
! Doubly-linked-list "object" type for the LASSI grid.
! Five kinds: Inlet, Outlet, Bubble, Section, Slug  (KIND_* in lassi_kinds).
!
! Layout along the pipe (left→right):
!   Inlet  ←→  Bubble ←→ Section ←→ Bubble ←→ Section ←→ Bubble
!         ←→ Slug  ←→ Bubble ←→ Section ←→ ... ←→ Outlet
!
! Convention for "ownership":
!   Each Section owns the Bubble immediately on its right (the right border).
!   A Slug stores only its mixture velocity Um; the Bubbles on its two sides
!   carry pressure / gas data (the Bubble at the slug's left = its tail; the
!   Bubble at the slug's right = its nose, after slug-status update).
!======================================================================
module lassi_objects
   use lassi_kinds
   implicit none
   private
   public :: object_t, ptr_wrap_t, link_left_right, count_objects, count_kind

   type :: object_t
      integer(ik) :: id   = 0
      integer(ik) :: kind = 0          ! KIND_INLET / KIND_OUTLET / KIND_BUBBLE / KIND_SECTION / KIND_SLUG

      ! Geometric position of the right border of this object [m]
      real(rk) :: zR = 0.0_rk
      ! Length of the object [m] (Bubble has L=0 after equilibration)
      real(rk) :: L  = 0.0_rk

      ! State variables (used selectively depending on kind)
      real(rk) :: beta = 0.0_rk        ! liquid holdup (Section: free; Slug: 1; Bubble: copy of section to its right)
      real(rk) :: Ul   = 0.0_rk        ! liquid velocity (Section); for Slug stored in Um
      real(rk) :: Um   = 0.0_rk        ! mixture velocity (Slug)
      real(rk) :: Um_n = 0.0_rk        ! Slug U_m at the start of the current PM step
                                       ! (frozen during the linear solve so all
                                       ! "U^n" terms in thesis Eq. 3.18 — U_F, U_b,
                                       ! |U_m|, RHS U_m^n — stay constant.  This is
                                       ! what makes the LASSI Eq. 3.18 row LINEAR in
                                       ! (p^{n+1}, U_m^{n+1}) so the simultaneous
                                       ! Thomas solve completes in a single pass.)
      real(rk) :: pres = 0.0_rk        ! pressure (Bubble)
      real(rk) :: rhog = 0.0_rk        ! gas density (Bubble)
      real(rk) :: Ugs  = 0.0_rk        ! gas superficial velocity α U_g (Bubble)
      real(rk) :: rhog_n = 0.0_rk      ! Bubble ρ_g at the start of the current PM step
      real(rk) :: Ugs_n  = 0.0_rk      ! Bubble U_g^S at the start of the current PM step
                                       ! (frozen during the linear solve so the
                                       ! interfacial-shear term in slug-momentum
                                       ! Eq. 3.18 stays at timestep n; without this
                                       ! the term D_slug0 becomes nonlinear in the
                                       ! Picard iterate via bubble%rhog/Ugs.)

      ! Auxiliary / runtime
      integer(ik) :: unit_id = 0       ! id of the Unit (between two slugs) used for gas-mass correction
      logical     :: is_nose = .false. ! true if this Bubble is a slug nose (not a front)
      logical     :: marked_for_remove = .false.
      logical     :: cfl_breach = .false.
      ! Previous-step β (Section) for the dα/dt term in pressure-CV g_J
      ! coefficient (thesis Eq. 3.13 + SPEC line 158 "geometric corrections").
      ! Updated at the END of each PM step to the current β.
      real(rk)    :: beta_prev = -1.0_rk    ! −1 sentinel ⇒ first step, no history

      ! Linked list
      type(object_t), pointer :: prev => null()
      type(object_t), pointer :: next => null()
   end type object_t

   ! Wrapper to allow arrays of pointers (Fortran forbids
   ! `type(object_t), pointer :: arr(:)`).
   type :: ptr_wrap_t
      type(object_t), pointer :: p => null()
   end type ptr_wrap_t

contains

   subroutine link_left_right(left, right)
      type(object_t), pointer :: left, right
      if (associated(left))  left%next  => right
      if (associated(right)) right%prev => left
   end subroutine link_left_right

   function count_objects(head) result(n)
      type(object_t), pointer, intent(in) :: head
      integer(ik) :: n
      type(object_t), pointer :: p
      n = 0
      p => head
      do while (associated(p))
         n = n + 1
         p => p%next
      end do
   end function count_objects

   function count_kind(head, kind) result(n)
      type(object_t), pointer, intent(in) :: head
      integer(ik),             intent(in) :: kind
      integer(ik) :: n
      type(object_t), pointer :: p
      n = 0
      p => head
      do while (associated(p))
         if (p%kind == kind) n = n + 1
         p => p%next
      end do
   end function count_kind

end module lassi_objects
