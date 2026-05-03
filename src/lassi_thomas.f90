!======================================================================
! lassi_thomas.f90
! Plain in-place Thomas algorithm for a tridiagonal system
!     a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i,   i = 1..n
! with a_1 = c_n = 0 implicitly (set by caller).
! Returns x in the d() array.
!======================================================================
module lassi_thomas
   use lassi_kinds
   implicit none
   private
   public :: thomas_solve

contains

   subroutine thomas_solve(a, b, c, d, n)
      integer(ik), intent(in)    :: n
      real(rk),    intent(in)    :: a(n), c(n)
      real(rk),    intent(inout) :: b(n)        ! diagonal — modified
      real(rk),    intent(inout) :: d(n)        ! rhs/result
      integer(ik) :: i
      real(rk)    :: w
      ! Forward elimination
      do i = 2, n
         if (abs(b(i-1)) < EPS_TINY) cycle      ! safety against pivot=0
         w = a(i)/b(i-1)
         b(i) = b(i) - w*c(i-1)
         d(i) = d(i) - w*d(i-1)
      end do
      ! Back substitution
      d(n) = d(n)/max(abs(b(n)), EPS_TINY)*sign(1.0_rk, b(n))
      do i = n-1, 1, -1
         d(i) = (d(i) - c(i)*d(i+1))/max(abs(b(i)), EPS_TINY)*sign(1.0_rk, b(i))
      end do
   end subroutine thomas_solve

end module lassi_thomas
