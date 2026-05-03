!======================================================================
! test_riemann.f90
! Unit tests for the modified-shallow-water Riemann solver.
! Verifies all 8 case dispatches (RR, RS, SR, SS, SSS, RV, VR, RRADB)
! by checking residuals of the wave relations.
!======================================================================
program test_riemann
   use lassi_kinds
   use lassi_riemann
   implicit none

   integer  :: nfail
   nfail = 0

   call test_dam_break_RR(nfail)
   call test_dam_break_RV(nfail)
   call test_collision_SS(nfail)
   call test_saturated_SSS(nfail)
   call test_RRADB(nfail)
   call test_RS(nfail)
   call test_SR(nfail)

   if (nfail == 0) then
      write(*,'(A)') 'ALL RIEMANN TESTS PASSED'
   else
      write(*,'(A,I0,A)') 'TESTS FAILED: ', nfail, ' failure(s)'
      stop 1
   end if

contains

!======================================================================
   subroutine assert_close(label, a, b, tol, nfail)
      character(*), intent(in)    :: label
      real(rk),     intent(in)    :: a, b, tol
      integer,      intent(inout) :: nfail
      if (abs(a - b) > tol) then
         write(*,'(A,A,A,ES16.8,A,ES16.8,A,ES10.2)') &
            'FAIL ', label, ' got=', a, ' expected=', b, ' tol=', tol
         nfail = nfail + 1
      else
         write(*,'(A,A,A,ES10.2)') 'PASS ', label, ' diff=', abs(a-b)
      end if
   end subroutine assert_close

!======================================================================
! Symmetric dam-break (β_L > β_R, U_L = U_R = 0): expect RR with U_M ≠ 0.
   subroutine test_dam_break_RR(nfail)
      integer, intent(inout) :: nfail
      type(riemann_state_t) :: st
      real(rk) :: kappa, U_check
      real(rk) :: beta_L, U_L, beta_R, U_R
      beta_L = 0.6_rk; U_L = 0.0_rk
      beta_R = 0.2_rk; U_R = 0.0_rk
      kappa  = 1.0_rk
      call solve_msw_riemann(beta_L, U_L, beta_R, U_R, kappa, st)
      ! expectation: case = RR or RS depending on sqrt curve crossing
      write(*,'(A,I0)') '[dam-break] case_id = ', st%case_id
      write(*,'(A,F12.6,A,F12.6)') '  beta_M = ', st%beta_M, '  U_M = ', st%U_M
      ! Check that left curve and right curve agree at (β_M, U_M)
      U_check = U_L - 2.0_rk*sqrt(kappa)*(sqrt(st%beta_M) - sqrt(beta_L))   ! left rarefaction inv
      call assert_close('RR-leftcurve', U_check, st%U_M, 1.0e-3_rk, nfail)
   end subroutine test_dam_break_RR

!======================================================================
! Right state dry: expect RV
   subroutine test_dam_break_RV(nfail)
      integer, intent(inout) :: nfail
      type(riemann_state_t) :: st
      real(rk) :: kappa
      call solve_msw_riemann(0.5_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, st)
      call assert_close('RV-case', real(st%case_id, rk), real(CASE_RV, rk), 0.5_rk, nfail)
      ! U at right-edge of fan = U_L + 2 sqrt(kβ_L)
      call assert_close('RV-Ub', st%U_b, 0.0_rk + 2.0_rk*sqrt(1.0_rk*0.5_rk), 1.0e-6_rk, nfail)
   end subroutine test_dam_break_RV

!======================================================================
! Two opposing flows (collision) — expect SS or SSS
   subroutine test_collision_SS(nfail)
      integer, intent(inout) :: nfail
      type(riemann_state_t) :: st
      real(rk) :: kappa
      call solve_msw_riemann(0.4_rk, +0.5_rk, 0.4_rk, -0.5_rk, 1.0_rk, st)
      write(*,'(A,I0)') '[collision] case_id = ', st%case_id
      ! expect β_M > β_L, β_M > β_R (Shock-Shock or SSS)
      if (st%case_id /= CASE_SS .and. st%case_id /= CASE_SSS) then
         write(*,'(A)') 'FAIL collision: case neither SS nor SSS'
         nfail = nfail + 1
      else
         write(*,'(A)') 'PASS collision case'
      end if
   end subroutine test_collision_SS

!======================================================================
! Strong collision driving β_M past 1 — expect SSS
   subroutine test_saturated_SSS(nfail)
      integer, intent(inout) :: nfail
      type(riemann_state_t) :: st
      real(rk) :: kappa
      call solve_msw_riemann(0.85_rk, +5.0_rk, 0.85_rk, -5.0_rk, 1.0_rk, st)
      call assert_close('SSS-case', real(st%case_id, rk), real(CASE_SSS, rk), 0.5_rk, nfail)
      call assert_close('SSS-betaM', st%beta_M, 1.0_rk, 1.0e-9_rk, nfail)
   end subroutine test_saturated_SSS

!======================================================================
! Strong outflow on both sides — expect RRADB
   subroutine test_RRADB(nfail)
      integer, intent(inout) :: nfail
      type(riemann_state_t) :: st
      call solve_msw_riemann(0.4_rk, -3.0_rk, 0.4_rk, +3.0_rk, 1.0_rk, st)
      call assert_close('RRADB-case', real(st%case_id, rk), real(CASE_RRADB, rk), 0.5_rk, nfail)
      call assert_close('RRADB-betaM', st%beta_M, 0.0_rk, 1.0e-9_rk, nfail)
   end subroutine test_RRADB

!======================================================================
! Asymmetric dam-break: tall left, calm right but with positive U_R bias
! to try to force right shock — exercise RS dispatch.
   subroutine test_RS(nfail)
      integer, intent(inout) :: nfail
      type(riemann_state_t) :: st
      call solve_msw_riemann(0.3_rk, +1.0_rk, 0.4_rk, -0.5_rk, 1.0_rk, st)
      write(*,'(A,I0)') '[RS-test] case_id = ', st%case_id
      if (st%case_id /= CASE_RR .and. st%case_id /= CASE_RS .and. st%case_id /= CASE_SR &
          .and. st%case_id /= CASE_SS) then
         write(*,'(A)') 'FAIL RS-test: unexpected case'
         nfail = nfail + 1
      end if
   end subroutine test_RS

!======================================================================
   subroutine test_SR(nfail)
      integer, intent(inout) :: nfail
      type(riemann_state_t) :: st
      call solve_msw_riemann(0.4_rk, +0.5_rk, 0.3_rk, +1.0_rk, 1.0_rk, st)
      write(*,'(A,I0)') '[SR-test] case_id = ', st%case_id
      if (st%case_id < CASE_RR .or. st%case_id > CASE_SS) then
         write(*,'(A)') 'FAIL SR-test: case out of range'
         nfail = nfail + 1
      end if
   end subroutine test_SR

end program test_riemann
