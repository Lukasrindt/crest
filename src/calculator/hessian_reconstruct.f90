module hessian_reconstruct
   ! use iso_fortran_env,only:wp => real64
   use hessupdate_module
   use optimize_maths
   use crest_parameters
   ! use hessian_quality
   ! use hr_utils
   implicit none
   private

   public cashed_hessian

   type :: cashed_hessian

      integer :: steps = 10
      real(wp), allocatable :: gradient(:, :, :)
      real(wp), allocatable :: coords(:, :, :)
      real(wp), allocatable :: energy(:)
      real(wp), allocatable :: H(:, :), S(:, :), Y(:, :)
      integer, allocatable :: order(:), natm
      integer :: stepcount = 0
      real(wp) :: hguess = 0.02_wp
      real(wp), allocatable ::hess(:)
      logical :: track_step = .true.
      integer :: initialize_type != 0
      integer :: hu_type != 0
      ! type(step_monitor) :: mon

   contains

      procedure :: alloc => cashed_hessian_allocate
      procedure :: dealloc => cashed_hessian_deallocate
      procedure :: update => update_cashed_hessian
      procedure :: construct_hessian
      procedure :: mode_quality_analysis

   end type cashed_hessian

contains

   subroutine cashed_hessian_allocate(self, N, steps, initialize_type, hu_type, hguess) !> maybe make keywords optional later
      integer, intent(in) :: N, steps, initialize_type, hu_type
      class(cashed_hessian), intent(inout) :: self
      real(wp), intent(in), optional :: hguess

      self%steps = steps
      if (present(hguess)) self%hguess = hguess
      self%natm = N
      self%initialize_type = initialize_type
      self%hu_type = hu_type
      allocate (self%gradient(steps, 3, N))
      allocate (self%coords(steps, 3, N))
      allocate (self%energy(steps))
      allocate (self%order(steps))
      allocate (self%hess((3*N*(3*N + 1))/2))
      allocate (self%H(3*N, 3*N))

      self%order(:) = 0

   end subroutine cashed_hessian_allocate

   subroutine cashed_hessian_deallocate(self)
      class(cashed_hessian), intent(inout) :: self

      if (allocated(self%gradient)) deallocate (self%gradient)
      if (allocated(self%coords)) deallocate (self%coords)
      if (allocated(self%energy)) deallocate (self%energy)
      if (allocated(self%order)) deallocate (self%order)

   end subroutine cashed_hessian_deallocate

   subroutine update_cashed_hessian(self, gradient, coords, energy)
      class(cashed_hessian), intent(inout) :: self
      real(wp), intent(in) :: gradient(:, :), coords(:, :)
      real(wp), intent(in), optional :: energy
      integer :: idx, i

      self%stepcount = self%stepcount + 1
      idx = minloc(self%order, 1)
      self%order(idx) = self%stepcount
      self%gradient(idx, :, :) = gradient
      if (present(energy)) self%energy(idx) = energy
      self%coords(idx, :, :) = coords
   end subroutine update_cashed_hessian

   subroutine construct_hessian(self)
      class(cashed_hessian), intent(inout) :: self
      integer :: i, j, k, nat3
      real(wp), allocatable :: tmp(:), tmp_coords(:, :), tmp_grads(:, :), dx(:)
      real(wp), allocatable :: S(:, :), Y(:, :)
      real(wp) :: gnorm
      integer :: unit, iter, made_iters

      nat3 = 3*self%natm
      if (self%hu_type < 5) then !Single secant updates!
         allocate (tmp_coords(self%steps, nat3))
         allocate (tmp_grads(self%steps, nat3))
         allocate (tmp(self%steps))
         allocate (dx(nat3))

         tmp = self%order

         tmp_coords = reshape(self%coords, [self%steps, nat3])
         tmp_grads = reshape(self%gradient, [self%steps, nat3])

         made_iters = self%steps

         !>Hessian guess is installed previously in optimize routine but could also be read in explicitly for better readability?

         if (minval(tmp) == 0) then !> Implement keyword like exact HU that kills the process
            made_iters = maxval(tmp) !> if made_iters<steps
            write (stdout, *) "Requsted Number of reconstruction steps is", self%steps, &
            & "but only", made_iters, "geometry optimization steps were made!"
            write (stdout, *) "Hessian is reconstructed with", made_iters, "update steps only!"

            do while (minval(tmp) == 0)
               j = minloc(tmp, 1)
               tmp(j) = HUGE(tmp(j))
            end do
         end if

         do i = 1, made_iters
            if (i == 1) then
               j = minloc(tmp, 1)
               tmp(j) = HUGE(tmp(j))
            else
               j = minloc(tmp, 1) !> This only happens if made_iters>steps
               if (j == 1) then  !> => Not affected if too many steps requested
                  dx = tmp_coords(j, :) - tmp_coords(self%steps, :)
                  call update_hessian(nat3, gnorm, tmp_grads(j, :), tmp_grads(self%steps, :), dx, self%hess(:), self%hu_type)
               else
                  dx = tmp_coords(j, :) - tmp_coords(j - 1, :)
                  call update_hessian(nat3, gnorm, tmp_grads(j, :), tmp_grads(j - 1, :), dx, self%hess(:), self%hu_type)
               end if
               tmp(j) = HUGE(tmp(j))
            end if
         end do

         call dhtosq(nat3, self%H(:, :), self%hess(:))
      else !Multisecnt updates
         call dhtosq(nat3, self%H(:, :), self%hess(:))

         allocate (S(nat3, self%steps - 1)) !CAUTION, steps needs to be initialized properly.
         allocate (self%S(nat3, self%steps - 1), self%Y(nat3, self%steps - 1))
         allocate (Y(nat3, self%steps - 1)) ! SHould be refactored to come from length of coord vector?

         do i = 1, self%steps-1
            S(:, i) = reshape(self%coords(self%steps, :, :) - self%coords(i-1, :, :), [nat3])
            Y(:, i) = reshape(self%gradient(self%steps, :, :) - self%gradient(i-1, :, :), [nat3])
         end do
         self%S = S
         self%Y = Y
         ! write (*, *) self%Y
         if (self%hu_type == 5) call ms_bfgs_update(nat3, self%steps - 1, self%H, S, Y) !For this, init need to be placed in self%H
         if (self%hu_type == 6) call ms_psb_update(nat3, self%steps - 1, self%H, S, Y)
         if (self%hu_type == 7) call ms_bfgs_polar_regularized(nat3, self%steps - 1, self%H, S, Y)
         if (self%hu_type == 8) call ms_rsr_inverse_update(nat3, self%steps - 1, self%H, S, Y)
         if (self%hu_type == 9) call ms_rsr_polar_update(nat3, self%steps - 1, self%H, S, Y)

      end if

   end subroutine construct_hessian

   subroutine update_hessian(nat3, gnorm, grd1, gold, dx, hess, hu_type)
      !==============================================
      !Wrapper for hessian update scheme selection
      !==============================================
      !class(cashed_hessian),intent(inout) :: self
      integer, intent(in) :: nat3
      real(wp), intent(in) :: dx(:), grd1(:), gold(:)
      real(wp), intent(in) :: gnorm
      real(wp), intent(inout) :: hess(:)
      integer, intent(in) :: hu_type

      select case (hu_type)
      case (0)
         call bfgs(nat3, gnorm, grd1, gold, dx, hess)
      case (1)
         call powell(nat3, gnorm, grd1, gold, dx, hess)
      case (2)
         call sr1(nat3, gnorm, grd1, gold, dx, hess)
      case (3)
         call bofill(nat3, gnorm, grd1, gold, dx, hess)
      case (4)
         call schlegel(nat3, gnorm, grd1, gold, dx, hess)
      case default
         write (*, *) 'invalid update selection for hessian reconstruction'
         stop
      end select

   end subroutine update_hessian

   subroutine mode_quality_analysis(self, type, init_hess, modes, n_modes, quality)
      !*******************************************************************
      !* Assess how well each Hessian eigenvector was sampled during reconstruction.
      !*
      !* Updated: uses QR decomposition to obtain an orthonormal basis Q
      !* so that projection is mathematically correct:
      !*   quality(i) = || Q Q^T v_i ||^2 = sum_j (q_j · v_i)^2
      !*
      !*******************************************************************
      implicit none

      real(wp), intent(in)  :: modes(:, :)   ! (n3, n_modes)
      real(wp), intent(inout) :: init_hess(:, :)
      integer, intent(in)  :: n_modes
      class(cashed_hessian) :: self
      real(wp), intent(out) :: quality(:)! (n_modes)
      integer, intent(in) :: type

      real(wp), allocatable :: Q(:, :), lambda(:), C(:, :), P(:, :)
      real(wp), allocatable :: tau(:), work(:), v_approx(:, :)
      real(wp) :: proj, temp, overlap
      integer  :: i, j, nat3, m, n, lwork, info, n_upmodes
      external :: dgeqrf
      external :: dorgqr
      real(wp), external :: ddot

      quality = 0.0_wp
      nat3 = self%natm*3

      select case (type)
      case (1)

         do i = 1, nat3
            do j = 1, nat3
               temp = ddot(nat3, init_hess(:, j), 1, modes(:, i), 1)
               if (temp > quality(i)) quality(i) = temp
            end do
         end do

         quality = abs(1 - quality)
         do i = 1, nat3
            write (*, *) "quality of the", i, "th mode", quality(i)
         end do
      case (2)
         quality = 0.0_wp
         allocate (lambda(self%steps - 1), C(self%steps - 1, self%steps - 1), v_approx(nat3, self%steps - 1))
         allocate (P(nat3, self%steps - 1), Q(nat3, self%steps - 1))
         P = self%Y
         Q = self%S
         call generalized_secant_eigen(Q, P, nat3, self%steps - 1, lambda, C)
         do i = 1, self%steps - 1
            v_approx(:, i) = matmul(self%Y, C(:, i))
            if (norm2(v_approx(:, i)) > 1e-12_wp) then
               v_approx(:, i) = v_approx(:, i)/norm2(v_approx(:, i))
            else
               v_approx(:, i) = 0.0_wp
            end if
         end do
         do j = 1, nat3
            do i = 1, self%steps - 1
               overlap = abs(ddot(nat3, v_approx(:, i), 1, modes(:, j), 1))
               quality(j) = quality(j) + overlap**2
            end do
         end do
         quality = quality/(maxval(quality) + 1.0e-12_wp)
         do i = 1, nat3
            write (*, *) "quality of the", i, "th mode", quality(i)
         end do
      end select

      !THIS IS A PROJECTION INTO THE S VECTORSPACE WITH PREVIOUS ORTHOGONALIZATION
      ! m = nat3
      ! n = self%steps - 1   ! number of step vectors
      !
      ! ! Copy S into Q (will be overwritten by QR)
      ! allocate(Q(m, n))
      ! Q = self%S(:, 1:n)
      !
      ! allocate(tau(min(m,n)))
      !
      ! ! Workspace query
      ! allocate(work(1))
      ! call dgeqrf(m, n, Q, m, tau, work, -1, info)
      ! lwork = int(work(1))
      ! deallocate(work)
      ! allocate(work(lwork))
      !
      ! ! Compute QR factorization
      ! call dgeqrf(m, n, Q, m, tau, work, lwork, info)
      !
      ! ! Generate explicit Q
      ! deallocate(work)
      ! allocate(work(1))
      ! call dorgqr(m, n, min(m,n), Q, m, tau, work, -1, info)
      ! lwork = int(work(1))
      ! deallocate(work)
      ! allocate(work(lwork))
      !
      ! call dorgqr(m, n, min(m,n), Q, m, tau, work, lwork, info)
      !
      ! ! ---- Compute quality using orthonormal Q ----
      ! do i = 1, n_modes
      !    quality(i) = 0.0_wp
      !    do j = 1, n
      !       proj = dot_product(Q(:, j), modes(:, i))
      !       quality(i) = quality(i) + proj**2
      !    end do
      ! end do
      ! n_upmodes = 0
      ! do i=1,nat3
      ! write(*,*) quality(i)
      ! if (quality(i)< 0.002_wp)n_upmodes = n_upmodes +1
      ! enddo
      ! write(*,*) n_upmodes
      !*******************************************************************
      ! Original implementation (incorrect if S not orthonormal)
      !*******************************************************************
      !
      ! real(wp), allocatable :: proj_vec(:)
      ! real(wp) :: proj
      ! integer  :: i, j, nat3
      !
      ! nat3 = self%natm*3
      ! allocate (proj_vec(nat3))
      !
      ! do i = 1, n_modes
      !    proj_vec = 0.0_wp
      !    do j = 1, self%steps-1
      !       proj = dot_product(self%S(:, j), modes(:, i))
      !       proj_vec = proj_vec + proj*self%S(:,j)
      !    end do
      !
      !    quality(i) = sum(proj_vec**2)
      !    quality(i) = max(0.0_wp, min(1.0_wp, quality(i)))
      ! end do
      !
      ! deallocate (proj_vec)
      !
      !*******************************************************************

      ! deallocate(Q, tau, work)

   end subroutine mode_quality_analysis

   subroutine generalized_secant_eigen(S, Y, n, msec, lambda, C)
      implicit none
      integer, intent(in) :: n, msec
      real(wp), intent(in) :: S(n, msec), Y(n, msec)
      real(wp), intent(out):: lambda(msec)
      real(wp), intent(out):: C(msec, msec)

      real(wp), allocatable :: A(:, :), B(:, :)
      real(wp), allocatable :: alphar(:), alphai(:), beta(:)
      real(wp), allocatable :: vl(:, :), vr(:, :), work(:)
      integer :: lwork, info, i

      !--------------------------------------------------------
      ! Build reduced secant matrices
      ! A = Y^T Y
      ! B = S^T Y
      !--------------------------------------------------------

      allocate (A(msec, msec))
      allocate (B(msec, msec))

      A = matmul(transpose(Y), Y)
      B = matmul(transpose(S), Y)

      allocate (alphar(msec), alphai(msec), beta(msec))
      allocate (vl(msec, msec), vr(msec, msec))

      !--------------------------------------------------------
      ! Workspace query
      !--------------------------------------------------------

      lwork = -1
      allocate (work(1))

      call dggev('N', 'V', msec, A, msec, B, msec, &
                 alphar, alphai, beta, &
                 vl, msec, vr, msec, work, lwork, info)

      lwork = int(work(1))
      ! if (info /= 0) then
      !    print *, "DGGEV workspace query failed, info=", info
      !    stop
      ! end if
      deallocate (work)
      allocate (work(lwork))
      !--------------------------------------------------------
      ! Solve generalized eigenproblem
      !--------------------------------------------------------

      call dggev('N', 'V', msec, A, msec, B, msec, &
                 alphar, alphai, beta, &
                 vl, msec, vr, msec, work, lwork, info)

      !--------------------------------------------------------
      ! Convert eigenvalues
      !--------------------------------------------------------

      do i = 1, msec
         if (abs(beta(i)) > 1.0e-12_wp) then
            lambda(i) = alphar(i)/beta(i)
         else
            lambda(i) = 0.0_wp
         end if
      end do

      C = vr

      deallocate (A, B, alphar, alphai, beta, vl, vr, work)

   end subroutine

end module hessian_reconstruct
