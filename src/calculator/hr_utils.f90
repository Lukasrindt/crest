module hr_utils
   use iso_fortran_env, only: wp => real64
   use crest_calculator
   use crest_parameters
   use optimize_maths
   use modelhessian_module
   use axis_module
   use strucrd
   use irmsd_module
   use, intrinsic :: ieee_arithmetic
   implicit none
   private

   public initialize_hessian
   public generate_chess_list
   public diagonalize_matrix
   public prj_hess

contains

   subroutine initialize_hessian(calc, type, xyz, nat, at, hess, hguess, pr) !>Matrix is forced to be positive definite
      type(calcdata), intent(inout) :: calc
      integer, intent(in) :: type
      real(wp), intent(in) :: xyz(3, nat)
      integer, intent(in) :: at(nat), nat
      real(wp), intent(inout) :: hess(:)
      real(wp), optional, intent(in) :: hguess
      logical, intent(in) :: pr
      type(calcdata), allocatable :: newcalc
      type(calculation_settings) :: clevel
      type(mhparam) :: mhset
      integer :: k, i, j, idx, io, nat3

      real(wp), allocatable :: hess_full(:, :)

      real(wp), allocatable :: pmode(:, :), grad(:, :)
      real(wp) :: rot(3), dumi
      logical :: linear
      type(coord) :: mol

      nat3 = 3*nat

    !!$omp critical
      !allocate (pmode(nat3,1)) ! dummy allocated
    !!$omp end critical

      !$omp critical
      allocate (newcalc)
      allocate (hess_full(nat3, nat3), source=0.0_wp)
      !$omp end critical

      select case (type)
      case (0) !>Initialize as a scaled identity
         if (present(hguess)) then
            k = 0
            do i = 1, nat3
               do j = 1, i
                  k = k + 1
                  if (i /= j) then
                     hess(k) = 0.0_wp
                  else
                     hess(k) = hguess
                  end if
               end do
            end do
         else
            write (stdout, *) "No hguess provided"
         end if
      case (1)
         !$omp critical
         !write(stdout,*) calc%calcs(1)%chrg
         call clevel%create('gfnff', chrg=calc%calcs(1)%chrg, uhf=calc%calcs(1)%uhf) !> Different levels?? and what happens to solvent??
         call newcalc%add(clevel)
         !$omp end critical
         call numhess1(nat, at, xyz, newcalc, hess_full(:, :), io)
         ! call dsqtoh(nat3, hess_full(:, :), hess(:)) !>Pack Hessian
      case (2)
         !$omp critical
         call clevel%create('gfn0', chrg=calc%calcs(1)%chrg, uhf=calc%calcs(1)%uhf) !> Different levels?? and what happens to solvent??
         call newcalc%add(clevel)
         !$omp end critical
         call numhess1(nat, at, xyz, newcalc, hess_full(:, :), io)
         ! call dsqtoh(nat3, hess_full(:, :), hess(:))
      case (3)
         !$omp critical
         call clevel%create('gfn1', chrg=calc%calcs(1)%chrg, uhf=calc%calcs(1)%uhf) !> Different levels?? and what happens to solvent??
         call newcalc%add(clevel)
         !$omp end critical
         call numhess1(nat, at, xyz, newcalc, hess_full(:, :), io)
         ! call dsqtoh(nat3, hess_full(:, :), hess(:))
      case (4)
         !$omp critical
         call clevel%create('gfn2', chrg=calc%calcs(1)%chrg, uhf=calc%calcs(1)%uhf) !> Different levels?? and what happens to solvent??
         call newcalc%add(clevel)
         !$omp end critical
         call numhess1(nat, at, xyz, newcalc, hess_full(:, :), io)
         ! call dsqtoh(nat3, hess_full(:, :), hess(:))
      case (5)
         !$omp critical
         mhset%model = calc%mh_type
         call modhes(calc, mhset, nat, xyz, at, hess(:), pr)
         !$omp end critical
      end select

      ! call force_psd_eig(hess_full, nat3)
      hess_full = 0.5*(hess_full + transpose(hess_full))
      call dsqtoh(nat3, hess_full(:, :), hess(:))
      ! call force_positive_definiteness(hess, nat3)

   end subroutine initialize_hessian

   subroutine force_psd_svd(hess, nat3)
      implicit none
      integer, intent(in) :: nat3
      real(wp), intent(inout) :: hess(nat3, nat3)

      real(wp), allocatable :: U(:, :), VT(:, :), S(:)
      real(wp), allocatable :: work(:)
      integer :: lwork, info
      integer :: i, j, l

      ! SVD storage
      allocate (U(nat3, nat3), VT(nat3, nat3), S(nat3))

      ! --- workspace query ---
      lwork = -1
      allocate (work(1))
      call dgesvd('A', 'A', nat3, nat3, hess, nat3, S, U, nat3, VT, nat3, work, lwork, info)
      lwork = int(work(1))
      deallocate (work)
      allocate (work(lwork))

      ! --- compute SVD: hess = U * diag(S) * VT ---
      call dgesvd('A', 'A', nat3, nat3, hess, nat3, S, U, nat3, VT, nat3, work, lwork, info)

      if (info /= 0) then
         write (*, *) "SVD failed, info=", info
         stop
      end if

      ! --- reconstruct PSD matrix: H_psd = V * S * V^T ---
      ! Note: VT = V^T → V(i,l) = VT(l,i)

      hess = 0.0_wp

      do l = 1, nat3
         do i = 1, nat3
            do j = 1, nat3
               hess(i, j) = hess(i, j) + S(l)*VT(l, i)*VT(l, j)
            end do
         end do
      end do

      deallocate (U, VT, S, work)

   end subroutine force_psd_svd

   subroutine force_psd_eig(hess, nat3)
      implicit none
      integer, intent(in) :: nat3
      real(wp), intent(inout) :: hess(nat3, nat3)

      real(wp), allocatable :: eigvec(:, :), eigval(:)
      real(wp), allocatable :: work(:)
      integer :: lwork, info
      integer :: i, j, l

      ! --- allocate ---
      allocate (eigvec(nat3, nat3), eigval(nat3))

      ! --- 1. symmetrize ---
      hess = 0.5_wp*(hess + transpose(hess))

      ! copy because dsyev overwrites input
      eigvec = hess

      ! --- workspace query ---
      lwork = -1
      allocate (work(1))
      call dsyev('V', 'U', nat3, eigvec, nat3, eigval, work, lwork, info)
      lwork = int(work(1))
      deallocate (work)
      allocate (work(lwork))

      ! --- 2. diagonalize ---
      call dsyev('V', 'U', nat3, eigvec, nat3, eigval, work, lwork, info)

      if (info /= 0) then
         write (*, *) "dsyev failed, info=", info
         stop
      end if

      ! --- 3. flip negative eigenvalues ---
      ! do i = 1, nat3
      !    if (eigval(i) < 0.0_wp) eigval(i) = -eigval(i)
      ! end do

      do i = 1, nat3
         if (eigval(i) < 0.0_wp) eigval(i) = 0
      end do
      ! --- 4. reconstruct H = Q Λ Q^T ---
      hess = 0.0_wp

      do l = 1, nat3
         do i = 1, nat3
            do j = 1, nat3
               hess(i, j) = hess(i, j) + eigval(l)*eigvec(i, l)*eigvec(j, l)
            end do
         end do
      end do

      deallocate (eigvec, eigval, work)

   end subroutine force_psd_eig

   subroutine force_positive_definiteness(hess, nat3)
      real(wp), intent(inout) :: hess(:)
      integer, intent(in) :: nat3
      real(wp), allocatable :: eigvec(:, :), eigval(:)
      real(wp), allocatable :: work(:)
      integer, allocatable :: iwork(:)
      integer :: lwork, liwork, info, i, j, k, l
      real(wp) :: elow, damp

      allocate (eigvec(nat3, nat3), eigval(nat3))
      lwork = 1 + 6*nat3 + 2*nat3*nat3
      liwork = 8*nat3
      allocate (work(lwork), iwork(liwork))

      call dspevd('V', 'U', nat3, hess(:), eigval, eigvec, nat3, &
                  work, lwork, iwork, liwork, info)

      if (info /= 0) then
         write (*, *) "dspevd failed, info = ", info
         stop
      end if

      elow = minval(eigval)
      damp = max(1.0e-4_wp - elow, 0.0_wp)
      eigval = eigval + damp

      hess(:) = 0.0_wp
      k = 0
      do j = 1, nat3
         do i = 1, j
            k = k + 1
            hess(k) = 0.0_wp
            do l = 1, nat3
               hess(k) = hess(k) + eigval(l)*eigvec(i, l)*eigvec(j, l)
            end do
         end do
      end do

      deallocate (eigvec, eigval, work, iwork)

   end subroutine force_positive_definiteness

   subroutine generate_chess_list(calc, nall, list, nstruc, structures, cos_thresh, S, Y, H_init)
      type(calcdata), intent(in) :: calc
      integer, intent(in) :: nall
      real(wp), intent(in) :: cos_thresh
      integer, intent(inout) :: nstruc
      logical, intent(inout):: list(:)
      real(wp), intent(in) :: H_init(:, :)
      type(coord), intent(in) :: structures(:)
      integer :: i, last, nat3, k, j, case
      real(wp) :: rmsdval
      integer :: n, m

      real(wp), intent(inout), allocatable :: S(:, :), Y(:, :)
      real(wp), allocatable :: R(:, :)
      real(wp) :: norm_s(nall - 1), sy(nall - 1), ci(nall - 1), norm_r(nall - 1)

      real(wp) :: cos_ij, max_cos, max_r
      real(wp) :: c_best, norm_1, norm_2

      real(wp), parameter :: eps = 1d-14

      case = 2

      nat3 = structures(1)%nat*3
      n = nat3
      list(:) = .false.
      ! list(nall) = .true.
      ! list(1:10) = .true.
      ! list(nall) = .true.
      ! list(nall-15:nall) = .true.
      ! write(*,*) list
      !
      !
      !
      !-------------------------------------------------------
      ! STEP 1: build secants
      !-------------------------------------------------------
      if (case == 1) then
         allocate (S(nat3, nall - 1), Y(nat3, nall - 1), R(nat3, nall - 1))
         do i = 1, nall - 1
            S(:, i) = reshape(structures(i + 1)%xyz - structures(i)%xyz, [nat3])
            Y(:, i) = reshape(structures(i + 1)%gradient - structures(i)%gradient, [nat3])
         end do

         m = nall

         !-------------------------------------------------------
         ! STEP 2: basic curvature screening
         !-------------------------------------------------------
         c_best = 0d0

         do k = 1, m - 1

            norm_s(k) = sqrt(dot_product(S(:, k), S(:, k)))
            sy(k) = dot_product(S(:, k), Y(:, k))

            ! if (norm_s(k) < eps) then
            !    list(k) = .false.
            !    cycle
            ! end if

            if (sy(k) <= 1d-10) then
               list(k) = .false.
               cycle
            end if

            ! norm_1 = norm2(S(:, k))
            ! norm_2 = norm2(Y(:, k))
            ! write (*, *) norm_2/norm_1
            ! if (norm_2/norm_1 < 0.0016_wp) then
            !    list(k) = .false.
            !    cycle
            ! end if

            ci(k) = sy(k)/(norm_s(k)**2 + eps)

            if (ci(k) > c_best) c_best = ci(k)

            list(k) = .true.

         end do
         write (*, *) "List after first filter:", list
         !-------------------------------------------------------
         ! STEP 3: greedy geometric filtering
         !-------------------------------------------------------

         do k = 1, m - 1

            if (.not. list(k)) cycle

            max_cos = 0d0

            do i = 1, k - 1
               if (.not. list(i)) cycle

               cos_ij = dot_product(S(:, i), S(:, k))/ &
                        (sqrt(dot_product(S(:, i), S(:, i)))* &
                         sqrt(dot_product(S(:, k), S(:, k))) + eps)

               if (abs(cos_ij) > max_cos) max_cos = abs(cos_ij)
            end do

            ! reject nearly collinear directions
            if (max_cos > cos_thresh) then
               list(k) = .false.
               cycle
            end if

         end do
      else
         allocate (S(nat3, nall - 1), Y(nat3, nall - 1))
         do i = 1, nall - 1
            S(:, i) = reshape(structures(nall)%xyz - structures(i)%xyz, [nat3])
            Y(:, i) = reshape(structures(nall)%gradient - structures(i)%gradient, [nat3])
         end do

         R = Y - matmul(H_init, S)
         m = nall

         list(:) = .true.
         do k = 1, m - 1

            if (.not. list(k)) cycle

            max_cos = 0d0

            do i = 1, k - 1
               if (.not. list(i)) cycle

               cos_ij = dot_product(R(:, i), R(:, k))/ &
                        (sqrt(dot_product(R(:, i), R(:, i)))* &
                         sqrt(dot_product(R(:, k), R(:, k))) + eps)

               if (abs(cos_ij) > max_cos) max_cos = abs(cos_ij)
            end do

            ! reject nearly collinear directions
            if (max_cos > cos_thresh) then
               list(k) = .false.
               cycle
            end if

         end do

         ! max_r = 0

         ! do i = 1, nall - 1
         !    norm_r(i) = norm2(R(:, i))
         !    write(*,*) norm_r(i)
         !    if (norm_r(i) > max_r) max_r = norm_r(i)
         ! end do

         ! do i = 1, nall - 1
         !    if (norm_r(i) > cos_thresh*max_r) list(i) = .true.
         ! end do
         ! do i = 1, nall - 1
         !    norm_1 = norm2(S(:, i))
         !    norm_2 = norm2(Y(:, i))
         !    write (*, *) norm_2/norm_1
         !    if (norm_2/norm_1 > cos_thresh) list(i) = .true.
         ! end do
         ! list(nall) = .true.
         ! nstruc = 0
         ! if (nint(cos_thresh*n) > 0) nstruc = max(nint(cos_thresh*n), 5)
         ! write (*, *) nstruc
         ! nstruc = min(nall, nstruc)
         ! write (*, *) "NSTRUC:", nstruc
         ! if (nstruc > 0) then
         !    list(nall-nstruc+1:nall) = .true.
         !    list(nall) = .true.
         ! else
         !    list(:) = .false.
         !    list(nall) = .true.
         ! end if
      end if

      nstruc = count(list)
      write(*,*) "nstruc", nstruc
      write (*, *) "Accepted Structures List:", list
   end subroutine generate_chess_list

   subroutine prj_hess(nat, nat3, xyz, hess, phess_ut)
!***************************************************************
!* Projection of the translational and rotational DOF out of
!* the numerical Hessian (wrapper)
!***************************************************************
      implicit none

      integer, intent(in) :: nat, nat3
      real(wp), intent(inout) :: hess(nat3, nat3)
      real(wp), intent(in) ::  xyz(3, nat)
      real(wp), intent(in), optional, target :: phess_ut(:)
      !real(wp) ::  hess_ut(nat3*(nat3+1)/2),pmode(nat3,1)
      real(wp), allocatable, target :: hess_ut(:)
      real(wp), allocatable :: pmode(:, :)
      real(wp), pointer :: phess(:)
      integer :: i

      if (present(phess_ut)) then
         phess => phess_ut
      else
         !$omp critical
         allocate (hess_ut(nat3*(nat3 + 1)/2), source=0.0_wp)
         phess => hess_ut
         !$omp end critical
      end if
      !$omp critical
      allocate (pmode(nat3, 1), source=0.0_wp)
      !$omp end critical

      !> Transforms matrix of the upper triangle vector
      call dsqtoh(nat3, hess, phess)

      !> Projection
      call trproj(nat, nat3, xyz, phess, .false., 0, pmode, 1)

      !> Transforms vector of the upper triangle into matrix
      call dhtosq(nat3, hess, phess)
   end subroutine prj_hess

   subroutine diagonalize_matrix(n, A, evals, info)
!*******************************************************************
!* LAPACK wrapper for symmetric matrix diagonalization
!*
!* Solves: A x = λ x
!* A is overwritten by eigenvectors
!*******************************************************************

      implicit none

      integer, intent(in)    :: n
      real(wp), intent(inout) :: A(n, n)
      real(wp), intent(out)  :: evals(n)
      integer, intent(out)   :: info

      real(wp), allocatable :: work(:), workspace
      integer, allocatable   :: iwork(:)

      integer :: lwork, liwork

      external :: dsyevd

      !---------------------------------------------------------------
      ! Query optimal workspace
      !---------------------------------------------------------------
      lwork = -1
      liwork = -1

      allocate (work(1), iwork(1))

      call dsyevd('V', 'U', n, A, n, evals, work, lwork, iwork, liwork, info)

      if (info /= 0) then
         print *, "Workspace query failed, info=", info
         return
      end if

      ! optimal sizes returned in work(1), iwork(1)
      lwork = int(work(1))
      liwork = iwork(1)

      deallocate (work, iwork)

      allocate (work(lwork), iwork(liwork))

      !---------------------------------------------------------------
      ! Actual diagonalization
      !---------------------------------------------------------------
      call dsyevd('V', 'U', n, A, n, evals, work, lwork, iwork, liwork, info)

      deallocate (work, iwork)

   end subroutine diagonalize_matrix

end module hr_utils
