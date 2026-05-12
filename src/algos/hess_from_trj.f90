!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2026 Lukas Rindt
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

subroutine hess_from_trj(env, tim)
   use crest_parameters
   use crest_data
   use crest_calculator
   use strucrd
   use hessian_reconstruct
   use hr_utils
   use thermochem_module
   use hessian_quality
   use optimize_maths
   use, intrinsic :: ieee_arithmetic
   implicit none
   type(systemdata), intent(inout) :: env
   type(timer), intent(inout)      :: tim
   type(coord) :: mol, molnew
   integer :: i, j, k, l, io, ich, idx, info
   logical :: pr, wr
!========================================================================================!
   ! type(calcdata) :: calc
   real(wp) :: accuracy, etemp

   integer :: V, maxgen
   integer, allocatable :: A(:, :)
   logical, allocatable :: rings(:, :)
   integer, allocatable :: tmp(:)
   logical :: connected, fail, doreturn

   real(wp) :: energy
   real(wp), allocatable :: grad(:, :), geo(:, :), csv(:, :), q(:)
   type(coord), allocatable :: structures(:)
   type(coord) :: init_mol
   integer :: nall, steps, nstruc, nat3, n_repaired
   real(wp) :: etot
   logical, allocatable :: list(:)
   real(wp), allocatable :: hess(:, :), freqs(:), quality(:), final_hess(:, :), init_hess(:, :), init_eigenvalues(:)
   real(wp), allocatable :: S(:,:), Y(:,:)

!========================================================================================!
   call tim%start(14, 'Test implementation')
!========================================================================================!
   !call system('figlet welcome')
   write (*, *) "              _                          "
   write (*, *) "__      _____| | ___ ___  _ __ ___   ___ "
   write (*, *) "\ \ /\ / / _ \ |/ __/ _ \| '_ ` _ \ / _ \"
   write (*, *) " \ V  V /  __/ | (_| (_) | | | | | |  __/"
   write (*, *) "  \_/\_/ \___|_|\___\___/|_| |_| |_|\___|"
   write (*, *)
!========================================================================================!
   call env%ref%to(mol)
   write (*, *)
   write (*, *) 'Input structure:'
   call mol%append(stdout)
   write (*, *)
!!========================================================================================!
!
   pr = .true.

   call rdensemble(env%inputcoords, nall, structures)
   allocate (env%calc%chess)
   nat3 = 3*structures(1)%nat

   ! nstruc = ceiling(real(nall)/real(env%calc%chess_space))
   ! write (*, *) "Number of Strucs:", nstruc
   allocate (list(nall-1))
   init_mol = structures(1)
  call env%calc%chess%alloc(structures(1)%nat,7,env%calc%initialize_hr_type,env%calc%hr_hu_type,hguess=env%calc%chess_id_guess)
   call initialize_hessian(env%calc, env%calc%chess%initialize_type, structures(1)%xyz, &
     & init_mol%nat, init_mol%at, env%calc%chess%hess(:), env%calc%chess%hguess, pr)
   call dhtosq(nat3, env%calc%chess%H, env%calc%chess%hess(:))
   call generate_chess_list(env%calc, nall, list, nstruc, structures, env%calc%cos_thresh,S,Y, env%calc%chess%H)
            write (*, *) "S has NaN =", any(ieee_is_nan(S))
            write (*, *) "S has Inf =", any(.not. ieee_is_finite(S))

   ! do i = 1, nall - 1, env%calc%chess_space
   !    call env%calc%chess%update(structures(i)%gradient, structures(i)%xyz)
   ! end do

   ! do i = 1, nall
   !    if (list(i)) then
   !       call env%calc%chess%update(structures(i)%gradient, structures(i)%xyz)
   !    end if
   ! end do

   ! if (nstruc > 1) then
      ! idx = minloc(env%calc%chess%order, 1)
      ! if (minval(env%calc%chess%order) .eq. 0) idx = 1
   ! else
   !    idx = nall
   ! end if
   write(*,*) nstruc, idx

   allocate (init_hess(nat3, nat3))
   if (nstruc > 1) call env%calc%chess%construct_hessian(S,Y,list,nstruc)
   if (nstruc == 1) call dhtosq(nat3, env%calc%chess%H,env%calc%chess%hess)
   ! call dhtosq(nat3,init_hess,env%calc%chess%hess)

   ! call prj_hess(structures(nall)%nat, nat3, structures(nall)%xyz, env%calc%chess%H)
   ! allocate(freqs(nat3), final_hess(nat3,nat3))
   ! call diagonalize_matrix(nat3, env%calc%chess%H,freqs,info)

   ! allocate(hess(nat3,nat3))
   ! hess = env%calc%chess%H !Will store modes
   ! allocate(quality(nat3))

   ! allocate(init_eigenvalues(nat3))
   ! call prj_hess(structures(idx)%nat, nat3, structures(idx)%xyz, init_hess)
   ! call diagonalize_matrix(nat3, init_hess,init_eigenvalues, info)

   ! call env%calc%chess%mode_quality_analysis(1,init_hess,hess,nat3,quality)
   ! quality = 0.000000001_wp
   ! quality(1:20) = 0.00001
   pr = .true.
   ! call selective_hessian_repair_v2(structures(nall),env%calc,hess,freqs,quality,nat3,0.2_wp,0.005_wp,n_repaired, final_hess,pr)
   !
   etot = structures(nall)%energy
   !  call calcthermo_from_modes(structures(nall),&
   !    & freqs, pr, env%calc%ithr, env%calc%fscal, env%calc%sthr, env%calc%nt, env%calc%temperatures, &
   !    & env%calc%et, env%calc%ht, env%calc%gt, env%calc%stot, etot, emodel=env%calc%emodel)

   call calc_thermo_from_hess(structures(nall), env%calc%chess%H, pr, &
   & env%calc%nt, env%calc%temperatures, env%calc%ithr, env%calc%fscal, env%calc%sthr, env%calc%et, &
   & env%calc%ht, env%calc%gt, env%calc%stot, etot, env%calc%emodel)
!========================================================================================!
   call tim%stop(14)
   return
end subroutine hess_from_trj
