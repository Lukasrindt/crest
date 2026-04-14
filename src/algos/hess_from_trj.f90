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
   integer :: nall, steps, nstruc,nat3, n_repaired
   real(wp) :: etot
   logical, allocatable :: list(:)
   real(wp), allocatable :: hess(:,:), freqs(:), quality(:), final_hess(:,:)
   

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
   allocate (list(nall))
   call generate_chess_list(env%calc, nall, list, nstruc, structures)
  call env%calc%chess%alloc(structures(1)%nat,nstruc,env%calc%initialize_hr_type,env%calc%hr_hu_type,hguess=env%calc%chess_id_guess)

   ! do i = 1, nall - 1, env%calc%chess_space
   !    call env%calc%chess%update(structures(i)%gradient, structures(i)%xyz)
   ! end do

   do i = 1, nall
      if (list(i)) then
         call env%calc%chess%update(structures(i)%gradient, structures(i)%xyz)
      end if
   end do

   ! call env%calc%chess%update(structures(nall)%gradient, structures(nall)%xyz)

   ! write (*, *) env%calc%chess%order

   idx = maxloc(env%calc%chess%order, 1)
   ! if (minval(env%calc%chess%order) .eq. 0) idx = 1

   init_mol = structures(idx)
   call initialize_hessian(env%calc, env%calc%chess%initialize_type, env%calc%chess%coords(idx, :, :), &
     & init_mol%nat, init_mol%at, env%calc%chess%hess(:), env%calc%chess%hguess, pr)

   call env%calc%chess%construct_hessian()

  call prj_hess(structures(nall)%nat, nat3, structures(nall)%xyz, env%calc%chess%H)
  allocate(freqs(nat3), final_hess(nat3,nat3))
  call diagonalize_matrix(nat3, env%calc%chess%H,freqs,info)

   ! call prj_mw_hess(structures(nall)%nat, structures(nall)%at, nat3, structures(nall)%xyz,env%calc%chess%H)
   allocate(hess(nat3,nat3))
   hess = env%calc%chess%H !Will store modes
   ! call frequencies(structures(nall)%nat, structures(nall)%at, structures(nall)%xyz,nat3,hess,freqs,io)
   allocate(quality(nat3))

   call env%calc%chess%mode_quality_analysis(hess,nat3,quality)
   call selective_hessian_repair(structures(nall),env%calc,hess,freqs,quality,nat3,0.005_wp,0.001_wp,n_repaired, final_hess)
  !
   etot = structures(nall)%energy
  !  call calcthermo_from_modes(structures(nall),&
  !    & freqs, pr, env%calc%ithr, env%calc%fscal, env%calc%sthr, env%calc%nt, env%calc%temperatures, &
  !    & env%calc%et, env%calc%ht, env%calc%gt, env%calc%stot, etot, emodel=env%calc%emodel)

   call calc_thermo_from_hess(structures(nall), final_hess, pr, &
   & env%calc%nt, env%calc%temperatures, env%calc%ithr, env%calc%fscal, env%calc%sthr, env%calc%et, &
   & env%calc%ht, env%calc%gt, env%calc%stot, etot, env%calc%emodel)
!========================================================================================!
   call tim%stop(14)
   return
end subroutine hess_from_trj
