!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> Implementation of whatever, for testing implementations
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> Input/Output:
!>  env  -  crest's systemdata object
!>  tim  -  timer object
!>-----------------------------------------------
subroutine crest_playground(env, tim)
   use crest_parameters
   use crest_data
   use crest_calculator
   use strucrd
   use hessian_reconstruct
   use hr_utils
   use irmsd_module
   use hessian_quality
   use thermochem_module
   use optimize_maths
   use hessupdate_module
   implicit none
   type(systemdata), intent(inout) :: env
   type(timer), intent(inout)      :: tim
   type(coord) :: mol, molnew
   integer :: i, j, k, l, io, ich
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
!  allocate (grad(3,mol%nat),source=0.0_wp)
   ! call env2env%calc(env,env%calc,mol)
   ! env%calc%env%calcs(1)%rdwbo = .true.
   ! call env%calc%info(stdout)
!
!  call engrad(mol,env%calc,energy,grad,io)
!  call env%calculation_summary(env%calc,mol,energy,grad)
!========================================================================================!

   ! allocate(mol%gradient(3,mol%nat), source=1.0_wp)
   ! call mol%write('dummy.extxyz')
   !
   !
   ! call molnew%open("dummy.extxyz")
   ! call molnew%write("dummy2.extxyz")

   pr = .true.

   block
      type(coord), allocatable :: structures(:), structures2(:)
      type(coord) :: init_mol
      type(step_monitor) :: mon
      ! type(cashed_hessian) :: chess
      integer :: i, nall, steps, nstruc, last_struc, nat3, idx
      real(wp) :: etot, rmsdval, gnorm
      real(wp), allocatable :: s_vec(:), temp_s(:, :), Y(:, :), S(:, :)
      real(wp) :: gain
      logical:: accepted
      logical, allocatable :: acc_list(:)

      call rdensemble(env%inputcoords, nall, structures)
      allocate (env%calc%chess)
   nat3 = 3*structures(1)%nat

   ! allocate (list(nall-1))
   init_mol = structures(nall)

  call env%calc%chess%alloc(structures(1)%nat,7,env%calc%initialize_hr_type,env%calc%hr_hu_type,hguess=env%calc%chess_id_guess)
   call initialize_hessian(env%calc, env%calc%chess%initialize_type, structures(nall)%xyz, &
     & init_mol%nat, init_mol%at, env%calc%chess%hess(:), env%calc%chess%hguess, pr)
  call dhtosq(nat3, env%calc%chess%H, env%calc%chess%hess(:))

  do i=1,nat3
    env%calc%chess%H(i,i) = env%calc%chess%H(i,i) - 1d-9
  enddo

   pr = .true.
   etot = structures(nall)%energy
   call calc_thermo_from_hess(structures(nall), env%calc%chess%H, pr, &
   & env%calc%nt, env%calc%temperatures, env%calc%ithr, env%calc%fscal, env%calc%sthr, env%calc%et, &
   & env%calc%ht, env%calc%gt, env%calc%stot, etot, env%calc%emodel)
   end block

end subroutine crest_playground
