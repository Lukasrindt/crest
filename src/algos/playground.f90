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
subroutine crest_playground(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use hessian_reconstruct
  use hr_utils
  use thermochem_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc
  real(wp) :: accuracy,etemp

  integer :: V,maxgen
  integer,allocatable :: A(:,:)
  logical,allocatable :: rings(:,:)
  integer,allocatable :: tmp(:)
  logical :: connected,fail,doreturn

  real(wp) :: energy
  real(wp),allocatable :: grad(:,:),geo(:,:),csv(:,:),q(:)

!========================================================================================!
  call tim%start(14,'Test implementation')
!========================================================================================!
  !call system('figlet welcome')
  write (*,*) "              _                          "
  write (*,*) "__      _____| | ___ ___  _ __ ___   ___ "
  write (*,*) "\ \ /\ / / _ \ |/ __/ _ \| '_ ` _ \ / _ \"
  write (*,*) " \ V  V /  __/ | (_| (_) | | | | | |  __/"
  write (*,*) "  \_/\_/ \___|_|\___\___/|_| |_| |_|\___|"
  write (*,*)
!========================================================================================!
  call env%ref%to(mol)
  write (*,*)
  write (*,*) 'Input structure:'
  call mol%append(stdout)
  write (*,*)
!!========================================================================================!
!
!  allocate (grad(3,mol%nat),source=0.0_wp)
  call env2calc(env,calc,mol)
  calc%calcs(1)%rdwbo = .true.
  ! call calc%info(stdout)
!
!  call engrad(mol,calc,energy,grad,io)
!  call calculation_summary(calc,mol,energy,grad)
!========================================================================================!

  ! allocate(mol%gradient(3,mol%nat), source=1.0_wp)
  ! call mol%write('dummy.extxyz')
  !
  !
  ! call molnew%open("dummy.extxyz")
  ! call molnew%write("dummy2.extxyz")
  
  pr = .true.

  block
    type(coord),allocatable :: structures(:)
    type(coord) :: init_mol
    ! type(cashed_hessian) :: chess 
    integer :: i, nall, steps
    real(wp) :: etot
     call rdensemble(env%inputcoords,nall,structures)
    ! write(*,*) nall,'structures read from ',env%inputcoords
    ! do i=1,5
    ! write(*,*) structures(1)%xyz(:,i)
    ! enddo
    ! do i=1,5
    ! write(*,*) structures(1)%gradient(:,i)
    ! enddo
    ! call wrensemble('dummyensemble.xyz',nall,structures)
    allocate(calc%chess)
    call calc%chess%alloc(structures(1)%nat,calc%hu_steps,calc%initialize_hr_type,calc%hr_hu_type)
    
    do i=1,nall
      call calc%chess%update(structures(i)%gradient,structures(i)%xyz)
    enddo

    init_mol = structures(1)
    call initialize_hessian(calc, calc%chess%initialize_type,calc%chess%coords(1,:,:), &
      &init_mol%nat,init_mol%at,calc%chess%hess(:),calc%chess%hguess,pr)

    call calc%chess%construct_hessian()
   
    write(stdout,*) calc%sthr, env%thermo%sthr
    call calc_thermo_from_hess(structures(nall), calc%chess%H,pr, &
    & calc%nt, calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
    & calc%ht, calc%gt, calc%stot, etot, calc%emodel)
  end block
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_playground
