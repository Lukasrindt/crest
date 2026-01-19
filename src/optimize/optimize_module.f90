!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
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
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

!> This module wrapps the different optimization algorithms,
!> i.e., this is what can be called for geometry opt.

module optimize_module
  use iso_fortran_env,only:wp => real64
  use crest_parameters
  use crest_calculator
  use strucrd
  use ancopt_module
  use gradientdescent_module
  use rfo_module
  use lbfgs_module
  use optimize_utils
  use thermochem_module
  use hessian_reconstruct
  use newton_raphson_module
  use hr_utils
  !use hessian_tools
  implicit none
  private

  public :: optimize_geometry
  public :: print_opt_data

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine optimize_geometry(mol,molnew,calc,etot,grd,pr,wr,iostatus)
    implicit none
    !> Input
    type(coord)    :: mol
    type(calcdata) :: calc
    logical,intent(in)        :: pr
    logical,intent(in)        :: wr
    !> Output
    type(coord)   :: molnew
    integer,intent(out)       :: iostatus
    real(wp),intent(inout)    :: etot
    real(wp),intent(inout)    :: grd(3,mol%nat)
    real(wp),allocatable :: H_init(:,:),freq(:)
    integer :: nat3
    integer :: io,idx
    integer :: i,j, stepno, step
    real(wp) :: identity(6), steps_incl(5)


    iostatus = -1
    !> do NOT overwrite original geometry
    !$omp critical
    molnew%at = mol%at
    molnew%xyz = mol%xyz
    molnew%nat = mol%nat
    !$omp end critical
    nat3 = 3*mol%nat

    !> Check for optimization-individual calculation setup
    if (calc%optnewinit) then
      !$omp critical
      call calc%dealloc_params()
      !$omp end critical
    end if

    !> Check if Hessian Reconstruct is called and initialize the type
    if (calc%do_HR) then
      allocate (calc%chess)
      allocate (H_init(nat3,nat3))
      call calc%chess%alloc(mol%nat,calc%hu_steps,calc%hguess,calc%initialize_hr_type, calc%hr_hu_type)
    end if

    !> initial singlepoint
    if (calc%do_HR) calc%chess%track_step = .false. !initial calc should not be tracked to avoid duplicate
    call engrad(molnew,calc,etot,grd,iostatus)
    if (calc%do_HR) calc%chess%track_step = .true.

    !> optimization
    select case (calc%opt_engine)
    case (0)
      call ancopt(molnew,calc,etot,grd,pr,wr,iostatus)
    case (1)
      !> l-bfgs goes here
      !write(stdout,'(a)') 'L-BFGS currently not implemented'
      !stop
      call lbfgs_optimize(molnew,calc,etot,grd,pr,iostatus)
    case (2)
      !> rfo goes here
      call rfopt(molnew,calc,etot,grd,pr,wr,iostatus)
    case (3)
      !> newton-raphson step goes here, this is a newton step with updated hessians, i.e. quasi Newton
      call newton_raphson(molnew,calc,etot,grd,pr,wr,iostatus)
    case (-1)
      call gradientdescent(molnew,calc,etot,grd,pr,wr,iostatus)
    case default
      write (stdout,'(a)') 'Unknown optimization engine!'
      stop
    end select
    molnew%energy = etot


    if (calc%do_HR .and. iostatus .eq. 0) then !> Hessian construction and post-processing happen here, only do it if geometry relaxation successful
    !identity = [0.001_wp,0.01_wp,0.02_wp,0.1_wp,0.5_wp,1.0_wp]
    !identity = [0.02_wp,0.04_wp,0.06_wp,0.08_wp,0.10_wp,0.12_wp,0.14_wp,0.16_wp,0.18_wp,0.20_wp]
    steps_incl = [0.1_wp,0.5_wp,1.0_wp,2.0_wp,100.0_wp]
    do i = 1,5 !> loop for hessian init
    if (i .ne. 4 .and. i .ne. 2 .and. i .ne. 3) then !> dont do gfn0 1 or 2 guesses duh
    do j = 1,5
    if (calc%do_HR) then !> Hessian construction and post-processing happen here

        write (stdout,*)
        write (stdout,*) "THERMO FROM RECONSTRUCTED HESSIAN with", i,steps_incl(j) 
        write (stdout,*)
      
        step = nint(steps_incl(j)*molnew%nat)
        stepno = max(5,step)
        idx = calc%chess%made_iters - stepno
        if (idx<0) then
        idx = 1
        endif
        call initialize_hessian(calc,i,calc%chess%coords(idx,:,:),molnew%nat,molnew%at,calc%chess%hess(:),calc%chess%hguess,pr)  !> This hguess is set through the hguess variable of the optimizer and needs to be hardcoded/set explicitly before initialization for benchmarking!!
        !call dhtosq(nat3,H_init,calc%chess%hess) !> maybe this should all be inside the construct bfgs function later? -> cannot due to circular import!!!
        !write(stdout,*)                                                                                                   !> Hessian type (gfnff,mod,identity) is set through input file and is already encoded into the calc object

        call calc%chess%construct_hessian_bfgs_stepsvar(stepno)

        call calc_thermo_from_hess(molnew,calc%chess%H(:,:),pr, &
        & calc%nt,calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
        & calc%ht,calc%gt,calc%stot,i,j,etot)
    
    end if
    enddo
    endif
    enddo
    endif

    if (allocated(calc%chess)) then
      call calc%chess%dealloc()
      deallocate (calc%chess)
    endif

    return
  end subroutine optimize_geometry

!========================================================================================!

  subroutine print_opt_data(calc,ich,natoms,tag)
    implicit none
    type(calcdata) :: calc
    integer,intent(in) :: ich
    integer,intent(in),optional :: natoms
    character(len=*),intent(in),optional :: tag
    integer :: tight,nat
    real(wp) :: ethr,gthr
    character(len=:),allocatable :: ttag
    if (present(tag)) then
      ttag = tag
    else
      ttag = ' '
    end if
    if (present(natoms)) then
      nat = natoms
    else
      nat = 0
    end if

    write (ich,'(a,a)',advance='no') ttag,'Optimization engine: '
    select case (calc%opt_engine)
    case (0)
      write (ich,'(a)') 'ANCOPT'
    case (1)
      write (ich,'(a)') 'L-BFGS'
    case (2)
      write (ich,'(a)') 'RFO'
    case (-1)
      write (ich,'(a)') 'Gradient Descent'
    case default
      write (ich,'(a)') 'Unknown'
    end select
    if (calc%opt_engine >= 0) then
      write (ich,'(a,a)',advance='no') ttag,'Hessian update type: '
      select case (calc%iupdat)
      case (0)
        write (ich,'(a)') 'BFGS'
      case (1)
        write (ich,'(a)') 'Powell'
      case (2)
        write (ich,'(a)') 'SR1'
      case (3)
        write (ich,'(a)') 'Bofill'
      case (4)
        write (ich,'(a)') 'Farkas-Schlegel'
      end select
    end if

    tight = calc%optlev
    call get_optthr(nat,tight,calc,ethr,gthr)
    write (ich,'(a,a,e10.3,a,e10.3,a)') ttag,'E/G convergence criteria: ',&
    & ethr,' Eh,',gthr,' Eh/a0'

    write (ich,'(a,a,i0)') ttag,'maximum optimization steps: ',calc%maxcycle

  end subroutine print_opt_data

!========================================================================================!
!========================================================================================!
end module optimize_module
