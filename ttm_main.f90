! Include all separate modules
include 'Reading.f90'
include 'Objects.f90'
include 'Constants.f90'
include 'Transfer_Matrix.f90'
include 'Output.f90'
include 'additional_procedures.f90'
include 'Solver.f90'
!****************************************************************************

program ttm_main
!!!!!!!!!!!!!!!!!!!!!!!!
! load modules
use Reading
use Objects
use Constants
use Transfer_Matrix
use Output
use IFPORT
use add_procedures
use Solver
!!!!!!!!!!!!!!!!!!!!!!!!
implicit none

! Variables
type(Multilayer) intarget
type(source) laser
type(TTM) params
character(100) output_path
! Body of ttm_main
output_path = 'Saved'

call Initialize_objects(laser, intarget, params)
call time_evolution(laser, intarget, params)
call Save_TTM(params, laser, output_path)
!call save_absorption(laser, intarget, output_path)
call Deallocate_all(InTarget, params)

print*, 'Calculation done!'
pause
end program ttm_main

