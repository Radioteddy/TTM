! Include all separate modules
include 'Objects.f90'
include 'Constants.f90'
include 'additional_procedures.f90'
include 'parameter_functions.f90'
include 'Transfer_Matrix.f90'
include 'Output.f90'
include 'Reading.f90'
include 'Solver.f90'
include 'Bspline.f90'
!****************************************************************************

program ttm_main

!!!!!!!!!!!!!!!!!!!!!!!!
    ! load modules
    !use bspline_module
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
    integer, dimension(8) :: curvalues

    ! Body of ttm_main
    call date_and_time(VALUES=curvalues)
    print'(a)', '***************************************************************'
    print'(a)', '       Two-temperature model written by F. Akhmetov            '
    print'(a)', '***************************************************************'
    output_path = 'Saved'
    call Initialize_objects(laser, intarget, params)
    call time_evolution(laser, intarget, params)
    call Save_TTM(params, laser, output_path)
    !call calculate_absorption(laser, intarget)
    !call save_absorption(laser, intarget, output_path)
    call Deallocate_all(InTarget, params)
    print'(a)', '***************************************************************'
    print 123, 'Calculation is started at ', curvalues(5),curvalues(6),curvalues(7),curvalues(8),curvalues(3),curvalues(2),curvalues(1)
    print '(a, f6.2, a)', 'Total calculation time is: ', params%calctime, ' seconds' 
    call date_and_time(VALUES=curvalues)
    print 123, 'Calculation is finished at ', curvalues(5),curvalues(6),curvalues(7),curvalues(8),curvalues(3),curvalues(2),curvalues(1)
    print'(a)', '***************************************************************'
    123    format (a, i2.2, ':', i2.2, ':', i2.2, ':', i3.3, ' on ', i2.2, '-', i2.2, '-', i4.4)
    !pause 'Calculation done! Press Enter to return...'

    end program ttm_main
    
    
    
    ! Desirable upgrades:
    ! - HDF5 file format for profiles
    ! - support of tabulated TTM parameters
    ! - addition of flags for various functional dependencies
    ! - pasrser of user-defined functions
    ! - add account of transmission to Beer-Lambert source and homogeneous heating