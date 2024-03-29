module add_procedures
! this module contains useful procedures for visualizing of calculations
! and working with arrays
    implicit none
    
    contains
    
    subroutine progress(string, ndone, ntotal)
    ! show progress bar of calculation
        character*(*) string
        character*255 prog, oldprog
        integer ndone, ntotal, i, z

        if (100.0*ndone/ntotal .GE. 100.0d0) then
            write(0,'(a,$)') '                                                                                   ',char(13)
        else
            write(prog,'(a20,1x,''['')') trim(string)
            do i=1,40
                prog(22+i:22+i)=' '
            enddo
            write(prog(39:45),'(f5.2,''%'')') 100.0*ndone/ntotal
            do i=1,40
                if ((1.0*ndone/ntotal).gt.(1.0*i/40)) then
                    if (prog(22+i:22+i).eq.' ') prog(22+i:22+i)='-'
                endif
            enddo
            prog(62:62)=']'
            write(0,'(a,a,$)') prog(1:77),char(13)
            return
        endif
    end subroutine progress
    
    subroutine get_path_separator(path_sep, read_well)
    ! check type of system and return path separator according to it
        CHARACTER(len=1), intent(out) :: path_sep   ! path separator
        logical, intent(inout) :: read_well ! did the data read well?
        CHARACTER(len = 100) :: path, Err_data
        CALL get_environment_variable("PATH",path)
        read_well = .true.
        if (path(1:1) .EQ. '/') then        !unix based OS
            path_sep = '/'
        else if (path(3:3) .EQ. '\') then   !Windows OS
            path_sep = '\'
        else
            Err_data = 'Path separator is not defined' ! unknown OS
            print*, '----------------'
            print*, Err_data
            print*, '----------------'
            read_well = .false. ! didn't read well this data
        endif
    end subroutine
    
    subroutine concatenate(array, subarray)
    ! analog of numpy.concatenate for 1d arrays in python
        integer :: isize, esize
        real(8), dimension(:), intent(in) :: subarray
        real(8), dimension(:), allocatable, intent(inout) :: array
        real(8), dimension(:), allocatable :: clist

        esize = size(subarray)
        if(allocated(array)) then
            isize = size(array)
            allocate(clist(isize+esize))
            clist(1:isize) = array(:)
            clist(isize+1:isize+esize) = subarray(:)
            deallocate(array)
            call move_alloc(clist, array)
        else
            allocate(array(esize))
            array(:) = subarray(:)
        end if
        return
    end subroutine concatenate
    
    subroutine append(array, element)
    ! analog of numpy.append for 1d arrays
        integer :: isize
        real(8), intent(in) :: element
        real(8), dimension(:), allocatable, intent(inout) :: array
        real(8), dimension(:), allocatable :: clist

        if(allocated(array)) then
            isize = size(array)
            allocate(clist(isize+1))
            clist(1:isize) = array(:)
            clist(isize+1) = element
            deallocate(array)
            call move_alloc(clist, array)
        else
            allocate(array(1))
            array(1) = element
        end if
        return
    end subroutine append
    
    subroutine column_stack(array, column)
    ! analog of numpy.column_stack
        integer :: size_1, size_2
        real(8), dimension(:), intent(in) :: column ! 1d array need to stack
        real(8), dimension(:,:), allocatable, intent(inout) :: array ! 2d array to which stack
        real(8), dimension(:,:), allocatable :: carr
        
        if(allocated(array)) then
            size_1 = size(array, 1)
            size_2 = size(array, 2)
            allocate(carr(size_1 + 1, size_2))
            carr(1:size_1, :) = array(:, :)
            carr(size_1+1, :) = column(:)
            deallocate(array)
            call move_alloc(carr, array)
        else
            allocate(array(1, size(column)))
            array(1, :) = column(:)
        end if            
        return
    end subroutine Column_stack
    
    pure subroutine Interpolate(Iflag, E1, E2, Sigma1, Sigma2, E_needed, OUT_value)
        ! Interpolation between two values with different methods:
        integer, intent(in) :: Iflag ! what kind of interpolation to use
        real(8), intent(in) :: E1, E2, Sigma1, Sigma2, E_needed  ! input data: X and Y points
        real(8), intent(out) :: OUT_value    ! interpolated value
        real(8) E2log, E1log, E_needed_log, Sigma1log, Sigma2log
        select case(Iflag) ! what interpolation to use:
            case(1) ! linear x and y
                OUT_value = Sigma1 + (Sigma2 - Sigma1)/(E2 - E1)*(E_needed - E1)
            case(3)	! logarithmic x, linear y
                E2log = log(E2)
                E1log = log(E1)
                E_needed_log = log(E_needed)
                OUT_value = Sigma1 + (Sigma2 - Sigma1)/(E2log - E1log)*(E_needed_log - E1log)
            case(4)	! linear x, logarithmic y
                Sigma1log = log(Sigma1)
                Sigma2log = log(Sigma2)
                OUT_value = Sigma1log + (Sigma2log - Sigma1log)/(E2 - E1)*(E_needed - E1)
                OUT_value = exp(OUT_value)
            case(5)	! logarithmic x and y
                E2log = log(E2)
                E1log = log(E1)
                E_needed_log = log(E_needed)
                Sigma1log = log(Sigma1)
                Sigma2log = log(Sigma2)
                OUT_value = Sigma1log + (Sigma2log - Sigma1log)/(E2log - E1log)*(E_needed_log - E1log)
                OUT_value = exp(OUT_value)
            case default ! linear x and y
                OUT_value = Sigma1 + (Sigma2 - Sigma1)/(E2 - E1)*(E_needed - E1) 
        end select
    end subroutine Interpolate

    subroutine linspace(sp, ep, Numpoints, arr)
    ! create evenly spaced grid for layers
    real(8), intent(in) :: sp   ! start point
    real(8), intent(in) :: ep  ! end point point
    integer, intent(in) :: Numpoints    !   points in grid
    real(8), dimension(:), intent(out) :: arr ! resulting array
    
    integer i
    
    do i = 1, Numpoints
        arr(i) = sp + (ep - sp) * (i-1) / (Numpoints - 1)
    enddo

    return    
    end subroutine linspace  
    
    pure subroutine Find_in_1D_array(Array, Value, Number)
        REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
        REAL(8), INTENT(in) :: Value   ! to be found in the array as near as possible
        integer, INTENT(out) :: Number ! number of the element which we are looking for 
        integer i
        i = 1
        do while (Array(i) .LT. Value)
            i = i + 1
        enddo
        Number = i
    end subroutine Find_in_1D_array

    pure subroutine Find_in_monotonous_1D_array(Array, Value0, Number)
       REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
       REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
       integer, INTENT(out) :: Number ! number of the element which we are looking for 
       integer i, N, i_cur, i_1, i_2, coun
       real(8) temp_val, val_1, val_2

           N = size(Array)
           if (Value0 .LE. Array(1)) then ! it's the first value, no need to search
               i_cur = 1
           else if (Value0 .GE. Array(N)) then ! it's the last value, no need to search
               i_cur = N
           else
               i_1 = 1
               i_2 = N
               i_cur = FLOOR((i_1+i_2)/2.0)
               temp_val = Array(i_cur)
               coun = 0
               do ! until the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                    if ((Value0 .GE. Array(i_cur)) .AND. (Value0 .LT. Array(i_cur+1))) exit ! when the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                    if (temp_val .LE. Value0) then
                       i_1 = i_cur
                    else
                       i_2 = i_cur
                    endif
                    i_cur = FLOOR((i_1+i_2)/2.0)
                    temp_val = Array(i_cur)
                    coun = coun + 1
               enddo
           endif
       Number = i_cur
    end subroutine Find_in_monotonous_1D_array

    
    
end module add_procedures