module Reading

    use Objects
    use Constants
    use add_procedures
    use bspline_sub_module, only: db1ink, get_status_message
    
    implicit none

    contains

    subroutine Read_layer_file(matnames, layer, nlayer, syscomp, npoints)
        ! material parameters of multilayer system and laser specifications
        character(20), dimension(:), allocatable, intent(out) :: matnames      ! names of every layer
        real(8), dimension(:), allocatable, intent(out) :: layer      ! width of layer [nm]
        complex(8), dimension(:), allocatable, intent(out) :: nlayer  ! optical constants per layer
        character(100), intent(out) :: syscomp  ! elements in system
        integer, intent(out) :: npoints
        
        ! internal variables
        integer FN, Reason, i, j, N_layers
        character(100) temp
        character(100) filename ! file with input parameters
        character(100) Error_descript  ! to write a description of an error, if any
        logical file_exist    ! to check where file to be open exists
        logical file_opened   ! to check if a file is still opened
        logical read_well     ! to check if we read the file without errors

        filename='layer_parameters.txt'
        FN=100
        inquire(file=trim(adjustl(filename)), exist=file_exist) ! check if file exists
        if (file_exist) then
            open(unit=FN, file=trim(adjustl(filename)), status='old', action='read') ! if yes, open and read
        else ! error message about it:
            Error_descript = 'File '//trim(adjustl(filename))//' is not found!' ! if no, print an error
            print*, trim(adjustl(Error_descript))
        endif

        i = 0 ! start reading      
        syscomp = ''
        READ(FN,*,IOSTAT=Reason) N_layers
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 1111
        if (.not. allocated(nlayer)) allocate(nlayer(N_layers))
        if (.not. allocated(layer)) allocate(layer(N_layers))
        if (.not. allocated(matnames)) allocate(matnames(N_layers))
        do j = 1, N_layers
            READ(FN,*,IOSTAT=Reason) temp, nlayer(j), layer(j) 
            call read_file(Reason, i, read_well) ! reports if everything read well
            if (.not. read_well) goto 1111
            write(*,'(a,i1,a,a)') 'Layer ', j, ' is: ', trim(adjustl(temp))   
            write(*,'(a, i1, a, f10.3, sp, f10.3, "i")') 'refractive index of ', j, ' layer is: ', nlayer(j)
            write(*,'(a, i1, a, e10.3)') 'width of ', j, ' layer is: ', layer(j)
            matnames(j) = trim(adjustl(temp))
            syscomp = trim(adjustl(syscomp))//'_'//trim(adjustl(temp))
        enddo
        write(*,'(a)') '---------------------------------------------------------------'
        READ(FN,*,IOSTAT=Reason) npoints
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 1111
        

        1111 if (.not. read_well) print*, 'Error in '//trim(adjustl(filename))//' file!' 
        inquire(unit=FN, opened=file_opened)    ! check if this file is opened
        if (file_opened) close(FN)             ! and if it is, close it
    end subroutine Read_layer_file
    
    subroutine read_TTM_file(laser, parameters, Npoints)
    ! read parameters for TTM calculation
        type(Source), intent(out) :: laser    
        type(TTM), intent(out) :: parameters ! object with parameters of TTM
        integer, intent(out) :: Npoints
        ! internal variables
        character(100) matname
        real(8) Te, Tl, Ce, Cl, ke, kl, Gel, dtt, tfin, dl, lb, ts, Cli, Hfus, Tm, A, B, D
        complex(8) nl
        integer FN, Reason, i, j, fl, fl_cl, fl_mel, fl_cel, fl_ke, fl_g
        character(100) temp
        character(100) filename ! file with input parameters
        character(100) Ce_filename ! file with Ce data
        character(100) Error_descript  ! to write a description of an error, if any
        logical file_exist    ! to check where file to be open exists
        logical file_opened   ! to check if a file is still opened
        logical read_well     ! to check if we read the file without errors
        ! B-spline interpolation variables
        real(8), dimension(:,:), allocatable :: Ce_tab ! tabulated values of electron heat capacity
        integer iflag
        
        filename = 'TTM_parameters.txt'
        FN = 200
        inquire(file=trim(adjustl(filename)), exist=file_exist) ! check if file exists
        if (file_exist) then
            open(unit=FN, file=trim(adjustl(filename)), status='old', action='read') ! if yes, open and read
        else ! error message about it:
            Error_descript = 'File '//trim(adjustl(filename))//' is not found!' ! if no, print an error
            print*, trim(adjustl(Error_descript))
        endif
        i = 0
        READ(FN,*,IOSTAT=Reason) matname 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) ! skip indicator of electron subsystem section
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222  
        READ(FN,*,IOSTAT=Reason) laser%lambda 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222   
        READ(FN,*,IOSTAT=Reason) laser%in_angle 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222  
        READ(FN,*,IOSTAT=Reason) temp   
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        laser%pol = trim(adjustl(temp))   
        READ(FN,*,IOSTAT=Reason) laser%Fluence
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222   
        READ(FN,*,IOSTAT=Reason) laser%tpulse 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) laser%tpeak 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) fl 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) nl, dl, lb
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) ! skip indicator of electron subsystem section
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222        
        READ(FN,*,IOSTAT=Reason) Te 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) fl_cel 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222 
        READ(FN,*,IOSTAT=Reason) Ce_filename 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222 
        READ(FN,*,IOSTAT=Reason) Ce 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) fl_ke 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222 
        READ(FN,*,IOSTAT=Reason) ke 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) A, B, D 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222 
        READ(FN,*,IOSTAT=Reason) fl_g 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222 
        READ(FN,*,IOSTAT=Reason) Gel 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222        
        READ(FN,*,IOSTAT=Reason) ! indicator of ion subsystem section
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222        
        READ(FN,*,IOSTAT=Reason) Tl 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) fl_cl 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) Cl 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) fl_mel 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) Tm 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) Hfus 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) Cli 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) kl 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) ! indicator of numercial parameters section
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222        
        READ(FN,*,IOSTAT=Reason) Npoints
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) dtt
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) ts
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222        
        READ(FN,*,IOSTAT=Reason) tfin
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        
        ! Print initial parameters onto terminal
        write(*,'(a)') 'Laser paramteres are:'
        select case(fl)
            case(0)
                write(*,'(a)') 'Homogeneous absorption profile'
            case(1)
                write(*,'(a)') 'Homogeneous absorption profile * (1-R)'
            case(2)
                write(*,'(a)') 'Beer-Lambert absorption profile'
            case(3)
                write(*,'(a)') 'Multilayer (TMM) absorption profile'
        end select 
        write(*,'(a,f8.3,a,f6.3,a,a,a)') 'lambda = ', laser%lambda, '; theta_0 = ', laser%in_angle, '; polarisation is ', laser%pol, ' type'
        write(*, '(a,f10.3,a,e10.3,a,e10.3)') 'F_inc = ', laser%fluence, '; FWHM = ', laser%tpulse, '; t_peak = ', laser%tpeak
        write(*,'(a)') '---------------------------------------------------------------'

        write(*,'(a)') 'Initial temperatures for TTM are:'
        write(*,'(a,f8.3,a,f8.3)') 'Te = ', Te, '; Tlat = ', Tl
        if (fl_mel .ne. 0) then
            write(*,'(a)') '---------------------------------------------------------------'
            write(*,'(a)') 'Melting parameters are:'
            write(*, '(a, 3e12.3)') 'Tmelt, Hf, Cliq: ', Tm, Hfus, Cli
        endif
        write(*,'(a)') '---------------------------------------------------------------'
        write(*,'(a)') 'Simulation parameters are:'
        write(*, '(a,1x,i5,1x,3e12.3)') 'Ncells, dt, tsave, tend:', Npoints, dtt, ts, tfin
        write(*,'(a)') '---------------------------------------------------------------'
        ! save to "parameters" object
        parameters%mat = trim(adjustl(matname))
        parameters%Clat = Cl
        parameters%Cliq = Cli
        parameters%Hf = Hfus
        parameters%Tmelt = Tm
        parameters%klat = kl
        parameters%Cel = Ce
        parameters%kel = ke
        parameters%G = Gel
        parameters%dt = dtt
        parameters%tsave = ts
        parameters%tend = tfin
        parameters%flag = fl
        parameters%cl_flag = fl_cl
        parameters%cel_flag = fl_cel
        parameters%kel_flag = fl_ke
        parameters%g_flag = fl_g
        parameters%melt_flag = fl_mel
        parameters%dlayer = dl*1.d-9 ! convert from nm to m
        parameters%l_bal = lb*1.d-9  ! convert from nm to m
        parameters%n = nl
        
        
        ! if Ce tabulated read the table and compute spline coeffs
        if (fl_cel .eq. 1) then
            call read_2Ddatafile(ce_filename, Ce_tab)
            parameters%ce_spline%nx = size(Ce_tab,1)
            parameters%ce_spline%kx = 4
            if (.not. allocated(parameters%ce_spline%tx)) allocate(parameters%ce_spline%tx(parameters%ce_spline%nx+parameters%ce_spline%kx))
            if (.not. allocated(parameters%ce_spline%bcoef)) allocate(parameters%ce_spline%bcoef(parameters%ce_spline%nx))
            call db1ink(Ce_tab(:,1), parameters%ce_spline%nx, Ce_tab(:,2), parameters%ce_spline%kx, 0, parameters%ce_spline%tx, parameters%ce_spline%bcoef, iflag)
            if (iflag .ne. 0) print*, 'Error initializing B-spline: '//get_status_message(iflag)
        endif
        
        if (.not. allocated(parameters%Tel)) allocate(parameters%Tel(Npoints))
        parameters%Tel(:) = Te
        
        if (.not. allocated(parameters%Tlat)) allocate(parameters%Tlat(Npoints))
        parameters%Tlat(:) = Tl
        
        ! check if all is read well
2222    if (.not. read_well) print*, 'Error in '//trim(adjustl(filename))//' file!' 
        inquire(unit=FN, opened=file_opened)    ! check if this file is opened
        if (file_opened) close(FN)             ! and if it is, close it
    end subroutine Read_TTM_file
        
    subroutine read_2Ddatafile(filename, array)
        character(100), intent(in) :: filename ! filename
        real(8), dimension(:,:), allocatable, intent(inout) :: array ! array where store dat
        integer :: i, Nlin, Ncol, Reason, count_lines, FN, counter
        character(100) Error_descript  ! to write a description of an error, if any
        logical file_exist    ! to check where file to be open exists
        logical file_opened   ! to check if a file is still opened
        logical read_well     ! to check if we read the file without errors  
        !
        FN=600
        inquire(file=trim(adjustl(filename)), exist=file_exist) ! check if file exists
        if (file_exist) then
            open(unit=FN, file=trim(adjustl(filename)), status='old', action='read') ! if yes, open and read
        else ! error message about it:
            Error_descript = 'File '//trim(adjustl(filename))//' is not found!' ! if no, print an error
            print*, trim(adjustl(Error_descript))
        endif
        !
        call Count_lines_in_file(FN, Nlin)
        call Count_columns_in_file(FN, Ncol)
        if (.not. allocated(array)) then 
            allocate(array(Nlin, Ncol))
        else
            deallocate(array)
            allocate(array(Nlin, Ncol))
        endif
        !
        counter = 0
        do i = 1, Nlin
            READ(FN,*,IOSTAT=Reason) array(i, :Ncol)
            call read_file(Reason, counter, read_well) ! reports if everything read well
            if (.not. read_well) goto 3333
        end do
        !    
3333    if (.not. read_well) print*, 'Error in '//trim(adjustl(filename))//' file!' 
        inquire(unit=FN, opened=file_opened)    ! check if this file is opened
        if (file_opened) close(FN)             ! and if it is, close it
    end subroutine read_2Ddatafile
    
    subroutine read_file(Reason, i, read_well)
        integer, intent(in) :: Reason    ! file number where to read from
        integer, intent(inout) :: i      ! line number
        logical, intent(inout) :: read_well  ! did we read ok?
        i = i + 1    ! it's next line
        IF (Reason .GT. 0)  THEN ! ... something wrong ...
            write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', wrong type of variable'
            read_well = .false.
        ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
            write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', unexpected END of file'
            read_well = .false.
        ELSE   ! normal reading
            read_well = .true.  ! it read well, nothing to report
        END IF
    end subroutine read_file         

    subroutine Initialize_objects(Laser, InTarget, parameters)
    !   initialize all atributes of used objects
        type(Source), intent(inout) :: Laser        !   object with parameters of laser
        type(Multilayer), intent(inout) :: InTarget !   object with parameters of the system
        type(TTM), intent(inout) :: parameters      ! object with TTM parameters
        integer j, Npoints, ind
        real(8) width
        
        ! Read initial parameters for TTM calculation
        call read_TTM_file(laser, parameters, Npoints)
        ! create spacegrid for TTM
        if (.not. allocated(parameters%X)) allocate(parameters%X(Npoints))
        select case(parameters%flag)
            case(0:2) ! Beer-Lambert law or homogeneous absorption
                call linspace(0.0d0, parameters%dlayer, Npoints, parameters%X)
            case(3) ! Trasfer-matrix
                ! Read initial parameters of multilayer from file
                call Read_layer_file(InTarget%materials, InTarget%layer, Intarget%nlayer, InTarget%syscompound, InTarget%points)
                ! create spacegrid
                ind = findloc(intarget%materials, parameters%mat, dim=1)
                parameters%dlayer = intarget%layer(ind)
                width = intarget%layer(ind) * 1.d-9 ! convert from [nm] to [m]     
                call linspace(0.0d0, width, Npoints, parameters%X)
                ! Allocate and initialize parameters for absorption calculation
                if (.not. allocated(InTarget%theta)) allocate(InTarget%theta(size(InTarget%layer)))
                InTarget%theta(1) = cmplx(Laser%in_angle * g_Pi / 180)
                InTarget%theta(2:) = (0.0d0, 0.0d0)
                if (.not. allocated(InTarget%phase)) allocate(InTarget%phase(size(InTarget%layer)))
                InTarget%phase = (0.0d0, 0.0d0)
                if (.not. allocated(InTarget%ref)) allocate(InTarget%ref(size(InTarget%layer-1)))
                InTarget%ref = (0.0d0, 0.0d0)
                if (.not. allocated(InTarget%trans)) allocate(InTarget%trans(size(InTarget%layer-1)))
                InTarget%trans = (0.0d0, 0.0d0)
                if (.not. allocated(InTarget%Mlayer)) allocate(InTarget%Mlayer(size(InTarget%layer), 2, 2))
                InTarget%Mlayer = (0.0d0, 0.0d0)
                if (.not. allocated(InTarget%M)) allocate(InTarget%M(2,2))
                InTarget%M = (0.0d0, 0.0d0)
                forall(j = 1:2) InTarget%M(j,j) = (1.0d0, 0.0d0) 
                if (.not. allocated(InTarget%fb)) allocate(InTarget%fb(size(InTarget%layer), 2))
                InTarget%fb = (0.0d0, 0.0d0)
            case default ! Beer-Lambert law
                call linspace(0.0d0, parameters%dlayer, Npoints, parameters%X)
        end select

        return
    end subroutine Initialize_objects 

    !!! The next three procedures are stealed from TREKIS 4 code by Nikita Medvedev !!!
    
    subroutine Count_lines_in_file(File_num, N, skip_lines)
        integer, INTENT(in) :: File_num     ! number of file to be opened
        integer, INTENT(out) :: N           ! number of lines in this file
        integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
        integer i
        if (present(skip_lines)) then	! If the user specified to start counting lines not from the first one:
            do i=1,skip_lines
                read(File_num,*, end=604) 
            enddo
            604 continue
        endif
        i = 0
        do	! count all the lines in the file from the specified one to the end:
            read(File_num,*, end=603)
            i = i + 1
        enddo
        603 continue
        rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
        N = i	! print out the number of lines from the specified one to the last one
    end subroutine Count_lines_in_file    
    
    subroutine Count_columns_in_file(File_num, N, skip_lines)
        integer, INTENT(in) :: File_num     ! number of file to be opened
        integer, INTENT(out) :: N           ! number of columns in this file
        integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
        real(8) temp
        character(1000) temp_ch
        integer i, Reason
        integer :: temp_i
        if (present(skip_lines)) then
            do i=1,skip_lines
                read(File_num,*, end=601) 
            enddo
            601 continue
        endif

        read(File_num,'(a)', IOSTAT=Reason) temp_ch ! count columns in this line
        N = number_of_columns(trim(adjustl(temp_ch))) ! see below

        rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
    end subroutine Count_columns_in_file
    
    pure function number_of_columns(line)
       integer :: number_of_columns
       character(*), intent(in) :: line
       integer i, n, toks
       logical :: same_space
       same_space = .false.
       i = 0
       n = len(line)
       number_of_columns = 0
       do while(i < n) ! scan through all the line
          i = i + 1
          selectcase (line(i:I))
          case (' ', '	') ! space or tab can be a separator between the columns
             if (.not.same_space) number_of_columns = number_of_columns + 1
             same_space = .true. ! in case columns are separated by more than one space or tab
          case default ! column data themselves, not a space inbetween
             same_space = .false.
          endselect
       enddo
       number_of_columns = number_of_columns + 1	! number of columns is by 1 more than number of spaces inbetween
    end function number_of_columns 

end module Reading