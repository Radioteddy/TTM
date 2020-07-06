module Reading
    
    use Objects
    use Constants
    use add_procedures
    
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
        READ(FN,*,IOSTAT=Reason) npoints
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 1111
        

        1111 if (.not. read_well) print*, 'Error in '//trim(adjustl(filename))//' file!' 
        inquire(unit=FN, opened=file_opened)    ! check if this file is opened
        if (file_opened) close(FN)             ! and if it is, close it
    end subroutine Read_layer_file

    subroutine read_laser_file(laser)
    ! read parameters of laser pulse
        type(Source), intent(out) :: laser  ! object with parameters of laser
        ! internal variables
        integer FN, Reason, i, j
        character(100) temp
        character(100) filename ! file with input parameters
        character(100) Error_descript  ! to write a description of an error, if any
        logical file_exist    ! to check where file to be open exists
        logical file_opened   ! to check if a file is still opened
        logical read_well     ! to check if we read the file without errors
        
        filename = 'laser_parameters.txt'
        FN = 300
        inquire(file=trim(adjustl(filename)), exist=file_exist) ! check if file exists
        if (file_exist) then
            open(unit=FN, file=trim(adjustl(filename)), status='old', action='read') ! if yes, open and read
        else ! error message about it:
            Error_descript = 'File '//trim(adjustl(filename))//' is not found!' ! if no, print an error
            print*, trim(adjustl(Error_descript))
        endif

        i = 0
        READ(FN,*,IOSTAT=Reason) laser%lambda 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 7777   
        
        READ(FN,*,IOSTAT=Reason) laser%in_angle 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 7777   

        READ(FN,*,IOSTAT=Reason) temp   
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 7777
        laser%pol = trim(adjustl(temp))   
        
        READ(FN,*,IOSTAT=Reason) laser%Fluence
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 7777   
        
        READ(FN,*,IOSTAT=Reason) laser%tpulse 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 7777

        READ(FN,*,IOSTAT=Reason) laser%tpeak 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 7777

        write(*,'(a,f8.3,a,f6.3,a,a,a)') 'Laser paramteres are: lambda = ', laser%lambda, '; theta_0 = ', laser%in_angle, '; polarisation is ', laser%pol, ' type'
        write(*, '(a,f10.3,a,e15.5,a,e15.5)') 'incident fluence F = ', laser%fluence, '; FWHM tpulse = ', laser%tpulse, '; center of peak tpeak = ', laser%tpeak
7777    if (.not. read_well) print*, 'Error in '//trim(adjustl(filename))//' file!' 
        inquire(unit=FN, opened=file_opened)    ! check if this file is opened
        if (file_opened) close(FN)             ! and if it is, close it
        
    end subroutine read_laser_file
    
    subroutine read_TTM_file(parameters, Npoints)
    ! read parameters for TTM calculation
        type(TTM), intent(out) :: parameters ! object with parameters of TTM
        integer, intent(out) :: Npoints
        ! internal variables
        character(100) matname
        real(8) Te, Tl, Ce, Cl, ke, kl, Gel, dtt, tfin, dl, lb, ts, Cli, Hfus, Tm
        complex(8) nl
        integer FN, Reason, i, j, fl
        character(100) temp
        character(100) filename ! file with input parameters
        character(100) Error_descript  ! to write a description of an error, if any
        logical file_exist    ! to check where file to be open exists
        logical file_opened   ! to check if a file is still opened
        logical read_well     ! to check if we read the file without errors
        
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
        READ(FN,*,IOSTAT=Reason) fl 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) nl, dl, lb
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222 
        READ(FN,*,IOSTAT=Reason) Te 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) Ce 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) ke 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) Tl 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) Cl, Cli 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) kl 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) Gel 
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2222
        READ(FN,*,IOSTAT=Reason) Hfus, Tm 
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
        write(*,'(a,f8.3,a,f8.3)') 'TTM parameters are: Te = ', Te, '; Tlat = ', Tl
        write(*, '(a, 6e12.3)') 'Ce, Cl, Cliq, ke, kl, G: ', Ce, Cl, Cli, ke, kl, Gel
        write(*, '(a, 2e12.3)') 'melting parameters Tmelt, Hf: ', Tm, Hfus        
        write(*, '(a, 3e12.3)') 'dt, tsave, tend: ', dtt, ts, tfin
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
        parameters%dlayer = dl*1.d-9 ! convert from nm to m
        parameters%l_bal = lb*1.d-9  ! convert from nm to m
        parameters%n = nl
        if (.not. allocated(parameters%Tel)) allocate(parameters%Tel(Npoints))
        parameters%Tel(:) = Te
        if (.not. allocated(parameters%Tlat)) allocate(parameters%Tlat(Npoints))
        parameters%Tlat(:) = Tl
        ! check if all is read well
2222    if (.not. read_well) print*, 'Error in '//trim(adjustl(filename))//' file!' 
        inquire(unit=FN, opened=file_opened)    ! check if this file is opened
        if (file_opened) close(FN)             ! and if it is, close it
    end subroutine Read_TTM_file
        
    
    
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
        
        ! Read parameters of laser
        call read_laser_file(Laser)
        ! Read initial parameters of electron and lattice systems from file
        call read_TTM_file(parameters, Npoints)
        ! create spacegrid for TTM
        if (.not. allocated(parameters%X)) allocate(parameters%X(Npoints))
        select case(parameters%flag)
            case(1) ! Beer-Lambert law
                call linspace(0.0d0, parameters%dlayer, Npoints, parameters%X)
            case(2) ! Trasfer-matrix
                ! Read initial parameters of multilayer from file
                call Read_layer_file(InTarget%materials, InTarget%layer, Intarget%nlayer, InTarget%syscompound, InTarget%points)
                ! create spacegrid
                ind = findloc(intarget%materials, parameters%mat, dim=1)
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
        
end module Reading