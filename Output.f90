module Output
    
    use Objects ! save data storred in derived types
    use add_procedures
    use IFPORT
    implicit none
    
    contains
    
    subroutine Save_absorption(Laser, InTarget, opath)
	! save absorption profile to the .dat file
    type(Multilayer), intent(in) :: InTarget    !   object with parameters of the system
    type(Source), intent(in) :: Laser    !   object with parameters of laser
    character(100), intent(in) :: opath     !   path for output data saving
    
    integer i, j, k, FN, Nloops, dd, ddd, errnum
    character(200) Filename, Filename2, command, header
    character(8) lamb, ind
    character(1) sep
    logical file_exist, file_opened, read_well
    
    call get_path_separator(sep, read_well)
    if (.not. read_well) goto 1234
    
    write(lamb, '(f6.2)') Laser%lambda
    write(Filename,'(10a)') trim(adjustl(opath)), sep, 'Absorption', sep, 'lambda_', trim(adjustl(lamb)), &
                                            '_nm_', Laser%pol, '_pol', trim(adjustl(InTarget%syscompound))
    inquire(DIRECTORY=trim(adjustl(Filename)),exist=file_exist)    ! check if input file excists
    if (.not. file_exist) then
        command='mkdir '//trim(adjustl(Filename)) ! to create a folder use this command
        i = system(trim(adjustl(command)))  ! create the folder
        if (i .eq. -1) then
            errnum = ierrno( )
            print *, 'Error ', errnum
        endif
    endif
    write(*,'(a)') '--------------------------------'
    write(*,'(a,a)') 'The outputs with absorption profile will be storred in the folder:', trim(adjustl(Filename))
    
    i = 0
    Filename2 = trim(adjustl(Filename))//sep//'absorption.dat'
    inquire(FILE=trim(adjustl(Filename2)),exist=file_exist)    ! check if input file excists
    do while (file_exist)
        i = i + 1
        write(ind,'(i8)') i
        write(Filename2,'(a,a,a,a,a)') trim(adjustl(Filename)), sep, 'absorption_', trim(adjustl(ind)), '.dat'
        inquire(FILE=trim(adjustl(Filename2)),exist=file_exist)    ! check if input file excists
    enddo
    
    FN = 200
    !Filename2 = 'test/lambda_800.00_nm_s_pol-Air-Ru-Si-Air/absorption.dat'
    open(unit=FN, file=trim(adjustl(Filename2)))
    write(FN, '(a,a,a)', advance='no') '# multilayer: ', trim(adjustl(InTarget%syscompound)), ' ('
    do i = 1, size(InTarget%layer)
        write(FN, '(a,e11.3)', advance='no') '_', InTarget%layer(i)
    enddo
    write(FN, '(a,f5.2,a)') '); theta_0 = ', Laser%in_angle, ' (deg)'
        
    header = '#depth (nm)       I_abs (1/nm)'
    write(FN, '(a)') trim(adjustl(header)) 

    Nloops = size(InTarget%Absorp)
    ddd = int(Nloops/100.0d0)
    dd = ddd  
    do i = 1, size(InTarget%Absorp)
        write(FN, '(e10.3,a,e10.3)') InTarget%spacegrid(i), '       ', InTarget%Absorp(i)
        if (i .EQ. dd) then
            call progress('Save to output:   ', dd, Nloops)
            dd = dd + ddd
        endif
    enddo
    
    1234 if (.not. read_well) print*, 'Unknown OS'  
    inquire(unit=FN,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN)             ! and if it is, close it
    
    end subroutine Save_absorption
    
    subroutine Save_TTM(parameters, laser, opath)
    ! save TTM calculation results
    type(Source), intent(in) :: laser       ! object with parameters of laser
    type(TTM), intent(in) :: parameters     ! object with TTM parameters
    character(100), intent(in) :: opath     ! path for output data saving
    ! internal variables
    integer i, j, k, FN, FNN, Nloops, dd, ddd, errnum
    character(200) Filename, Filename2, Filename3, command, header
    character(20) flu, tpu, lam, ind
    character(1) sep
    logical file_exist, file_opened, read_well    
    
    call get_path_separator(sep, read_well)
    if (.not. read_well) goto 12345    
    
    write(lam, '(f20.2)') laser%lambda
    write(flu, '(f20.2)') laser%Fluence
    write(tpu, '(f20.2)') laser%tpulse*1.d15
    write(Filename,'(12a)') trim(adjustl(opath)), sep, 'TTM', sep, trim(adjustl(parameters%mat)), sep, &
                                                trim(adjustl(flu)), '_J_m_2_', trim(adjustl(tpu)), '_fs_', trim(adjustl(lam)), '_nm'
    inquire(DIRECTORY=trim(adjustl(Filename)),exist=file_exist)    ! check if input file excists
    if (.not. file_exist) then
        command='mkdir '//trim(adjustl(Filename)) ! to create a folder use this command
        i = system(trim(adjustl(command)))  ! create the folder
        if (i .eq. -1) then
            errnum = ierrno( )
            print *, 'Error ', errnum
        endif
    endif
    write(*,'(a)') '--------------------------------'
    write(*,'(a,a)') 'The outputs with TTM profiles will be storred in the folder: ', trim(adjustl(Filename))
    
    ! check if profile.dat file exists
    i = 0
    Filename2 = trim(adjustl(Filename))//sep//'profiles.dat'
    inquire(FILE=trim(adjustl(Filename2)),exist=file_exist)    ! check if input file excists
    do while (file_exist)
        i = i + 1
        write(ind,'(i8)') i
        write(Filename2,'(a,a,a,a,a)') trim(adjustl(Filename)), sep, 'profiles_', trim(adjustl(ind)), '.dat'
        inquire(FILE=trim(adjustl(Filename2)),exist=file_exist)    ! check if input file excists
    enddo
    ! check if conservation.dat file exists
    i = 0
    Filename3 = trim(adjustl(Filename))//sep//'conservation.dat'
    inquire(FILE=trim(adjustl(Filename3)),exist=file_exist)    ! check if input file excists
    do while (file_exist)
        i = i + 1
        write(ind,'(i8)') i
        write(Filename3,'(a,a,a,a,a)') trim(adjustl(Filename)), sep, 'conservation_', trim(adjustl(ind)), '.dat'
        inquire(FILE=trim(adjustl(Filename3)),exist=file_exist)    ! check if input file excists
    enddo
    
    FN = 300
    open(unit=FN, file=trim(adjustl(Filename2)))
    header = '# time (ps)       depth (nm)       Te (K)      Tl(K)'
    write(FN, '(a)') trim(adjustl(header)) 
    
    FNN = 400
    open(unit=FNN, file=trim(adjustl(Filename3)))
    header = '# time (ps)       Fabs (J/m^2)       E_el (J/m^2)      E_lat (J/m^2)       Delta_E (%)'
    write(FNN, '(a)') trim(adjustl(header)) 
    
    Nloops = size(parameters%res_Tel)
    ddd = int(Nloops/100.0d0)
    dd = ddd  
    do i = 1, size(parameters%res_Tel, 1)
        write(FNN, '(4(f18.8,3x),f10.5)') parameters%res_time(i)*1.d12, parameters%res_Fabs(i), parameters%res_Eel(i), parameters%res_Elat(i), parameters%res_dE(i)
        do j = 1, size(parameters%res_Tel, 2)
            write(FN, '(4(e18.10,3x))') parameters%res_time(i)*1.d12, parameters%X(j)*1.d9, parameters%res_Tel(i,j), parameters%res_Tlat(i,j)
            if (int(i*j) .EQ. dd) then
                call progress('Save to output:   ', dd, Nloops)
                dd = dd + ddd
            endif
        enddo
        write(FN, '(a)') ''
    enddo
    
    
    12345 if (.not. read_well) print*, 'Unknown OS'  
    inquire(unit=FN,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN)             ! and if it is, close it
    inquire(unit=FNN,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FNN)             ! and if it is, close it    
    
    end subroutine Save_TTM
    
    
    subroutine Deallocate_all(InTarget, parameters)
        ! Deallocate all arrays after calculation
        type(Multilayer), intent(inout) :: InTarget !   object with parameters of the system
        type(TTM), intent(inout) :: parameters      ! object with TTM parameters
        ! deallocate multilayer parameters
        if (allocated(InTarget%nlayer)) deallocate(InTarget%nlayer)
        if (allocated(InTarget%layer)) deallocate(InTarget%layer)
        if (allocated(InTarget%theta)) deallocate(InTarget%theta)
        if (allocated(InTarget%phase)) deallocate(InTarget%phase)
        if (allocated(InTarget%ref)) deallocate(InTarget%ref)
        if (allocated(InTarget%trans)) deallocate(InTarget%trans)
        if (allocated(InTarget%Mlayer)) deallocate(InTarget%Mlayer)
        if (allocated(InTarget%M)) deallocate(InTarget%M) 
        if (allocated(InTarget%fb)) deallocate(InTarget%fb)
        if (allocated(InTarget%Absorp)) deallocate(InTarget%Absorp)
        if (allocated(InTarget%spacegrid)) deallocate(InTarget%spacegrid)
        ! deallocate 2T model parameters
        if (allocated(parameters%X)) deallocate(parameters%X)
        if (allocated(parameters%Tel)) deallocate(parameters%Tel)
        if (allocated(parameters%Tlat)) deallocate(parameters%Tlat)
        if (allocated(parameters%res_time)) deallocate(parameters%res_time)
        if (allocated(parameters%res_Tel)) deallocate(parameters%res_Tel)
        if (allocated(parameters%res_Tlat)) deallocate(parameters%res_Tlat)
        
    end subroutine Deallocate_all
    
end module Output