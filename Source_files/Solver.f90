Module Solver
    
    use Constants, only: g_Pi
    use Objects
    use Transfer_Matrix, only: Fresnel, calculate_absorption
    use add_procedures, only: progress, append, column_stack, interpolate, Find_in_1D_array
    use param_funcs
    
    implicit none
    
    contains
    
    function T_Source(t, Laser)
    ! time shape of laser source
    ! space shape is calculated according to transfer-matrix algortihm
        real(8) t            ! current timestep [s]
        type(Source), intent(in) :: Laser   ! object with parameters of laser
        real(8) T_Source                    ! fucntion itself
        real(8) sigma, tpos                 ! width of gaussian     
        
        sigma  = 4.d0*log(2.0d0) / Laser%tpulse**2
        T_Source = sqrt(sigma/g_pi) * Laser%Fluence * exp(-sigma*(t-Laser%tpeak)**2.d0)
    end function T_Source
    
    subroutine Beer_Lambert(laser, parameters, space_source)
    ! space part of source according Beer-Lambert law
    type(Source), intent(in) :: laser   ! object with parameters of laser
    type(TTM), intent(in) :: parameters ! object with TTM parameters
    real(8), dimension(:), intent(out) :: space_source    ! space part of source [1/m]
    real(8) delta, R, T
    complex(8) rr, tt, theta_t, theta_i, n_i
    
    theta_i = cmplx(laser%in_angle)   ! complex incident angle
    n_i = (1.0d0, 0.0d0) ! complex index of refraction of vacuum ambient
    call Fresnel(theta_i, n_i, parameters%n, laser, theta_t, rr, tt) ! calculate Fresnel equations
    R = abs(rr)**2
    !select case(laser%pol)
    !    case('s')
    !        T = abs(tt)**2 * real(parameters%n*cos(theta_t)) / real(n_i*cos(theta_i))
    !    case('p')
    !        T = abs(tt)**2 * real(parameters%n*cos(conjg(theta_t))) / real(n_i*cos(conjg(theta_i)))
    !    case default 
    !        T = abs(tt)**2 * real(parameters%n*cos(theta_t)) / real(n_i*cos(theta_i))
    !    end select
    delta = (laser%lambda * 1.d-9) / (4*g_pi*aimag(parameters%n)) + parameters%l_bal ! penetration depth
    space_source(:) = (1.0d0 - R) * exp(-parameters%X(:)/delta) / delta !(delta*(1 - exp(-parameters%dlayer/delta)))
    !space_source(:) =exp(-parameters%X(:)/delta) / delta !(delta*(1 - exp(-parameters%dlayer/delta)))
    return
    end subroutine Beer_Lambert
    
    subroutine Homogeneous(laser, parameters, space_source)
    ! space part of source according Beer-Lambert law
    type(Source), intent(in) :: laser   ! object with parameters of laser
    type(TTM), intent(in) :: parameters ! object with TTM parameters
    real(8), dimension(:), intent(out) :: space_source    ! space part of source [1/m]
    real(8) delta, R, T
    complex(8) rr, tt, theta_t, theta_i, n_i
    
    theta_i = cmplx(laser%in_angle)   ! complex incident angle
    n_i = (1.0d0, 0.0d0) ! complex index of refraction of vacuum ambient
    call Fresnel(theta_i, n_i, parameters%n, laser, theta_t, rr, tt) ! calculate Fresnel equations
    R = abs(rr)**2
    !select case(laser%pol)
    !    case('s')
    !        T = abs(tt)**2 * real(parameters%n*cos(theta_t)) / real(n_i*cos(theta_i))
    !    case('p')
    !        T = abs(tt)**2 * real(parameters%n*cos(conjg(theta_t))) / real(n_i*cos(conjg(theta_i)))
    !    case default 
    !        T = abs(tt)**2 * real(parameters%n*cos(theta_t)) / real(n_i*cos(theta_i))
    !    end select
    space_source(:) = (1.0d0 - R) / parameters%dlayer
    return
    end subroutine Homogeneous
    
    subroutine Explicit_lattice(dt, dz, Tlat, Clat, klat, Tel, G)
    ! Explicit scheme for lattice temperature, initial predictions
    ! index 2 in arrays corresponds to the current and next timesteps
        real(8), dimension(:,:), intent(inout) :: Tlat  ! Lattice temperature [K]
        real(8), dimension(:), intent(in) :: Clat  ! lattice heat capacity [J/(K m^3)] on current timestep
        real(8), dimension(:), intent(in) :: klat  ! lattice thermal conductvivty [W/(K m)] on current timestep
        real(8), dimension(:), intent(in) :: Tel	    ! Electron temperature  [K] on current timestep
        real(8), dimension(:), intent(in) :: G        ! electron-phonon coupling [W/(K m^3)] on current timestep
        real(8), intent(in) :: dt   ! timestep [s]
        real(8), intent(in) :: dz   ! sapcestep [m] ! currently isn't used
        integer Narr    ! size of arrays
        
        Narr = size(Tlat,2)
        ! left boundary
        Tlat(2,1) = dt/(dz*dz) * (klat(1) + klat(2))*(Tlat(1,2) - Tlat(1,1))/Clat(1) &
                    + Tel(1) * dt*G(1)/Clat(1) + Tlat(1,1)*(1.0d0 - dt*G(1)/Clat(1)) 
        ! right boundary
        Tlat(2,Narr) = -dt/(dz*dz) * (klat(Narr-1) + klat(Narr))*(Tlat(1,Narr) - Tlat(1,Narr-1))/Clat(Narr) &
                        + Tel(Narr) * dt*G(Narr)/Clat(Narr) + Tlat(1,Narr)*(1.0d0 - dt*G(Narr)/Clat(Narr)) 
        ! intermediate points
        Tlat(2,2:Narr-1) = 0.5d0*dt/Clat(2:Narr-1)/(dz*dz) * ( (klat(2:Narr-1) + klat(3:Narr))*(Tlat(1,3:Narr)-Tlat(1,2:Narr-1)) &
                                            + (klat(2:Narr-1) + klat(1:Narr-2))*(Tlat(1,2:Narr-1)-Tlat(1,1:Narr-2)) ) &
                            + Tel(2:Narr-1) * dt*G(2:Narr-1)/Clat(2:Narr-1) &
                            + Tlat(1,2:Narr-1) * (1.0d0 - dt*G(2:Narr-1)/Clat(2:Narr-1))
        return
    end subroutine Explicit_lattice
    
    subroutine corrector_lattice(dt, dz, Tlat, Clat, klat, Tel, G, dTlat)
    ! Crank-Nicolson scheme for lattice temperature, which is used as corrector
    ! index 2 in arrays corresponds to the current and next timesteps
        real(8), dimension(:,:), intent(inout) :: Tlat  ! Lattice temperature [K]
        real(8), dimension(:,:), intent(in) :: Clat  ! lattice heat capacity [J/(K m^3)] on current timestep
        real(8), dimension(:,:), intent(in) :: klat  ! lattice thermal conductvivty [W/(K m)] on current timestep
        real(8), dimension(:,:), intent(in) :: Tel	    ! Electron temperature  [K] on current timestep
        real(8), dimension(:,:), intent(in) :: G        ! electron-phonon coupling [W/(K m^3)] on current timestep
        real(8), intent(in) :: dt   ! timestep [s]
        real(8), intent(in) :: dz   ! sapcestep [m] ! currently isn't used
        real(8), dimension(:,:), intent(inout) :: dTlat  ! lattice temperature derivative
        integer Narr    ! size of arrays
        
        Narr = size(Tlat,2)
        ! left boundary
        dTlat(:,1) = (klat(:,1) + klat(:,2)) * (Tlat(:,2) - Tlat(:,1))/Clat(:,1)/(dz*dz) &
                    + (Tel(:,1) - Tlat(:,1)) * G(:,1)/Clat(:,1) 
        ! right boundary
        dTlat(:,Narr) = -(klat(:,Narr-1) + klat(:,Narr)) * (Tlat(:,Narr) - Tlat(:,Narr-1))/Clat(:,Narr)/(dz*dz) &
                        + (Tel(:,Narr) - Tlat(:,Narr)) * G(:,Narr)/Clat(:,Narr)
        ! intermediate points
        dTlat(:,2:Narr-1) = 0.5d0/Clat(:,2:Narr-1)/(dz*dz) * ( (klat(:,2:Narr-1) + klat(:,3:Narr))*(Tlat(:,3:Narr)-Tlat(:,2:Narr-1)) &
                                            - (klat(:,2:Narr-1) + klat(:,1:Narr-2))*(Tlat(:,2:Narr-1)-Tlat(:,1:Narr-2)) ) &
                            + (Tel(:,2:Narr-1) - Tlat(:,2:Narr-1)) * G(:,2:Narr-1)/Clat(:,2:Narr-1)
        ! Crank-Nicolson scheme
        Tlat(2,:) = Tlat(1,:) + 0.5d0*dt * (dTlat(1,:) + dTlat(2,:))
        !
!        deallocate(dTlat)
        return
    end subroutine corrector_lattice
    
    subroutine CN_electrons(tcurr, dt, dz, Tel, Tlat, Cel, kel, G, laser, S_source, dTel)
    ! Crank-Nikolson scheme for lattice temperature
    ! index 2 in arrays corresponds to the current and next timesteps
        real(8), dimension(:,:), intent(in) :: Tlat     ! Lattice temperature [K]
        real(8), dimension(:,:), intent(in) :: Cel   ! electron heat capacity [J/(K m^3)]
        real(8), dimension(:,:), intent(inout) :: Tel	! Electron temperature  [K]
        real(8), dimension(:,:), intent(in) :: G     ! electron-phonon coupling [W/(K m^3)]
        real(8), dimension(:,:), intent(in) :: kel   ! electron heat conductivity [W/(K m^3)]
        type(Source), intent(in) :: laser               ! object with parameters of laser
        real(8), dimension(:), intent(in) :: S_Source   ! space part of profile
        real(8), dimension(:,:), intent(inout) :: dTel  ! electron temperature derivative
        real(8), intent(in) :: tcurr    ! current timepoint
        real(8), intent(in) :: dt       ! timestep [s]
        real(8), intent(in) :: dz       ! sapcestep [m]
        real(8), dimension(:), allocatable :: an   ! coefficient for previous spacestep 
        real(8), dimension(:), allocatable :: bn   ! coefficient for current spacestep
        real(8), dimension(:), allocatable :: cn   ! coefficient for next spacestep
        real(8), dimension(:), allocatable :: rn   ! right side
        integer Narr   ! size of temporary arrays
        
        Narr = size(Tel,2)
        
        ! Time derivative of electron temperature
        ! left boundary
        dTel(:,1) = (kel(:,1) + kel(:,2)) * (Tel(:,2) - Tel(:,1))/Cel(:,1)/(dz*dz) &
                    + (Tlat(:,1) - Tel(:,1)) * G(:,1)/Cel(:,1) 
        ! right boundary
        dTel(:,Narr) = -(kel(:,Narr-1) + kel(:,Narr)) * (Tel(:,Narr) - Tel(:,Narr-1))/Cel(:,Narr)/(dz*dz) &
                        + (Tlat(:,Narr) - Tel(:,Narr)) * G(:,Narr)/Cel(:,Narr)
        ! intermediate points
        dTel(:,2:Narr-1) = 0.5d0/Cel(:,2:Narr-1)/(dz*dz) * ( (kel(:,2:Narr-1) + kel(:,3:Narr))*(Tel(:,3:Narr)-Tel(:,2:Narr-1)) &
                                            - (kel(:,2:Narr-1) + kel(:,1:Narr-2))*(Tel(:,2:Narr-1)-Tel(:,1:Narr-2)) ) &
                            + (Tlat(:,2:Narr-1) - Tel(:,2:Narr-1)) * G(:,2:Narr-1)/Cel(:,2:Narr-1)
        
        dTel(1,:) = dTel(1,:) + T_Source(tcurr, laser)*S_Source(:)/Cel(1,:)
        dTel(2,:) = dTel(2,:) + T_Source(tcurr+dt, laser)*S_Source(:)/Cel(2,:)
        ! Triagonal system solving
        ! allocate all arrays of coefficients
        allocate(an(Narr))
        allocate(bn(Narr))
        allocate(cn(Narr))
        allocate(rn(Narr))
        ! left boundary
        an(1) = 0.0d0
        bn(1) = 1.0d0 + dt/(dz*dz)*(kel(2,1)+kel(2,2))/Cel(2,1) + 0.5d0*dt*G(2,1)/Cel(2,1)
        cn(1) = - dt/(dz*dz)*(kel(2,1) + kel(2,2))/Cel(2,1) 
        rn(1) = Tel(1,1) * (1.0d0 - dt/(dz*dz)*(kel(1,1) + kel(1,2))/Cel(1,1) - 0.5d0*dt*G(1,1)/Cel(1,1)) &
                + Tel(1,2) * dt/(dz*dz)*(kel(1,1) + kel(1,2))/Cel(1,1)  &
                + Tlat(2,1) * 0.5d0*dt*G(2,1)/Cel(2,1) + Tlat(1,1) * 0.5d0*dt*G(1,1)/Cel(1,1) &
                + 0.5d0*dt*T_Source(tcurr, laser)*S_Source(1)/Cel(1,1) + 0.5d0*dt*T_Source(tcurr+dt, laser)*S_Source(1)/Cel(2,1)
        ! right boundary
        an(Narr) = - dt/(dz*dz)*(kel(2,Narr-1) + kel(2,Narr))/Cel(2,Narr)
        bn(Narr) = 1.0d0 + dt/(dz*dz)*(kel(2,Narr-1)+kel(2,Narr))/Cel(2,Narr) + 0.5d0*dt*G(2,Narr)/Cel(2,Narr)
        cn(Narr) = 0.0d0  
        rn(Narr) = Tel(1,Narr) * (1.0d0 - dt/(dz*dz)*(kel(1,Narr-1) + kel(1,Narr))/Cel(1,Narr) - 0.5d0*dt*G(1,Narr)/Cel(1,Narr)) &
                + Tel(1,Narr-1) * dt/(dz*dz)*(kel(1,Narr-1) + kel(1,Narr))/Cel(1,Narr)  &
                + Tlat(2,Narr) * 0.5d0*dt*G(2,Narr)/Cel(2,Narr) + Tlat(1,Narr) * 0.5d0*dt*G(1,Narr)/Cel(1,Narr) &
                + 0.5d0*dt*T_Source(tcurr, laser)*S_Source(Narr)/Cel(1,Narr) + 0.5d0*dt*T_Source(tcurr+dt, laser)*S_Source(Narr)/Cel(2,Narr)        
        ! define coefficients for internal gridpoints
        an(2:Narr-1) = -0.5d0*dt/(dz*dz) * 0.5d0*(kel(2,1:Narr-2)+kel(2,2:Narr-1))/Cel(2,2:Narr-1)
        bn(2:Narr-1) =  1.0d0 + 0.5d0*dt/(dz*dz) * 0.5d0*((kel(2,1:Narr-2)+kel(2,2:Narr-1))+(kel(2,2:Narr-1)+kel(2,3:Narr)))/Cel(2,2:Narr-1) &
                            + 0.5d0*dt*G(2,2:Narr-1)/Cel(2,2:Narr-1) 
        cn(2:Narr-1) = -0.5d0*dt/(dz*dz) * 0.5d0*(kel(2,2:Narr-1)+kel(2,3:Narr))/Cel(2,2:Narr-1)
        rn(2:Narr-1) = Tel(1,2:Narr-1) * ( 1.0d0 - 0.5d0*dt/(dz*dz) * 0.5d0*((kel(1,1:Narr-2)+kel(1,2:Narr-1))+(kel(1,2:Narr-1)+kel(1,3:Narr)))/Cel(1,2:Narr-1) &
                                                        - 0.5d0*dt*G(1,2:Narr-1)/Cel(1,2:Narr-1) ) &
                           + Tel(1,1:Narr-2) * 0.5d0*dt/(dz*dz) * 0.5d0*(kel(1,1:Narr-2)+kel(1,2:Narr-1))/Cel(1,2:Narr-1) &
                           + Tel(1,3:Narr) * 0.5d0*dt/(dz*dz) * 0.5d0*(kel(1,2:Narr-1)+kel(1,3:Narr))/Cel(1,2:Narr-1) &
                           + 0.5d0*dt*G(1,2:Narr-1)/Cel(1,2:Narr-1)*Tlat(1,2:Narr-1) &
                           + 0.5d0*dt*G(2,2:Narr-1)/Cel(2,2:Narr-1)*Tlat(2,2:Narr-1) &
                           + 0.5d0*dt*T_Source(tcurr, laser)*S_Source(2:Narr-1)/Cel(1,2:Narr-1) &
                           + 0.5d0*dt*T_Source(tcurr+dt, laser)*S_Source(2:Narr-1)/Cel(2,2:Narr-1)              
        call tridiag(an, bn, cn, rn, Tel(2,:))
        ! deallocate coefficients
        deallocate(an)
        deallocate(bn)
        deallocate(cn)
        deallocate(rn)
        return
    end subroutine CN_electrons
    
    subroutine time_evolution(laser, intarget, TTM_parameters)
        ! time evolution of electron and lattice temperatures within TTM approach
        type(Source), intent(in) :: laser           !   object with parameters of laser
        type(Multilayer), intent(inout) :: intarget    !   object with parameters of the system    
        type(TTM), intent(inout) :: TTM_parameters  !   object with TTM functions (Te, Tl, ke, kl, Ce, Cl, G)
         
        real(8) dz, t_curr, summ, t1, t2
        real(8), dimension(:,:), allocatable :: Tel_temp, dTel_temp, Tlat_temp, dTlat_temp          !
        real(8), dimension(:,:), allocatable :: kel_temp, Cel_temp, Clat_temp, G_temp, klat_temp    ! temporary arrays of 2T system
        real(8), dimension(:), allocatable :: dTlat_dz ! array for spatial gradient of lattice temperature
        real(8), dimension(:), allocatable :: Tel_corrector ! array for saving Tel for predictor-corrector scheme
        real(8), dimension(:), allocatable :: s_source      ! space part of source
        real(8), parameter :: eps = 1.d-4   ! predictor-corrector error
        integer corrector_step, Nstep
        integer q   ! every q steps save results
        integer dd, ddd, nloops ! for progress bar
        real(8), dimension(:), allocatable :: temp_space
        real(8) width   !   width of layers before layer used in TTM 
        integer ind ! index of layer where TTM will be calculated
        integer i, j    ! loop variables
        real(8) out_absorp  ! absorption for given X grid point
        
        real(8) Fabs, Elat, Eel, DeltaE, res ! variables for energy conservation checking
        real(8), dimension(:), allocatable :: Tint, Cint    ! temporary arrays of temperature and heat capacity
        integer, parameter :: Nint = 20     ! number of integration steps in Simpson integration routine
        
        ! prepare space part of source !
        allocate(s_source(size(TTM_parameters%X)))
        select case(TTM_parameters%flag)
            case(0) ! Homogeneous heating with energy density equal to incident
                s_source(:) = 1.0d0/TTM_parameters%dlayer
            case(1) ! Homogeneous heating accounting reflection (NOT FINISHED YET!)
                call Homogeneous(laser, TTM_parameters, s_source)
            case(2) ! Beer-Lambert law
                call Beer_Lambert(laser, TTM_parameters, s_source) !(NOT FINISHED YET!)
            case(3) ! Transfer-matrix algorithm
                ind = findloc(intarget%materials, TTM_parameters%mat, dim=1)
                width = 0.0d0
                if (ind .gt. 3) then    
                    width = sum(intarget%layer(2:ind-1))    ! here and below we don't account ambient layer
                else if (ind .gt. 2) then
                    width = intarget%layer(2)
                end if
                call calculate_absorption(laser, intarget)
                allocate(temp_space(size(intarget%spacegrid)))
                temp_space(:) = (intarget%spacegrid(:) - width) * 1.d-9 ! convert from [nm] to [m]
                do i = 1, size(TTM_parameters%X)
                    call find_in_1D_array(temp_space(:), TTM_parameters%X(i), j)
                    if (temp_space(j) .eq. TTM_parameters%X(i)) then
                        s_source(i) = intarget%absorp(j) * 1.d9 ! convert form [1/nm] to [1/m]
                    else
                        call interpolate(1, temp_space(j-1), temp_space(j), intarget%absorp(j-1), intarget%absorp(j), TTM_parameters%X(i), out_absorp)
                        s_source(i) = out_absorp * 1.d9 ! convert form [1/nm] to [1/m]
                    end if
                end do
            case default
                call Beer_Lambert(laser, TTM_parameters, s_source)
        end select
        ! end of preparation !
        q = int(TTM_parameters%tsave/TTM_parameters%dt) !
        if (q.lt.1) q = 1
        Nstep = 0   ! timestep for saving results
        
        nloops = int(TTM_parameters%tend/TTM_parameters%dt)   ! total amount of timesteps 
        ddd = int(nloops/100.0d0)   ! nloops in percents
        dd = ddd
        
        dz = TTM_parameters%X(2)-TTM_parameters%X(1)    ! space step
        
        ! allocation of all intermediate arrays
        allocate(Tel_temp(2, size(TTM_parameters%Tel)))
        Tel_temp(1,:) = TTM_parameters%Tel(:)
        Tel_temp(2,:) = 0.0d0
        allocate(Tlat_temp(2, size(TTM_parameters%Tel)))
        Tlat_temp(1,:) = TTM_parameters%Tlat(:)
        Tlat_temp(2,:) = 0.0d0
        allocate(dTlat_temp(2, size(TTM_parameters%Tel)))
        dTlat_temp(:,:) = 0.0d0
        allocate(dTel_temp(2, size(TTM_parameters%Tel)))
        dTel_temp(:,:) = 0.0d0
        allocate(dTlat_dz(size(TTM_parameters%Tel))) 
        ! later don't forget to add functions or reading from file for parameters below
        allocate(kel_temp(2, size(TTM_parameters%Tel)))
        allocate(Cel_temp(2, size(TTM_parameters%Tel)))
        allocate(Clat_temp(2, size(TTM_parameters%Tel)))
        allocate(klat_temp(2, size(TTM_parameters%Tel)))
        allocate(G_temp(2, size(TTM_parameters%Tel)))
        allocate(Tel_corrector(size(TTM_parameters%Tel)))   ! for corrections
        ! initialize internal variables
        Fabs = 0.0d0
        Elat = 0.0d0
        Eel = 0.0d0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        t_curr = 0.0d0
        corrector_step = 0  ! how many steps to achieve required precision
        call CPU_TIME(t1)	! show how long calculation is doing
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !write(*,'(a)') 't (ps)     Te_surf (K)     Tl_surf (K)     Fabs (J/m^2)    Ee (J/m^2)      El(J/m^2)     Delta_E (%)'
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Timeloop: do
        Nstep = Nstep + 1
            ! Predictor part
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Initial values for lattice temperature are calculated explicitly to use them in Crank-Nikolson scheme for electron temperature         
            ! calculate spatial gradient of Tlat
            dTlat_dz(1) = (Tlat_temp(1,2) - Tlat_temp(1,1))
            dTlat_dz(2:size(dTlat_dz)-1) = (Tlat_temp(1,3:) - Tlat_temp(1,:size(dTlat_dz)-2))
            dTlat_dz(size(dTlat_dz)) = (Tlat_temp(1,size(dTlat_dz)) - Tlat_temp(1,size(dTlat_dz)-1))
            ! calculate parameters
            Clat_temp(1,:) = lattice_cap(Tlat_temp(1,:), dTlat_dz(:), TTM_parameters)
            klat_temp(1,:) = lattice_conduct(Tlat_temp(1,:), TTM_parameters)
            G_temp(1,:) = coupling(Tel_temp(1,:), Tlat_temp(1,:), TTM_parameters)
            Cel_temp(1,:) = electron_cap(Tel_temp(1,:), TTM_parameters)
            kel_temp(1,:) = electron_conduct(Tel_temp(1,:), Tlat_temp(1,:), TTM_parameters)
            call Explicit_lattice(TTM_parameters%dt, dz, Tlat_temp, Clat_temp(1,:), klat_temp(1,:), Tel_temp(1,:), G_temp(1,:))
            !call CN_lattice(TTM_parameters%dt, dz, Tlat_temp, Clat_temp, Tel_temp, G_temp)
            !
            ! predictor parameters for electron temperature equation
            ! calculate spatial gradient of Tlat
            dTlat_dz(1) = (Tlat_temp(2,2) - Tlat_temp(2,1))
            dTlat_dz(2:size(dTlat_dz)-1) = (Tlat_temp(2,3:) - Tlat_temp(2,:size(dTlat_dz)-2))
            dTlat_dz(size(dTlat_dz)) = (Tlat_temp(2,size(dTlat_dz)) - Tlat_temp(2,size(dTlat_dz)-1))
            ! calculate parameters
            Cel_temp(2,:) = electron_cap(Tel_temp(1,:), TTM_parameters)
            kel_temp(2,:) = electron_conduct(Tel_temp(1,:), Tlat_temp(2,:), TTM_parameters)
            G_temp(2,:) = coupling(Tel_temp(1,:), Tlat_temp(1,:), TTM_parameters)
            Clat_temp(2,:) = lattice_cap(Tlat_temp(2,:), dTlat_dz(:), TTM_parameters)
            klat_temp(2,:) = lattice_conduct(Tlat_temp(2,:), TTM_parameters)
            call CN_electrons(t_curr, TTM_parameters%dt, dz, Tel_temp, Tlat_temp, Cel_temp, kel_temp, G_temp, laser, s_source, dTel_temp)
            ! End of predictor part
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Corrector part. New values for temperatures calculated semi-explicitly, and correct parameters of system ke, Ce, Cl, G
            Cel_temp(2,:) = electron_cap(Tel_temp(2,:), TTM_parameters) ! corrector parameters for electron heat capacity
            kel_temp(2,:) = electron_conduct(Tel_temp(2,:), Tlat_temp(2,:), TTM_parameters)
            G_temp(2,:) = coupling(Tel_temp(2,:), Tlat_temp(2,:), TTM_parameters)
            Tel_corrector(:) = Tel_temp(2,:)
            summ = 1.0d0
            do while (summ .ge. eps)
                call corrector_lattice(TTM_parameters%dt, dz, Tlat_temp, Clat_temp, klat_temp, Tel_temp, G_temp, dTlat_temp)
                call CN_electrons(t_curr, TTM_parameters%dt, dz, Tel_temp, Tlat_temp, Cel_temp, kel_temp, G_temp, laser, s_source, dTel_temp)
                summ = sum(Tel_temp(2,:) - Tel_corrector(:))
                corrector_step = corrector_step + 1
			    if (corrector_step .ge. 5000000) then
					    write (*,*) t_curr*1.d12, Tel_corrector(1), Tel_temp(2,1), summ, corrector_step
					    pause "Could not reach the neccessary precision, 5 000 000 corrections attempted"
                endif
                Tel_corrector(:) = Tel_temp(2,:)
                ! calculate spatial gradient of Tlat
                dTlat_dz(1) = (Tlat_temp(2,2) - Tlat_temp(2,1))
                dTlat_dz(2:size(dTlat_dz)-1) = (Tlat_temp(2,3:) - Tlat_temp(2,:size(dTlat_dz)-2))
                dTlat_dz(size(dTlat_dz)) = (Tlat_temp(2,size(dTlat_dz)) - Tlat_temp(2,size(dTlat_dz)-1))
                ! calculate parameters
                Cel_temp(2,:) = electron_cap(Tel_temp(2,:), TTM_parameters)
                Clat_temp(2,:) = lattice_cap(Tlat_temp(2,:), dTlat_dz(:), TTM_parameters)
                klat_temp(2,:) = lattice_conduct(Tlat_temp(2,:), TTM_parameters)
                G_temp(2,:) = coupling(Tel_temp(2,:), Tlat_temp(2,:), TTM_parameters)
            enddo
            ! end of corrector part
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !energy conservation checking
            Fabs = Fabs + sum(s_source(:))*dz * T_source(t_curr, laser)*TTM_parameters%dt
            Eel = Eel + sum(Cel_temp(2,:)*dTel_temp(2,:))*dz*TTM_parameters%dt
            Elat = Elat + sum(Clat_temp(2,:)*dTlat_temp(2,:))*dz*TTM_parameters%dt
            DeltaE = abs(Eel+Elat-Fabs)/Fabs
            ! new initial values for next timestep
            t_curr = t_curr + TTM_parameters%dt
            Tlat_temp(1,:) = Tlat_temp(2,:)
            Tel_temp(1,:) = Tel_temp(2,:)
            if (Nstep .EQ. dd) then
                call progress('Calculation:   ', dd, Nloops)
                !write(*,'(7(f12.4,2x))') t_curr*1.0d12, Tel_temp(2,1), Tlat_temp(2,1), Fabs, Elat, Eel, DeltaE
                dd = dd + ddd
            endif
            if (mod(Nstep, q) .eq. 0 .or. Nstep .eq. 1) then
                ! save data every tsave/dt steps 
                call append(ttm_parameters%res_fabs, fabs)
                call append(ttm_parameters%res_eel, eel)
                call append(ttm_parameters%res_elat, elat)
                call append(ttm_parameters%res_de, deltae)
                call append(TTM_parameters%res_time, t_curr)
                call column_stack(TTM_parameters%res_Tel, Tel_temp(1,:))
                call column_stack(TTM_parameters%res_Tlat, Tlat_temp(1,:))
			endif
            if (t_curr .ge.(TTM_parameters%tend)) exit Timeloop        ! checking if simulation time has ended            
        
        enddo Timeloop
        call CPU_TIME(t2)
        ttm_parameters%calctime = t2 - t1
        ! save absoprtion profile
        call move_alloc(s_source, ttm_parameters%res_abs)
        ! deallocate all temporary arrays
        deallocate(Tel_temp)
        deallocate(Tel_corrector)
        deallocate(Tlat_temp)
        deallocate(dTlat_temp)
        deallocate(dTlat_dz)
        deallocate(kel_temp)
        deallocate(Cel_temp)
        deallocate(Clat_temp)
        deallocate(G_temp)
        if (allocated(temp_space)) deallocate(temp_space)
            
        return
    end subroutine time_evolution       

    !!!!!!!!!!!!!!!!!!!!!!!!!TRIDAG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Press, W. H., Teukolsky, S. A., Flannery, B. P., & Vetterling, W. T. (1992). 
    !Numerical recipes in Fortran 77: volume 1, volume 1 of Fortran numerical recipes: the art of scientific computing. Cambridge university press.
    subroutine tridiag(a, b, c, r, T)
    !Solves tridiagonal matrix equation for T(1:n)
        real(8), dimension(:), intent(in) :: a  ! coefficients of previous spacestep
        real(8), dimension(:), intent(in) :: b  ! coefficients of current spacestep
        real(8), dimension(:), intent(in) :: c  ! coefficietns of next spacestep
        real(8), dimension(:), intent(in) :: r  ! right side of equation
        real(8), dimension(:), intent(inout) :: T ! temperature
        integer, parameter :: NMAX = 500
        integer n, j
        real(8) ki, temp(NMAX)

        n = size(T)
        if (b(1) .lt. 1.d-100) pause "tridag: equations must be rewritten"
        ki=b(1)
        T(1)=r(1)/ki
        do j=2,n ! Decomposition and forward substitution
            temp(j)=c(j-1)/ki
            ki=b(j)-a(j)*temp(j)
            if (ki.lt.1.d-100) pause "tridag: equations must be rewrittenz`"
            T(j)=(r(j)-a(j)*T(j-1))/ki
        enddo
        do j=n-1,1,-1 ! Backsubstitution
            T(j)=T(j)-temp(j+1)*T(j+1)
        enddo
        return
    end subroutine tridiag
    
    subroutine CN_lattice(dt, dz, Tlat, Clat, Tel, G)
    ! Crank-Nikolson scheme for lattice temperature
    ! This scheme supposes that lattice heat capacity is excluded from equation for lattice temperature
    ! currently useless
    ! index 2 in arrays corresponds to the current and next timesteps
        real(8), dimension(:,:), intent(inout) :: Tlat  ! Lattice temperature [K]
        real(8), dimension(:,:), intent(in) :: Clat  ! lattice heat capacity [J/(K m^3)]
        real(8), dimension(:,:), intent(in) :: Tel	    ! Electron temperature  [K]
        real(8), dimension(:,:), intent(in) :: G        ! electron-phonon coupling [W/(K m^3)]
        real(8), intent(in) :: dt   ! timestep [s]
        real(8), intent(in) :: dz   ! sapcestep [m]
        real(8), dimension(:), allocatable :: old_lat, new_lat, old_el, new_el    ! electron and lattice temperature coefficients for current and next timesteps
        
        ! allocate all arrays of coefficients
        allocate(old_lat(size(Tlat,2)))
        allocate(new_lat(size(Tlat,2)))
        allocate(old_el(size(Tlat,2)))
        allocate(new_el(size(Tlat,2)))
        ! don't forget to add diffusion part and boundaries
        ! lattice temperature on next timestep
        old_lat(:) = 1.0d0 - 0.5d0*dt*G(1,:)/Clat(1,:)
        new_lat(:) = 1.0d0 + 0.5d0*dt*G(2,:)/Clat(2,:)
        old_el(:) = 0.5d0*dt*G(1,:)/Clat(1,:)
        new_el(:) = 0.5d0*dt*G(2,:)/Clat(2,:)
        Tlat(2,:) = ( Tlat(1,:)*old_lat(:) + Tel(1,:)*old_el(:) + Tel(2,:)*new_el(:) ) / new_lat(:)
        ! deallocate coefficients
        deallocate(old_lat)
        deallocate(new_lat)
        deallocate(old_el)
        deallocate(new_el)
        return
    end subroutine CN_lattice  
    
end module Solver