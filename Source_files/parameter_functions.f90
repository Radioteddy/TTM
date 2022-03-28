module param_funcs
! this module contains functions for parameters 
! of electron and lattice subsystems
! atribute pure elemental is required for applying functions to both scalars and arrays
! functions must be not recursive. If you want use recursions, change atribute to pure recursive 

! [1] Bulgakova et al. App.Phys.A 81(2), 345-356, 2005 
    
    !use bspline_module
    use Objects
    use Constants
    use add_procedures, only: interpolate, find_in_1D_array, find_in_Monotonous_1D_array
    use bspline_sub_module, only: db1val, get_status_message
    
    implicit none
    
    contains
    
    pure elemental function electron_conduct(Te, Tl, params) result (ke)
    ! function for electron thermal conductivity
    real(8), intent(in) :: Te ! electron temperature [K]
    real(8), intent(in) :: Tl ! lattice temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8) ke ! electron heat thermal conductivity [W/m K]
    real(8), parameter :: a3 = 0.647e-4
    ! here we define function
    !ke = params%kel
    !a3 = params%kel !* Te / Tl ! for example. Or specify any shape of dependence
    select case(params%kel_flag)
        case (0) ! constant value
            ke = params%kel
        case (1) ! linear dependence
            ke = params%kel * Te / Tl
        case (2) ! Drude model
            ke = params%D*electron_cap(Te, params)/(params%A * Te**2 + params%B * Tl)
        case (3) ! Rethfeld model NOT READY YET
            ke = params%kel 
        case (4) ! Petrov parametrization for Ru
            ke = ( ( ( 1.0d3 * (1.0d0 + 0.4017d0*a3*Te + 1.7877d0*(a3*Te)**2.0d0 + 0.3725d0*(a3*Te)**3.0d0) & 
            / (a3*Te*(25.123d0 + 0.2524d0*a3*Te)) )**(-1.0d0) + ( 1.8d6*(a3*Te*(1.0d0 + 0.2704d0*(a3*Te)**2.0d0) & 
            / (1.0d0 + 0.1991d0*(a3*Te)**1.9371d0))/Tl )**(-1.0d0) )**(-1.0d0) )
        case default ! constant value
            ke = params%kel
    end select
    end function electron_conduct

    pure elemental function lattice_conduct(Tl, params) result (kl)
    ! function for lattice thermal conductivity
    real(8), intent(in) :: Tl ! lattice temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8) kl ! electron heat thermal conductivity [W/m K]
    ! here we define function
    kl = params%klat ! * Tl ! for example. Or specify any shape of dependence
    !
    end function lattice_conduct
       
    pure elemental function electron_cap(Te, params) result (Ce)
    ! function for electron heat capacity
    real(8), intent(in) :: Te ! electron temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    integer iflag
    logical extrap
    real(8) Ce, f ! electron heat capacity [J/m^3 K]
    
    
    ! fitted Te-dependence for Ru
    real(8), dimension(10), parameter :: A = (/141530.53695, -186.1678, 0.30496, -6.83914E-5, 7.37374E-9, -3.8321E-13, &
            2.53974E-18, 7.16624E-22, -3.25301E-26, 4.51366E-31/)
    
    ! here we define function
    select case(params%cel_flag) ! what model of electron heat capacity is used
        case (0) ! linear dependence
            Ce = params%Cel * Te
        case (1) ! use B-splines for interpolation/extrapolation over tabulated data
            extrap = .true. 
            call db1val(Te,0,params%ce_spline%tx,params%ce_spline%nx,params%ce_spline%kx,params%ce_spline%bcoef,f,iflag,params%ce_spline%inbvx,extrap) 
            Ce = f
        case (2) ! fitted values, needs to be redefined and recompiled for every material
        ! fitted XTANT data for Ru. Or specify any shape of dependence
            Ce = A(1) + A(2)*Te + A(3)*Te**2.0d0 + A(4)*Te**3.0d0 + A(5)*Te**4.0d0 + &
                A(6)*Te**5.0d0 + A(7)*Te**6.0d0 + A(8)*Te**7.0d0 + A(9)*Te**8.0d0 + A(10)*Te**9.0d0
        case default
            Ce = params%Cel * Te
    end select
    !
    end function electron_cap
    
    pure elemental function lattice_cap(Tl, dTl, params) result (Cl)
    ! function for lattice heat capacity
    real(8), intent(in) :: Tl ! lattice temperature [K]
    real(8), intent(in) :: dTl ! gradient of lattice temperature in free points [K]; used for estimation of gamma. NOT READY YET
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8) Cl ! Lattice heat capacity [J/m^3 K]
    real(8) gamma ! [K]
    integer power
    real(8), dimension(6), parameter :: A = (/2053288.1421, 3886.559, -3.9569, 0.002105, -4.5995e-07, 4.0864e-11/) ! coefficients of polynomial fitting

    gamma = 1.0d0
    ! here we define function
    select case(params%cl_flag) ! what model of heat capacity will be used
        case (0) ! Constant value
            Cl = params%Clat !for example
        case (1) ! functional dependence
            Cl = a(1) + a(2)*tl + a(3)*tl**2.0d0 + a(4)*tl**3.0d0 + a(5)*tl**4.0d0 + a(6)*tl**5.0d0  !  c_ru parametrization from igor
        case default
            Cl = params%Clat
        end select
    if (params%melt_flag .ne. 0) then
        if (tl .gt. params%tmelt) cl = params%cliq ! constant value for liquid phase
    ! here we define effective heat capacity with accounting of latent heat of melting, see [1]
        !if (gamma .lt. dtl) then
        !    power = nint(log10(abs(dtl))) ! set width of delta-function no more than 3 spatial points
        !    gamma = 10**power ! as we use central derivative it will be automatically satisfied just setting gamma = 10^(order of dtl)
        !end if
        ! sign "-" before hf means that during the heating and melting temperature propagates from surface in depth, so temperature gradient is negative
        !cl = cl - sign(params%hf, dtl)  * 1.0d0 / sqrt(2.0d0 * g_pi) / gamma * dexp (-(tl - params%tmelt)**2 /2.0d0/gamma**2 )
        cl = cl + params%hf  * 1.0d0 / sqrt(2.0d0 * g_pi) / gamma * dexp (-(tl - params%tmelt)**2 /2.0d0/gamma**2 )
    end if
    
    end function lattice_cap
    
    pure elemental function coupling(Te, Tl, params) result (G)
    ! function for electron-phonon coupling
    real(8), intent(in) :: Te ! electron temperature [K]
    real(8), intent(in) :: Tl ! lattice temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8) G ! electron-phonon coupling [W/m^3 K]
    real(8) G0 ! electron-phonon coupling [W/m^3 K] at Tl = 300 K
    real(8), parameter :: alpha = 0.55 ! parameters of lattice temperature dependence
    real(8), parameter :: Tl0 = 300  !
    real(8), parameter :: a = 4.23d17, xc = 3943.2898, k = 4.696d-4  ! coefficients of sigmoidal fitting   
    ! here we define function
    ! fitted XTANT data for Ru. Or specify any shape of dependence
    G0 = a / (1.0d0 + exp(-k*(Te - xc)))
    !G = G0
    G = G0 * (1.0d0 + alpha*(Tl/Tl0 - 1.0d0)) 
    !
    ! Petrov parametrization
    !G = (18 - 12.5*Te/(50.0d3 + Te))*1.0d17
    end function coupling    
    
end module param_funcs