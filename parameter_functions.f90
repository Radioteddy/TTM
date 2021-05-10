module param_funcs
! this module contains functions for parameters 
! of electron and lattice subsystems
! atribute pure elemental is required for applying functions to both scalars and arrays
! functions must be not recursive. If you want use recursions, change atribute to pure recursive 

! [1] Bulgakova et al. App.Phys.A 81(2), 345-356, 2005 
    
    use Objects
    use Constants
    
    implicit none
    
    contains
    
    pure elemental function electron_conduct(Te, Tl, params) result (ke)
    ! function for electron thermal conductivity
    real(8), intent(in) :: Te ! electron temperature [K]
    real(8), intent(in) :: Tl ! lattice temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8) ke ! electron heat thermal conductivity [W/m K]
    real(8) a3
    ! here we define function
    a3 = params%kel !* Te / Tl ! for example. Or specify any shape of dependence
    ke = ( ( ( 1.0d3 * (1.0d0 + 0.4017d0*a3*Te + 1.7877d0*(a3*Te)**2.0d0 + 0.3725d0*(a3*Te)**3.0d0) & 
            / (a3*Te*(25.123d0 + 0.2524d0*a3*Te)) )**(-1.0d0) + ( 1.8d6*(a3*Te*(1.0d0 + 0.2704d0*(a3*Te)**2.0d0) & 
            / (1.0d0 + 0.1991d0*(a3*Te)**1.9371d0))/Tl )**(-1.0d0) )**(-1.0d0) )
    !
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
    real(8) Ce ! electron heat capacity [J/m^3 K]
    real(8), dimension(10), parameter :: A = (/112688.21596, -188.5018, 0.32547, -8.0733E-5, 1.0725E-8, &
                        -8.82236E-13, 4.59064E-17, -1.46521E-21, 2.61331E-26, -1.99355E-31/) ! coefficients of polynomial fitting 
    ! here we define function
    ! fitted XTANT data for Ru. Or specify any shape of dependence
    Ce = A(1) + A(2)*Te + A(3)*Te**2.0d0 + A(4)*Te**3.0d0 + A(5)*Te**4.0d0 + &
        A(6)*Te**5.0d0 + A(7)*Te**6.0d0 + A(8)*Te**7.0d0 + A(9)*Te**8.0d0 + A(10)*Te**9.0d0
    !
    end function electron_cap
    
    pure elemental function lattice_cap(Tl, params) result (Cl)
    ! function for lattice heat capacity
    real(8), intent(in) :: Tl ! lattice temperature [K]
    ! real(8), intent(in) :: dTl ! gradient of lattice temperature in free points [K]; used for estimation of gamma. NOT READY YET
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8) Cl ! Lattice heat capacity [J/m^3 K]
    real(8), parameter :: gamma = 10.0d0 ! [K]
    real(8), dimension(5), parameter :: A = (/2165402.11596, 3308.94109, -2.92427, 0.00128, -1.59295E-7/) ! coefficients of polynomial fitting
    ! here we define function
    !Cl = params%Clat
    if (Tl .le. params%Tmelt) then 
        Cl = A(1) + A(2)*Tl + A(3)*Tl**2.0d0 + A(4)*Tl**3.0d0 + A(5)*Tl**4.0d0  ! !for example. Or specify any shape of dependence
    else
        Cl = params%Cliq  !
    end if
    if (params%Tmelt .gt. 0.0d0) then
    ! here we define effective heat capacity with accounting of latent heat of melting, see [1]
    Cl = Cl + params%Hf * 1.0d0 / sqrt(2.0d0 * g_Pi) / gamma * dexp (-(Tl - params%Tmelt)**2 /2.0d0/gamma**2 )
    end if
    !
    end function lattice_cap
    
    pure elemental function coupling(Te, Tl, params) result (G)
    ! function for electron-phonon coupling
    real(8), intent(in) :: Te ! electron temperature [K]
    real(8), intent(in) :: Tl ! lattice temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8) G ! electron-phonon coupling [W/m^3 K]
    real(8) G0 ! electron-phonon coupling [W/m^3 K] at Tl = 300 K
    real(8), parameter :: alpha = 0.45 ! parameters of lattice temperature dependence
    real(8), parameter :: Tl0 = 300  !
    real(8), dimension(10), parameter :: A = (/4.1746E16, 7.1666E13, -2.47402E10, 1.03739E7, &
                                -2121.74454, 0.23247, -1.47187E-5, 5.42251E-10, -1.08188E-14, 9.04948E-20/) ! coefficients of polynomial fitting   
    ! here we define function
    ! fitted XTANT data for Ru. Or specify any shape of dependence
    G0 = A(1) + A(2)*Te + A(3)*Te**2.0d0 + A(4)*Te**3.0d0 + A(5)*Te**4.0d0 + &
        A(6)*Te**5.0d0 + A(7)*Te**6.0d0 + A(8)*Te**7.0d0 + A(9)*Te**8.0d0 + A(10)*Te**9.0d0
    G = G0 * (1 + alpha*(Tl/Tl0 - 1)) 
    !
    end function coupling    
    
end module param_funcs