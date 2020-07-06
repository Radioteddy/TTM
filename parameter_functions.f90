module param_funcs
! this module contains functions for parameters 
! of electron and lattice subsystems
! atribute pure elemental is required for applying functions to both scalars and arrays
! functions must be not recursive. If you want use recursions, change atribute to pure recursive 
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
    ! here we define function
    ke = params%kel * Te / Tl ! for example. Or specify any shape of dependence
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
    ! here we define function
    Ce = params%Cel * Te ! for example. Or specify any shape of dependence
    !
    end function electron_cap
    
    pure elemental function lattice_cap(Tl, params) result (Cl)
    ! function for lattice heat capacity
    real(8), intent(in) :: Tl ! lattice temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8) Cl ! Lattice heat capacity [J/m^3 K]
    real(8), parameter :: gamma = 1.0d-3
    ! here we define function
    if (Tl .le. params%Tmelt) then 
        Cl = params%Clat  ! * Tl !for example. Or specify any shape of dependence
    else
        Cl = params%Cliq  !
    end if
    Cl = Cl + params%Hf/params%Tmelt * gamma/g_pi / ((Tl/params%Tmelt - 1.0d0)**2 + gamma**2 )! here we define effective heat capacity with accounting of latent heat of melting
    !
    end function lattice_cap
    
    pure elemental function coupling(Te, params) result (G)
    ! function for electron-phonon coupling
    real(8), intent(in) :: Te ! electron temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8) G ! electron heat capacity [J/m^3 K]
    ! here we define function
    G = params%G ! * Te ! for example. Or specify any shape of dependence
    !
    end function coupling    
    
end module param_funcs