module param_funcs
! this module contains functions for parameters 
! of electron and lattice subsystems
    use Objects
    use Constants
    
    implicit none
    
    contains
    
    function electron_conduct(Te, Tl, params) result (ke)
    ! function for electron thermal conductivity
    real(8), dimension(:), intent(in) :: Te ! electron temperature [K]
    real(8), dimension(:), intent(in) :: Tl ! lattice temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8), dimension(size(Te)) :: ke ! electron heat thermal conductivity [W/m K]
    ! here we define function
    ke(:) = params%kel ! * Te(:) * Tl(:) ! for example. Or specify any shape of dependence
    !
    end function electron_conduct

    function lattice_conduct(Tl, params) result (kl)
    ! function for lattice thermal conductivity
    real(8), dimension(:), intent(in) :: Tl ! lattice temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8), dimension(size(Tl)) :: kl ! electron heat thermal conductivity [W/m K]
    ! here we define function
    kl(:) = params%klat ! * Tl(:) ! for example. Or specify any shape of dependence
    !
    end function lattice_conduct
       
    function electron_cap(Te, params) result (Ce)
    ! function for electron heat capacity
    real(8), dimension(:), intent(in) :: Te ! electron temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8), dimension(size(Te)) :: Ce ! electron heat capacity [J/m^3 K]
    ! here we define function
    Ce(:) = params%Cel * Te(:) ! for example. Or specify any shape of dependence
    !
    end function electron_cap
    
    function lattice_cap(Tl, params) result (Cl)
    ! function for lattice heat capacity
    real(8), dimension(:), intent(in) :: Tl ! lattice temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8), dimension(size(Tl)) :: Cl ! Lattice heat capacity [J/m^3 K]
    ! here we define function
    Cl(:) = params%Clat ! * Tl(:) !for example. Or specify any shape of dependence
    !
    end function lattice_cap
    
    function coupling(Te, params) result (G)
    ! function for electron-phonon coupling
    real(8), dimension(:), intent(in) :: Te ! electron temperature [K]
    type(TTM), intent(in) :: params ! object with TTM parameters
    real(8), dimension(size(Te)) :: G ! electron heat capacity [J/m^3 K]
    ! here we define function
    G(:) = params%G ! * Te(:) ! for example. Or specify any shape of dependence
    !
    end function coupling    
    
end module param_funcs