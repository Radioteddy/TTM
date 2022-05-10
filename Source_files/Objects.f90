module Objects
    
    !use bspline_module
    
    implicit none
    
    type :: Multilayer 
    ! Multilayer system
        character(100) syscompound
        character(20), dimension(:), allocatable :: materials ! array of material names
        real(8), dimension(:), allocatable :: layer     ! array of layer widths
        real(8), dimension(:), allocatable :: spacegrid ! array of points in WHOLE multilayer system
        real(8), dimension(:), allocatable :: Absorp ! array of points in WHOLE multilayer system
        complex(8), dimension(:), allocatable :: nlayer ! array of complex indices of refraction
        complex(8), dimension(:), allocatable :: phase  ! phase of filed in layer
        complex(8), dimension(:), allocatable :: theta  ! incident angle in layer
        complex(8), dimension(:), allocatable :: ref    ! reflectance in layer 
        complex(8), dimension(:), allocatable :: trans  ! transmittance in layer
        complex(8), dimension(:,:), allocatable :: M    ! TMM matrix
        complex(8), dimension(:,:,:), allocatable :: Mlayer ! matrix of layer in TMM
        complex(8), dimension(:,:), allocatable :: fb   ! vector of wave amplitudes
        complex(8) tt ! total transmittance
        complex(8) rr ! total reflectance
        real(8) T ! Transmittivity ~|tt|^2
        real(8) R ! Reflectivity ~|rr|^2
        real(8) A ! Absorption
        integer points ! grid points for absorption calculation
    end type Multilayer
    
    type :: Source
    ! Laser source
        character(1) pol ! type of polarisation (s or p)
        real(8) in_angle ! angle of incidence [deg]
        real(8) lambda   ! wavelength of laser [nm]
        real(8) tpulse	 ! pulse duration [s]
        real(8) tpeak    ! center of peak time [s]
        real(8) Fluence	 ! incident fluence [J/m^2]
    end type Source

    type :: spln
        integer nx       !! the number of interpolation points in \(x\).
        integer kx       !! order of polynomial pieces in \(x\).
        integer inbvx    !! initialization parameter which must be set to 1 the first time this routine is called,
        real(8), dimension(:), allocatable :: tx    !! sequence of knots defining the piecewise polynomial in the \(x\) direction.
        real(8), dimension(:), allocatable :: bcoef !! the b-spline coefficients computed by [[db1ink]].       
    end type spln
        
    type :: TTM
    ! two-temperature model parameters
        ! initial parameters
        character(20) mat    ! name of material
        integer flag    ! flag of space part of source. 0 is homogeneous heating, 1 is Beer-Lambert law, 2 is Transfer-Matrix method
        real(8), dimension(:), allocatable :: X     ! spacegrid [m]
        real(8), dimension(:), allocatable :: Tel   ! initial electron tempearture [K]
        real(8), dimension(:), allocatable :: Tlat  ! initial lattice temperature [K]
        real(8) dlayer  ! thickness of layer [m]. Read only for Iso and BL law. In the case of Transfer-Matrix filled by value from multilayer.
        real(8) l_bal   ! ballistic range [m]. Otpional, used in the case of BL law
        complex(8) n   ! complex refractive index of layer. Used only in the case of BL law
        ! simulation parameteres
        real(8) dt      ! timestep [s]
        real(8) tsave   ! time to save data [s]
        real(8) tend    ! end of calculation time [s]      
        ! electron subsystem parameters
        integer cel_flag ! flag of which model for electron capacity is used
        real(8) Cel     ! electron heat capacity constant [J/m^3]: Cel -> Cel*Te
        type(spln) :: ce_spline ! parameters of B-spline interpolation of tabulated data
        integer kel_flag ! flag of which model for electron conductivity is used
        real(8) kel     ! electron thermal conductivity constant: [W/m K]
        real(8) A, B, D    ! constants for thermal conductivity in Drude model: ke = D*Ce(Te)/(A*Te^2 + B*Tl)
        ! coupling parameters
        integer g_flag ! flag of which model for e-ph coupling is used. 0 = constant G, 1 = G(Te), 2 = G(Te,Ti)
        real(8) G       ! electron-phonon coupling constant [W/m^3 K]
        type(spln) :: g_spline ! parameters of B-spline interpolation of tabulated data
        real(8) g_alpha ! linear coefficient for G(Ti) dependence
        ! lattice subsystem parameters
        real(8) klat      ! lattice thermal conductivity [W/m K]
        integer cl_flag ! flag of which model for lattice capacity is used
        integer melt_flag ! flag if melting is turned on
        real(8) Clat    ! lattice heat capacity [J/m^3 K] in solid phase
        real(8) Cliq    ! lattice heat capacity [J/m^3 K] in liquid phase
        real(8) Hf    ! latent heat of melting [J/m^3]
        real(8) Tmelt    ! melting temperature [K]
        ! results of calculation
        real(8), dimension(:), allocatable :: res_time      ! timegrid [s]
        real(8), dimension(:), allocatable :: res_abs       ! absorption profile [1/nm]
        real(8), dimension(:,:), allocatable :: res_Tel     ! final electron temperature (res_time,X)
        real(8), dimension(:,:), allocatable :: res_Tlat    ! final lattice temperature (res_time,X)
        real(8), dimension(:), allocatable :: res_Fabs      ! absorbed fluence (res_time)
        real(8), dimension(:), allocatable :: res_Eel       ! electron energy denisty (res_time)
        real(8), dimension(:), allocatable :: res_Elat      ! lattice energy density (res_time)
        real(8), dimension(:), allocatable :: res_dE        ! energy deviation in percents (res_time)
        real(8) calctime ! time of TTM evolution execution
    end type TTM
     
    !
    !type :: lattice_cap
    !    integer model ! flag of model is used: 0=constant value, 1=analytical function (e.g. linear), 2=tabulated form
    !    integer melting ! flag of enthalpy of melting account: 0=no, 1=yes
    !    real(8) Clat    ! lattice heat capacity [J/m^3 K] in solid phase, constant value
    !    real(8) Tmelt    ! melting temperature [K]
    !    real(8) Hf    ! latent heat of melting [J/m^3] 
    !    real(8) Cliq    ! lattice heat capacity [J/m^3 K] in liquid phase, constant value with melting=yes
    !    character(200) Clat_filename ! filename where heat capacity is storred. If melting=yes, for Cliq will be considered as values T > Tm
    !end type lattice_cap
    !
    !type :: electron_conduct
    !    integer model ! flag of model is used: 0=constant value, 1=Rethfeld-1, 2=Rehfeld-2, 3=tabulated form
    !end type electron_conduct
    
end module Objects