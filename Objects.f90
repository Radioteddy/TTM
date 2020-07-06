module Objects
    
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
    
    type :: TTM
    ! two-temperature model parameters
        real(8), dimension(:), allocatable :: X     ! spacegrid [m]
        real(8), dimension(:), allocatable :: Tel   ! initial electron tempearture [K]
        real(8), dimension(:), allocatable :: Tlat  ! initial lattice temperature [K]
        real(8) Cel     ! electron heat capacity [J/m^3 K]
        real(8) kel     ! electron thermal conductivity [W/m K]
        real(8) klat      ! lattice thermal conductivity [W/m K]
        real(8) Clat    ! lattice heat capacity [J/m^3 K] in solid phase
        real(8) Cliq    ! lattice heat capacity [J/m^3 K] in liquid phase
        real(8) Hf    ! latent heat of melting [J/m^3]
        real(8) Tmelt    ! melting temperature [K]
        real(8) G       ! electron-phonon coupling [W/m^3 K]
        real(8) dt      ! timestep [s]
        real(8) tsave   ! time to save data [s]
        real(8) tend    ! end of calculation time [s]
        real(8), dimension(:), allocatable :: res_time      ! timegrid [s]
        real(8), dimension(:,:), allocatable :: res_Tel     ! final electron temperature (res_time,X)
        real(8), dimension(:,:), allocatable :: res_Tlat    ! final lattice temperature (res_time,X)
        real(8), dimension(:), allocatable :: res_Fabs      ! absorbed fluence (res_time)
        real(8), dimension(:), allocatable :: res_Eel       ! electron energy denisty (res_time)
        real(8), dimension(:), allocatable :: res_Elat      ! lattice energy density (res_time)
        real(8), dimension(:), allocatable :: res_dE        ! energy deviation in percents (res_time)
        character(20) mat    ! name of material
        integer flag    ! flag of space part of source. 1 is Beer-Lambert law, 2 is Transfer-Matrix method
        real(8) dlayer  ! thickness of layer [m]. Used only in the case of BL law
        real(8) l_bal   ! ballistic range [m]. Otpional, used in the case of BL law
        complex(8) n   ! complex refractive index of layer. Used only in the case of BL law
    end type TTM
    
end module Objects