module Transfer_Matrix
!   This module contains the subroutines for calculation 
!   of absorption profile into layers using transfer-matrix method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module is a copy of python code, developed in UDCM Group Stockholm University, http://udcm.fysik.su.se   !
! https://github.com/udcm-su/AbsorptionTMM                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    use Constants ! load module with universal constants
    use Objects   ! Objects with initial paramters
    use add_procedures
    
    implicit none
    
    contains
    
    subroutine Fresnel(theta_in, n_in, n_out, Laser, theta_out, r, t) 
    !   calculation of reflectivity and transmittivity per layer using Fresnel equations     
        complex(8), intent(in) :: theta_in, n_in, n_out !   incident angle of layer i, refraction coefficients of layers i and i+1 
        type(Source), intent(in) :: Laser !   object with parameters of laser
        complex(8), intent(out) :: theta_out, r, t !        incedent angle of layer i+1, reflection and transmission coefficients
        integer polarization
        
        theta_out = asin(n_in*sin(theta_in)/n_out)
        select case(Laser%pol)
            case('s')
                r = (n_in*cos(theta_in) - n_out*cos(theta_out))/(n_in*cos(theta_in) + n_out*cos(theta_out))
                t = 2*n_in*cos(theta_in)/(n_in*cos(theta_in) + n_out*cos(theta_out))
            case('p')
                r = (n_out*cos(theta_in) - n_in*cos(theta_out))/(n_out*cos(theta_in) + n_in*cos(theta_out))
                t = 2*n_in*cos(theta_in)/(n_out*cos(theta_in) + n_in*cos(theta_out))
            case default ! default is p polarization
                r = (n_out*cos(theta_in) - n_in*cos(theta_out))/(n_out*cos(theta_in) + n_in*cos(theta_out))
                t = 2*n_in*cos(theta_in)/(n_out*cos(theta_in) + n_in*cos(theta_out))
        end select
        return
    end subroutine Fresnel

    subroutine TransferMatrix(Laser, InTarget)
    !   calculation of matrix of Reflectance and transmittance
        type(Source), intent(in) :: Laser !   object with parameters of laser
        type(Multilayer), intent(inout) :: InTarget ! object with parameters of the system
        complex(8) Tlayer(2,2), Rlayer(2,2), Tlayer0(2,2)   ! matrices for intermediate calculations
!         complex(8) t, r
        integer i, k
        
        
        do i = 1, size(InTarget%nlayer)-1
            call Fresnel(InTarget%theta(i), InTArget%nlayer(i), InTArget%nlayer(i+1), Laser, InTarget%theta(i+1), InTarget%ref(i), InTarget%trans(i))
        enddo
        !   M = M0*M1*M2*M4*....
        do k = 2, size(InTArget%nlayer)-1
            InTarget%phase(k) = 2*g_Pi*InTArget%nlayer(k)*cos(InTarget%theta(k))*InTarget%layer(k) / Laser%lambda
            Tlayer(1,1) = exp((0.0d0,-1.0d0)*InTarget%phase(k))/InTarget%trans(k)
            Tlayer(1,2) = (0.0d0, 0.0d0)
            Tlayer(2,1) = (0.0d0, 0.0d0)
            Tlayer(2,2) = exp((0.0d0,1.0d0)*InTarget%phase(k))/InTarget%trans(k)
            Rlayer(1,1) = (1.0d0, 0.0d0)
            Rlayer(1,2) = InTarget%ref(k)
            Rlayer(2,1) = InTarget%ref(k)
            Rlayer(2,2) = (1.0d0, 0.0d0)            
            InTarget%Mlayer(k, :, :) = matmul(Tlayer, Rlayer)
            InTarget%M = matmul(InTarget%M, InTarget%Mlayer(k, :, :))
        enddo          
        !   compute parameters for the first interface
        Tlayer0(1,1) = 1.0d0/InTarget%trans(1)
        Tlayer0(1,2) = InTarget%ref(1)/InTarget%trans(1)
        Tlayer0(2,1) = InTarget%ref(1)/InTarget%trans(1)
        Tlayer0(2,2) = 1.0d0/InTarget%trans(1)
        InTarget%M = matmul(Tlayer0, InTarget%M)
        !   compute transmission and reflection amplitudes
        InTarget%tt = 1.0d0/InTarget%M(1,1)
        InTarget%rr = InTarget%M(2,1)/InTarget%M(1,1)
        !   relative transmission, reflection and absorption for both type of polarization
        select case(Laser%pol)
            case('s')
                InTarget%T = abs(InTarget%tt)**2 *real(InTArget%nlayer(size(InTArget%nlayer))*cos(InTarget%theta(size(InTarget%theta)))) / &
                dreal(InTArget%nlayer(1)*cos(InTarget%theta(1)))
            case('p')
                InTarget%T = abs(InTarget%tt)**2 *real(InTArget%nlayer(size(InTArget%nlayer))*cos(conjg(InTarget%theta(size(InTarget%theta))))) / &
                real(InTArget%nlayer(1)*cos(conjg(InTarget%theta(1))))
            case default ! default is p polarization
                InTarget%T = abs(InTarget%tt)**2 *real(InTArget%nlayer(size(InTArget%nlayer))*cos(conjg(InTarget%theta(size(InTarget%theta))))) / &
                real(InTArget%nlayer(1)*cos(conjg(InTarget%theta(1))))
        end select
        InTarget%R = abs(InTarget%rr)**2
        InTarget%A = 1.0d0 - InTarget%T - InTarget%R
        return
    end subroutine TransferMatrix
    

    subroutine Ampl(Laser, InTarget)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   After r & t have been calculated and all the respective matrices M_n
    !   for each layer are known, we can go 'backwards', i.e. from the last to the
    !   first layer, and compute all the amplituedes for the forward f_n and 
    !   backward b_n traveling wave. -> [f_n,b_n].T = M_n @ [f_{n+1},b_{n+1}].T
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        type(Source), intent(in) :: Laser !   object with parameters of laser
        type(Multilayer), intent(inout) :: InTarget ! object with parameters of the system
        complex(8) fb_layer(2) !    forward and backward wave amplitudes per layer
        integer i
        
        call TransferMatrix(Laser, InTarget)
        fb_layer(1) = InTarget%tt
        fb_layer(2) = 0
        InTarget%fb(size(InTarget%fb, 1), :) = fb_layer !transpose(fb_layer)
        do i = size(InTArget%nlayer)-1, 2, -1
            fb_layer = matmul(InTarget%Mlayer(i, :, :), fb_layer)
            InTarget%fb(i, :) = fb_layer !transpose(fb_layer)
        enddo
        return
    end subroutine Ampl
    
    subroutine Calculate_Absorption(Laser, InTarget)
    !   absorption profile calculation in whole system layer by layer
        type(Source), intent(in) :: Laser                   !   object with parameters of laser
        type(Multilayer), intent(inout) :: InTarget             ! object with parameters of the system
        complex(8) kz, Af, Ab                               !   wavenumber and field amplitudes for i-th layer
        complex(8), dimension(:), allocatable :: Ef         !   forward traveling field for i-th layer
        complex(8), dimension(:), allocatable :: Eb         !   backward traveling field for i-th layer
        complex(8), dimension(:), allocatable :: Absorption !   Absorption profile for i-th layer
        real(8), dimension(:), allocatable :: grid          !   Spacegrid per layer
        integer i, total, layerpoints
        real(8) total_len
        
        !   reload the forward and backward wave coefficients for every layer
        total_len = sum(InTarget%layer(2:size(InTarget%layer)-1))
        total = 0
        call Ampl(Laser, InTarget)
        do i = 2, size(InTarget%nlayer)-1
            layerpoints = nint(InTarget%points * InTarget%layer(i)/total_len)
            if(.not. allocated(grid)) allocate(grid(layerpoints))
            grid = 0.0d0
            call linspace(0.0d0, InTarget%layer(i), layerpoints, grid)
            if (.not. allocated(Ef)) allocate(Ef(size(grid)))
            if (.not. allocated(Eb)) allocate(Eb(size(grid)))
            if (.not. allocated(Absorption)) allocate(Absorption(size(grid)))
            kz = 2*g_pi*InTarget%nlayer(i)*cos(InTarget%theta(i))/Laser%lambda
            Af = InTarget%fb(i, 1)
            Ab = InTarget%fb(i, 2)
            Ef(:) = Af*exp((0.0d0,1.0d0)*kz*grid(:))
            Eb(:) = Ab*exp((0.0d0,-1.0d0)*kz*grid(:))
            select case(Laser%pol)
                case('s')
                    Absorption(:) = aimag( kz * InTArget%nlayer(i) * cos(InTarget%theta(i))) * abs(Ef(:)+Eb(:))**2 / &
                        real( InTarget%nlayer(1) * cos(InTarget%theta(1)))                     
                case('p')
                    Absorption(:) = aimag( InTArget%nlayer(i) * conjg(cos(InTarget%theta(i))) * (kz * abs(Ef(:)-Eb(:))**2 - conjg(kz) * abs(Ef(:)+Eb(:))**2) ) / &
                        real( InTarget%nlayer(1) * conjg(cos(InTarget%theta(1))) )
                case default ! default is p polarization
                     Absorption(:) = aimag( InTArget%nlayer(i) * conjg(cos(InTarget%theta(i))) * (kz * abs(Ef(:)-Eb(:))**2 - conjg(kz) * abs(Ef(:)+Eb(:))**2) ) / &
                        real( InTarget%nlayer(1) * conjg(cos(InTarget%theta(1))) ) 
            end select
            call concatenate(InTarget%Absorp, real(Absorption))            
            deallocate(Ef)
            deallocate(Eb)
            deallocate(Absorption)
            deallocate(grid)
            total = total + layerpoints
        enddo
        if (total .eq. size(InTarget%Absorp)) then
            if (.not. allocated(InTarget%spacegrid)) allocate(InTarget%spacegrid(total))
            call linspace(0.0d0, total_len, total, InTarget%spacegrid)
        else
            print*, 'Something went wrong!'
        endif
        return
    end subroutine Calculate_Absorption
            


end module Transfer_Matrix