Module functions

    use common
    use iterativeMethods
    use integration

    implicit None

    Contains
    
    ! --- Function that computes the cross product 
    !     multiplication of two vectors. 
    !     Input : a, b (must be of rank one and of size 3)
    !     Output: cross (of dimension 3)
    Function cross_product(a, b)
        Real(wp), Dimension(3), Intent(in)  :: a, b
        Real(wp), Dimension(3)              :: cross_product

        cross_product(1) = a(2) * b(3) - a(3) * b(2)
        cross_product(2) = a(3) * b(1) - a(1) * b(3)
        cross_product(3) = a(1) * b(2) - a(2) * b(1)

    end Function cross_product
    
    ! --- Function that returns the magnetic field at (x,y,z).
    !     Input : position (x,y,z)
    !     Output: magnetic field strength (Bx, By, Bz)
    Function get_Bfield(xyz) result (Bfield)
        Real(wp), Intent(in)    :: xyz(:)
        Real(wp), Dimension(3)  :: Bfield
        Real(wp)                :: rPos, zPos, theta, radialFieldStrength, axialFieldStrength, orthoradialFieldStrength
        
        if (type_of_field == 'constant') then
            Bfield = [ 0.0_wp, 0.0_wp, Bz ]

        else if (type_of_field == 'gradient') then
            Bfield = [ 0.0_wp, 0.0_wp , Bz + 0.1_wp * xyz(2) ]

        else if (type_of_field == 'curved') then
            if (xyz(1) > 0) then
                theta = atan(xyz(2)/xyz(1))
            else
                theta = atan(xyz(2)/xyz(1)) + pi
            end if
            orthoradialFieldStrength = 1.0_wp
            Bfield = [ -orthoradialFieldStrength*sin(theta), orthoradialFieldStrength*cos(theta), 0.0_wp ]

        else if (type_of_field == 'helmholtz') then
            rPos = sqrt(xyz(1)**2 + xyz(2)**2)
            zPos = xyz(3)
            radialFieldStrength = integrate(r_integrand, 0._wp, 2._wp*pi, 100, method='simpson')
            axialFieldStrength = integrate(z_integrand, 0._wp, 2._wp*pi, 100, method='simpson')
            Bfield = [ radialFieldStrength*cos(theta), radialFieldStrength*sin(theta), axialFieldStrength ]

        end if

    end Function get_Bfield

    function r_integrand(theta)
        real(wp)             :: r_integrand
        real(wp), intent(in) :: theta
        real(wp)             :: cTheta, sTheta
        
        cTheta = cos(theta)
        sTheta = sin(theta)

        r_integrand = (zPos-1._wp) * ( (rPos - coilRadius*cTheta)**2 + (coilRadius*sTheta)**2+(zPos-1._wp)**2)**(-1.5_wp) 
        r_integrand = r_integrand+((zPos+1._wp)*((rPos - coilRadius*cTheta)**2+(coilRadius*sTheta)**2+(zPos+1._wp)**2)**(-1.5_wp))
        r_integrand = r_integrand * (1._wp+coilRadius**2)**1.5_wp / (4._wp*pi*coilRadius) * cTheta
    
    end function r_integrand

    function z_integrand(theta)
        real(wp)             :: z_integrand
        real(wp), intent(in) :: theta
        real(wp)             :: cTheta, sTheta
        
        cTheta = cos(theta)
        sTheta = sin(theta)
 
        z_integrand = ((rPos - coilRadius*cTheta)**2 + (coilRadius*sTheta)**2 + (zPos-1.0_wp)**2)**(-1.5_wp) 
        z_integrand = z_integrand + ((rPos - coilRadius*cTheta)**2 + (coilRadius*sTheta)**2 + (zPos+1.0_wp)**2)**(-1.5_wp)
        z_integrand = z_integrand * (1._wp+coilRadius**2)**1.5_wp / (4._wp*pi) * (1._wp - cTheta * rPos / coilRadius)
    
    end function z_integrand

    ! --- Function that returns the gradient of the magnetic field at (x,y,z).
    !     Input : position (x,y,z)
    !     Output: gradient of magnetic field strength grad( [Bx, By, Bz] )
    Function get_grad_Bfield(xyz) result (Bfield)
        Real(wp), Intent(in)    :: xyz(:)
        Real(wp), Dimension(3)  :: Bfield
        
        Bfield = [ 0.0_wp, 0.0_wp , 0.0_wp ]

    end Function get_grad_Bfield
    
    function dydt(pos_vel, unused) result(derivative)
        Real(wp), Intent(in)    :: pos_vel(:), unused
        Real(wp), dimension(6)  :: derivative
        
        ! --- Evaluate the slope of the position.
        derivative(1:3) = pos_vel(4:6)

        ! --- Evaluate the slope of the velocity.
        derivative(4:6) = sgn * ( Efield + cross_product(pos_vel(4:6), get_Bfield(pos_vel(1:3))) )

    end function dydt

    ! --- Subroutine that runs the mid-point scheme for one particle.
    !     Input : none
    !     Output: none but position and velocity array is updated
    Subroutine iterate(pos_vel, dt, method)
        Implicit None
        Real(wp), Intent(in)    :: dt
        Real(wp)                :: pos_vel(:,:)
        Integer                 :: i, j
        character(len=*), intent(in), optional :: method
         
        ! --- Set the starting position for the particle.
        pos_vel(1,1:3) = [x0, y0, z0]
        
        ! --- Set the starting velocity for the particle.
        pos_vel(1,4:6) = [u0, v0, w0]

        ! --- Iterate through all the time steps using your method of choice.
        if ( present(method) ) then
            select case (method)
            case ('euler')
                write(*,*) 'Iterate through time using the Euler scheme.'
                do j = 2, nSteps
                    pos_vel(j,:) = forward_euler(dydt, pos_vel(j-1,:), j*dt, dt)
                end do
            case ('midpoint')
                write(*,*) 'Iterate through time using the midpoint scheme.'
                do j = 2, nSteps
                    pos_vel(j,:) = midpoint(dydt, pos_vel(j-1,:), j*dt, dt)
                end do
            case default
                write(*,*) 'Iterate through time using the 4th order Runge-Kutta scheme.'
                do j = 2, nSteps
                    pos_vel(j,:) = rk4(dydt, pos_vel(j-1,:), j*dt, dt)
                end do
            end select
        else
            do j = 2, nSteps
                pos_vel(j,:) = rk4(dydt, pos_vel(j-1,:), j*dt, dt)
            end do
        end if

    end Subroutine iterate

    ! ----------------------------------------------------
    ! --- STUFF THAT HAS TO DO WITH THE EXACT SOLUTION --- 
    ! ----------------------------------------------------   

    ! --- Function that returns the ExB drift velocity.
    !     Input : Ey, Bz
    !     Output: updated position
    Function calculate_perp_drift(E, B) result(vd)
        Implicit None
        Real(wp), Dimension(3), Intent(in)  :: E, B
        Real(wp), Dimension(3)              :: vd

        vd = 0.0_wp !cross_product(E, B) / norm2(B)**2
        

    end Function calculate_perp_drift

    ! --- Function that returns the magnetic curvature drift velocity.
    !     Input : Ey, Bz
    !     Output: updated position
    Function calculate_curv_drift(B, divB) result(vd)
        Implicit None
        Real(wp), Dimension(3), Intent(in)  :: B, divB
        Real(wp), Dimension(3)              :: vd

        vd = 0.0_wp!sgn * 0.5_wp * v_perp**2 * cross_product(B, divB) / norm2(B)**2

    end Function calculate_curv_drift

    Function calculate_exact_trajectory() result(pos)
        real(wp), allocatable   :: pos(:,:)
        real(wp)                :: t
        real(wp), dimension(3)  :: vd_curv, vd_perp
        integer                 :: n

        allocate(pos(Nsteps/dSteps, 6))

        vd_perp = calculate_perp_drift(Efield, Bfield)

        vd_curv = calculate_curv_drift(Bfield, get_grad_Bfield([0.0_wp,0.0_wp,0.0_wp]))

        write(*,*) 'ExB drift: ', vd_perp, 'Curvature drift: ', vd_curv
        
        ! --- only accounts for vd_perp and vd_curv in x-direction
        do n = 1, Nsteps/dSteps
            t = n*dSteps*timestep
            pos(n,1) = (vd_perp(1) + vd_curv(1))*t + (1.0_wp - vd_perp(1) - vd_curv(1)) * sin(t) + x0
            pos(n,2) = sgn * (1.0_wp - vd_perp(1) - vd_curv(1)) * ( cos(t) - 1 ) + y0
            pos(n,3) = v_para/v_perp * t + z0

            pos(n, 4) = 0.0_wp
            pos(n, 5) = 0.0_wp
            pos(n, 6) = 0.0_wp
        end do

    end Function calculate_exact_trajectory

    ! -------------------------------------------------
    ! --- STUFF THAT HAS TO DO WITH WRITING TO FILE --- 
    ! -------------------------------------------------
    
    ! --- Subroutine that append a vector to a file.
    !     If the file does not exist it will be made.
    !     Input : vector and filename
    !     Output: none but the file is updated
    Subroutine append_vector_to_file(vector, filename)
        Implicit None
        Real(wp), Intent(in)            :: vector(:)
        Character(len=*), Intent(in)    :: filename
        Logical                         :: exist   

        ! --- Open output file named filename.
        inquire(file=filename, exist=exist)
        if (exist) then
            open(12, file=filename, status="old", position="append", action="write")
        else
            open(12, file=filename, status="new", action="write")
        end if
        
        ! --- Append vector to file.
        write(12,*) vector
        
        ! --- Close output file.
        close(unit=12)

    end Subroutine append_vector_to_file 

    ! --- Subroutine that write vectors to a file.
    !     Input : vectors and filename
    !     Output: none but the file is updated
    Subroutine write_to_file(time, matrix, filename)
        Implicit None
        Real(wp), Intent(in)            :: time(:), matrix(:,:)
        Character(len=*), Intent(in)    :: filename
        
        ! --- Open output file named <filename>.
        open(unit=1, file=filename, form="formatted", status="replace", action="write")
        
        ! --- Write positions and velocities to file.
        write(1,*) time(:)
        write(1,*) matrix(:, 1)
        write(1,*) matrix(:, 2)
        write(1,*) matrix(:, 3)
        write(1,*) matrix(:, 4)
        write(1,*) matrix(:, 5)
        write(1,*) matrix(:, 6)
        
        ! --- Close output file.
        close(unit=1)

    end Subroutine write_to_file  

        

    ! ------------------------------------------------
    ! --- STUFF THAT HAS TO DO WITH RANDOM NUMBERS --- 
    ! ------------------------------------------------

    ! --- Subroutine that initializes the random seed for the 
    !     random number generator random_number()
    !     Input : none
    !     Output: none but the random seed is updated
    Subroutine initialize_random_seed()
        Implicit None
        integer                             :: i, n, clock
        integer, dimension(:), allocatable  :: seed
            
        ! --- Allocate seed array.
        call random_seed(size = n)
        Allocate(seed(n))
        
        ! --- Set each element of the seed array 
        !     to a function of the current time.
        call system_clock(count=clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)
        
        ! --- Dealloacate seed array.
        Deallocate(seed)

    end Subroutine

end Module functions
