Module functions

    ! --- The modules we need to use in the module.
    Use common
    Implicit None

    Contains

    ! --- Function that compute the potential U(x, t)
    !     Input : the position 'x' and the time 't'
    !     Output: the potential at position 'x' and time 't'
    Function calculate_potential(x, t) result(potential)
        Implicit None
        Real(wp)                    :: potential
        Real(wp), Intent(in)        :: x, t
        Real(wp)                    :: x_temp, t_temp
        
        ! --- Relocate to the corresponing position in the interval 
        !     [0, 1] and time in the interval [0, tau*omega]
        x_temp = x - floor(x)
        t_temp = t - floor(t/(tau*omega))*tau*omega
        
        ! --- Calculate what the potential is at that location and time
        if (dU == 0.0_wp .or. (tau /= 0.0_wp .and. t_temp >= 0.0_wp .and. t_temp < fraction_off*tau*omega)) then
            potential = 0.0_wp
        else
            if (x_temp >= 0.0_wp .and. x_temp < alpha) then
                potential = x_temp/alpha
            else if (x_temp >= alpha .and. x_temp < 1.0_wp) then
                potential = (1.0_wp-x_temp)/(1.0_wp-alpha)
            end if
        end if

    end Function calculate_potential

    ! --- Function that compute the force F(x, t) due to the 
    !     potential U(x,t) (implicit)
    !     Input : the position 'x' and the time 't'
    !     Output: the force at position 'x' and time 't'
    Function calculate_force2(x, t) result(force)
        Implicit None
        Real(wp)                    :: force
        Real(wp), Intent(in)        :: x, t
        Real(wp)                    :: dx = 1.0e-3_wp
        
        ! --- Calculate what the potential is at that location and time
        if (calculate_potential(x,t) == 0.0_wp) then
            force = 0.0_wp
        else
            force = - ( calculate_potential(x + dx, t) - calculate_potential(x - dx, t) ) / (2.0_wp*dx)
        end if

    end Function calculate_force2

    ! --- Function that compute the force F(x, t) due to the 
    !     potential U(x,t) (explicit)
    !     Input : the position 'x' and the time 't'
    !     Output: the force at position 'x' and time 't'
    Function calculate_force(x, t) result(force)
        Implicit None
        Real(wp)                    :: force
        Real(wp), Intent(in)        :: x, t
        Real(wp)                    :: x_temp, t_temp, dx = 1.0e-3_wp
    
        ! --- Relocate to the position in the corresponding interval [0, 1]
        x_temp = x - floor(x)

        ! --- Relocate and scale the time to the corresponding interval [0,1]
        t_temp = t/(tau*omega) - floor(t/(tau*omega))
        
        ! --- Calculate what the force is at that location and time
        if (dU == 0.0_wp .or. (tau /= 0.0_wp .and. t_temp >= 0.0_wp .and. t_temp < fraction_off)) then
            force = 0.0_wp
        else
            if (x_temp >= 0.0_wp .and. x_temp < alpha) then
                force = -1.0_wp/alpha
            else if (x_temp >= alpha .and. x_temp < 1.0_wp) then
                force = 1.0_wp/(1.0_wp-alpha)
            end if
        end if

    end Function calculate_force

    ! --- Function that returns the updated position (using the Euler scheme).
    !     Input : position, time, timestep
    !     Output: updated position
    Function update_position(x, t, dt) result(new_x)
        Implicit None
        Real(wp), Intent(in)    :: x, t, dt
        Real(wp)                :: new_x       

        new_x = x + calculate_force(x, t) * dt + sqrt(2*kT/dU*dt)*get_random_gauss()

    end Function update_position

    ! --- Subroutine that runs the Euler scheme for one particle.
    !     Input : none
    !     Output: none but position array is updated
    Subroutine euler_scheme(positions, dt)
        Implicit None
        Real(wp), Intent(in)    :: dt
        Real(wp)                :: positions(:)
        Integer                 :: i, j
                
            ! --- Set the starting position for the particle.
            positions(1) = x0

            ! --- Iterate through all the time steps.
            do j = 2, nSteps
                positions(j) = update_position(positions(j-1), j*dt, dt)
                !write(*,*) 't: ', j*dt/omega, 'F: ', calculate_force(positions(j-1), j*dt), 'x: ', positions(j)*L*1e6
            end do

    end Subroutine euler_scheme

    ! --- Function that checks the time criterion.
    !     Input : none
    !     Output: logical true or false
    Function check_time_criterion() result(is_time_criterion)
        Implicit None
        Logical     :: is_time_criterion
        Real(wp)    :: lhs, f1, f2

        f1 = abs(calculate_force(0.5_wp*alpha, (fraction_off + 1.0_wp)/2.0_wp*tau*omega))
        f2 = abs(calculate_force(alpha, (fraction_off + 1.0_wp)/2.0_wp*tau*omega))

        lhs = (max(f1, f2)*timestep + 4.0_wp*sqrt(2.0_wp*kT/dU)*timestep)/alpha

        write(*,*) 'Time criterion: ', lhs*100, ' % << 100 %.'

        if (lhs <= 0.1_wp) then
            is_time_criterion = .True.
        else
            is_time_criterion = .False.
        end if

    end Function check_time_criterion

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

    ! ------------------------------------------------
    ! --- STUFF THAT HAS TO DO WITH RANDOM NUMBERS --- 
    ! ------------------------------------------------
    
    ! --- Function that returns a Gaussian distributed random number
    !     with mean 0 and unit standard deviation.
    !     Input : none
    !     Output: a Gaussian random number
    Function get_random_gauss() result(random_gauss)
        Implicit None
        Real(wp)            :: random_gauss
        Real(wp)            :: rand1, rand2, w
        Integer             :: i

        ! --- Normal distribution is implemented by the polar 
        !     Box-Müller algorithm.
        w = 0.0_wp
        do while(w == 0.0_wp .OR. w >= 1.0_wp)

            ! --- Generate two random numbers from the uniform 
            !     distribution U(0,1).
            call random_number(rand1)
            call random_number(rand2)
            
            rand1 = rand1 * 2.0_wp - 1.0_wp
            rand2 = rand2 * 2.0_wp - 1.0_wp
            w = rand1**2 + rand2**2
        end do

        ! --- Calculate an uncorrelated Gaussian (or normal) 
        !     deviate of mean 0 and standard deviation 1 N(0,1).
        random_gauss = rand1 * Sqrt( - 2.0_wp * Log(w) / w)
    end Function get_random_gauss
    
    ! --- Subroutine that checks if the numbers from get_random_gauss() 
    !     are Gaussian distributed with mean 0 and unit standard deviation.
    !     Input : none
    !     Output: none
    Subroutine check_random_gauss(N)
        Implicit None
        Integer, Intent(In)     :: N
        Integer                 :: i, res
        
        ! --- Open file.
        open(unit=1,file="check_gaussian.txt", form="formatted", status="replace", action="write", iostat = res)
        
        ! --- Check if file is okay.
        if (res /= 0) then
            write(*,*) "Error in opening file, status: ", res
            stop
        end if
        
        ! --- Write to console. (?)
        write(*,*) get_random_gauss()
        write(*,*) get_random_gauss()
        
        ! --- Draw N numbers with the get_random_gauss() function and write 
        !     them to file.
        do i = 0, N
            write(1,*) get_random_gauss()
        end do
        
        ! --- Close file.
        close(unit=1)
    
    end Subroutine

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
