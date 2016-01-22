Module functions

    ! --- The modules we need to use in the module.
    Use common
    Implicit None

    Contains

    ! --- Function that compute the potential U(x, t)
    !     Input: the position 'x' and the time 't'
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
        if (tau /= 0.0_wp .and. t_temp >= 0.0_wp .and. t_temp < fraction_off*tau*omega) then
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
    !     Input: the position 'x' and the time 't'
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
    !     Input: the position 'x' and the time 't'
    !     Output: the force at position 'x' and time 't'
    Function calculate_force(x, t) result(force)
        Implicit None
        Real(wp)                    :: force
        Real(wp), Intent(in)        :: x, t
        Real(wp)                    :: x_temp, t_temp, dx = 1.0e-3_wp
    
        ! --- Relocate to the corresponing position in the interval
        !     [0, 1] and time in the interval [0, tau*omega]
        x_temp = x - floor(x)
        t_temp = t - floor(t/(tau*omega))*tau*omega
        
        ! --- Calculate what the force is at that location and time
        if (tau /= 0.0_wp .and. t_temp >= 0.0_wp .and. t_temp < fraction_off*tau*omega) then
            force = 0.0_wp
        else
            if (x_temp >= 0.0_wp .and. x_temp < alpha) then
                force = -1.0_wp/alpha
            else if (x_temp >= alpha .and. x_temp < 1.0_wp) then
                force = 1.0_wp/(1.0_wp-alpha)
            end if
        end if

    end Function calculate_force
    
    ! --- Function that returns a Gaussian distributed random number
    !     with mean 0 and unit standard deviation.
    !     Input: none
    !     Output: a Gaussian random number
    Function get_random_gauss() result(random_gauss)
        Real(wp)            :: random_gauss
        Real(wp)            :: rand1, rand2, w
        Integer             :: i

        ! --- Normal distribution is implemented by the polar 
        !     Box-Müller algorithm.
        w = 0.0_wp
        do while(w == 0.0_wp .OR. w >= 1.0_wp)
            ! --- Assign a random number from the uniform distribution
            !     using random_number() to rand1 and rand2.
            call random_number(rand1)
            call random_number(rand2)
            rand1 = rand1 * 2.0_wp - 1.0_wp
            rand2 = rand2 * 2.0_wp - 1.0_wp
            w = rand1**2 + rand2**2
        end do
        random_gauss = rand1 * Sqrt( - 2.0_wp * Log(w) / w)
    end Function get_random_gauss
    
    ! --- Subroutine that checks if the numbers from get_random_gauss() 
    !     are Gaussian distributed with mean 0 and unit standard deviation.
    !     Input: none
    !     Output: none
    Subroutine check_random_gauss()
        Integer     :: i, res
        
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
        do i = 0, 10000
            write(1,*) get_random_gauss()
        end do
        
        ! --- Close file.
        close(unit=1)
    
    end Subroutine

    ! --- Subroutine that initializes the random seed for the 
    !     random number generator random_number()
    !     Input: none
    !     Output: none but the random seed is updated
    Subroutine initialize_random_seed()
        integer                             :: i, n, clock
        integer, dimension(:), allocatable  :: seed
            
        ! --- (?)
        call random_seed(size = n)
        allocate(seed(n))
        
        ! --- (?)
        call system_clock(count=clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)
        
        ! --- (?)
        deallocate(seed)

    end Subroutine

end Module functions
