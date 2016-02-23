Module functions

    ! --- The modules we need to use in the module.
    Use common
    Implicit None

    Contains

    
    ! --- Function that returns the ExB drift velocity.
    !     Input : Efield, Bfield
    !     Output: updated position
    Function calculate_drift_velocity(Efield, Bfield) result(v_ExB)
        Implicit None
        Real(wp), Intent(in)    :: Efield, Bfield
        Real(wp)                :: v_ExB    

        v_Exb = (Efield*Bfield) / (Bfield**2)

    end Function calculate_drift_velocity

    ! --- Subroutine that runs the Euler scheme for one particle.
    !     Input : none
    !     Output: none but position and velocity array is updated
    Subroutine euler_scheme(xs, ys, zs, us, vs, ws, dt)
        Implicit None
        Real(wp), Intent(in)    :: dt
        Real(wp)                :: xs(:), ys(:), zs(:), us(:), vs(:), ws(:)
        Integer                 :: i, j
                
            ! --- Set the starting position for the particle.
            xs(1) = x0
            ys(1) = y0
            zs(1) = z0
            
            ! --- Set the starting velocity for the particle.
            us(1) = u0
            vs(1) = v0
            ws(1) = w0

            ! --- Iterate through all the time steps.
            do j = 2, nSteps
                us(j) = us(j-1) + dt * sgn * (Efield + Bfield * vs(j-1))
                vs(j) = vs(j-1) - dt * sgn * Bfield * us(j-1)
                ws(j) = ws(j-1)
                
                xs(j) = xs(j-1) + dt * 0.5_wp*(us(j) + us(j-1))
                ys(j) = ys(j-1) + dt * 0.5_wp*(vs(j) + vs(j-1))
                zs(j) = zs(j-1) + dt * 0.5_wp*(ws(j) + ws(j-1))
            end do

    end Subroutine euler_scheme

    ! --- Subroutine that runs the mid-point scheme for one particle.
    !     Input : none
    !     Output: none but position and velocity array is updated
    Subroutine midpoint_scheme(xs, ys, zs, us, vs, ws, dt)
        Implicit None
        Real(wp), Intent(in)    :: dt
        Real(wp)                :: xs(:), ys(:), zs(:), us(:), vs(:), ws(:)
        Integer                 :: i, j
                
            ! --- Set the starting position for the particle.
            xs(1) = x0
            ys(1) = y0
            zs(1) = z0
            
            ! --- Set the starting velocity for the particle.
            us(1) = u0
            vs(1) = v0
            ws(1) = w0

            ! --- Iterate through all the time steps.
            do j = 2, nSteps
                us(j) = us(j-1) + 0.5_wp*dt * sgn * (Efield + Bfield * vs(j-1))
                vs(j) = vs(j-1) - 0.5_wp*dt * sgn * Bfield * us(j-1)
                us(j) = us(j-1) + dt * sgn * (Efield + Bfield * vs(j))
                vs(j) = vs(j-1) - dt * sgn * Bfield * us(j)
                ws(j) = ws(j-1)
                
                xs(j) = xs(j-1) + dt * 0.5_wp*(us(j) + us(j-1))
                ys(j) = ys(j-1) + dt * 0.5_wp*(vs(j) + vs(j-1))
                zs(j) = zs(j-1) + dt * 0.5_wp*(ws(j) + ws(j-1))
            end do

    end Subroutine midpoint_scheme

     ! --- Subroutine that runs the Runge-Kutta scheme for one particle.
    !     Input : none
    !     Output: none but position and velocity array is updated
    Subroutine rk4_scheme(xs, ys, zs, us, vs, ws, dt)
        Implicit None
        Real(wp), Intent(in)    :: dt
        Real(wp)                :: xs(:), ys(:), zs(:), us(:), vs(:), ws(:)
        Real(wp)                :: a1, a2, a3, a4, b1, b2, b3, b4
        Integer                 :: i, j
                
            ! --- Set the starting position for the particle.
            xs(1) = x0
            ys(1) = y0
            zs(1) = z0
            
            ! --- Set the starting velocity for the particle.
            us(1) = u0
            vs(1) = v0
            ws(1) = w0

            ! --- Iterate through all the time steps.
            do j = 2, nSteps
                a1 = + sgn * (Efield + Bfield * vs(j-1))
                b1 = - sgn * Bfield * us(j-1)
                a2 = + sgn * (Efield + Bfield * (vs(j-1) + dt/2*b1))
                b2 = - sgn * Bfield * (us(j-1) + dt/2*a1)
                a3 = + sgn * (Efield + Bfield * (vs(j-1) + dt/2*b2))
                b3 = - sgn * Bfield * (us(j-1) + dt/2*a2)
                a4 = + sgn * (Efield + Bfield * (vs(j-1) + dt*b3))
                b4 = - sgn * Bfield * (us(j-1) + dt*a3)
                
                us(j) = us(j-1) + dt/6 * (a1 + 2*a2 + 2*a3 + a4)
                vs(j) = vs(j-1) + dt/6 * (b1 + 2*b2 + 2*b3 + b4)
                ws(j) = ws(j-1)
                
                xs(j) = xs(j-1) + dt * 0.5_wp*(us(j) + us(j-1))
                ys(j) = ys(j-1) + dt * 0.5_wp*(vs(j) + vs(j-1))
                zs(j) = zs(j-1) + dt * 0.5_wp*(ws(j) + ws(j-1))
            end do

    end Subroutine rk4_scheme

   

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
