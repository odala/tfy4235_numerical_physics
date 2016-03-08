program main

    ! --- The modules we need to use in the program.
    use common
    use readinput
    use functions
    use integration
    use functionHolder

    implicit none

    ! --- Define variables.
    integer                 :: i, j
    real(wp)                :: start, middle1, middle2, finish
    real(wp), allocatable   :: pos_vel(:,:), pos_vel_exact(:,:)

    ! --- Define variables for checking integration schemes.
    real(wp)    :: correct_answer, startPoint, endPoint
    integer     :: Nin
    
    ! --- Define name of files.
    character(len=*), parameter :: input_file = 'input/input.txt'
    character(len=*), parameter :: e_file     = 'output/e_output.txt'
    character(len=*), parameter :: mp_file    = 'output/mp_output.txt'
    character(len=*), parameter :: rk4_file   = 'output/rk4_output.txt'
    character(len=*), parameter :: exact_file = 'output/exact_output.txt'
    
    ! --- Get initial time and say hello.
    write(*,*) 'Program starting!'
    write(*,*)
    call cpu_time(start)
    
    ! --- Read input file.
    !     Assign values to all common parameters as given 
    !     in the input file.
    call read_input_file(input_file)

    ! --- Write to command line a few words to check
    !     the parameters have read correctly.
    write(*,*) 'Check parameters.'
    write(*,'(A15, I14)') 'Nparticles: '  , Nparticles
    write(*,'(A15, I14)') 'Nsteps: '      , Nsteps
    write(*,'(A15, E14.4)') 'timestep (s): ', timestep
    write(*,'(A15, E14.4, A2, E14.4, A2, E14.4, A2)') '(x0, y0, z0) (m): ( ', x0, ', ', y0, ', ', z0, ' )'
    write(*,'(A15, E14.4, A2, E14.4, A2, E14.4, A2)') '(v0, u0, w0) (m/s): ( ', u0, ', ', v0, ', ', w0, ' )'
    !write(*,format_1)  'Mass: (kg)'   , mass
    write(*,'(A15, I14)') 'dSteps: '      , dSteps
    write(*,*) 'B = ', Bfield 
    write(*,*) 'E = ', Efield
    write(*,*)    

    ! --- Introduce reduced units.
    !x0 = x0/r !y0 = y0/r !z0 = z0/r!u0 = u0/vPerp!v0 = v0/vPerp!w0 = w0/!vPerpvPerp = vPerp/vPerp !vPara = vPara/vPerp !Bz = Bz*e/(mass*omega) !Ey = Ey*e/(mass*r*omega**2)
    
    ! --- Allocate space for the position- and velocity arrays.
    allocate(pos_vel(Nsteps, 6), pos_vel_exact(Nsteps/dSteps, 6))
    
    ! --- Initialize random number generator.
    !call initialize_random_seed()

    ! --- Initialize the magnetic and electric fields.
    write(*,*) 'Type of field: ', type_of_field
    
    ! --- Make the particles evolve during Nsteps.
    call cpu_time(middle1)
    call iterate(pos_vel, timestep, 'euler')
    call cpu_time(middle2)
    write(*,*) 'Euler: Time =', middle2 - middle1, 'seconds.'
    call write_to_file([ ((1._wp*i*dSteps*timestep),i=0,Nsteps/dSteps-1)  ], pos_vel(1:Nsteps:dSteps,:), e_file)
    call cpu_time(middle1)
    call iterate(pos_vel, timestep, 'midpoint')
    call cpu_time(middle2)
    write(*,*) 'Midpoint: Time =', middle2 - middle1, 'seconds.'
    call write_to_file([ ((1._wp*i*dSteps*timestep),i=0,Nsteps/dSteps-1)  ], pos_vel(1:Nsteps:dSteps,:), mp_file)
    call cpu_time(middle1)
    call iterate(pos_vel, timestep, 'rk4')
    call cpu_time(middle2)
    write(*,*) 'RK4: Time =', middle2 - middle1, 'seconds.'
    call write_to_file([ ((1._wp*i*dSteps*timestep),i=0,Nsteps/dSteps-1)  ], pos_vel(1:Nsteps:dSteps,:), rk4_file)

    ! --- Get exact solution to the trajectory.
    pos_vel_exact = calculate_exact_trajectory()
    call write_to_file([ ((1._wp*i*dSteps*timestep),i=0,Nsteps/dSteps-1) ], pos_vel_exact, exact_file)

    ! --- Deallocate arrays.
    deallocate(pos_vel, pos_vel_exact)

    ! --- Check for Helmholtz field.
    if (type_of_field == 'helmholtz') then
        rPos = 0._wp
        zPos = 0._wp
        write(*,*) 'Radial field strength in (0,0,0): ', integrate(r_integrand, 0._wp, 2._wp*pi, 100, method='simpson')
        write(*,*) 'Axial field strength in (0,0,0): ' , integrate(z_integrand, 0._wp, 2._wp*pi, 100, method='simpson')
        zPos = 1.0_wp
        write(*,*) 'Radial field strength in (0,0,1): ', integrate(r_integrand, 0._wp, 2._wp*pi, 100, method='simpson')
        write(*,*) 'Axial field strength in (0,0,1): ' , integrate(z_integrand, 0._wp, 2._wp*pi, 100, method='simpson')
        correct_answer = (coilRadius**2 + (zPos - 1._wp)**2)**(-1.5_wp) + (coilRadius**2 + (zPos + 1._wp)**2)**(-1.5_wp)
        correct_answer = correct_answer * 0.5_wp*(coilRadius**2 + 1._wp)**1.5_wp 
        write(*,*) 'Should be: ', correct_answer
    end if

    correct_answer = 20.035749854820_wp
    startPoint = 0._wp
    endPoint = 10._wp
    Nin = 100
    
    write(*,*) 'simpson: ', (integrate(afun, startPoint, endPoint, Nin, method='simpson') )
    write(*,*) 'leftrect: ', (integrate(afun, startPoint, endPoint, Nin, method='leftrect'))
    write(*,*) 'midrect: ', (integrate(afun, startPoint, endPoint, Nin, method='midrect'))
    write(*,*) 'rightrect: ', (integrate(afun, startPoint, endPoint, Nin, method='rightrect'))
    write(*,*) 'trapezoid: ', (integrate(afun, startPoint, endPoint, Nin, method='trapezoid'))

    ! --- Get final time and say goodbye.
    call cpu_time(finish)
    write(*,*) 'The run is done!'
    write(*,*)
    write(*,*) 'Time =', finish - start, 'seconds.'
    write(*,*)
	
end program main
