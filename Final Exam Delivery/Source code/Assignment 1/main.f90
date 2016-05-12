Program main

    ! --- The modules we need to use in the program.
    Use common
    Use readinput
    Use functions

    ! --- Define variables.
    Implicit None
    Integer                 :: i
    Real(wp)                :: start, finish
    Real(wp), allocatable   :: positions(:)
    Real(wp), allocatable   :: drift_velocities(:)
    Real(wp)                :: avg_vd
    Real(wp)                :: std_dev
    Character(len=*), parameter :: input_file = 'input.txt'
    Character(len=*), parameter :: trajectory_file = 'output.txt'
    Character(len=*), parameter :: drift_velocity_file = 'drift_velocities.txt'
    Character(len=*), parameter :: format_0 = '(A15, I14)'
    Character(len=*), parameter :: format_1 = '(A15, E14.4)'
    
    
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
    write(*,format_0) 'Nparticles: '  , Nparticles
    write(*,format_0) 'Nsteps: '      , Nsteps
    write(*,format_1) 'timestep (s): ', timestep
    write(*,format_1) 'radius (m): '  , radius
    write(*,format_1) 'x0 (m): '      , x0
    write(*,format_1) 'tau: (s)'      , tau
    write(*,format_1) 'dU: (eV)'      , dU/1.60217657E-19_wp
    write(*,format_1) 'L: (m)'        , L             
    write(*,format_1) 'alpha: '       , alpha
    write(*,format_1) 'zeta: (Pa*s)'  , zeta
    write(*,format_1) 'kT: (eV)'      , kT/1.60217657E-19_wp
    write(*,format_1) 'fraction_off: ', fraction_off
    write(*,format_0) 'dSteps: '      , dSteps
    write(*,*)
    
    ! --- Introduce reduced units.
    timestep    = timestep*omega

    ! --- Check time criterion.
    if (check_time_criterion()) then
        write(*,*) 'Time criterion is fulfilled.'
    else 
        write(*,*) 'Time criterion is not fulfilled.'
        write(*,*) 'This program is going to exit.'
        !call exit(0)
    end if

    ! --- Write N numbers from the
    !     get_random_gauss() to file.
    !call check_random_gauss(1000000)

    ! --- Allocate space for the position- and drift velocity arrays.
    Allocate(positions(Nsteps), drift_velocities(Nparticles))

    ! --- Initialize random number generator.
    call initialize_random_seed()

    ! --- Open output file named 'output.txt'.
    !     It will contain the positions of 
    !     the particles through the simulation.
    open(unit=1, file=trajectory_file, form="formatted", status="replace", action="write")
    
    ! --- Make the particles evolve during Nsteps.
    !     Using the euler scheme on the N particles.
    do i = 1, Nparticles
        write(*,*) i
        call euler_scheme(positions, timestep)
        write(1,*) positions(1:Nsteps:dSteps)
    end do

    ! --- Close output file.
    close(unit=1)

    ! --- Repeat the simulation for different values of tau.
    do while (tau <= -1.0_wp)
        write(*,*) 'tau = ', tau, ': '

        ! --- Calculate the individual drift velocites in real units.
        do i = 1, Nparticles
            call euler_scheme(positions, timestep)
            drift_velocities(i) = (positions(Nsteps) - positions(1))*L/(Nsteps*timestep/omega)
        end do

        write(*,*) 'Total time: ', Nsteps*timestep/omega
        
        ! --- Calculate average drift velocity.
        avg_vd = sum(drift_velocities)/Nparticles
        write(*,*) 'Average v: ', avg_vd

        ! --- Calculate standard deviation.
        std_dev = 0.0_wp
        do i = 1, Nparticles
            std_dev = std_dev + (avg_vd - drift_velocities(i))**2
        end do
        std_dev = sqrt(std_dev/(Nparticles-1))
        
        ! --- Write to the file 'drift_velocities.txt'.
        !     It will contain the average drift velocities 
        !     with standard deviation for different flashing periods.
        call append_vector_to_file([tau, avg_vd, std_dev], drift_velocity_file)

        tau = tau + 0.01_wp
    end do

    ! --- Deallocate arrays.
    Deallocate(positions, drift_velocities)

    ! --- Get final time and say goodbye.
    call cpu_time(finish)
    write(*,*) 'The run is done!'
    write(*,*)
    write(*,*) 'Time =', finish - start, 'seconds.'
    write(*,*)
	
end Program main
