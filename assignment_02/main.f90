Program main

    ! --- The modules we need to use in the program.
    Use common
    Use readinput
    Use functions

    ! --- Define variables.
    Implicit None
    Integer                 :: i
    Real(wp)                :: start, middle1, middle2, finish
    Real(wp), allocatable   :: xs(:), ys(:), zs(:), us(:), vs(:), ws(:)
    
    ! --- Define name of files.
    Character(len=*), parameter :: input_file = 'input/input.txt'
    Character(len=*), parameter :: trajectory_file = 'output/e_output.txt'
    Character(len=*), parameter :: trajectory_file_2 = 'output/mp_output.txt'
    Character(len=*), parameter :: trajectory_file_3 = 'output/rk4_output.txt'
    Character(len=*), parameter :: drift_velocity_file = 'drift_velocities.txt'
    Character(len=*), parameter :: format_0 = '(A15, I14)'
    Character(len=*), parameter :: format_1 = '(A15, E14.4)'
    Character(len=*), parameter :: fmt_xyz = '(A15, E14.4, A2, E14.4, A2, E14.4, A2)'
    
    
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
    write(*,fmt_xyz) '(x0, y0, z0) (m): ( ', x0, ', ', y0, ', ', z0, ' )'
    write(*,fmt_xyz) '(v0, u0, w0) (m/s): ( ', u0, ', ', v0, ', ', w0, ' )'
    write(*,format_1)  'Mass: (kg)'   , mass
    write(*,format_0) 'dSteps: '      , dSteps
    write(*,format_1) 'B (T): '       , Bfield
    write(*,format_1) 'E (T): '       , Efield
    write(*,*)

    ! --- Introduce reduced units.
    x0 = x0/r
    y0 = y0/r
    z0 = z0/r
    u0 = u0/vPerp
    v0 = v0/vPerp
    w0 = w0/vPerp
    vPerp = vPerp/vPerp
    vPara = vPara/vPerp
    write(*,*) 'omega: ', omega
    !timestep    = timestep*omega
    Bfield = Bfield*e/(mass*omega)
    Efield = Efield*e/(mass*r*omega**2)
    write(*,*) 'Bfield: ', Bfield
    write(*,*) 'Efield: ', Efield

    write(*,*) 'Drift velocity = ', calculate_drift_velocity(Efield, Bfield)

    ! --- Allocate space for the position- and velocity arrays.
    Allocate(xs(Nsteps), ys(Nsteps), zs(Nsteps), us(Nsteps), vs(Nsteps), ws(Nsteps))

    ! --- Initialize random number generator.
    !call initialize_random_seed()

    ! --- Open output file named 'output.txt'.
    !     It will contain the positions of 
    !     the particles through the simulation.
    open(unit=1, file=trajectory_file, form="formatted", status="replace", action="write")
    open(unit=2, file=trajectory_file_2, form="formatted", status="replace", action="write")
    open(unit=3, file=trajectory_file_3, form="formatted", status="replace", action="write")
    
    ! --- Make the particles evolve during Nsteps.
    !     Using the euler scheme on the N particles.
    do i = 1, Nparticles
        write(*,*) i
        call cpu_time(middle1)
        call euler_scheme(xs, ys, zs, us, vs, ws, timestep)
        call cpu_time(middle2)
        write(*,*) 'Euler: Time =', middle2 - middle1, 'seconds.'
        write(1,*) xs(1:Nsteps:dSteps)
        write(1,*) ys(1:Nsteps:dSteps)
        write(1,*) zs(1:Nsteps:dSteps)
        call cpu_time(middle1)
        call midpoint_scheme(xs, ys, zs, us, vs, ws, timestep)
        call cpu_time(middle2)
        write(*,*) 'Midpoint: Time =', middle2 - middle1, 'seconds.'
        write(2,*) xs(1:Nsteps:dSteps)
        write(2,*) ys(1:Nsteps:dSteps)
        write(2,*) zs(1:Nsteps:dSteps)
        call cpu_time(middle1)
        call rk4_scheme(xs, ys, zs, us, vs, ws, timestep)
        call cpu_time(middle2)
        write(*,*) 'RK4: Time =', middle2 - middle1, 'seconds.'
        write(3,*) xs(1:Nsteps:dSteps)
        write(3,*) ys(1:Nsteps:dSteps)
        write(3,*) zs(1:Nsteps:dSteps)
        
    end do

    ! --- Close output file.
    close(unit=1)
    close(unit=2)
    close(unit=3)

    ! --- Deallocate arrays.
    Deallocate(xs, ys, zs, us, vs, ws)

    ! --- Get final time and say goodbye.
    call cpu_time(finish)
    write(*,*) 'The run is done!'
    write(*,*)
    write(*,*) 'Time =', finish - start, 'seconds.'
    write(*,*)
	
end Program main
