Program main

    ! --- The modules we need to use in the program.
    Use common
    Use readinput
    !Use functions

    ! --- Define variables.
    Implicit None
    Real(wp)                :: start, finish
    !Real(wp), allocatable   :: x(:)
    
    ! --- Get initial time and say hello.
    write(*,*) 'Program starting!'
    write(*,*)
    call cpu_time(start)
    
    ! --- Read input file.
    call read_input_file()

    ! --- Now param%Nspecies and param%Nsteps have values
    !     given in the input file.

    ! --- Write to command line a few words to check
    !     the parameters have read correctly.
    write(*,*) 'Check parameters.'
    write(*,*) param%Nparticles, param%Nsteps, param%gama, param%omega


    ! --- Allocate space for the ecosystem array.
    !Allocate(ecosystem(param%Nspecies))

    ! --- Initialize random number generator.
    !call initialize_random_seed()


    ! --- Open output file named 'minloc.dat.
    !     It will contain the location of 
    !     all minima through the simulation.
    !open(unit=1,file='minloc.dat')


    ! --- Initialize ecosystem.
    !call initialize_ecosystem(ecosystem)


    ! --- Make the ecosystem evolve during Nsteps.
    !call evolution(ecosystem,param%Nsteps)

    ! --- Close output file.
    !close(1)

    ! --- Deallocate array.
    !Deallocate(ecosystem)

    ! --- Get final time and say goodbye.
    call cpu_time(finish)
    write(*,*) 'The run is done!'
    write(*,*) 'Time =', finish - start, 'seconds.'
	
end Program main
