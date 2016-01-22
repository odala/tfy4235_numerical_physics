Module readinput

    ! --- The modules we need to use in the module.
    Use common
    Contains

    ! --- Subroutine that reads an input file
    !     and assign values to parameters.
    Subroutine read_input_file()
        Implicit none
        Integer             :: Nparticles, Nsteps
        Real(wp)            :: radius, x0, tau, dU, L, alpha, zeta, kT, fraction_off

        ! --- List of names to look for in the input file.
        !     These are local variables.
        Namelist / parameters / & 
        Nparticles,             &
        Nsteps,                 &
        radius,                 &
        x0,                     &
        tau,                    &
        dU,                     &
        L,                      &
        alpha,                   &
        zeta,                   &
        kT,                     &
        fraction_off            

        ! --- Read the input file and assign corresponding values
        !     to named variables.
        open(unit=1, file='input.txt', status='old')
        read(1,nml=parameters)
        close(1)

        ! --- Assign the values read to the parameters (global variables).
        param%Nparticles = Nparticles
        param%Nsteps   = Nsteps
        param%radius = radius
        param%x0 = x0
        param%tau = tau
        param%dU = dU
        param%L = L
        param%alpha = alpha
        param%zeta = zeta
        param%kT = kT
        param%fraction_off = fraction_off
        
        ! --- Calculate values from these parameters
        param%gama = 6.0_wp*pi*zeta*radius
        param%omega = dU/(param%gama*L**2)

    end subroutine read_input_file

end Module readinput
