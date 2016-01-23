Module readinput

    ! --- The modules we need to use in the module.
    Use common
    Contains

    ! --- Subroutine that reads an input file
    !     and assign values to parameters.
    Subroutine read_input_file()
        Implicit none
        Integer             :: iNparticles, iNsteps
        Real(wp)            :: iradius, ix0, itau, idU, iL, ialpha, izeta, ikT, ifraction_off

        ! --- List of names to look for in the input file.
        !     These are local variables.
        Namelist / parameters / & 
        iNparticles,            &
        iNsteps,                &
        iradius,                &
        ix0,                    &
        itau,                   &
        idU,                    &
        iL,                     &
        ialpha,                 &
        izeta,                  &
        ikT,                    &
        ifraction_off            

        ! --- Read the input file and assign corresponding values
        !     to named variables.
        open(unit=1, file='input.txt', status='old')
        read(1,nml=parameters)
        close(1)

        ! --- Assign the values read to the parameters (global variables).
        Nparticles      = iNparticles
        Nsteps          = iNsteps
        radius          = iradius
        x0              = ix0
        tau             = itau
        dU              = idU
        L               = iL
        alpha           = ialpha
        zeta            = izeta
        kT              = ikT
        fraction_off    = ifraction_off
        
        ! --- Calculate values from these parameters
        gama    = 6.0_wp*pi*zeta*radius
        omega   = dU/(gama*L**2)

    end subroutine read_input_file

end Module readinput
