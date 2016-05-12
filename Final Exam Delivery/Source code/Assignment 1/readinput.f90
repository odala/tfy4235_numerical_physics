Module readinput

    ! --- The modules we need to use in the module.
    Use common
    Contains

    ! --- Subroutine that reads an input file
    !     and assign values to parameters.
    Subroutine read_input_file(filename)
        Implicit none
        Character(len=*), Intent(in) :: filename
        Integer                      :: iNparticles, iNsteps, idSteps
        Real(wp)                     :: itimestep, iradius, ix0, itau, idU, iL, ialpha, izeta, ikT, ifraction_off

        ! --- List of names to look for in the input file.
        !     These are local variables.
        Namelist / parameters / & 
        iNparticles,            &
        iNsteps,                &
        itimestep,              &        
        iradius,                &
        ix0,                    &
        itau,                   &
        idU,                    &
        iL,                     &
        ialpha,                 &
        izeta,                  &
        ikT,                    &
        ifraction_off,          &
        idSteps                  

        ! --- Read the input file and assign corresponding values
        !     to named variables.
        open(unit=1, file=filename, status='old')
        read(1,nml=parameters)
        close(1)

        ! --- Assign the values read to the parameters (global variables).
        Nparticles      = iNparticles
        Nsteps          = iNsteps
        timestep        = itimestep
        radius          = iradius
        x0              = ix0
        tau             = itau
        dU              = idU*1.60217657E-19_wp     ! ! Scaling factor from eV to J
        L               = iL
        alpha           = ialpha
        zeta            = izeta
        kT              = ikT*1.60217657E-19_wp     ! Scaling factor from eV to J
        fraction_off    = ifraction_off
        dSteps          = idSteps
        
        ! --- Calculate values from these parameters
        gama    = 6.0_wp*pi*zeta*radius
        omega   = dU/(gama*L**2)

    end subroutine read_input_file

end Module readinput
