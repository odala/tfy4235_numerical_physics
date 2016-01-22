Module readinput

    ! --- The modules we need to use in the module.
    Use common
    Contains

    ! --- Subroutine that reads an input file
    !     and assign values to parameters.
    Subroutine read_input_file()
        Implicit none
        Integer             :: Nspecies, Nsteps

        ! --- List of names to look for in the input file.
        !     These are local variables.
        Namelist / parameters / & 
        Nspecies,           &
        Nsteps

        ! --- Read the input file and assign corresponding values
        !     to named variables.
        open(unit=1, file='input.txt', status='old')
        read(1,nml=parameters)
        close(1)

        ! --- Assign the values read to the parameters (global variables).
        param%Nspecies = Nspecies
        param%Nsteps   = Nsteps

    end subroutine read_input_file

end Module readinput
