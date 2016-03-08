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
        Real(wp)                     :: itimestep, isgn, ix0, iy0, iz0, iu0, iv0, iw0, iBz, iEx, iEy, iEz

        ! --- List of names to look for in the input file.
        !     These are local variables.
        Namelist / parameters / & 
        iNparticles,            &
        iNsteps,                &
        idSteps,                &
        itimestep,              &        
        iBz,                    &
        iEx,                    &
        iEy,                    &
        iEz,                    &
        isgn,                   &
        ix0,                    &
        iy0,                    &
        iz0,                    &
        iu0,                    &
        iv0,                    &
        iw0  

        ! --- Read the input file and assign corresponding values
        !     to named variables.
        open(unit=1, file=filename, status='old')
        read(1, nml=parameters)
        close(1)

        ! --- Assign the values read to the parameters (global variables).
        Nparticles      = iNparticles
        Nsteps          = iNsteps
        dSteps          = idSteps
        timestep        = itimestep
        sgn             = isgn
        x0              = ix0
        y0              = iy0
        z0              = iz0
        u0              = iu0
        v0              = iv0
        w0              = iw0
        Bz              = iBz
        Ex              = iEx
        Ey              = iEy
        Ez              = iEz
                
        ! --- Calculate values from these parameters
        !if (sgn > 0) then
        !    mass = m_p
        !else
        !    mass = m_e
        !end if
        !omega = e*Bz/mass
        v_perp = sqrt((u0)**2+v0**2)
        v_para = w0
        !r = vPerp/omega

        ! --- Initialize electric and magnetic field arrays.
        Efield(1) = Ex
        Efield(2) = Ey
        Efield(3) = Ez
        Bfield(1) = 0.0_wp
        Bfield(2) = 0.0_wp
        Bfield(3) = Bz
        
    end subroutine read_input_file

end Module readinput
