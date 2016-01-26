Module common

    Implicit None
    ! --- Define default precision.
    Integer, Parameter      :: sp  = KIND(1.0)
    Integer, Parameter      :: dp  = KIND(1.0D0)
    Integer, Parameter      :: wp  = dp

    ! --- Define global constants
    Real(wp), Parameter     :: pi = 4.0_wp*atan(1.0_wp)

    ! --- Define global parameters
    Integer     :: Nparticles, Nsteps, dSteps
    Real(wp)    :: timestep, radius, x0, tau, dU, L, alpha, zeta, kT, fraction_off, gama, omega
        
    !Type :: Parameter_Container
    ! -----------------------------------------------------------
    ! --- This is the container for the simulation parameters ---
    ! -----------------------------------------------------------

    !Integer :: Nparticles, Nsteps
    !Real    :: radius, x0, tau, dU, L, alpha, zeta, kT, fraction_off, gama, omega

    !end Type Parameter_Container

    ! --- The shared container named param.
    !Type(Parameter_Container), Save :: param

end Module common
