Module common

    Implicit None
    ! --- Define default precision.
    Integer, Parameter      :: sp  = KIND(1.0)
    Integer, Parameter      :: dp  = KIND(1.0D0)
    Integer, Parameter      :: wp  = dp

    ! --- Define global constants
    Real(wp), Parameter     :: pi   = 4.0_wp*atan(1.0_wp)
    Real(wp), Parameter     :: e   = 1.602176565E-19 ! elementary charge (C)
    Real(wp), Parameter     :: m_e = 9.10938291E-31  ! electron mass (kg)
    Real(wp), Parameter     :: m_p = 1.672621777E-27 ! proton mass (kg)

    ! --- Define global parameters
    Integer     :: Nparticles, Nsteps, dSteps, sgn
    Real(wp)    :: timestep, x0, y0, z0, u0, v0, w0, omega, r, mass, vPerp, vPara, Bfield, Efield
    
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
