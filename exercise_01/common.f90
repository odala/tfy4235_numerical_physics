Module common

    Implicit None
    ! --- Define default precision.
    Integer, Parameter :: sp  = KIND(1.0)
    Integer, Parameter :: dp  = KIND(1.0D0)
    Integer, Parameter :: wp  = sp
        
    Type :: Parameter_Container
    ! -----------------------------------------------------------
    ! --- This is the container for the simulation parameters ---
    ! -----------------------------------------------------------

    Integer :: Nspecies, Nsteps

    end Type Parameter_Container

    ! --- The shared container named param.
    Type(Parameter_Container), Save :: param

end Module common
