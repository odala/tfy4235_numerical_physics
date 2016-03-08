module common

    implicit none
    ! --- Define default precision.
    integer, parameter      :: sp  = KIND(1.0)
    integer, parameter      :: dp  = KIND(1.0D0)
    integer, parameter      :: wp  = dp

    ! --- Define global constants
    real(wp), parameter     :: pi  = 4.0_wp*atan(1.0_wp)
    real(wp), parameter     :: e   = 1.602176565E-19 ! elementary charge (C)
    real(wp), parameter     :: m_e = 9.10938291E-31  ! electron mass (kg)
    real(wp), parameter     :: m_p = 1.672621777E-27 ! proton mass (kg)
    character(len=*), parameter :: type_of_field = 'curved' !
    real(wp), parameter     :: coilRadius = 0.2

    ! --- Define global parameters
    integer     :: Nparticles, Nsteps, dSteps, sgn
    real(wp)    :: timestep, x0, y0, z0, u0, v0, w0, mass, v_perp, v_para, Bz, Ex, Ey, Ez, zPos, rPos
    
    ! --- Define electric and magnetic field arrays.
    real(wp), dimension(1:3)    :: Efield, Bfield    

end module common
