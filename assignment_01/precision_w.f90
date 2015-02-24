module precision_w

    implicit none

    private         ! Default to private

    public :: sp, dp, wp


    ! ######## Precision #########
    
    ! -- single
    integer, parameter :: sp = kind(0.0)

    ! -- double
    integer, parameter :: dp = kind(0.0d0)

    ! -- working
    integer, parameter :: wp = dp

end module
