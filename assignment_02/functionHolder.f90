module functionHolder

    use common

    implicit none

    contains

    function afun(x)
        real(wp)             :: afun
        real(wp), intent(in) :: x
 
        afun = x - x**3
    end function afun

end module functionHolder
