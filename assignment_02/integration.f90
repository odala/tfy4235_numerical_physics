module integration

    use common
    
    implicit none
 
    contains
 
    ! --- Function that evaluates an integral.
    !     input : function, lower limit, upper limit, number of steps, method:
    !             > 'leftrect'
    !             > 'midrect'
    !             > 'rightrect'
    !             > 'trapezoid'
    !             > 'simpson' 
    !     output: the solution to the integral
    function integrate(f, a, b, n, method)
        real(wp) :: integrate
        real(wp), intent(in) :: a, b
        integer, intent(in)  :: n
        character(len=*), intent(in), optional :: method
        interface
            function f(x)
                use common
                real(wp)             :: f
                real(wp), intent(in) :: x
            end function f
        end interface
 
        integer  :: i, m
        real(wp) :: h
        real(wp), dimension(:), allocatable :: xpoints
        real(wp), dimension(:), target, allocatable :: fpoints
        real(wp), dimension(:), pointer :: fleft, fmid, fright
 
        h = (b - a) / n
 
        allocate(xpoints(0:2*n), fpoints(0:2*n))
 
        xpoints = (/ (a + h*i/2, i = 0,2*n) /)
 
        do i = 0,2*n      
            fpoints(i) = f(xpoints(i))
        end do
        fleft  => fpoints(0 : 2*n-2 : 2)
        fmid   => fpoints(1 : 2*n-1 : 2)
        fright => fpoints(2 : 2*n   : 2)

        if ( present(method) ) then
            select case (method)
            case ('leftrect')
                m = 1
            case ('midrect')
                m = 2
            case ('rightrect')
                m = 3
            case ( 'trapezoid' )
                m = 4
            case default
                m = 0
            end select
        else
            m = 0
        end if
 
        select case (m)
        case (0) ! simpson
            integrate = h / 6.0_wp * sum(fleft + 4.0_wp*fmid + fright)
        case (1) ! leftrect
            integrate = h * sum(fleft)
        case (2) ! midrect
            integrate = h * sum(fmid)
        case (3) ! rightrect
            integrate = h * sum(fright)
        case (4) ! trapezoid
            integrate = h / 2.0_wp * sum(fleft + fright)
        end select
 
        deallocate(xpoints, fpoints)

    end function integrate
 
end module integration
