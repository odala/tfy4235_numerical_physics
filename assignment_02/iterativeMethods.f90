module iterativeMethods

    use common
    
    implicit none

    abstract interface
        function derivative(y, t)
            use common
            real(wp), intent(in)   :: y(:)
            real(wp), intent(in)   :: t
            real(wp), dimension(6) :: derivative
        end function
    end interface
 
    contains    
    ! ---------------------------
    ! --- RUNGE-KUTTA SCHEMES --- 
    ! ---------------------------
    
    ! --- Function that returns the updated value 
    !     of y using the forward Euler scheme.
    !     Input : old y (1D-array), old_t, h
    !     Output: none but new_y is updated
    function forward_euler(f, old_y, old_t, h) result(new_y)
        procedure(derivative)   :: f
        real(wp), intent(in)    :: old_y(:), old_t, h
        real(wp), dimension(size(old_y)):: new_y

        if (h <= 0) stop "negative step size"

        new_y = old_y + h * f(old_y, old_t)

    end function forward_euler

    ! --- Function that returns the updated value 
    !     of y using the forward midpoint scheme.
    !     Input : old y (1D-array), old_t, h
    !     Output: none but new_y is updated
    function midpoint(f, old_y, old_t, h) result(new_y)
        procedure(derivative)           :: f
        real(wp), intent(in)            :: old_y(:), old_t, h
        real(wp), dimension(size(old_y)):: new_y, k1

        if (h <= 0) stop "negative step size"
        
        ! --- Evaluate the slope of y at the midpoint of the interval.
        k1 = f(old_y, old_t)

        ! --- Evaluate y at the end of the interval.
        new_y = old_y + h * f(old_y + 0.5_wp * h * k1, old_t + 0.5_wp*h)

    end function midpoint

    ! --- Function that returns the updated value 
    !     of y using the 4th order Runge-Kutta scheme.
    !     Input : old y (1D-array), old_t, h
    !     Output: none but new_y is updated
    function rk4(f, old_y, old_t, h) result(new_y)
        procedure(derivative)   :: f
        real(wp), intent(in)    :: old_y(:), old_t, h
        real(wp), dimension(size(old_y)):: new_y, k1, k2, k3, k4

        if (h <= 0) stop "negative step size"
        
        ! --- Evaluate the slope of y at the beginning of the interval.
        k1 = f(old_y, old_t)

        ! --- Evaluate the slope of y at the midpoint of the interval.
        k2 = f(old_y + 0.5_wp*h*k1, old_t+0.5_wp*h)     

        ! --- Evaluate the slope of y at the midpoint of the interval.
        k3 = f(old_y + 0.5_wp*h*k2, old_t+0.5_wp*h) 

        ! --- Evaluate the slope of y at the end of the interval.
        k4 = f(old_y + h*k3, old_t+h)  

        ! --- Evaluate y at the end of the interval.
        new_y = old_y + h/6.0_wp * ( k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4 )

    end function rk4

end module iterativeMethods
