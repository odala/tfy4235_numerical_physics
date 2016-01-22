Module functions

    ! --- The modules we need to use in the module.
    Use common
    Implicit None

    Contains

    ! --- Subroutine that initializes the ecosystem
    !     Input: the array 'ecosystem'
    !     Output: none but the array is updated
    Subroutine initialize_ecosystem(ecosystem)
        Implicit None
        Integer  :: i, N
        Real(wp) :: r, ecosystem(:)

        ! --- Get size of array
        N = size(ecosystem,1)

        ! --- Loop over the array
        do i = 1,N
            call random_number(r)
            ecosystem(i) = r
        end do

    end Subroutine initialize_ecosystem

    ! --- Function that locates the minimum in a array
    !     Input: the array 'ecosystem'
    !     Output: the location of the minimum in the array imin
    Function locate_minimum(ecosystem) result(imin)
        Implicit None
        Integer  :: imin, i, N
        Real(wp) :: mini, ecosystem(:)

        ! --- Get the size
        N = size(ecosystem,1)

        mini = 1._wp
        imin = 1

        ! --- Locate the minimum by looping through the array and
        !     keep track of the minimum and its position
        do i = 1, N
        
        if (mini > ecosystem(i)) then
            mini = ecosystem(i)
            imin = i
        end if

        end do

    end Function locate_minimum

    ! --- Subroutine that updates the ecosystem
    !     Input: the array 'ecosystem'
    !     Ouput: none but array updated
    Subroutine update_ecosystem(ecosystem)
        Implicit None
        Integer  :: imin, N
        Real(wp) :: r, ecosystem(:)

        ! --- Get size of the array
        N = size(ecosystem,1)

        ! --- Locate minimum
        imin = locate_minimum(ecosystem)

        ! --- Write location of the minimum to outpute file
        write(1,*) imin

        ! --- Assign a random fitness to the minimum
        call random_number(r)
        ecosystem(imin) = r

        ! --- Assign a random fitness to the left neighbor of the minimum
        if (imin > 1) then
            call random_number(r)
            ecosystem(imin-1) = r
        else
            call random_number(r)
            ecosystem(N) = r
        end if

        ! --- Assign a random fitness to the right neighbor of the minimum
        if (imin < N) then
        call random_number(r)
        ecosystem(imin+1) = r
        else
        call random_number(r)
        ecosystem(1) = r
        end if

    end Subroutine update_ecosystem
        
    ! --- Subroutine that makes the ecosystem evolve
    !     Input: the array 'ecosystem' and number of steps Nsteps
    !     Output: none but the ecosystem evolve several times
    Subroutine evolution(ecosystem,Nsteps)
        Implicit None
        Integer :: i, Nsteps
        Real(wp) :: ecosystem(:)

        ! --- Loop over number of steps
        do i = 1, Nsteps

            !write(*,*) ecosystem(:)
            ! --- Update the ecosystem (i.e. locate min and update fitness)
            call update_ecosystem(ecosystem)

        end do

    end Subroutine evolution

end Module functions
