!
!// The code is in reduced units: new t = old t * omega, new x = old x/L where L is the period of saw-tooth potential, new U = old U/dU where dU is the amplitude of the potential.

! diffusjonsligning i bevegelig koordinatsystem

program biased_brownian_motion
    
    implicit none    !// Prevent the use of implicit declaration
    
    ! ######## Precision #########    
    ! -- single
    integer, parameter :: sp = kind(0.0)

    ! -- double
    integer, parameter :: dp = kind(0.0d0)

    ! -- working
    integer, parameter :: wp = dp
    
    ! ######## Variables #########    
    integer, parameter                  :: N = 10
    integer, parameter                  :: nSteps =  3000000!1100000
    real(wp)                            :: r                    !// m
    real(wp)                            :: x0
    real(wp)                            :: kT                   !// eV
    real(wp)                            :: dU                   !// eV
    real(wp)                            :: L                    !// m
    real(wp)                            :: alfa                 !// asymmetry factor
    real(wp)                            :: tau                  !// the period of the flash
    real(wp)                            :: zeta, gama, omega    !// dynamic viscosity, friction constant, kT/dU, dU/gama/L**2
    real(wp)                            :: dt
    real(wp), dimension(0:nSteps-1)     :: x
    real(wp), dimension(0:N-1)          :: vd
    real(wp)                            :: F_ext
    integer                             :: i, j, res
    integer, parameter                  :: measurement_period = 1000
    real(wp), parameter                 :: pi = 4.0_wp* atan(1.0_wp)
    
    ! -- standard values
    real(wp), parameter     :: std_r        = 12e-9  !12E-9_wp  ! m
    real(wp), parameter     :: std_L        = 20E-6_wp  ! m
    real(wp), parameter     :: std_dU       = 80.0_wp   ! eV
    real(wp), parameter     :: std_alfa     = 0.2_wp
    real(wp), parameter     :: std_tau      = 0.5_wp     
    real(wp), parameter     :: std_zeta     = 1E-3_wp   ! Pa*s
    real(wp), parameter     :: std_kT       = 0.026_wp  ! Termisk energi i eV
    real(wp), parameter     :: std_x0       = 0.0_wp
    real(wp), parameter     :: fraction_off = 3.0_wp / 4.0_wp

    ! -- files
    character(len=1024)                 :: file_trajectory
    character(len=1024)                 :: file_constants
    character(len=1024)                 :: file_vd
    character(len=*), parameter         :: file_fmt = '(A,I0,A,F0.1,A,F0.2,A,F0.4,A)'
    
    !// Initialise the random seed
    call init_random_seed    
    
    !// Read in values from the user
    call input_potential()    
    call input_particle(-1)        ! place inside loop to get individual r/x0
    gama = 6.0_wp*pi*zeta*r
    omega = dU/(gama*L**2)

    !// Criterion for the time step
    call get_dt()
    write(*,*) 'Timestep [s]: ', dt/omega
    write(*, fmt='(A,F0.1)') 'Total time [s]: ', dt*nSteps/omega
    
    !// Write to file
    write(file_trajectory, fmt=file_fmt) 'data/trajectory_N', N, '_rnm', r/1.0e-9, '_tau', tau, '_dU', dU*6.24150934e18, '.txt'
    write(file_constants,  fmt=file_fmt) 'data/constants_N', N, '_rnm', r/1.0e-9, '_tau', tau, '_dU', dU*6.24150934e18, '.txt'
    write(file_vd, fmt=file_fmt) 'data/drift_velocity_N', N, '_rnm', r/1.0e-9, '_tau', tau, '_dU', dU*6.24150934e18, '.txt'
    file_trajectory = trim(file_trajectory)
    file_constants = trim(file_constants)
    file_vd = trim(file_vd)
    
    call write_constants_to_file(file_constants)
    call create_empty_trajectory_file(file_trajectory)

    !// Euler scheme
    do while (tau <= 0.0)
        write(*,*) tau
        !// Start the diffusion process for the N particles
        do i = 0, (N-1)
            !if (mod(i, 10) == 0) then
            !    write(*,fmt='(I0, A)', advance='no') i, ' '
            !end if
            
            !// Set the starting position for the particles and then update x
            x(0) = x0
            do j = 1, nSteps - 1
                x(j) = updateX(x(j-1), j*dt, dt)
            end do
            
            call append_trajectory_to_file(file_trajectory)
            !vd(i) = (x(nSteps-1) - x(0))*L / (dt*nSteps/omega)
        end do
        !call append_drift_velocity_to_file(file_vd)
        tau = tau + 0.01_wp
    end do
    
    !call system ('python make_plot.py')
    
contains    

!
! Update the position of the particle using the Euler scheme
! using the explicit expression for the force
function updateX(x, t, dt) result(new_x)
    real(wp), intent(in)    :: x, t, dt
    real(wp)                :: new_x       

    new_x = x + F_e(x, t) * dt + sqrt(2*kT/dU*dt)*random_gauss()
    
end function


! ADDED
!Compute the potential U(x, t)
function U(xPos, time) result(pot)
    real(wp)                    :: pot
    real(wp), intent(in)        :: time, xpos
    real(wp)                    :: t, x
    
    x = xPos - floor(xPos)
    t = time - floor(time/(tau*omega))*tau*omega
    
    if (tau /= 0.0_wp .and. t >= 0.0_wp .and. t < fraction_off*tau*omega) then
        pot = 0.0_wp
    else
        if (x >= 0.0_wp .and. x < alfa) then
            pot = x/alfa
        else if (x >= alfa .and. x < 1.0_wp) then
            pot = (1.0_wp-x)/(1.0_wp-alfa)
        end if
    end if

end function


! ADDED
! Compute the force due to the potential U(x, t) (implicit)
function F_i(x, t) result(force)
    real(wp)                    :: force
    real(wp), intent(in)        :: x, t
    real(wp)                    :: dx = 1.0e-3_wp
    
    if (U(x,t) == 0.0_wp) then
        force = 0.0_wp
    else
        force = - ( U(x + dx, t) - U(x - dx, t) ) / (2.0_wp*dx)
    end if
    
end function

! ADDED
! Compute the force due to the potential U(x, t) (explicit)
function F_e(x, t) result(force)
    real(wp)                    :: force
    real(wp), intent(in)        :: x, t
    real(wp)                    :: x_temp
    real(wp)                    :: dx = 1.0e-3_wp
    
    x_temp = x - floor(x)

    if (U(x_temp,t) == 0.0_wp) then
        force = 0.0_wp
    else
        if (x_temp >= 0.0_wp .and. x_temp < alfa) then
            force = -1.0_wp/alfa
        else if (x_temp >= alfa .and. x_temp < 1.0_wp) then
            force = 1.0_wp/(1.0_wp-alfa)
        end if
    end if
    
end function

!
! Implement the time criterion
subroutine get_dt()
    dt = 0.1_wp
    
    if ( 1/alfa >= 1/(1-alfa) ) then
        do while ( 1/alfa*dt + 4*sqrt(2*kT/dU*dt) > 0.1_wp*alfa )
            dt = dt/2
        end do
    else 
        do while ( 1/(1-alfa)*dt + 4*sqrt(2*kT/dU*dt) > 0.1_wp*alfa )
            dt = dt/2
        end do
    end if
    
    if ( 1/alfa >= 1/(1-alfa) ) then
        write(*,*) (1/alfa*dt + 4*sqrt(2*kT/dU*dt)) / alfa * 100
    else
        write(*,*) (1/(1-alfa)*dt + 4*sqrt(2*kT/dU*dt)) / alfa * 100
    end if
    
end subroutine



!
! Append the drift velocitys to file (diffusion length per s) )
subroutine append_drift_velocity_to_file(filename)
    character(len=*), intent(in)    :: filename
    real(wp)                        :: avg_vd
    real                            :: std_dev
    logical                         :: exist
    integer                         :: i

    if (N >= 50) then
        inquire(file=filename, exist=exist)
        if (exist) then
            open(12, file=filename, status="old", position="append", action="write")
        else
            open(12, file=filename, status="new", action="write")
        end if
    
        avg_vd = sum(vd)/N

        std_dev = 0.0_wp
        do i = 0, N-1
            std_dev = std_dev + (avg_vd - vd(i))**2
        end do
        std_dev = sqrt(std_dev/(N-1))

        write(12, *) tau, avg_vd, std_dev
        close(12)
    else
        write(*,*) 'The number of particles is too small to get an accurate measurment of the average drift velocity.'
        write(*,*) 'You should use N >= 50. You have N = ', N
    end if
end subroutine

subroutine create_empty_trajectory_file(file_trajectory)
    character(len=*), intent(in)        :: file_trajectory
    open(unit=13,file=file_trajectory, form="formatted", status="replace", action="write")
    close(unit=13) 
end subroutine
    
!
! Writes the trajectory of the particle(s) to file
subroutine append_trajectory_to_file(filename)
    character(len=*), intent(in)   :: filename
    open(unit=1, file=filename, status="old", position="append", action="write")
    write(1,*) x(0:nSteps-1:measurement_period)
    close(unit=1)
end subroutine


!
! Writes constants to file
subroutine write_constants_to_file(filename)
    character(len=*), intent(in)        :: filename
    
    open(unit=13,file=filename, form="formatted", status="replace", action="write", iostat = res)
    if (res /= 0) then
        write(*,*) "Error in opening file, status", res
        stop
    end if

    write(13,*) N, nSteps/measurement_period   !// #of particles, #of time steps
    write(13,*) dt * measurement_period              !// delta_t
    write(13,*) alfa                                !// alfa
    write(13,*) tau                                 !// tau
    write(13,*) zeta                                !// zeta
    write(13,*) r                                   !// radius of particle
    write(13,*) dU                                  !// dU
    write(13,*) L                                   !// L
    write(13,*) kT                                  !// kT
    
    close(unit=13)
end subroutine



!
! Takes input from the user
subroutine input_potential()
    character(len=*), parameter :: format_1 = "(A, T30, A, F14.10, A)"
    
    write(*, *) 'INPUT PARAMETERS:'
    write(*, format_1, advance='no') 'Period of flash, tau ','(', std_tau, '):'
    call input_w_default(tau, std_tau)
    write(*, format_1, advance='no') 'Potential strength, dU [eV] ','(', std_dU, '):'
    call input_w_default(dU, std_dU)
    dU = dU * 1.60217657E-19_wp
    write(*, format_1, advance='no') 'Period of saw-tooth potential, L ','(', std_L, '):'
    call input_w_default(L, std_L)
    write(*, format_1, advance='no') 'Asymmetric factor, alfa ','(', std_alfa, '):'
    call input_w_default(alfa, std_alfa)
    write(*, format_1, advance='no') 'Dynamic viscosity, zeta ','(', std_zeta, '):'
    call input_w_default(zeta, std_zeta)
    write(*, format_1, advance='no') 'k_B * T [eV] ','(', std_kT, '):'
    call input_w_default(kT, std_kT)
    kT = kT * 1.60217657E-19_wp
    
    !// Validation
    if (dU == 0 .OR. L == 0 .OR. zeta == 0 .OR. kT == 0 .OR. alfa == 0) then
        write(*, *) "Neither dU, L, zeta, kT nor alfa are allowed to be 0!"
        write(*, *) "Exiting."
        call EXIT(1)
    end if
end subroutine



!
! Takes input from the user
subroutine input_particle(i)
    integer, intent(in)                :: i
    
    character(len=*), parameter :: format_2 = "(A, I3.0, T30, A, F14.10, A)"

    write(*, format_2, advance='no') 'Radius r of particle ', i+1, '(', std_r, '):'
    call input_w_default(r, std_r)
    write(*, format_2, advance='no') 'Starting position x0 of particle ', i+1, '(', std_x0, '):'
    call input_w_default(x0, std_x0)
    
    !// Validation
    if (r == 0) then
        write(*, *) "r are not allowed to be 0!"
        write(*, *) "Exiting."
        call EXIT(1)
    end if
end subroutine


!
! If no input is given the parameter is set to the default parameter
subroutine input_w_default(param, default_param)
    real(wp), intent(out) :: param
    real(wp), intent(in)  :: default_param

    integer, parameter              :: MAX_LENGTH   = 100
    character(len=MAX_LENGTH)       :: input_str
    real(wp)                        :: input_real
    integer                         :: flag

    read(*, "(A)") input_str

    read(input_str, *, iostat=flag) input_real
    if (flag == 0) then
        param = input_real
    else
        param = default_param
    end if
end subroutine



! ADDED
! Gaussian random number
! Draw Gaussian distributed random numbers. Normal distribution is 
! implemented by the polar Box-Müller algorithm using the uniform 
! distribution from random_number().
function random_gauss() result(rand_gauss)
    real(wp)            :: rand_gauss
    real(wp)            :: rand1, rand2, w
    integer             :: i
    
    w = 0.0_wp
    do while(w == 0.0_wp .OR. w >= 1)
        call random_number(rand1)
        call random_number(rand2)
        rand1 = rand1 * 2.0_wp - 1.0_wp
        rand2 = rand2 * 2.0_wp - 1.0_wp
        w = rand1**2 + rand2**2
    end do
    rand_gauss = rand1 * SQRT( - 2.0_wp * LOG(w) / w)
end function



! ADDED
! Draw N numbers with the random_gauss() function. Writes them to file to 
! check if it is a normal distribution
subroutine check_gaussian()
    open(unit=1,file="check_gaussian.txt", form="formatted", status="replace", action="write", iostat = res)
    
    if (res /= 0) then
        write(*,*) "Error in opening file, status", res
        stop
    end if
    
    write(*,*) random_gauss()
    write(*,*) random_gauss()
    
    do i = 0, 10000
        write(1,*) random_gauss()
    end do
    close(unit=1)
    
end subroutine



! ADDED
! Initialise the random seed for the random number generator random_number()
subroutine init_random_seed()
  integer                                 :: i, n, clock
  integer, dimension(:), allocatable     :: seed

  call random_seed(size = n)
  allocate(seed(n))

  call system_clock(count=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)

  deallocate(seed)
end subroutine

end program biased_brownian_motion

