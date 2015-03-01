program biased_brownian_motion
	
	implicit none	!// Prevent the use of implicit declaration
	
	! ######## Precision #########    
    ! -- single
    integer, parameter :: sp = kind(0.0)

    ! -- double
    integer, parameter :: dp = kind(0.0d0)

    ! -- working
    integer, parameter :: wp = dp
    
	! ######## Variables #########    
	integer, parameter					:: N = 50
	real(wp)							:: r				!// m
	real(wp)							:: x0
	real(wp)							:: kT				!// eV
	real(wp)							:: dU				!// eV
	real(wp)							:: L				!// m
	real(wp)							:: alfa				!// asymmetry factor
	real(wp)							:: tau				!// the period of the flash
	real(wp) 							:: zeta, gama, omega !// dynamic viscosity, friction constant, kT/dU, dU/gama/L**2
	real(wp)							:: dt
	real(wp), dimension(0:N-1,0:999999)	:: x
	real(wp)							:: pi = 4.0_wp* atan(1.0_wp)
	integer								:: i, j, res
	integer 							:: measurementPeriod = 10
	real								:: fractionOff = 0.5_wp !3.0_wp / 4.0_wp
	real(wp)			 				:: F_ext
	
	! -- standard values
	real(wp), parameter     :: std_r       	= 12E-9_wp	! m
	real(wp), parameter     :: std_L        = 20E-6_wp	! m
	real(wp), parameter     :: std_dU       = 80.0_wp 	! eV
    real(wp), parameter     :: std_alfa    	= 0.2_wp
    real(wp), parameter     :: std_tau		= 4.0_wp 	
    real(wp), parameter     :: std_zeta     = 1E-3_wp	! Pa*s
    real(wp), parameter     :: std_kT	    = 0.026_wp	! Termisk energi i eV
    real(wp), parameter     :: std_x0       = 0.0_wp
    
    ! -- standard values for the experiment
	!real(wp), parameter     :: std_r       	= 2E-6_wp	! m
	!real(wp), parameter     :: std_L        = 10E-6_wp	! m
	!real(wp), parameter     :: std_dU       = 1250.0_wp	! eV
    !real(wp), parameter     :: std_alfa    	= 0.1_wp
    !real(wp), parameter     :: std_tau		= 4.0_wp 	
    !real(wp), parameter     :: std_zeta     = 1E-3_wp	! Pa*s
    !real(wp), parameter     :: std_kT	    = 0.026		! Termisk energi i eV
    !real(wp), parameter     :: std_x0       = 0.0_wp
    
    
    !// Initialise the random seed
	call init_random_seed	
	
	!// Read in values from the user
	call input_potential()	
	call input_particle(-1)		!// place this inside a loop to get individual radiuses and starting positions for the different particles; set -1 to i.

	!// Converts to reduced units: new t = old t * omega, introduce new x = old x/L where L is  the period of saw-tooth potential, introduce new U = old U/dU where dU is the amplitude of the potential

	!// Criterion for the time step
	call get_dt()
	write(*,*) "The timestep is ", dt
	
	!// The external force due to gravity
	F_ext = 0.0_wp 	!- 4.0_wp/3.0_wp * pi * r**3 * (1.05e3 - 0.9982071e3_wp)*9.81*sin(pi/4.0_wp)
	
		
	!// Start the diffusion process for the N particles
	do i = 0, (N-1)
		gama = 6.0_wp*pi*zeta*r
		omega = dU/(gama*L**2)
		
		!// Set the starting position for the particles
		x(i, 0) = x0

		do j = 1, size(x, dim=2) - 1
			x(i, j) = updateX(x(i, j-1), j*dt, dt)
		end do
		
	end do
	
	!call append_drift_velocity_to_file()
	call write_constants_to_file()
	call write_trajectory_to_file()
	
	call system ('python make_plot.py')
	
contains	

!
! Update the position of the particle using the Euler scheme
real(wp) function updateX(x, t, dt)
	real(wp), intent(in)		:: x, t, dt
	
	updateX = x + F(x, t) * dt + sqrt(2*kT/dU*dt)*random_gauss()
	return
	
end function


!
!Compute the potential U(x, t)
real(wp) function U(xPos, time)
	real(wp), intent(in)		:: time
	real(wp), intent(in)		:: xPos
	real(wp)					:: x
	real(wp)					:: t
	
	x = xPos - floor(xPos)
	t = time - floor(time/(tau*omega))*tau*omega
	
	if (tau /= 0.0_wp .and. t >= 0.0_wp .and. t < fractionOff*tau*omega) then
		U = 0.0_wp
	else
		if (x >= 0.0_wp .and. x < alfa) then
			U = x/alfa
		else if (x >= alfa .and. x < 1.0_wp) then
			U = (1.0_wp-x)/(1.0_wp-alfa)
		end if
	end if
	
	return
end function


!
! Compute the force due to the potential U(x, t)
real(wp) function F(x, t)
	real(wp), intent(in) 		:: x, t
	real(wp)					:: dx = 1.0e-3_wp
	
	if (U(x,t) == 0.0_wp) then
		F = F_ext
	else
		F = - ( U(x + dx, t) - U(x - dx, t) ) / (2.0_wp*dx) + F_ext
	end if
	
	return
end function

!
! Implement the time criterion
subroutine get_dt()
	dt = 1.0e-3_wp
	
	if ( 1/alfa >= 1/(1-alfa) ) then
		do while ( 1/alfa*dt + 4*sqrt(2*kT/dU*dt) > 0.05_wp*alfa )
			dt = dt/2
		end do
	else 
		do while ( 1/(1-alfa)*dt + 4*sqrt(2*kT/dU*dt) > 0.05_wp*alfa )
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
	integer, intent(in)				:: i
	
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
! Draw N numbers with the random_gauss() function. Writes them to file to check
! if it is a normal distribution
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

! Gaussian random number
! Draw Gaussian distributed random numbers. Normal distribution is implemented 
! by the polar Box-Müller algorithm using the uniform distribution from
! random_number().
real(wp) function random_gauss()
	real(wp)		:: rand1, rand2, w
	integer			:: i
	
    w = 0.0_wp
    do while(w == 0.0_wp .OR. w >= 1)
        call random_number(rand1)
        call random_number(rand2)
        rand1 = rand1 * 2.0_wp - 1.0_wp
        rand2 = rand2 * 2.0_wp - 1.0_wp
        w = rand1**2 + rand2**2
    end do

    random_gauss = rand1 * SQRT( - 2.0_wp * LOG(w) / w)
end function

!
! If no input is given the parameter is set to the default parameter
subroutine input_w_default(param, default_param)
	real(wp), intent(out) :: param
	real(wp), intent(in)  :: default_param

	integer, parameter              :: MAX_LENGTH   = 100
	!character(len=*), parameter     :: format_0     = "(A)"
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


!
! Append the drift velocitys to file
subroutine append_drift_velocity_to_file()
	real(wp)	:: drift_velocity
	real		:: std_dev
	logical 	:: exist
	integer		:: i
	if (N >= 25) then 
		inquire(file="drift_velocity.txt", exist=exist)
		if (exist) then
			open(12, file="drift_velocity.txt", status="old", position="append", action="write")
		else
			open(12, file="drift_velocity.txt", status="new", action="write")
		end if
	
		drift_velocity = 0.0_wp
	
		do i = 0, N - 1
			drift_velocity = drift_velocity + (x(i, size(x,dim=2)-1) - x(i, 0)) / size(x,dim=2)
		end do
	
		drift_velocity = drift_velocity / N
		
		std_dev = 0.0_wp
		do i = 0, N-1
			std_dev = std_dev + (drift_velocity - (x(i, size(x,dim=2)-1) - x(i, 0)) / size(x,dim=2))**2
		end do
			std_dev = sqrt(std_dev/(N-1))

		write(12, *) tau, drift_velocity, std_dev
		close(12)
	else
		write(*,*) 'The number of particles is too small to get an accurate measurment of the average drift velocity.'
		write(*,*) 'You should use N >= 25. You have N = ', N
	end if
end subroutine
    
    
!
! Writes the trajectory of the particle(s) to the file "trajectory.txt"
subroutine write_trajectory_to_file()
	open(unit=1,file="trajectory.txt", form="formatted", status="replace", action="write", iostat = res)
	if (res /= 0) then
		write(*,*) "Error in opening file, status", res
		stop
	end if
	
	do i = 0, size(x, dim=1) - 1
		write(1,*) x(i,0:size(x, dim=2)-1:measurementPeriod)		!// x = [[]]
	end do
	close(unit=1)
	
end subroutine

!
! Writes most of the constants to the file "constants.txt"
subroutine write_constants_to_file()
	open(unit=13,file="constants.txt", form="formatted", status="replace", action="write", iostat = res)
	if (res /= 0) then
		write(*,*) "Error in opening file, status", res
		stop
	end if

	write(13,*) size(x, dim=1), size(x, dim=2)/10	!// #of particles, #of time steps
	write(13,*) dt									!// dt
	write(13,*) alfa								!// alfa
	write(13,*) tau									!// tau
	write(13,*) omega								!// omega
	write(13,*) dU									!// dU
	write(13,*) L									!// L
	write(13,*) kT									!// kT
	write(13,*) measurementPeriod
	
	close(unit=13)
end subroutine

!
! Initialise the random seed for the random number generator random_number()
subroutine init_random_seed()
  integer 								:: i, n, clock
  integer, dimension(:), allocatable 	:: seed

  call random_seed(size = n)
  allocate(seed(n))

  call system_clock(count=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)

  deallocate(seed)
end subroutine

end program biased_brownian_motion

