program biased_brownian_motion

	implicit none	!// Prevent the use of implicit declaration
	
	integer							:: nParticle1, nParticle2
	real							:: r1, r2, m1, m2
	real 							:: gama					!// friction constant; 6 pi ny ri
	real							:: kB = 1.3806488e-23 	!// m^2 kg s^-2 K^-1
	real							:: rand
	
	!// Initialiserer et random seed
	call init_random_seed
	
	rand = gaussianRand(0.0, 1.0)
	
	write(*,*) rand

contains


!// introduce new t = old t * omega and new tau = old tau * omega
real function flash(t)
	real, intent(in)		:: t
	real					:: tau 		!// the period of the flash
	
	if (t < 3*tau/4 .and. t >= 0) then
		flash = 0.0
	else if (t < tau .and. t >= 3*tau/4) then
		flash = 1.0
	end if
	
	return
end function


!// introduce new x = old x/L where L is  the period of saw-tooth potential
!// introduce new U = old U/dU where dU is the amplitude of the potential
real function Ur(x)
	real, intent(in)		:: x
	real					:: alfa		!// asymmetry factor
	
	if (x >= 0 .and. x < alfa) then
		Ur = x/alfa
	else if (x >= alfa .and. x < 1) then
		Ur = (1-x)/((1-alfa))
	end if
		
	return
end function


real function U(x, t)
	real, intent(in)		:: x, t
	
	U = Ur(x)*flash(t)
		
	return
end function


real function gaussianRand(mu, sigma)
	real, intent(in)	:: mu
	real, intent(in)	:: sigma
	logical 			:: haveSpare =.false.
	real				:: rand1, rand2, pi = 4.* atan(1.)
 
	if (haveSpare) then
		haveSpare = .false.
		gaussianRand = sigma * sqrt(rand1) * sin(rand2) + mu
		return
	end if
 
	haveSpare = .true.
 
	call random_number(rand1)
	
	if (rand1 < 1e-10) then 
		rand1 = 1e-10
	end if
		
	rand1 = -2 * log(rand1);
	call random_number(rand2)
	rand2 = rand2 * 2 * pi
 	
 	gaussianRand = sigma * sqrt(rand1) * cos(rand2) + mu
 	
	return
end function


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
