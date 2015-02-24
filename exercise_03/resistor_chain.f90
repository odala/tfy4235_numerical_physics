program resistor_chain

	implicit none	!// Prevent the use of implicit declaration
	
	integer							:: nResistors = 1000
	real, dimension(0:999) 			:: c
	real, dimension(1:999)			:: b
	real, dimension(0:1000)			:: v
	real, dimension(0:999)			:: vnew
	real							:: temp, error = 1.0, pres = 2.e-4
	integer							:: i, j, k, n = 0
	
	!// Genererer 1000 konduktanser
	call init_random_seed
	call random_number(c)
		
	!// Initialiserer b
	b(1) = c(0)*v(0)
	do i = 2, 998
		b(i) = 0
	end do
	b(999) = c(999)*v(1000)
	
	write(*,*) 'JACOBI RELAXATION METHOD'
	call jacobi_relaxation(c, b)
	
	write(*,*) 'GAUSS-SEIDEL RELAXATION METHOD'
	call gauss_seidel_relaxation(c, b)
	
	write(*,*) 'SUCCESSIVE OVER RELAXATION'
	call sor(c,b)

end program resistor_chain

subroutine sor(c, b)
	real, dimension(0:999) 			:: c
	real, dimension(1:999)			:: b
	real, dimension(0:1000)			:: v
	real, dimension(0:999)			:: vnew
	real							:: temp, omega, error = 1.0, pres = 2.e-4
	integer							:: i, j, k, n = 0
	
	!// Leser inn parameteren omega
	write(*,*) 'Omega = '
	read(*,*) omega
	
	!// Initialiserer v
	do i = 0, 1000
		v(i) = float(i)/1000
	end do	
	vnew(0) = 0
	
	!// Oppdaterer v
	do k = 1, 100000
		do i = 1, 999
			vnew(i) = (b(i) - (-c(i-1)*vnew(i-1) - c(i)*v(i+1))) / (c(i-1) + c(i))
		end do
		
		do i = 1, 999
			vnew(i) = omega*vnew(i) + (1 - omega)*v(i)
		end do
		
		do i = 1, 999
		 v(i) = vnew(i)
		end do
	
		error = 0.0
		do i = 1, 999
			error = error + (c(i-1)*v(i-1) + c(i)*v(i+1) - (c(i-1) + c(i))*v(i))**2
		end do
		
		error = sqrt(error)

		n = n + 1
		
		if (error < pres) then
			exit
		end if
	end do
	
	write(*,*) n,error
	
	
end subroutine


subroutine gauss_seidel_relaxation(c, b)
	real, dimension(0:999) 			:: c
	real, dimension(1:999)			:: b
	real, dimension(0:1000)			:: v
	real, dimension(0:999)			:: vnew
	real							:: temp, error = 1.0, pres = 2.e-4
	integer							:: i, j, k, n = 0
	
	!// Initialiserer v
	do i = 0, 1000
		v(i) = float(i)/1000
	end do	
	vnew(0) = 0
	
	!// Oppdaterer v
	do k = 1, 100000
		do i = 1, 999
			vnew(i) = (b(i) - (-c(i-1)*vnew(i-1) - c(i)*v(i+1))) / (c(i-1) + c(i))
		end do
		
		do i = 1, 999
		 v(i) = vnew(i)
		end do
	
		error = 0.0
		do i = 1, 999
			error = error + (c(i-1)*v(i-1) + c(i)*v(i+1) - (c(i-1) + c(i))*v(i))**2
		end do
		
		error = sqrt(error)

		n = n + 1
		
		if (error < pres) then
			exit
		end if
	end do
	
	write(*,*) n,error

end subroutine


subroutine jacobi_relaxation(c, b)
	real, dimension(0:999) 			:: c
	real, dimension(1:999)			:: b
	real, dimension(0:1000)			:: v
	real, dimension(0:999)			:: vnew
	real							:: temp, error = 1.0, pres = 2.e-4
	integer							:: i, j, k, n = 0
	
	!// Initialiserer v
	do i = 0, 1000
		v(i) = float(i)/1000
	end do
	
	!// Oppdaterer v
	do k = 1, 150000
		do i = 1, 999
			vnew(i) = (b(i) - (-c(i-1)*v(i-1) - c(i)*v(i+1))) / (c(i-1) + c(i))
		end do
		
		do i = 1, 999
		 v(i) = vnew(i)
		end do
	
		error = 0.0
		do i = 1, 999
			error = error + (c(i-1)*v(i-1) + c(i)*v(i+1) - (c(i-1) + c(i))*v(i))**2
		end do
		
		error = sqrt(error)

		n = n + 1
		
		if (error < pres) then
			exit
		end if
	end do
	
	write(*,*) n,error

end subroutine


subroutine init_random_seed()
  integer 							:: i, n, clock
  integer, dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))

  call system_clock(count=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)

  deallocate(seed)
end subroutine
