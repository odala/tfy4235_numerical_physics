program to_bit
! Writes the bit contents of integers, single precision real (or float) numbers, double precision real numbers as they are represented on the computer.
	
	!dimension in(32)
	write(*,*) "integer ="
	read(*,*) ia
	write(*,*) ia
	!do i = 1, 32
	!	in(i) = iand(1,ishft(ia, 1-i)
	!end do
	!write(*,1) (in(i), i=32, 1, -1)
	!format(1x, 32(i1))

end program to_bit

