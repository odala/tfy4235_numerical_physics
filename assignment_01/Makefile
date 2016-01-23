# The compiler
FC = gfortran -O3#-fopenmp -O3

# Executable to built
PROGRAM = biased_brownian 

# Source files
SRCS =	common.f90	\
	readinput.f90	\
	main.f90

# Object files (can make a easier version where alle srcs are made :.f90-.o)
OBJS =	common.o	\
	readinput.o	\
	main.o

# 'make' builds all

$(PROGRAM): $(OBJS)
		$(FC) -o $(PROGRAM) $(OBJS)

%.mod: %.o %.f90
		$(FC) -c $<

%.o: %.f90
		$(FC) -c $<
%.o: %.f
		$(FC) -c $<

clean:
		rm -f *.o *.mod *~
