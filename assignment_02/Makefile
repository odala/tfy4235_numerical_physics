# The compiler
FC = gfortran
FFLAGS =  -O3#-fopenmp -O3

# Executable to build
PROGRAM = aurora 

# Linked libraries -lacml -lacml_mv -lblas -ltmglib - lbrefblas -llapack
LIBSPATH = /usr/bin/
LIBS	 = -L${LIBPATH} -llapack

# Source files
SRCS =	common.f90	\
	readinput.f90	\
	functionHolder.f90\
	integration.f90 \
	iterativeMethods.f90\
	functions.f90	\
	main.f90

# Object files
OBJS =	common.o	\
	readinput.o	\
	functionHolder.o\
	iterativeMethods.o\
	integration.o 	\
	functions.o	\
	main.o

# 'make' builds all

$(PROGRAM): $(OBJS)
		$(FC) $(FFLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)

%.mod: %.o %.f90
		$(FC) $(FFLAGS) -c $< $(LIBS)

%.o: %.f90
		$(FC) $(FFLAGS) -c $< $(LIBS)
%.o: %.f
		$(FC) $(FFLAGS) -c $< $(LIBS)

clean:
		rm -f *.o *.mod *~
