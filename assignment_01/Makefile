PROG =	{a1.f90

SRCS =	a1.f90

OBJS =	a1.o

C           = gcc
CFLAGS      = -O
CLIBS       = 
CLDFLAGS    = -s -L$(HOME)/C/Lib/

CPP         = g++
CPPFLAGS    = -O -I. -I$(HOME)/CC/include/
CPPLIBS     = 
CPPLDFLAGS  = -s -L$(HOME)/CC/Lib/

F77         = pgf77
F77LAGS     = -fast
F77LIBS     = -llapack -lpgplot -lU77 
F77LDFLAGS  =  -L$(HOME)/Fortran/Lib/

F90         = gfortran
F90FLAGS    = -Wall -Wextra
F90LIBS     =  
#F90LIBS     = -lmod -lnr -llapack90 -llapack -lpgplot -lU77 
F90LDFLAGS  = -s -L$(HOME)/Fortran/Lib/ 

LIBS        =  $(F90LIBS) 	

LDFLAGS     =  $(F90LDFLAGS) 	

ETAGS       = etags
DEPEND      = makedepf90
SVNVERSION  = svnversion



all: etags depend $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod  TAGS VERSION .depend

etags:
	$(ETAGS) $(SRCS) *.f90

depend .depend:
	rm -f .depend
	$(DEPEND) $(SRCS) > .depend 

version:
	$(SVNVERSION) > VERSION 


.SUFFIXES: $(SUFFIXES) .f90 .F90 .mod

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.f.o:
	$(F90) $(F90FLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c  $<

.cc.o:
	$(CPP) $(CPPFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<

# Override the modula2 default built-in rule! 
#   (without it, make will try to execute m2c .....) 
%.o:%.mod 


# ... Dependencies
# ......................

# Include the dependency list (created by $(DEPEND))
include .depend

