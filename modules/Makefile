# Path to the moduales
subroutines=./subroutines

LDFLAGS	= -llapack -lblas -lm -lc

FC=gfortran -fbounds-check -march=native -cpp -ffixed-line-length-none #-Wall -Wextra -Wconversion

all:master
master:master_obj
	f2py -c --fcompiler=gfortran $(LDFLAGS) *.o -m master master.f90
	make clean

master_obj:
	$(FC) -c $(F90flags) $(subroutines)/wigner/*.f 	

clean:
	$(RM) -r *.o *.mod *~ *dSYM *.pyc

cleanall:
	$(RM) *.o *.mod *.so *.pyf *~
