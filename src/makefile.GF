#-----------------------------------------------------------------
#
# LEDdiag
#
#-----------------------------------------------------------------
# exact diagonalization of 2D and 3D extended Lieb models
# see https://github.com/DisQS/LiebExactDiag
#---------------------------------------------------------------
# SCRTP desktop:
# module load GCC/11.3.0 OpenMPI/4.1.4 SuiteSparse/5.13.0-METIS-5.1.0

F77 = gfortran -ffixed-form -std=legacy -fimplicit-none
F90 = gfortran -std=gnu
#F90 = gfortran

F77FLAGS = -fbounds-check -Warray-bounds  # For debugging
#F77FLAGS = -fno-range-check -Warray-bounds  # For optimizing
#F77FLAGS = # For profiling

#F90FLAGS = -g -std=legacy -fbounds-check -fno-range-check -fimplicit-none -Warray-bounds -fallow-invalid-boz -D__UNDERSCORE__ -fcheck=all -march=core-avx2 #-Wall # For debugging #-Wall # For debugging
F90FLAGS = -ffree-line-length-none -fbounds-check -mcmodel=large -std=legacy -fno-range-check -fimplicit-none -fcray-pointer -fexternal-blas -D__UNDERSCORE__ -march=core-avx2 # For optimizing

CC = gcc
CCFLAGS = -O2

LIBDIR = $(HOME)/f77/lib
LIBFLAGS  = -mkl #-lblas -llapack #-lg2c
#LIBFLAGS  = -llapack -lblas 

PRINT1 = a2ps -g -o-
PRINT2 = > LEDdiag.ps # | fleet -odouble

LEDdiag.GF:    CommonModules.o LEDModules.o random.o main_tmp.o util.o inout.o etime.o makefile.GF
	$(F90) $(F90FLAGS) -o LEDdiag.GF CommonModules.o LEDModules.o\
                util.o main_tmp.o inout.o random.o etime.o $(LIBFLAGS)

CommonModules.o:      CommonModules.f90
	$(F90) $(F90FLAGS) -c CommonModules.f90

LEDModules.o:     LEDModules.f90
	$(F90) $(F90FLAGS) -c LEDModules.f90

random.o:     random.f90
	$(F90) $(F90FLAGS) -c random.f90

main_tmp.o:     main.f90
	sed "s/GITVERSION/`git describe --tags --long`/g" main.f90 | \
	sed "s/GITBRANCH/`git branch --show-current`/g" | \
	sed "s/COMPILED/`date`/g" >main_tmp.f90
	$(F90) $(F90FLAGS) -c main_tmp.f90 -cpp -Dgit=true

util.o:     util.f90
	$(F90) $(F90FLAGS) -c util.f90

inout.o:     inout.f90
	$(F90) $(F90FLAGS) -c inout.f90 -cpp

etime.o:     etime.f90
	$(F90) $(F90FLAGS) -c etime.f90

tar:	LEDdiag.tar	
LEDdiag.tar:   makefile.* main.f90 util.f90 inout.f90 \
                random.f90 LEDdiag.inp
	tar -cvf LEDdiag.tar makefile.* main.f90 util.f90 \
                inout.f90 random.f90 LEDdiag.inp

print:  
	$(PRINT1) LEDdiag.inp main.f90 util.f90 inout.f90 random.f90 $(PRINT2)

clean:
	rm -f core *.mod *.o
