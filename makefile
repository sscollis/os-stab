#==============================================================
#
#  Makefile for OS solvers (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 11-10-97
#
#==============================================================
#
# If you have a license to use Numerical Recipes 
# you can activate it with this define.  
#
ifdef USE_NR_ALL
	USE_NR = 1
	USE_NR_ODEINT = 1
	USE_NR_HUNT = 1
endif
ifdef USE_NR
ifdef USE_NR_ODEINT
  DEFINES += -DUSE_NR_ODEINT
endif
ifdef USE_NR_HUNT
  DEFINES += -DUSE_NR_HUNT
endif
	NR_OBJ = nr.o
endif

ifdef USE_RKF45
	DEFINES += -DUSE_RKF45
	RKF45_OBJ = rkf45.o
endif

ifdef USE_LSRK14
	DEFINES += -DUSE_LSRK14
endif

ifdef USE_RKCK45
	DEFINES += -DUSE_RKCK45
endif

ifdef USE_VODE
	DEFINES += -DUSE_VODE
	VODE_OBJ = zvode.o
endif
#
# Compilers
#
CC = gcc-9
CXX = gxx-9
FC = gfortran
F77 = gfortran
#
# Define all objects
#
ALL = bl conte contebl conteucbl orrncbl orrcolchan orrsom orrspace \
			orrfdbl orrfdchan orrucbfs orrucbl orrwong
#
OBJS = $(foreach module, $(ALL), $(module).o) $(NR_OBJ)
#
# Optimization/debug flags
#
OPT = -O2 -ffpe-trap=invalid,zero,overflow
#DEBUG = -g -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal
#
# Compiler flags
#
FFLAGS = -cpp -ffixed-line-length-120 -fdefault-integer-8 -fdefault-real-8 \
         -fdefault-double-8 -std=legacy $(DEFINES) $(OPT) $(DEBUG)
F90FLAGS = -cpp -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 \
					 $(DEFINES) $(OPT) $(DEBUG)
#
# External Libraries
#
LIB = -L$(HOME)/local/OpenBLAS/lib -lopenblas $(NR_OBJ) $(RKF45_OBJ) $(VODE_OBJ)
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 

all:
	$(MAKE) bl 
	$(MAKE) conte 
	$(MAKE) contebl 
	$(MAKE) conteucbl
	$(MAKE) orrncbl
	$(MAKE) orrcolchan
	$(MAKE) orrsom
	$(MAKE) orrspace
	$(MAKE) orrfdbl
	$(MAKE) orrfdchan
	$(MAKE) orrucbfs
	$(MAKE) orrucbl
	$(MAKE) orrwong

contebl: contebl.o $(NR_OBJ) $(RKF45_OBJ) $(VODE_OBJ)
	$(FC) $(LIB) contebl.o -o contebl

conteucbl: conteucbl.o $(NR_OBJ)
	$(FC) $(LIB) conteucbl.o -o conteucbl

orrsom: orrsom.o $(NR_OBJ)
	$(FC) $(LIB) orrsom.o -o orrsom

orrspace: orrspace.o $(NR_OBJ)
	$(FC) $(LIB) orrspace.o -o orrspace

orrncbl: orrncbl.o
	$(FC) $(LIB) orrncbl.o -o orrncbl

orrcolchan: orrcolchan.o
	$(FC) $(LIB) orrcolchan.o -o orrcolchan

orrfdbl: orrfdbl.o
	$(FC) $(LIB) orrfdbl.o -o orrfdbl

orrfdchan: orrfdchan.o
	$(FC) $(LIB) orrfdchan.o -o orrfdchan

orrucbfs: orrucbfs.o
	$(FC) $(LIB) orrucbfs.o -o orrucbfs

orrucbl: orrucbl.o
	$(FC) $(LIB) orrucbl.o -o orrucbl

orrwong: orrwong.o
	$(FC) $(LIB) orrwong.o -o orrwong

clean:
	/bin/rm -fr $(ALL) *.o *.mod *.dSYM

.f90.o:
	 $(FC) $(F90FLAGS) -c $*.f90

.f.o:
	 $(F77) $(FFLAGS) -c $*.f
