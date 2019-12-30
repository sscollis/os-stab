#==============================================================
#
#  Makefile for OS solvers (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 11-10-97
#
#==============================================================
DEBUG = -O2 -ffpe-trap=invalid,zero,overflow
#DEBUG = -g -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow
#,zero,overflow,underflow,denormal
#DEBUG = -g -fbacktrace -ffpe-trap=invalid
FFLAGS = -cpp -ffixed-line-length-120 -fdefault-integer-8 -fdefault-real-8 -std=legacy $(DEBUG)
F90FLAGS = -cpp -fdefault-integer-8 -fdefault-real-8 $(DEBUG)
LIB = -L$(HOME)/local/OpenBLAS/lib -lopenblas
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
ALL = conte contebl orrncbl orrcolchan orrsom bl orrspace conteucbl orrfdbl orrfdchan orrucbfs orrucbl orrwong
#
OBJS = $(foreach module, $(ALL), $(module).o)
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
	/bin/rm -fr $(ALL) $(OBJS) *.mod *.dSYM

.f90.o:
	 $(FC) $(F90FLAGS) -c $*.f90

.f.o:
	 $(F77) $(FFLAGS) -c $*.f
