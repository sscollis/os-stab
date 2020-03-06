## OS-Stab

Performs stability analsysis for incompressible flows by solving the Orr-Sommerfeld equation 
for planar channel flow and the Blasius boundary layer.  Also includes `bl.f` which is a 
spectral collocation Blasius flow solver and there are also shooting-based Blasius solvers
embedded within the stability solvers (see below).

### Building

Should build on platforms with Gfortran.  Start by typing:

    make help

Default build is:

    make 

Recommended build is:

    make USE_RKCK45=1

You can make the Doxygen documentation using:

    make docs

### Running

The simplest case is temporal stability analysis for a planar channel using
Conte's method and shooting:

    ./conte

Enter "d" for default values 

A similar analysis for the Blasius boundary layer is done using:

    ./contebl

Enter "d" for default values 

or the similar which takes input fron stdin

    ./orrsom < test.inp

Just computing the Blasisus boundary layer alone (no eigenanalysis):

    ./bl < bl.inp

A spatial analysis for BL is done using:

    ./orrsspace 

And then type `space.inp` for input file.

You can do a sweep through the stability curve using:

    ./orrspace

And enter `sweep.inp` as the input file.

#### Additional examples

These examples require that one currently builds with `USE_ALL_NR=1` that 
requires the user to supply commercially licenses code from Numerical Recipes
in a file called `nr.f`.  These routines are currently only used sparingly and
can/should be replaced with open source code as time permits. 

Nevertheless, you must have a valid license to use Numerical Recipes code
and you must now check such code into this repository!

To do a sweep in alpha using finite-differences:

    orrfdchan < fd-128.inp

or a single alpha and printing the eigenfunction:

    orrfdchan < fdchan.inp

NOTE:  There currently are not example inputs for all codes!

### Updates as of 12-30-2019

Updated Orr-Sommerfeld solvers for gfortran 

  Code      |  Description
------------|----------------------------------------------------------
  conte.f   | Solves OS for channel using Godunov-Conte shooting
  contebl.f | Solve Blasius and then OS for boundary layer using Conte
  bl.f      | Spectral collocation solver for Blasius equation
  orrsom.f  | Similar to contebl.f
  shoot.f   | Collection of general routines
  orrspace.f| Solves spatial problem for boundary layers

Numerous updates in this version include:

  1. Fixed bug in eigenfunction values at nstep
  2. Improved orthonormalization criteria
  3. Switch from IMSL to LAPACK (need to test still)
  4. Updated to support gfortran in makefile
  5. Added scripts for gnuplot output
  6. Improved some inline documentation
  7. Use OpenBLAS for numerical linear algebra
  8. Works on Mac OS-X (Darwin)

### Updates as of 7-13-97

Incompressible Orr--Sommerfeld Solvers for Channel flow and Blasius Boundary
Layer.  Also bl.f is a spectral blasius flow solver.

### Updates as of 3-07-2001

I have updated the following codes to use LAPack routines

`orrncbl.f`	This finds the neutral curve for Boundary layer and
channel profiles

`orrcolchan.f`	Solve Orr-Sommerfeld for channel profile

S. Scott Collis
