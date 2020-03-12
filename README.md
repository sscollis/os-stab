## OS-Stab

### Description

Performs stability analsysis for incompressible flows by solving the Orr-Sommerfeld equation 
for planar channel flow, the Blasius boundary layer, or for a given arbitary mean velocity
profile.  Note that in all cases, the parallel-flow approximation is made. 

### Background

Several different solvers are provided using different numerical methods including stablized
shooting, Chebyschev spectral collocation, and finite-difference discretization. Solvers
are also provided for both temporal analysis (real wavenumber, complex wave speed) and
the more realistic spatial analysis (real frequency, complex wavenumber).

Also includes `bl.f` which is a spectral collocation Blasius flow solver and there are 
also shooting-based Blasius solvers embedded within the stability solvers (see below).

As an example, the mean Blasius velocity profile and second derivative of that profile computed
using `contebl.f` is shown below.

![Mean Profile](https://github.com/sscollis/os-stab/blob/master/images/mean.png)

The eigenfunction for the most unstable Tollmien-Schlichting wave (TS-wave) at `Re=5800`, 
`alpha=0.179` is shown in the following figure with an complex valued wave speed 
(eigenvalue) of `c = 0.36412287E+00 + i 0.79597206E-02`.

![Eigenfunction](https://github.com/sscollis/os-stab/blob/master/images/phi.png)

These results are from running `contebl.f` using default input parameters as 
described below.

Note that there are example scripts for making plots of results using 
`gnuplot`.  These scripts have files `*.com` such as `mean.com` for 
plotting the mean flow profile and `phi.com` for plotting the eigenfunction, 
phi.

There are also example scripts for plotting these same results using Python 
and matplotlib.   See `mean.py` and `phi.py` and these are used by typing

    python mean.py

which generates the file `mean.png`.  Similarly for `phi.py`.

### Building

Should build on platforms with Gfortran.  Start by typing:

    make help

Default build is:

    make 

Recommended build is:

    make USE_RKCK45=1

You can make the Doxygen documentation using:

    make docs

and then you can view the documentation locally by pointing your browser at
`html/index.html`. For example, on a Mac you simply type

    open html/index.html

from within the main `os-stab` directory.

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

  Code          |  Description
----------------|---------------------------------------------------------------------------
  `conte.f`     | Solves OS for channel using Godunov-Conte shooting
  `contebl.f`   | Solve Blasius and then OS for boundary layer using Conte
  `bl.f`        | Spectral collocation solver for Blasius equation
  `orrsom.f`    | Similar to `contebl.f`
  `shoot.f`     | Collection of general routines for shooting
  `orrspace.f`  | Solves spatial problem for boundary layers can sweep in frequency
  `orrfdchan.f` | Solves temporal channel problem using finite-differences including sweeps 

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
Layer.  Also `bl.f` is a spectral blasius flow solver.

### Updates as of 3-07-2001

I have updated the following codes to use LAPack routines

`orrncbl.f`	This finds the neutral curve for Boundary layer and
channel profiles

`orrcolchan.f`	Solve Orr-Sommerfeld for channel profile

S. Scott Collis
