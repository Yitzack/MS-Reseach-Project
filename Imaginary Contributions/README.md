Imaginary Contributions
=======================

The files contained here calculate the imaginary contributions to the spectral functions.

Compile Directions
------------------

g++ Spectral.cpp \[-O3\] \[-o program name\] -D \<BB|CC\>= -D ORDER=\<37|97\> -D VERSION=\<22|24|42|EXP\> \[-D HALF=\]

The program will compile correctly at all common levels of optimization, which are optional. So -O, -O1, -O2, -O3 all work and are optional but highly recommended due the large number of loops.

Given the compile options, I recommend using a program name to match the selected options. I would recommend you use a naming convention that matches the output file names found on Spectral.cpp, lines 17-42.

You must selet bottomium or charmonium with the -D BB or -D CC macro.

You must select 37th order or 97th order with the -D ORDER= macro. You must use 97th order for time-like mesons or your results will not be accurate enough for the spatial correlation function. If you don't need the spatial correlation function, 37th order may be good enough. You may use 37th order for space-like mesons.

You must select the potential version with the -D VERSION= macro. I do recommend that you correct the zero/cusp/maximum of the potiential on Spectral.h lines 657, 662, and 722-725 to match the selected potential. I was not exceedingly hindered by them being wrong, but it might be just detectable in the relative error of the spectral function and spatial correlation.

You may select half of the in-medium self-energy with -D HALF=.

Physics
-------

See paper and thesis repositories stored elsewhere.
