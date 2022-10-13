# Real Contributions

The files contained here calculate the real contributions to the spectral functions.

## Compile Directions

g++ ReSpectral.cpp \[-O3\] \[-o program_name\] -D \<BB|CC\>= -D VERSION=\<22|24|42|EXP\> \[-D HALF=\]

The program will compile correctly at all common levels of optimization, which are optional. So -O, -O1, -O2, -O3 all work and are optional but highly recommended due to the large number of loops.

Given the compile options, I recommend using a program name to match the selected options. I would recommend you use a naming convention that matches the output file names found on Spectral.cpp, lines 17-36.

You must selet bottomium or charmonium with the -D BB or -D CC macro.

You must select the potential version with the -D VERSION= macro.

You may select half of the in-medium self-energy with -D HALF=.

## Execute Directions

./program\_name ProcessID Number\_of\_Threads Temprature Fraction\_of\_Coupling\_Constant Debye\_Mass Quark\_Mass Starting\_Point Ending\_Point Momentum\_Scale Fraction\_to\_Vacuum

The output will land in ./data/Spectral\*.Temprature.ProcessID as space seperated values. The output includes error estimates after the mode value, which are comma seperated. To collect all of the output I recommend using bash command `sort -unk2 -unk1 data/Spectral\*.Temprature.\* | sed 's/ /,/g' > data/Spectral\*.Temrature.csv`. This will put mode values and error values in seperate fields.

The number of threads used should be between 1 and 575. I don't recommend going higher as that won't help anything even if you're lying to it. There are only 575 columns and you can use the starting and stopping points to break up the columns. Mod(ProcessID,Number\_of\_Threads) will tell you which processes are working on the same set of columns.

Starting\_Point and Ending\_Point are see tin, the starting and ending points. What that exactly means in terms of invariant mass and center of mass momentum will have to be determined from the code or an appendix in the thesis.

Temprature is an integer between 0 and 4 inclusive that represents Vacuum and T=194, 258, 320, and 400 MeV.

Fraction of coupling constant will tell it that the coupling constant should be the vacuum value times that number for P=0. It will remain there if the fraction to vacuum is 0.

The Debye\_Mass is that number times the actual temprature, not the number to call for that temprature.

Quark\_Mass is see tin, the quark mass.

Fraction\_to\_Vacuum is how far to vacuum the system goes at the infinite momentum. 0 is no change. 1 means the quark mass goes to 1.8 GeV, fraction of coupling constant goes to 1, and the Debye mass becomes 0.

Momentum\_Scale is how fast the system goes the infinite momentum limit. At the momentum of the momentum scale, it is half way there.

## Around.h
This is the uncertain value object. It should be able to substitute for long double in most situtations. You may find that it lacks support for most floating point functions. But it does cover the basic four functions, abs(), and comparison. Not all constant methods are declared as such.

## Elements.h
Basically an array of 4 quantities to be integrated together. Element<> to Element<> comparison is done to the sum of elements. Element<> to scalar comarison requires the set to pass as it is used in integration and all elements need to pass the test. The basic 4 functions and abs() are included for integration. Not all constant methods are declared as such.

## Interpolation.h
This is a 2D interpolation object. It uses a dumbed down version of the NURBS (Non-Uniform Ratio B-Spline) algorithm. It needs an array of control points that are precalculated either by you or Mathematica. Mathematica will provide them with `Interpolation[table,Method->"Spline"][[4,1,5,1]]`. It expects the control points to be evenly spaced in both dimensions. It is a bicubic spline. It does not track the valid range for you. It will report the number of control points. It is up to user to ensure they remain in bounds and how the index value maps to the function input.
