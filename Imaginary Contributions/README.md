# Imaginary Contributions

The files contained here calculate the imaginary contributions to the spectral functions.

## Compile Directions

g++ Spectral.cpp \[-O3\] \[-o program_name\] -D \<BB|CC\>= -D ORDER=\<37|97\> -D VERSION=\<22|24|42|EXP\> \[-D HALF=\]

The program will compile correctly at all common levels of optimization, which are optional. So -O, -O1, -O2, -O3 all work and are optional but highly recommended due the large number of loops.

Given the compile options, I recommend using a program name to match the selected options. I would recommend you use a naming convention that matches the output file names found on Spectral.cpp, lines 17-42.

You must selet bottomium or charmonium with the -D BB or -D CC macro.

You must select 37th order or 97th order with the -D ORDER= macro. You must use 97th order for time-like mesons or your results will not be accurate enough for the spatial correlation function. If you don't need the spatial correlation function, 37th order may be good enough. You may use 37th order for space-like mesons.

You must select the potential version with the -D VERSION= macro. I do recommend that you correct the zero/cusp/maximum of the potiential on Spectral.h lines 657, 662, and 722-725 to match the selected potential. I was not exceedingly hindered by them being wrong, but it is probably detectable in the relative error of the spectral function and spatial correlation. These are set for the 42 VERSION. It is the only potential that is impecable up to 1.2 fm, the others start fail around 1-1.1 fm. The 22 VERSION also does fairly well without correct examination.

You may select half of the in-medium self-energy with -D HALF=.

## Execute Directions

./program\_name ProcessID Number\_of\_Threads Temprature Fraction\_of\_Coupling\_Constant Debye\_Mass Quark\_Mass Starting\_Point Ending\_Point Momentum\_Scale Fraction\_to\_Vacuum

The output will land in ./data/Spectral\*.Temprature.ProcessID as space seperated values. To collect all of the output I recommend using bash command `sort -unk2 -unk1 data/Spectral\*.Temprature.\* | sed 's/ /,/g' > data/Spectral\*.Temrature.csv`.

The number of threads used should be between 1 and 616. I don't recommend going higher as that won't help anything even if you're lying to it. There are only 616 columns and you can use the starting and stopping points to break up the columns. Mod(ProcessID,Number\_of\_Threads) will tell you which processes are working on the same set of columns.

Starting\_Point and Ending\_Point are see tin, the starting and ending points. What the exactly means in terms of invariant mass and center of mass momentum will have to reference from the code or an appendix in the thesis.

Temprature is an integer between 1 and 4 inclusive that represents T=194, 258, 320, and 400 MeV.

Fraction of coupling constant will tell it that the coupling constant should be the vacuum value times that number for P=0. It will remain there if the fraction to vacuum is 0.

The Debye\_Mass is that number times the actual temprature, not the number to call for that temprature.

Quark\_Mass is see tin, the quark mass.

Fraction\_to\_Vacuum is how far to vacuum the system goes at infinite momentum. 0 is no change. 1 means the quark mass goes to 1.8 GeV, fraction of coupling constant goes to 1, and the Debye mass becomes 0.

Momentum\_Scale is how fast the system goes the infinite momentum limit. At the momentum of the momentum scale, it is half way there.

## Physics

See paper and thesis repositories stored elsewhere.
