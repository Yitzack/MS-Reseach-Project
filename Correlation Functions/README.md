# Correlation Functions

I have here the original pretty printing Mathematica notebooks and an extremely simiplified Mathematica *.m file. The *.m file is here so that you can read the *.nb should you not have access to Mathematica. Mathematica notebooks are impossible to read in plain text.

## Data Collector.*
Because the data collector is recycled for every dataset and used with the notebook, the datafile sources and destinations are hard coded. There's also a point where the vacuum coupling constant get selected and that will have to be changes for each potential. It isn't too far off to change the *.m version to be a usable script on the command line called by `math -script Data\ Collector.m FileName`.

Data Collector collects the imaginary and real contributions to the spectral function and the parameters and outputs all of the data to a Mathematica *.m file. The data exported are spline interpolations. It may also collect previously calculated correlation functions and include those tables in the output.

## Correlation Function.*
Same story as Data Collector.\*. The Mathematica notebook is what I used. The Mathematica \*.m files are for *possible* reading. Mathematica is known as a write-only language and it has hit here particularlly hard on the spatial correaltion functions. Again, the \*.m isn't far from a useable script from the command line. Watch out for backup files writting over each other when multiple instance of the script are ran. Watch out for incomplete output of the writing of files as that will corrupt the files on disk.

Correlation Function.\* takes the file outputted by Data Collector.\* and calculates the momentum-dependant Euclidean correlation function and the spatial correlation function. It periodically saves the data to a backup so that progress isn't lost. When done with the Euclidean correlation functions it will output results to the input file. It does it again when it finishes with the spatial correlation function.

The input file, output file, and Euclidean correlation function tempratures are hard coded. As written, the Euclidean correlation function will calcluate all tempratures simultanously.

## Spatial Correlation.cpp
Mathematica was adding hours to days to the execute time for the space-like contributions to the spatial correlation function. So I wrote this to do it myself and faster.

### Compile directions
g++ Spatial Correlation.cpp \<-O3\> \<-o program_name\>

Optimisation is recommended. Program name of your choice. Nothing fancy.

### Execute directions
program\_name File\_name

The file name should point to a space seperated file that was exported from Mathematica as a CSV file. It contains, in this order, the coupling counstant fraction for the in-medium, the vacuum coupling constant, and 9 wide table of 39,109 control points. The colmns of the table go imaginary part of the continuum, the imaginary part of the constant numerator, the imaginary part of the linear numerator, the imaginary part of the quadratic numerator, the imaginary part of the denominator numerator, the real part of the constant numerator, the real part of the linear numerator, the real part of the quadratic numerator, and the real part of the denominator numerator. It assumes the table is row major order.
