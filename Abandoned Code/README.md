# Abandoned Code

"All hope abandon ye who enter here." - Inferno by Dante Alighieri, translated by Henry Francis Cary

All right, it isn't quite that bleak in here. I just haven't touched most of this code in something like 5 years as of 5/26/2022 and some of it is extremely over thought for its purpose. The XMLScipt? series converts space seperated value files into XML files for Mathematica to import. I could have converted space seperated value files to comma seperated value files and Mathematica would have imported those without a second thought.

The source code that produced something useful for this has been altered several times over. Basically, it no long has any valid input files, thus abandoned code. However, the bicubic and bilinear interpolation algorthms in here might be useful for recycling at some point.

## Euclidean1D.cpp

As I read through the code, it doesn't appear to be lost like the other code in here. That is why I rescued it from the dust bin. It appears to be here to do parameter searches for good potentials that corrispond to good P=0 Euclidean correlation functions.

It compiles with `g++ Euclidean1D.cpp`. I recommend a program name option and optimiser option. It takes in data from "data/Tables/SpectralccArg1.Arg2.Arg3.csv". It sends output to "data/Tables/SpectralccArg1.Arg2.Arg3.m". It writes out a Mathematica file where the first element is a table of the data for the various contributions to the spectral functions and the next five elements are tables of the Euclidean correlation function for each part of the spectral function.

Call it as `./program_name Arg1 Arg2 Arg3`.

## XMLScript

This is the key to most of these source codes in here as it calls most of them.

sh XMLScript FileName

XMLScript is a bash script. It does several steps:
1. Calls XMLScript1-3 to convert a space seperated value (SSV) file to an XML file for Mathematica using awk.
2. Call Mathematica script Derive.m to calulate derivatives of the functions and store them in a CSV.
3. Call sed to strip out the commas from the CSV to make an SSV.
4. Call Dispersion.cpp and calculate the real contributions and storing them in an SSV.
5. Use paste and cut to get everyting into one file.
6. Call on XMLscript4-6 to convert the SSV files to an XML file for Mathematica using awk.
7. Call Derive2.m repeats Derive.m except it is for the real contributions instead.
8. Use sed to turn CSV into SSV
9. Call Fourier.cpp that calculates the spatial correlation function for the provisioned spectral function.

## XMLScript?

Awk scripts to turn SSV files into XML files for Mathematica to import.

## Derive.m and Derive2.m

Mathematica scripts to calculate the first derivative with respect to s and P or i and j and the second derivative with respect to s and P or i and j. Store the numbers in a CSV. These derivatives and the value of the functions themselves are useful to my C++ program so I don't have have to calculate the derivatives for a spline interpolation.

## ReSpectral.m

Mathematica used to be fast enough on the Cyclotron Institute's cluster to run this. The College of Science only has 40 licenses to share now. That's not enough to keep the cluster going. You'll probably want to convert this to C++. Also C++ usually runs circles around Mathmatica anyways. Its a question of if you want to spend your time developing correct code or do you want spend your time waiting on correct code to complete.

This is not the proper way calculate the real contributions of anything. You'll need evaluate the dispersion relation of the imaginary propagator. That will be the first integral you do.

### Execute Directions

math -script ProcessID Temprature Version

It will then import "~/data/Tables/Spectralcc24.97.P\_Dep4.Temprature\_Version.csv". It appears to be configured to examine spectral functions with momentum dependence at the moment. This behavior can be altered on line 4.

ProcessID can be any number for which there exists a row of data in the first column of the input file.

Once done, it will export the data to "~/data/ReSpectralcc24.97.P\_Dep4.Temprature\_Version.ProcessID.csv". This behavior can be altered on line 86.

## Dispersion.cpp

Calculates the real contributions by dispersion relation. There should be a bicubic and bilinear interpolation given f, df/dx, df/dy, and d^2f/dxdy in here.

## Fourier.cpp

Calculates the spatial correlation function. There should be a bicubic and bilinear interpolation given f, df/dx, df/dy, and d^2f/dxdy in here too. It problably uses my oldest technique for calculating the spatial correlation: keep going until it completes at something like 2001*pi/z (1-25 TeV) and recycle the last wavelength of the spectral function.

## Epsilon.m

I have no clue. It has not been altered since it was created on 8/27/2014, 8 years ago. I think it had something to do with an analtyic representation of the interacting spectral function and I thought I could get the width out of something in some way or another and make something useful out of it. It was probably already dead when it got committed.

## Spatial Correlation.cpp
Mathematica was adding hours to days to the execute time for the space-like contributions to the spatial correlation function. So I wrote this to do it myself and faster.

It is faster order for order, but it is not better. As far as I can tell, there is nothing interesting in the space-like region. Interesting meaning that extra attention is called for in one area or another. But I hit it with my biggest hammer on the smallest reasonable intervals and I still get garbage. I hit it with my 97th order integral on regions that are identical to interpolation cells and I got answers very different from Mathematica. What will probably be called for an adaptive alogrithm similar to the one used by Mathematica. You have the quadrature rule it uses here (16th order). What you need is a tree of sub-intervals with a priority queue of pointers pointing at the leaves sorted by absolute estimated error. Sub-divide the biggest contributor to the error. Each sub-interval member of the tree should have the value and error estimate when it was evaluated, the value and error estimate total of its leaves, and pointers to its children sub-intervals or NULL. This should make it fairly easy to find "NIntegrate::slwcon" conditions where subdividing increases the estimated error instead of reducing it. It should make it fairly easy to find the other error states that Mathematica reports.

For your consideration of runtime:
- 7th order, ~6 minutes
- 16th order, ~30 minutes
- 37th order, ~2.75 hours
- 64th order, ~8.9 hours
- 97th order, ~20.2 hours

### Compile directions
g++ Spatial Correlation.cpp \[-O3\] \[-o program_name\] -D ORDER=\<7|16|37|64|97\>

Optimisation is recommended. Program name of your choice.

The order of integration has a compile time option of 7th, 16th, 37th, 64th, or 97th order. You must pick one.

### Execute directions
program\_name File\_name

The file name should point to a space seperated file that was exported from Mathematica as a CSV file. It contains, in this order, the coupling counstant fraction for the in-medium, the vacuum coupling constant, and 9 wide table of 39,109 control points. The colmns of the table go imaginary part of the continuum, the imaginary part of the constant numerator, the imaginary part of the linear numerator, the imaginary part of the quadratic numerator, the imaginary part of the denominator numerator, the real part of the constant numerator, the real part of the linear numerator, the real part of the quadratic numerator, and the real part of the denominator numerator. It assumes the table is row major order.
