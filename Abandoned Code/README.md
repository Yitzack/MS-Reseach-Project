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

## Dispersion.cpp

Calculates the real contributions by dispersion relation. There should be a bicubic and bilinear interpolation given f, df/dx, df/dy, and d^2f/dxdy in here.

## Fourier.cpp

Calculates the spatial correlation function. There should be a bicubic and bilinear interpolation given f, df/dx, df/dy, and d^2f/dxdy in here too. It problably uses my oldest technique for calculating the spatial correlation: keep going until it completes at something like 2001*pi/z (1-25 TeV) and recycle the last wavelength of the spectral function.

## Epsilon.m

I have no clue. It has not been altered since it was created on 8/27/2014, 8 years ago. I think it had something to do with an analtyic representation of the interacting spectral function and I thought I could get the width out of something in some way or another and make something useful out of it. It was probably already dead when it got committed.
