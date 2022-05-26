Schematic Spectral Function Search
==================================

The files contained here attempt to locate the best fit schematic function to achieve the lQCD spatial correlation function and Euclidean correlation function difference.

Optimiser.cpp is looking for a best fit that also has all of its parameters on a linear relation with temprature. Optimiser T.cpp is looking for the best for each temprature individually.

Compile Directions
------------------

g++ Optimiser.cpp \[-O3\] \[-o program_name\] -D \<Ps|sP\>= -D ORDER=\<16|37|64|97\>
g++ Optimiser\ T.cpp \[-O3\] \[-o program_name\] -D \<Ps|sP\>= -D ORDER=\<16|37|64|97\>

The program will compile correctly at all common levels of optimization, which are optional. So -O, -O1, -O2, -O3 all work and are optional and highly recommended.

You must select the compile the order of integration for the full spatial correlation function of either P then s with -D Ps= or s then P with -D sP=

You must select 16th, 37th, 64th, or 97th order Gauss-Kronrod integration with the -D ORDER= macro. I do recommend that you modify the program so that either you have one option permantly selected or you break up the order option for integral function. Before put the options in all of the integrating functions, I had 37th order for first level of the spatial (PInt, sInt), Lorentz, and P0_Int; 64th order for second level of the spatial (Spatial); and 97th order for Euclidean.

Execute Directions
------------------

Optimiser T.cpp has 5 modes.

Optimiser.cpp has 2 modes. The first the default behavior of searching for the best on a linear relation with the default bounding box. The second has the bounding box overwritten by the call from the command line.

First mode: ./Optimiser ProcessID

Second mode: ./Optimiser ProcessID \<32 floating point numbers\>

ProcessID is an integer. Output is sent to ./data/Optimiser_Output.ProcessID.csv

The 32 numbers are 16 pairs of lower boundary and upper boundary. You might need to futz with the schematic spectral function to understand what I mean by some of the names.
1. Deviation from J/Psi amplitude boundary for all T
2. Momentum scale of 1) for all T
3. Deviation from T=194 J/Psi mass
4. Deviation from T=400 J/Psi mass
5. Momentum scale of 3) and 4) for all T
6. Boundary for J/Psi width for all T
7. Boundary of momentum scale for 6)
8. Boundary for low energy tail width rate for all T
9. Momentum scale for 8)
10. Boundary for transition interval from low energy tail to constant high energy width for all T
11. Momentum scale for 10)
12. Deviation from T=194 quark mass
13. Deviation from T=400 quark mass
14. Momentum scale for 12) and 13)
15. Boundary for continuum power for all T
16. Momentum scale for 15)
