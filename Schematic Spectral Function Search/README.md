# Schematic Spectral Function Search

The files contained here attempt to locate the best fit schematic function to achieve the lQCD spatial correlation function and Euclidean correlation function difference.

Optimiser.cpp is looking for a best fit that also has all of its parameters on a linear relation with temprature. Optimiser T.cpp is looking for the best for each temprature individually.

## Compile Directions

g++ Optimiser.cpp \[-O3\] \[-o program_name\] -D \<Ps|sP\>= -D ORDER=\<16|37|64|97\>

g++ Optimiser\ T.cpp \[-O3\] \[-o program\_name\] -D \<Ps|sP\>= -D ORDER=\<16|37|64|97\>

The program will compile correctly at all common levels of optimization, which are optional. So -O, -O1, -O2, -O3 all work and are optional and highly recommended.

You must select the compile the order of integration for the full spatial correlation function of either P then s with -D Ps= or s then P with -D sP=

You must select 16th, 37th, 64th, or 97th order Gauss-Kronrod integration with the -D ORDER= macro. I do recommend that you modify the program so that either you have one option permantly selected or you break up the order option for each integral function. Before I put the ORDER option in all of the integrating functions, I had 37th order for the first level of the spatial (PInt, sInt), Lorentz, and P0_Int; 64th order for the second level of the spatial (Spatial); and 97th order for the Euclidean.

## Execute Directions

### Optimiser T.cpp

Optimiser T.cpp has 5 modes.
1. Default behavior
2. Evaluate list of parameters in a file
3. API hook
4. Default behavior but in different bounding box
5. Calculated list of parameters

First mode: ./program\_name Temprature ProcessID

Second mode: ./program\_name Temprature ProcessID InputFile

Third mode: ./program\_name Temprature \<14 floating point numbers\>

Fourth mode: ./program\_name Temprature ProcessID \<28 floating point numbers\>

Fifth mode: ./program\_name Temprature ProcessID Number\_of\_Threads Thread_number \<30 floating point numbers\>

Temprature is an integer between 1 and 4 inclusive that represents T=194, 258, 320, and 400 MeV

The first and fourth modes are best fit searches. The first mode uses a default bounding box for the search. The fourth mode overwrites that bounding box. The output is sent to ./data/Optimiser_Output.ProcessID.Temprature.csv. If ProcessID is low enough, it will start a random walk search from the last best known position which it has hardcoded. If ProcessID is high enough, it will start at a random position somewhere in the bounding box. It needs 14 pairs of upper and lower bounds.

The second and third modes are roughly the same except the API hook can only do one at a time and is recorded in ./data/Optimiser\_Output.API.Temprature.csv and the file input can do it in bulk and records its result to ./data/Optimiser\_Output.ProcessID.Temprature.csv.

The 14 parameters for the second and third modes or 14 pairs of bounds in the forth mode are:
1. J/Psi amplitude at infinte momentum
2. Momentum scale for 1)
3. J/Psi mass at infinte momentum
4. Momentum scale for 3)
5. J/Psi width at infinte momentum
6. Momentum scale for 5)
7. Low energy width rate at infinte momentum
8. Momentum scale for 7)
9. Transition interval in energy from low energy exponential to high energy constant at infinte momentum
10. Momentum scale for 9)
11. Quark mass at infinte momentum
12. Momentum scale for 11)
13. Continuum power at infinte momentum
14. Momentum scale for 13)

You might need to futz with the schematic spectral function to understand what I mean by some of the names.

The fifth mode is looking for 10 sets of starting points, stopping points, and steps. The data is recorded in ./data/Optimiser_Output.ProcessID.Temprature.csv. I was doing an array search so I could do an interpolation. I'm not going to say what the 10 sets of numbers are as you might want to examine some other set parameters than the 10 I had set out.

### Optimiser.cpp

Optimiser.cpp has 2 modes. The first the default behavior of searching for the best on a linear relation with the default bounding box. The second has the bounding box overwritten by the call from the command line.

First mode: ./program\_name ProcessID

Second mode: ./program\_name ProcessID \<32 floating point numbers\>

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
