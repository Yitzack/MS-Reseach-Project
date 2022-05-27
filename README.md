# Momentum Dependence of Charmonium Spectral Function in Quark-Gluon Plasma

This is the code base for my MS and incomplete PhD in Theortical Nuclear Physics. It was on the momentum dependence of the charmoinum spectral function in the QGP. It uses _T_-matrix formulaism to calculate the spectral function. The potential is seperable, which means it can analytically rearanged to be numerically tractible with a solved expression. The potential is also Lorentz invariant which means that the spectral function has no momentum dependence in vacuum. In the QGP, the medium will define a prefered rest frame and break the Lorentz invariance through the single quark self-energy and possibly the potential.

## Abandoned Code

Most of the code here no longer has a valid input and so can no longer perform its function. One script is completely unknown to me after 8 years. Euclidean1D.cpp still works but is only useful for parameter searches at P=0.

## Correlation Functions

These are Mathematica scipts that collect up the data of the spectral function and calculate the momentum-dependent Euclidean corralation function and the spatial correlation function

## Imaginary Contributions

This is where most of the compute time is used. This program will calculate the imaginary contributions to the spectral function. It is usually a 3-d integral. It could be 15 minutes to an hour to calculate a point. Best to have a cluster or super computer. I used about 8 million CPU-hours running this program.

## Real Contributions

The Mathematica script in here will take the imaginary contributions and calculate the real contributions useing a dispersion relation.

## Schematic Spectral Function Search

The programs in here attempt to find a best fit solution to the lattice QCD calcluated spatial correlation ratio and momentum-depenant Euclidean correlation ratio. It kept pointing to physcially unreasonable results and had difficulty reproducing existing results from the above microscopic calculations.
