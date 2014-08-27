MS-Reseach-Project
==================

This is the associated code for my MS reseach. It is on the momentum dependance of the T-matrix and spectral function of massive quarkonia.

The T-matrix of a system can be used to find the cross-section of an interaction. The interaction can then be used in simulations to test the predictions of the input against experimental results derived from RHIC, LHC, and other high energy, heavy ion collosions. I'm not looking into simulations. I'm comparing my results against the cream of computational results: lattice QCD.

Theory
------
T(P,E,k,k')=V(P,k,k')+int V(P,k,q)G(P,E,q)T(P,E,q,k')d^4q (Bethe-Salpeter Equation)
V=g(Lambda/(Lambda^2+k_mu k^mu))(Lambda/(Lambda^2+k'_mu k'^mu)) (Seperable potiential)
G(P,E,q) (Eq. 16, Physics Review 82, Quarkonia and heavy-quark relaxation times in the quark-gluon plasma", BbS propagator)
Because the BbS propagator contains a delta function, I can analytically evaluate the zeroth (time) dimention and turn it into the Lippman-Schwinger Equation. The seperable potiential allows me to rearrange the equation and explicityly solve the equation:
T(P,E,k,k')=V(P,k,k')/(1-int V(P,q,q)G(P,E,q)d^3q)


Now that I have the T-matrix, I can find the spectral function:
sigma(E,P)=int G(P,E,q)d^3q+int G(P,E,q) d^3q int G(P,E,q')T(P,E,q,q') d^3q'
Here, again the seperable potiential can be used to simplify things.
sigma(E,P)=int G(P,E,q)d^3q+(int G(P,E,q)V(P,k) d^3q)^2/(1-int V(P,q,q)G(P,E,q)d^3q)

This is all well and good, but now I need to compare it lattice. To do that, I need to calculate the spatial correlator. I could caluate the Euclidiean-time correlator, but that would not examine the momentum depenance of the system.
G(z)=int domega/omega int_-infty^infty dP e^(iPz)sigma(omega,P)
I don't have sigma(omega,P), I have sigma(E,P). Sigma is even symmetric in P. With these properties and E^2=s=omega^2-P^2, I can rewrite the correlator:
G(z)=int dE int_0^infty dP 4Ecos(Pz)/(E^2+P^2)sigma(E,P)

Current State
-------------

I've written a more or less satisfactory program that will calculate the T-matrix and spectral function and save them to a file. I've written another program with will readin from file a Spectral function and the first derivatives with respect to the E and P and the second derivative with respect to E and P. The derivatives are calcuated by Mathematica.
Spectral.cpp is the program that calculates the T-matrix and spectral function.
Fourier.cpp is the program that calculates the correlator.
Spectral.h is a common header file. The reason for this is because the Fourier.cpp has in past needed approximations that are held there. It may need access to them in future.
Derive.m is the Mathematica script that calculates the derivitives that Fourier.cpp needs.
Epsilon.m calculates the parameters of the approximation that Fourier.cpp may use.
