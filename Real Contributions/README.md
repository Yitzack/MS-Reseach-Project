# Real Contribtions

Mathematica used to be fast enough on the Cyclotron Institute's cluster to run this. The College of Science only has 40 licenses to share now. That's not enough to keep the cluster going. You'll probably want to convert this to C++. Also C++ usually runs circles around Mathmatica anyways. Its a question of if you want to spend your time developing correct code or do you want spend your time waiting on correct code to complete.

This is not the proper way calculate the real contributions of anything. You'll need evaluate the dispersion relation of the imaginary propagator. That will be the first integral you do.

## Execute Directions

math -script ProcessID Temprature Version

It will then import "~/data/Tables/Spectralcc24.97.P\_Dep4.Temprature\_Version.csv". It appears to be configured to examine spectral functions with momentum dependence at the moment. This behavior can be altered on line 4.

ProcessID can be any number for which there exists a row of data in the first column of the input file.

Once done, it will export the data to "~/data/ReSpectralcc24.97.P\_Dep4.Temprature\_Version.ProcessID.csv". This behavior can be altered on line 86.
