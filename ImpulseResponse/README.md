# Readme

This directory contains files to generate an impulse response using FDTD algorithm.

[MagFilter_0p01.m](MagFilter_0p01.m)
The script contains filter coefficients to model viscothermal losses using the Foster network representation. The filter coefficients are set up so they can be used for different tube radii.

[ThermConstants.m](ThermConstants.m)
This function returns the thermodynamic constants for a given temperature.

[TubeImpulseResponse.m](TubeImpulseResponse.m)
This function generates an impulse response for a given acoustic tube profile. The tube is assumed to be closed at the input and open at the opposite end.
The tube is excited using a volume source, where the first element has a value of 1 and the rest is 0.
The function returns the pressure value at the input of the instrument. The input impedance can be calculated by taking the fft of the output of this function.
