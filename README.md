# AcousticTube

This repository contains code to simulate wave propagation in acoustic tubes, based on my [PhD Thesis](https://www.researchgate.net/publication/325957710_Physical_Modelling_of_Brass_Instruments_using_Finite-Difference_Time-Domain_Methods). The model encorporates:
* time-domain wave propgation
* variable bore profile
* viscothermal losses using a low-order Foster network approximation to the Zwikker and Kosten model
* radiation losses using low order RLC network approximation to Levine and Schwinger unflanged open tube

## Impulse response
The [ImpulseResponse](ImpulseResponse) directory contains files to generate an impulse response. Example usage is shown in [ExampleUsage.m](ExampleUsage.m)

## Testing
The [Testing](Testing) directory contains code to compare agains the FDTD simulations.

## Development
The source code has been developed in Octave on MacOS, but should work on in Matlab and on other operating systems (although this has not been tested).

## Acknowledgement

If you wish to use this code in a project, please add some acknowledgement.
