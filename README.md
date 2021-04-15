This MATLAB code can be used to model the spin dynamics of spin-3/2 nuclei using irreducible spherical tensor operators (ISTOs).

Code for simulation of 23Na Magnetic Resonance Fingerprinting (MRF) sequences is presented. Both a CPU and a GPU-accelerated version are available.


Thanks to Arthur Magill, on whose work this code is based on: https://github.com/arthurmagill/isto_sim


The model used is described in:

A model for the dynamics of spins 3/2 in biological media: signal loss during radiofrequency excitation in triple-quantum-filtered sodium MRI
I. Hancua, J. R. C. van der Maarelb and F. E. Boada
JMR 147, 2000, p179-191
http://www.ncbi.nlm.nih.gov/pubmed/11097809

The longitudinal recovery modelled is described here:

Sodium Imaging Optimization Under Specific Absorption Rate Constraint
Robert Stobbe and Christian Beaulieu
Magnetic Resonance in Medicine 59:345–355 (2008)
http://onlinelibrary.wiley.com/doi/10.1002/mrm.21468/abstract

Further reading:

Thermal relaxation and coherence dynamics of spin 3/2
I. Static and fluctuating quadrupolar interactions in the multipole basis
J. R. C. van der Maarel
Concepts Magn Resn A 19A(2), 2003, pp97-116
http://onlinelibrary.wiley.com/doi/10.1002/cmr.a.10087/abstract

Fabian Kratzer, March 2021
f.kratzer@dkfz-heidelberg.de