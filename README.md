# FEMRI

This library is an extension of the original FEMRI library by Luca Antiga and David Steinman.

It is able to simulate MRI image acquisition by computing the Fourier transform of a volume bounded by linear or quadratic triangular surface elements using an adaptive numerical quadrature. The number of quadrature points for each of the elements for a desired error is computed using an a priori error estimate. This library can also generate any number of quadrature points by using the [Golub-Welsch](http://web.stanford.edu/class/cme335/spr11/S0025-5718-69-99647-1.pdf) algorithm implemented using VTK functions.

How the library was tested and how it can be used is explained in detail in my Master's [Thesis](https://tspace.library.utoronto.ca/handle/1807/67855).
