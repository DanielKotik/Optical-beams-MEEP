# -*- coding: utf-8 -*-
"""
file:    beamprofile.pxd
brief:   ...
author:  Daniel Kotik
version: 1.5-beta
release date: xx.xx.2020
creation date: 22.02.2020
"""
import cython

# function type declaration for spectrum amplitudes
ctypedef double complex (*f_spectrum_type)(double theta, double phi)

#def (double complex, double, double) complex_dblquad(f_spectrum_type func, double a, double b, double gfun, double hfun)

cdef double f_Gauss_spherical(double sin_theta, double theta, double phi, 
                              double W_y, double k)

cdef double complex f_Laguerre_Gauss_spherical(double sin_theta, double theta, 
                                               double phi, double W_y, int m, 
                                               double k)

#cpdef double complex psi_spherical(r, f_spectrum_type f, x, k)
