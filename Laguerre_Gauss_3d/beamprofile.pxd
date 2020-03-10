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

# declare C functions as "cpdef" to export them to the module
cdef extern from "stdlib.h":
    cpdef int abs (int n)

cdef extern from "math.h":
    cpdef double sin(double x)
    cpdef double cos(double x)
    cpdef double exp(double x)

cdef extern from "complex.h":
    cpdef double complex cexp(double complex z)

# function type declaration for spectrum amplitudes
ctypedef double complex (*f_spectrum_type)(double sin_theta, double theta, double phi, dict params)
ctypedef double complex (*integrand_type)(double theta, double phi)

# function prototypes
cdef double real_func(double x, double y, func)
cdef double imag_func(double x, double y, func)

@cython.locals(real=cython.double, imag=cython.double, real_tol=cython.double, imag_tol=cython.double)
cdef (double complex, double, double) complex_dblquad(func, double a, double b, double gfun, double hfun)

@cython.locals(W_y=cython.double, k=cython.double)
cdef double complex f_Gauss_spherical(double sin_theta, double theta, double phi,
                                      dict params)

@cython.locals(m=cython.int)
cdef double complex f_Laguerre_Gauss_spherical(double sin_theta, double theta,
                                               double phi, dict params)

cdef class PsiSpherical:
    cdef public int m
    cdef public double x, k
    cdef public dict params
    cdef f_spectrum_type f
    cdef double ry, rz
    
    cdef double phase(self, double theta, double phi, double x, double y, double z)
    cpdef double complex integrand(self, double theta, double phi)
