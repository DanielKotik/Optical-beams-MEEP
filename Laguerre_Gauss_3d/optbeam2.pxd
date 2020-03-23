# -*- coding: utf-8 -*-
"""
file:    optbeam2.pxd
brief:   ...
author:  Daniel Kotik
version: 1.5-beta
release date: xx.xx.2020
creation date: 23.02.2020
"""
cimport cython
from cpython.pycapsule cimport PyCapsule_New
from cpython cimport bool

# declare C functions as "cpdef" to export them to the module
cdef extern from "stdlib.h":
    cpdef int _abs "abs" (int n) nogil

cdef extern from "math.h":
    cpdef double _sin "sin" (double x) nogil
    cpdef double _cos "cos" (double x) nogil
    cpdef double _exp "exp" (double x) nogil
    cpdef double _acos "acos" (double x) nogil
    cpdef double _atan2 "atan2" (double y, double x) nogil

cdef extern from "complex.h":
    cpdef double complex _cexp "cexp" (double complex z) nogil
    cpdef double complex _csqrt "csqrt" (double complex z) nogil


# function prototypes
cdef double _imag_2d_func_c(int n, double *arr, void *func_ptr)
cdef double _real_2d_func_c(int n, double *arr, void *func_ptr)

@cython.locals(real=cython.double, imag=cython.double, real_tol=cython.double,
               imag_tol=cython.double)
cdef (double complex, double, double) _complex_dblquad(PsiCartesian func,
                                                       double a, double b,
                                                       double gfun, double hfun)

cdef double _phi(double k_y, double k_z) nogil
cdef double _theta(double k_y, double k_z, double k) nogil

@cython.locals(W_y=cython.double)
cdef double complex f_Gauss_cartesian(double k_y, double k_z, double W_y) nogil

@cython.locals(m=cython.int)
cdef double complex f_Laguerre_Gauss_cartesian(double k_y, double k_z,
                                               double W_y, double k, int m) nogil

cdef class PsiCartesian:
    cdef:
        readonly dict params
        int m
        double x, k, W_y
        double ry, rz
        bool called

    cdef double phase(self, double k_y, double k_z, double x, double y, double z) nogil
    cdef double complex integrand(self, double k_y, double k_z) nogil
    cdef double complex f_spectrum(self, double k_y, double k_z) nogil
