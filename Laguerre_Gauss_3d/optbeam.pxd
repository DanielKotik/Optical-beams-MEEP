# -*- coding: utf-8 -*-
"""
file:    optbeam.pxd
brief:   ...
author:  Daniel Kotik
version: 1.5-beta
release date: xx.xx.2020
creation date: 22.02.2020
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


# function declarations
cdef double _imag_2d_func_c(int n, double *arr, void *func_ptr)
cdef double _real_2d_func_c(int n, double *arr, void *func_ptr)

@cython.locals(real=cython.double, imag=cython.double, real_tol=cython.double,
               imag_tol=cython.double)
cdef (double complex, double, double) _complex_dblquad(Beam3d func,
                                                       double a, double b,
                                                       double gfun, double hfun)

# class declarations
cdef class Beam3d:
    cdef:
        readonly dict params
        double x, k
        bool called

    cdef double complex _integrand(self, double x, double y) nogil

cdef class Beam3dSpherical(Beam3d):
    cdef:
        double ry, rz

    cdef double phase(self, double sin_theta, double cos_theta, double phi,
                      double x, double y, double z) nogil
    cdef double complex spectrum(self, double sin_theta, double theta, double phi) nogil
    cdef double complex _integrand(self, double theta, double phi) nogil

cdef class Beam3dCartesian(Beam3d):
    cdef:
        double ry, rz

    cdef double phase(self, double k_y, double k_z, double x, double y, double z) nogil
    cdef double complex spectrum(self, double k_y, double k_z) nogil
    cdef double complex _integrand(self, double k_y, double k_z) nogil

cdef class LaguerreGauss3d(Beam3dSpherical):
    cdef:
      int m
      double W_y

    #@classmethod
    #@cython.binding(True)
    cdef double complex _f_Gauss_spherical(self, double sin_theta, double W_y, double k) nogil
    cdef double complex _f_Laguerre_Gauss_spherical(self, double sin_theta, double theta, double phi,
                                                    double W_y, double k, int m) nogil
    cdef double complex spectrum(self, double sin_theta, double theta, double phi) nogil

cdef class LaguerreGauss3d_(Beam3dCartesian):
    cdef:
      int m
      double W_y

    cdef double _phi(self, double k_y, double k_z) nogil
    cdef double _theta(self, double k_y, double k_z, double k) nogil
    cdef double complex _f_Gauss_cartesian(self, double k_y, double k_z, double W_y) nogil
    cdef double complex _f_Laguerre_Gauss_cartesian(self, double k_y, double k_z,
                                                    double W_y, double k, int m) nogil
    cdef double complex spectrum(self, double k_y, double k_z) nogil
