# -*- coding: utf-8 -*-
"""
file:    optbeam.pxd
brief:   ...
author:  Daniel Kotik
version: 1.5-beta
release date: xx.xx.2020
creation date: 03.04.2020
"""
cimport cython
from cpython.pycapsule cimport PyCapsule_New
from cpython cimport bool

cdef extern from "math.h":
    cpdef double _exp "exp" (double x) nogil
    cpdef double _sqrt "sqrt" (double x) nogil

cdef extern from "complex.h":
    cpdef double complex _cexp "cexp" (double complex z) nogil
