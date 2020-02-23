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


cdef double f_Gauss_spherical(double sin_theta, double theta, double phi, 
                              double W_y, double k)

cdef double complex f_Laguerre_Gauss_spherical(double sin_theta, double theta, 
                                               double phi, double W_y, int m, 
                                               double k)

