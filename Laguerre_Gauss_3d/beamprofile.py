# -*- coding: utf-8 -*-
"""
file:    beamprofile.py
brief:   ...
author:  Daniel Kotik
version: 1.5-beta
release date: xx.xx.2020
creation date: 22.02.2020
"""
import math
import numpy as np
import sys

from scipy.integrate import dblquad

def complex_dblquad(func, a, b, gfun, hfun, **kwargs):
    """Integrate real and imaginary part of the given function."""
    def real_func(x, y):
        return np.real(func(x, y))

    def imag_func(x, y):
        return np.imag(func(x, y))

    def real_integral():
        return dblquad(real_func, a, b, gfun, hfun, **kwargs)

    def imag_integral():
        return dblquad(imag_func, a, b, gfun, hfun, **kwargs)

    result = real_integral()[0] + 1j * imag_integral()[0]
    real_tol = real_integral()[1]
    imag_tol = imag_integral()[1]

    return result, real_tol, imag_tol

def f_Gauss_spherical(sin_theta, theta, phi, W_y, k):
    """2d-Gaussian spectrum amplitude.

    Impementation for spherical coordinates.
    """
    return math.exp(-(k*W_y*sin_theta/2)**2)

def f_Laguerre_Gauss_spherical(sin_theta, theta, phi, W_y, m, k):
    """Laguerre-Gaussian spectrum amplitude.

    Impementation for spherical coordinates.
    """
    return f_Gauss_spherical(sin_theta, theta, phi, W_y, k) * theta**abs(m) * \
        np.exp(1j*m*phi)

def psi_spherical(r, f, x, k):
        """Field amplitude function.

        Integration in spherical coordinates.
        """
        try:
            getattr(psi_spherical, "called")
        except AttributeError:
            psi_spherical.called = True
            print("Calculating inital field configuration. "
                  "This will take some time...")

        def phase(theta, phi, x, y, z):
            """Phase function."""
            sin_theta, sin_phi = math.sin(theta), math.sin(phi)
            cos_theta, cos_phi = math.cos(theta), math.cos(phi)

            return k*(sin_theta*(y*sin_phi - z*cos_phi) + cos_theta*x)

        try:
            (result,
             real_tol,
             imag_tol) = complex_dblquad(lambda theta, phi:
                                         math.sin(theta) * math.cos(theta) * \
                                         f(math.sin(theta), theta, phi) * \
                                         np.exp(1j*phase(theta, phi, x, r.y, r.z)),
                                         0, 2*math.pi, 0, math.pi/2)
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return k**2 * result