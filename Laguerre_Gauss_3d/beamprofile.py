# -*- coding: utf-8 -*-
"""
file:    beamprofile.py
brief:   ...
author:  Daniel Kotik
version: 1.5-beta
release date: xx.xx.2020
creation date: 22.02.2020
"""
import cython
import math
import sys

from scipy.integrate import dblquad

if not cython.compiled:
    from math import sin, cos, exp
    from cmath import exp as cexp
    print("Please consider compiling `beamprofile.py` via Cython:\n\n"
          "     `$ cythonize -3 -i beamprofile.py`")


def f_Gauss_spherical(sin_theta, theta, phi, params):
    """2d-Gaussian spectrum amplitude.

    Impementation for spherical coordinates.
    """
    W_y, k = params['W_y'], params['k']

    return exp(-(k*W_y*sin_theta/2)**2)


def f_Laguerre_Gauss_spherical(sin_theta, theta, phi, params):
    """Laguerre-Gaussian spectrum amplitude.

    Impementation for spherical coordinates.
    """
    m = params['m']

    return f_Gauss_spherical(sin_theta, theta, phi, params) * theta**abs(m) * \
        cexp(1j*m*phi)


class PsiSpherical:
    """Field amplitude class.

    Integration in spherical coordinates.

    Usage:
     psi_spherical = PsiSpherical(x=shift, params=params)
     psi_spherical(r)
    """

    def __init__(self, x, params):
        """..."""
        self.x = x
        self.params = params
        self.k = params['k']
        self.m = params['m']

        if self.m == 0:
            self.f = f_Gauss_spherical
        else:
            self.f = f_Laguerre_Gauss_spherical

    def __call__(self, r):
        """..."""
        self.ry = r.y
        self.rz = r.z
        try:
            (result,
             real_tol,
             imag_tol) = self.complex_dblquad(
                                         0, 2*math.pi, 0, math.pi/2,
                                         self.ry, self.rz)
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return self.k**2 * result

    def phase(self, theta, phi, x, y, z):
        """Phase function."""
        sin_theta, sin_phi = sin(theta), sin(phi)
        cos_theta, cos_phi = cos(theta), cos(phi)

        return self.k*(sin_theta*(y*sin_phi - z*cos_phi) + cos_theta*x)

    def integrand(self, theta, phi):
        """..."""
        return sin(theta) * cos(theta) * \
            self.f(sin(theta), theta, phi, self.params) * \
            cexp(1j*self.phase(theta, phi, self.x, self.ry, self.rz))

    def real_func(self, x, y):
        """Return real part of function."""
        return self.integrand(x, y).real

    def imag_func(self, x, y):
        """Return imag part of function."""
        return self.integrand(x, y).imag

    def complex_dblquad(self, a, b, gfun, hfun, ry, rz):
        """Integrate real and imaginary part of the given function."""
        real, real_tol = dblquad(self.real_func, a, b, gfun, hfun)
        imag, imag_tol = dblquad(self.imag_func, a, b, gfun, hfun)

        return real + 1j*imag, real_tol, imag_tol
    
    #try:
    #    getattr(psi_spherical, "called")
    #except AttributeError:
    #    psi_spherical.called = True
    #    print("Calculating inital field configuration. "
    #          "This will take some time...")


def main():

    import meep as mp

    x, y, z = -2.15, 0.3, 0.5
    r = mp.Vector3(0, y, z)

    k1 = 31.41592653589793
    w_0 = 0.25464790894703254
    m_charge = 2

    params = dict(W_y=w_0, m=m_charge, k=k1)

    psi_spherical = PsiSpherical(x=x, params=params)

    return lambda: psi_spherical(r)


if __name__ == '__main__':
    main()
