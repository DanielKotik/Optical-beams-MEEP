# -*- coding: utf-8 -*-
"""
file:    optbeam.py
brief:   ...
author:  Daniel Kotik
version: 1.5-beta
release date: xx.xx.2020
creation date: 03.04.2020
"""
import cmath
import math
import sys

from scipy.integrate import quad


def complex_quad(func, a, b, **kwargs):
    """Integrate real and imaginary part of the given function."""
    def real_func(x):
        return func(x).real

    def imag_func(x):
        return func(x).imag

    real, real_tol = quad(real_func, a, b, **kwargs)
    imag, imag_tol = quad(imag_func, a, b, **kwargs)

    return real + 1j*imag, real_tol, imag_tol


class Beam2d:
    """..."""

    def __init__(self, x, params, called=False):
        """..."""
        self.x = x
        self._k = params['k']
        self._params = params
        self.called = called

        self._W_y = params['W_y']

    def profile(self, r):
        """Field amplitude function."""
        self._ry = r.y
        self._rz = r.z

        if not self.called:
            print("Calculating inital field configuration. "
                  "This will take some time...")
            self.called = True

        try:
            (result,
             real_tol,
             imag_tol) = complex_quad(self._integrand, -self._k, self._k)
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return result

    def spectrum(self, k_y):
        """Spectrum amplitude function, f."""
        return self._f_Gauss(k_y, self._W_y)

    def _phase(self, k_y, x, y):
        """Phase function."""
        return x*math.sqrt(self._k**2 - k_y**2) + k_y*y

    def _integrand(self, k_y):
        """Integrand function."""
        return self.spectrum(k_y) * cmath.exp(1j*self._phase(k_y, self.x, self._ry))

    def _f_Gauss(self, k_y, W_y):
        """Gaussian spectrum amplitude."""
        return math.exp(-(k_y*W_y/2)**2)


def main():
    class Vector3:
        """Simple 3d vector class."""

        def __init__(self, x, y, z):
            self.x = x
            self.y = y
            self.z = z

    x, y, z = -2.15, 0.3, 0.5

    r = Vector3(0, y, z)
    # alternative:
    # import meep; r = meep.Vector3(0, y, z)

    k1 = 116.11326447667875
    w_0 = 0.1061032953945969
    params = dict(W_y=w_0, k=k1)

    beam = Beam2d(x=x, params=params)

    return (beam, r)


if __name__ == '__main__':
    main()
