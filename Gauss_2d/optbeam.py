# -*- coding: utf-8 -*-
"""
file:    optbeam.py
brief:   ...
author:  Daniel Kotik
version: 1.5-beta
release date: xx.xx.2020
creation date: 03.04.2020
"""
import cython
import cmath
import math
import sys

from scipy.integrate import quad
from types import MappingProxyType


if not cython.compiled:
    from math import (exp as _exp,
                      sqrt as _sqrt)
    from cmath import exp as _cexp
    print("\nPlease consider compiling `%s.py` via Cython: "
          "`$ cythonize -3 -i %s.py`\n" % (__name__, __name__))
else:
    from scipy import LowLevelCallable


def _real_1d_func(x, func):
    """Return real part of a 1d function."""
    return func(x).real


def _imag_1d_func(x, func):
    """Return imag part of a 1d function."""
    return func(x).imag


def _imag_1d_func_c(n, arr, func_ptr):
    """Return imag part of a 1d function.

    Cython implementation.
    """
    # pure python formulation of:
    # return (<Beam2dCartesian>func_ptr)(arr[0], arr[1]).imag
    return cython.cast(Beam2dCartesian, func_ptr)._integrand(arr[0]).imag


def _real_1d_func_c(n, arr, func_ptr):
    """Return real part of a 1d function.

    Cython implementation.
    """
    # pure python formulation of:
    # return (<Beam2dCartesian>func_ptr)(arr[0], arr[1]).real
    return cython.cast(Beam2dCartesian, func_ptr)._integrand(arr[0]).real


def complex_quad(func, a, b):
    """Integrate real and imaginary part of the given function."""
    real, real_tol = quad(_real_1d_func, a, b, (func,))
    imag, imag_tol = quad(_imag_1d_func, a, b, (func,))

    return real + 1j*imag, real_tol, imag_tol


class Beam2dCartesian:
    """..."""

    def __init__(self, x, params, called=False):
        """..."""
        self.x = x
        self._k = params['k']
        self._params = MappingProxyType(params)  # read-only view of a dict
        self.called = called

    @property
    def params(self):
        """Beam specific parameters.

        This is a read-only property.
        """
        return self._params

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
        raise NotImplementedError

    def _phase(self, k_y, x, y):
        """Phase function."""
        return x*_sqrt(self._k**2 - k_y**2) + k_y*y

    def _integrand(self, k_y):
        """Integrand function."""
        return self.spectrum(k_y) * _cexp(1j*self._phase(k_y, self.x, self._ry))


# -----------------------------------------------------------------------------
# specific beam implementations based on Beam2dCartesian base class
# -----------------------------------------------------------------------------
class Gauss2d(Beam2dCartesian):
    """2d Gauss beam."""

    def __init__(self, x, params, called=False):
        """..."""
        self._W_y = params['W_y']
        super().__init__(x, params, called)

    def profile(self, r):
        """..."""
        if self.x == 0:
            return NotImplemented
        else:
            return super().profile(r)

    def spectrum(self, k_y):
        """Spectrum amplitude function, f."""
        return self._f_Gauss(k_y, self._W_y)

    def _f_Gauss(self, k_y, W_y):
        """Gaussian spectrum amplitude."""
        return _exp(-(k_y*W_y/2)**2)


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

    beam = Gauss2d(x=x, params=params)

    return (beam, r)


if __name__ == '__main__':
    main()
