# -*- coding: utf-8 -*-
"""
file:    optbeam.py
brief:   ...
author:  Daniel Kotik
version: 1.5-beta
release date: xx.xx.2020
creation date: 12.04.2020
"""
import cython
import math
import sys

from scipy.integrate import quad
from types import MappingProxyType


if not cython.compiled:
    from math import sqrt as _sqrt
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
    # return (<Beam2dCartesian>func_ptr)(arr[0]).imag
    return cython.cast(Beam2dCartesian, func_ptr)._integrand(arr[0]).imag


def _real_1d_func_c(n, arr, func_ptr):
    """Return real part of a 1d function.

    Cython implementation.
    """
    # pure python formulation of:
    # return (<Beam2dCartesian>func_ptr)(arr[0]).real
    return cython.cast(Beam2dCartesian, func_ptr)._integrand(arr[0]).real


def _complex_quad(func, a, b, kwargs={}):
    """Integrate real and imaginary part of the given function."""
    if cython.compiled:
        # pure python formulation of: cdef void *f_ptr = <void*>func
        f_ptr = cython.declare(cython.p_void, cython.cast(cython.p_void, func))

        func_capsule = PyCapsule_New(f_ptr, cython.NULL, cython.NULL)

        current_module = sys.modules[__name__]

        ll_real_1d_func_c = LowLevelCallable.from_cython(current_module,
                                                         '_real_1d_func_c',
                                                         func_capsule)
        ll_imag_1d_func_c = LowLevelCallable.from_cython(current_module,
                                                         '_imag_1d_func_c',
                                                         func_capsule)
        real, real_tol = quad(ll_real_1d_func_c, a, b, **kwargs)
        imag, imag_tol = quad(ll_imag_1d_func_c, a, b, **kwargs)
    else:
        real, real_tol = quad(_real_1d_func, a, b, (func,), **kwargs)
        imag, imag_tol = quad(_imag_1d_func, a, b, (func,), **kwargs)

    return real + 1j*imag, real_tol, imag_tol


def critical(n1, n2):
    """Calculate critical angle in degrees."""
    assert n1 > n2, "\nWarning: Critical angle is not defined, since n1 <= n2!"
    return math.degrees(math.asin(n2/n1))


def brewster(n1, n2):
    """Calculate Brewster angle in degrees."""
    return math.degrees(math.atan(n2/n1))


class Beam2dCartesian:
    """..."""

    def __init__(self, x, params, called=False):
        """..."""
        self.x = x
        self._k = params['k']
        self._params = MappingProxyType(params)  # read-only view of a dict
        self.called = called

        # integration boundaries
        self._a = -self._k
        self._b = self._k

    @property
    def params(self):
        """Beam specific parameters.

        This is a read-only property.
        """
        return self._params

    def profile(self, r):
        """Field amplitude function.

        Plane wave decomposition: calculate field amplitude at light source
        position if not coinciding with beam waist.
        """
        self._ry = r.y
        self._rz = r.z

        if not self.called:
            print("Calculating inital field configuration. "
                  "This will take some time...")
            self.called = True

        try:
            (result,
             real_tol,
             imag_tol) = _complex_quad(self if cython.compiled else self._integrand,
                                       self._a, self._b, dict(limit=100))
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return result

    def spectrum(self, k_y):
        """Spectrum amplitude distribution function, f."""
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
class IncAiry2d(Beam2dCartesian):
    """2d incomplete Airy beam."""

    def __init__(self, x, params, called=False):
        """..."""
        self._W_y = params['W_y']
        self._M = params['M']
        self._W = params['W']
        super().__init__(x, params, called)

    def profile(self, r):
        """..."""
        if self.x == 0:
            # adjust integration boundaries
            self._a = self._M-self._W
            self._b = self._M+self._W

        return super().profile(r)

    def spectrum(self, k_y):
        """..."""
        return self._f_Airy(k_y, self._W_y, self._M, self._W)

    def _f_Airy(self, k_y, W_y, M, W):
        """Airy spectrum amplitude."""
        return W_y * _cexp(1.0j*(-1/3)*(k_y*W_y)**3) \
            * self._heaviside(W_y * k_y - (M - W)) \
            * self._heaviside((M + W) - W_y * k_y)

    def _heaviside(self, x):
        """Theta (Heaviside step) function."""
        return 0 if x < 0 else 1

    def _integrand(self, k_y):
        """..."""
        if self.x == 0:
            xi = k_y
            return _cexp(1.0j*(-(xi**3)/3 + (xi * self._ry/self._W_y)))
        else:
            # next line needs _integrand declared cpdef without nogil attribute,
            # and will execute slower than repeating the super class integration
            # routine function here
            #return super(IncAiry2d, self)._integrand(k_y)
            return self.spectrum(k_y) * _cexp(1j*self._phase(k_y, self.x, self._ry))


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

    k1 = 75.39822368615503
    w_0 = 0.09284038347027228
    M = 0
    W = 4
    params = dict(W_y=w_0, k=k1, M=M, W=W)

    beam = IncAiry2d(x=x, params=params)

    return (beam, r)


if __name__ == '__main__':
    main()
