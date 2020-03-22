# -*- coding: utf-8 -*-
"""
file:    optbeam2.py
brief:   ...
author:  Daniel Kotik
version: 1.5-beta
release date: xx.xx.2020
creation date: 20.03.2020
"""
import cython
import math
import sys

from scipy.integrate import dblquad


if not cython.compiled:
    from math import (sin as _sin,
                      cos as _cos,
                      exp as _exp,
                      acos as _acos,
                      atan2 as _atan2)
    from cmath import (exp as _cexp,
                       sqrt as _csqrt)
    from builtins import abs as _abs
    print("\nPlease consider compiling `%s.py` via Cython: "
          "`$ cythonize -3 -i %s.py`\n" % (__name__, __name__))
else:
    from scipy import LowLevelCallable


def _real_2d_func(x, y, func):
    """Return real part of a 2d function."""
    return func(x, y).real


def _imag_2d_func(x, y, func):
    """Return imag part of a 2d function."""
    return func(x, y).imag


def _imag_2d_func_c(n, arr, func_ptr):
    """Return imag part of a 2d function.

    Cython implementation.
    """
    # pure python formulation of:
    # return (<PsiCartesian>func_ptr)(arr[0], arr[1]).imag
    return cython.cast(PsiCartesian, func_ptr).integrand(arr[0], arr[1]).imag


def _real_2d_func_c(n, arr, func_ptr):
    """Return real part of a 2d function.

    Cython implementation.
    """
    # pure python formulation of:
    # return (<PsiCartesian>func_ptr)(arr[0], arr[1]).real
    return cython.cast(PsiCartesian, func_ptr).integrand(arr[0], arr[1]).real


def _complex_dblquad(func, a, b, gfun, hfun):
    """Integrate real and imaginary part of the given function."""
    if cython.compiled:
        # pure python formulation of: cdef void *f_ptr = <void*>func
        f_ptr = cython.declare(cython.p_void, cython.cast(cython.p_void, func))

        func_capsule = PyCapsule_New(f_ptr, cython.NULL, cython.NULL)

        current_module = sys.modules[__name__]

        ll_real_2d_func_c = LowLevelCallable.from_cython(current_module,
                                                         '_real_2d_func_c',
                                                         func_capsule)
        ll_imag_2d_func_c = LowLevelCallable.from_cython(current_module,
                                                         '_imag_2d_func_c',
                                                         func_capsule)
        real, real_tol = dblquad(ll_real_2d_func_c, a, b, gfun, hfun)
        imag, imag_tol = dblquad(ll_imag_2d_func_c, a, b, gfun, hfun)
    else:
        real, real_tol = dblquad(_real_2d_func, a, b, gfun, hfun, (func,))
        imag, imag_tol = dblquad(_imag_2d_func, a, b, gfun, hfun, (func,))

    return real + 1j*imag, real_tol, imag_tol


def _phi(k_y, k_z):
    """Azimuthal angle.

    Part of coordinate transformation from k-space to (theta, phi)-space.
    """
    return _atan2(k_y, -k_z)


def _theta(k_y, k_z, k):
    """Polar angle.

    Part of coordinate transformation from k-space to (theta, phi)-space.
    """
    return _acos(_csqrt(k**2 - k_y**2 - k_z**2).real / k)


def f_Gauss_cartesian(k_y, k_z, W_y):
    """2d-Gaussian spectrum amplitude.

    Impementation for Cartesian coordinates.
    """
    return _exp(-W_y**2 * (k_y**2 + k_z**2)/4)


def f_Laguerre_Gauss_cartesian(k_y, k_z, W_y, k, m):
    """Laguerre-Gaussian spectrum amplitude.

    Impementation for Cartesian coordinates.
    """
    return f_Gauss_cartesian(k_y, k_z, W_y) * \
        _cexp(1j*m*_phi(k_y, k_z)) * _theta(k_y, k_z, k)**_abs(m)


class PsiCartesian:
    """Field amplitude class.

    Integration in cartesian coordinates.

    Usage:
     psi_cartesian = PsiCartesian(x=shift, params=params)
     psi_cartesian(r)  # r...vector-like object with scalar y and z attributes
    """

    def __init__(self, x, params, called=False):
        """..."""
        self.x = x
        self.W_y = params['W_y']
        self.k = params['k']
        self.m = params['m']
        self.called = called

    def __call__(self, r):
        """..."""
        self.ry = r.y
        self.rz = r.z

        if not self.called:
            print("Calculating inital field configuration. "
                  "This will take some time...")
            self.called = True

        try:
            (result,
             real_tol,
             imag_tol) = _complex_dblquad(self if cython.compiled else self.integrand,
                                          0, 2*math.pi, 0, math.pi/2)
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return self.k**2 * result

    def phase(k, k_y, k_z, x, y, z):
        """Phase function."""
        return x*_csqrt(k**2 - k_y**2 - k_z**2).real + y*k_y + z*k_z

    def f_spectrum(self, sin_theta, theta, phi):
        """Spectrum amplitude function."""
        if self.m == 0:
            return f_Gauss_spherical(sin_theta, self.W_y, self.k)
        else:
            return f_Laguerre_Gauss_spherical(sin_theta, theta, phi,
                                              self.W_y, self.k, self.m)

    def integrand(self, theta, phi):
        """Integrand function."""
        sin_theta = _sin(theta)
        cos_theta = _cos(theta)

        return sin_theta * cos_theta * self.f_spectrum(sin_theta, theta, phi) * \
            _cexp(1j*self.phase(sin_theta, cos_theta, phi, self.x, self.ry, self.rz))


def main():

    import meep as mp

    x, y, z = -2.15, 0.3, 0.5
    r = mp.Vector3(0, y, z)

    k1 = 31.41592653589793
    w_0 = 0.25464790894703254
    m_charge = 2

    params = dict(W_y=w_0, m=m_charge, k=k1)

    psi_cartesian = PsiCartesian(x=x, params=params)

    return lambda: psi_cartesian(r)


if __name__ == '__main__':
    main()
