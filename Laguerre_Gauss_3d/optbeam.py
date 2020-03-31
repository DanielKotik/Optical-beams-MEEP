# -*- coding: utf-8 -*-
"""
file:    optbeam.py
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
    # return (<Beam3d>func_ptr)(arr[0], arr[1]).imag
    return cython.cast(Beam3d, func_ptr).integrand(arr[0], arr[1]).imag


def _real_2d_func_c(n, arr, func_ptr):
    """Return real part of a 2d function.

    Cython implementation.
    """
    # pure python formulation of:
    # return (<Beam3d>func_ptr)(arr[0], arr[1]).real
    return cython.cast(Beam3d, func_ptr).integrand(arr[0], arr[1]).real


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


class Beam3d:
    """Abstract base class."""

    def __init__(self, x, params, called=False):
        self.x = x   # TODO: rename x to x_shift
        self.k = params['k']
        self.params = params
        self.called = called

    def integrand(self, x, y):
        """Integrand function over two coordinates x and y."""
        raise NotImplementedError


class Beam3dSpherical(Beam3d):
    """Reference implementaion of a 3d beam in spherical coordinates."""

    def profile(self, r):
        """Beam profile function, psi."""
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

    def spectrum(self, sin_theta, theta, phi):
        """Spectrum amplitude function, f."""
        raise NotImplementedError

    def phase(self, sin_theta, cos_theta, phi, x, y, z):
        """Phase function."""
        sin_phi = _sin(phi)
        cos_phi = _cos(phi)

        return self.k*(sin_theta*(y*sin_phi - z*cos_phi) + cos_theta*x)

    def integrand(self, theta, phi):
        """Integrand function."""
        sin_theta = _sin(theta)
        cos_theta = _cos(theta)

        return sin_theta * cos_theta * self.spectrum(sin_theta, theta, phi) * \
            _cexp(1j*self.phase(sin_theta, cos_theta, phi, self.x, self.ry, self.rz))


class Beam3dCartesian(Beam3d):
    """Reference implementaion of a 3d beam in Cartesian coordinates."""

    def profile(self, r):
        """Beam profile function, psi."""
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
                                          -self.k, self.k, -self.k, self.k)
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return result

    def spectrum(self, k_y, k_z):
        """Spectrum amplitude function, f."""
        raise NotImplementedError

    def phase(self, k_y, k_z, x, y, z):
        """Phase function."""
        return x*_csqrt(self.k**2 - k_y**2 - k_z**2).real + y*k_y + z*k_z

    def integrand(self, k_y, k_z):
        """Integrand function."""
        return self.spectrum(k_y, k_z) * \
            _cexp(1j*self.phase(k_y, k_z, self.x, self.ry, self.rz))


class LaguerreGauss3d_(Beam3dCartesian):
    """This class serves only as an example for a Cartesian implementaion and
    will be removed in future versions.
    """

    def __init__(self, x, params, called=False):
        """Laguerre-Gauss beam specifc parameters."""
        super().__init__(x, params, called)
        self.W_y = params['W_y']
        self.m = params['m']

    def profile(self, r):
        """..."""
        if self.x == 0:
            # TODO: Consider calling simple Gauss for special case x=0
            return NotImplemented
        else:
            return super().profile(r)

    def spectrum(self, k_y, k_z):
        """Spectrum amplitude function, f."""
        if self.m == 0:
            return self._f_Gauss_cartesian(k_y, k_z, self.W_y)
        else:
            return self._f_Laguerre_Gauss_cartesian(k_y, k_z, self.W_y, self.k, self.m)

    def _phi(self, k_y, k_z):
        """Azimuthal angle.

        Part of coordinate transformation from k-space to (theta, phi)-space.
        """
        return _atan2(k_y, -k_z)

    def _theta(self, k_y, k_z, k):
        """Polar angle.

        Part of coordinate transformation from k-space to (theta, phi)-space.
        """
        return _acos(_csqrt(k**2 - k_y**2 - k_z**2).real / k)

    def _f_Gauss_cartesian(self, k_y, k_z, W_y):
        """2d-Gaussian spectrum amplitude.

        Impementation for Cartesian coordinates.
        """
        return _exp(-W_y**2 * (k_y**2 + k_z**2)/4)

    def _f_Laguerre_Gauss_cartesian(self, k_y, k_z, W_y, k, m):
        """Laguerre-Gaussian spectrum amplitude.

        Impementation for Cartesian coordinates.
        """
        return self._f_Gauss_cartesian(k_y, k_z, W_y) * \
            _cexp(1j*m*self._phi(k_y, k_z)) * self._theta(k_y, k_z, k)**_abs(m)


class LaguerreGauss3d(Beam3dSpherical):
    """3d Laguerre-Gauss beam.

    Usage:
     LGbeam = LaguerreGauss3d(x=x, params=params)
     LGbeam.profile(r)
     LGbeam.spectrum(sin_theta, theta, phi)
    """

    def __init__(self, x, params, called=False):
        """Laguerre-Gauss beam specifc parameters."""
        super().__init__(x, params, called)
        self.W_y = params['W_y']
        self.m = params['m']

    def profile(self, r):
        """..."""
        if self.x == 0:
            # TODO: Consider calling simple Gauss for special case x=0
            return NotImplemented
        else:
            return super().profile(r)

    def spectrum(self, sin_theta, theta, phi):
        """Spectrum amplitude function, f."""
        if self.m == 0:
            return self._f_Gauss_spherical(sin_theta, self.W_y, self.k)
        else:
            return self._f_Laguerre_Gauss_spherical(sin_theta, theta, phi,
                                                    self.W_y, self.k, self.m)

    #@classmethod
    #@cython.binding(True)
    def _f_Gauss_spherical(self, sin_theta, W_y, k):
        """2d-Gaussian spectrum amplitude.

        Implementation for spherical coordinates.
        """
        return _exp(-(k*W_y*sin_theta/2)**2)

    #@classmethod
    #@cython.binding(True)
    def _f_Laguerre_Gauss_spherical(self, sin_theta, theta, phi, W_y, k, m):
        """Laguerre-Gaussian spectrum amplitude.

        Implementation for spherical coordinates.
        """
        return self._f_Gauss_spherical(sin_theta, W_y, k) * theta**_abs(m) * \
            _cexp(1j*m*phi)


def main():
    class Vector3:
        """Simple vector class."""

        def __init__(self, x, y, z):
            self.x = x
            self.y = y
            self.z = z

    x, y, z = -2.15, 0.3, 0.5

    # import meep; r = meep.Vector3(0, y, z)
    r = Vector3(0, y, z)

    k1 = 31.41592653589793
    w_0 = 0.25464790894703254
    m_charge = 2

    params = dict(W_y=w_0, m=m_charge, k=k1)

    #psi_spherical = PsiSpherical(x=x, params=params)
    #psi_cartesian = PsiCartesian(x=x, params=params)

    #return (psi_spherical, psi_cartesian, r)

    LGbeam = LaguerreGauss3d(x=x, params=params)
    LGbeam_ = LaguerreGauss3d_(x=x, params=params)

    return (LGbeam, LGbeam_, r)


if __name__ == '__main__':
    main()
