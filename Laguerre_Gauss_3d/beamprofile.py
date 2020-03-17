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
    from math import (sin as _sin,
                      cos as _cos,
                      exp as _exp)
    from cmath import exp as _cexp
    from builtins import abs as _abs
    print("\nPlease consider compiling `%s.py` via Cython: "
          "`$ cythonize -3 -i %s.py`\n" % (__name__, __name__))
else:
    from scipy import LowLevelCallable


def _real_func(x, y, func):
    """Return real part of function."""
    return func(x, y).real


def _imag_func(x, y, func):
    """Return imag part of function."""
    return func(x, y).imag


def __imag_func(n, arr, func_ptr):
    """Return imag part of function."""   
    # pure python formulation of: 
    # return (<PsiSpherical>func_ptr)(arr[0], arr[1]).imag
    return cython.cast(PsiSpherical, func_ptr).integrand(arr[0], arr[1]).imag


def __real_func(n, arr, func_ptr):
    """Return real part of function."""
    # pure python formulation of: 
    # return (<PsiSpherical>func_ptr)(arr[0], arr[1]).real
    return cython.cast(PsiSpherical, func_ptr).integrand(arr[0], arr[1]).real


def _complex_dblquad(func, a, b, gfun, hfun):
    """Integrate real and imaginary part of the given function."""
    if cython.compiled:
        # pure python formulation of: cdef void *f_ptr = <void*>func
        f_ptr = cython.declare(cython.p_void, cython.cast(cython.p_void, func))
        
        func_capsule = PyCapsule_New(f_ptr, cython.NULL, cython.NULL)
        
        current_module = sys.modules[__name__]

        ll_real_func = LowLevelCallable.from_cython(current_module, 
                                                    '__real_func', func_capsule)
        ll_imag_func = LowLevelCallable.from_cython(current_module, 
                                                    '__imag_func', func_capsule)
        real, real_tol = dblquad(ll_real_func, a, b, gfun, hfun)
        imag, imag_tol = dblquad(ll_imag_func, a, b, gfun, hfun)       
    else:
        real, real_tol = dblquad(_real_func, a, b, gfun, hfun, (func,))
        imag, imag_tol = dblquad(_imag_func, a, b, gfun, hfun, (func,))

    return real + 1j*imag, real_tol, imag_tol


def f_Gauss_spherical(sin_theta, W_y, k):
    """2d-Gaussian spectrum amplitude.

    Implementation for spherical coordinates.
    """
    return _exp(-(k*W_y*sin_theta/2)**2)


def f_Laguerre_Gauss_spherical(sin_theta, theta, phi, W_y, k, m):
    """Laguerre-Gaussian spectrum amplitude.

    Implementation for spherical coordinates.
    """
    return f_Gauss_spherical(sin_theta, W_y, k) * theta**_abs(m) * \
        _cexp(1j*m*phi)


class PsiSpherical:
    """Field amplitude class.

    Integration in spherical coordinates.

    Usage:
     psi_spherical = PsiSpherical(x=shift, params=params)
     psi_spherical(r)
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

    
    def phase(self, sin_theta, cos_theta, phi, x, y, z):
        """Phase function."""
        sin_phi = _sin(phi)
        cos_phi = _cos(phi)

        return self.k*(sin_theta*(y*sin_phi - z*cos_phi) + cos_theta*x)
    
    
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

    psi_spherical = PsiSpherical(x=x, params=params)

    return lambda: psi_spherical(r)


if __name__ == '__main__':
    main()
