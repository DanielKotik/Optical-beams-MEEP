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
import numpy as np
import sys

from scipy.integrate import dblquad

if not cython.compiled:
    print("Please consider compiling `beamprofile.py` via Cython:\n\n"
          "     `$ cythonize -3 -i beamprofile.py`")

def complex_dblquad(func, a, b, gfun, hfun):
    """Integrate real and imaginary part of the given function."""
    def real_func(x, y):
        return np.real(func(x, y))

    def imag_func(x, y):
        return np.imag(func(x, y))
    
    real, real_tol = dblquad(real_func, a, b, gfun, hfun)
    imag, imag_tol = dblquad(imag_func, a, b, gfun, hfun)

    return real + 1j*imag, real_tol, imag_tol

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

        #@cython.cfunc
        #@cython.returns(cython.double)
        #@cython.locals(x=cython.double, y=cython.double, z=cython.double,
        #               theta=cython.double, phi=cython.double)
        def phase(theta, phi, x, y, z):
            """Phase function."""
            sin_theta, sin_phi = math.sin(theta), math.sin(phi)
            cos_theta, cos_phi = math.cos(theta), math.cos(phi)

            return k*(sin_theta*(y*sin_phi - z*cos_phi) + cos_theta*x)
        
        def integrand(theta, phi):
            """..."""
            return math.sin(theta) * math.cos(theta) * \
                f(math.sin(theta), theta, phi) * \
                np.exp(1j*phase(theta, phi, x, r.y, r.z))

        try:
            (result,
             real_tol,
             imag_tol) = complex_dblquad(integrand, 0, 2*math.pi, 0, math.pi/2)
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return k**2 * result
    
def main():
    
    import meep as mp
    
    k_y, k_z = 1.0, 5.2
    x, y, z = -2.15, 0.3, 0.5
    r = mp.Vector3(0, y, z)
    
    k1 = 31.41592653589793
    w_0 = 0.25464790894703254
    m_charge = 2
    
    f = lambda sin_theta, theta, phi: f_Laguerre_Gauss_spherical(sin_theta, theta, phi, W_y=w_0, m=m_charge, k=k1)
    
    return lambda: psi_spherical(r, f, x, k1)
    
    
if __name__ == '__main__':
    main()
    
