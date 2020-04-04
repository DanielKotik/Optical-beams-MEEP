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


def f_Gauss(k_y, params):
    """Gaussian spectrum amplitude."""
    W_y = params['W_y']

    return math.exp(-(k_y*W_y/2)**2)


def psi(r, x, params):
    """Field amplitude function."""
    k1 = params["k"]

    try:
        getattr(psi, "called")
    except AttributeError:
        psi.called = True
        print("Calculating inital field configuration. "
              "This will take some time...")

    def phase(k_y, x, y):
        """Phase function."""
        return x*math.sqrt(k1**2 - k_y**2) + k_y*y

    try:
        (result,
         real_tol,
         imag_tol) = complex_quad(lambda k_y:
                                  f_Gauss(k_y, params)
                                  * cmath.exp(1j*phase(k_y, x, r.y)),
                                  -k1, k1)
    except Exception as e:
        print(type(e).__name__ + ":", e)
        sys.exit()

    return result


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

    return (lambda r: psi(r, x, params), r)


if __name__ == '__main__':
    main()
