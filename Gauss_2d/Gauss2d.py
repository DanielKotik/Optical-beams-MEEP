#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file:    Gauss2d.py
brief:   ...
author:  Daniel Kotik
version: X.X.X
date:    17.12.2019
"""
import argparse
import math
import meep as mp
import scipy as sp
import sys

from datetime import datetime
from scipy.integrate import quad

print("Meep version:", mp.__version__, end="\n\n")


def interfaceType(string):
    """..."""
    value = string
    if (value != "planar" and
        value != "concave" and
        value != "convex"):
        raise argparse.ArgumentTypeError('Value has to be either concave, '
                                         'convex or planar (but %s is provided)'
                                         % value)
    return value


def complex_quad(func, a, b, **kwargs):
    """Integrate real and imaginary part of the given function."""
    def real_integral():
        return quad(lambda x: sp.real(func(x)), a, b, **kwargs)

    def imag_integral():
        return quad(lambda x: sp.imag(func(x)), a, b, **kwargs)

    return (real_integral()[0] + 1j * imag_integral()[0],
            real_integral()[1], imag_integral()[1])


def Critical(n1, n2):
    """Calculate critical angle in degrees."""
    assert n1 > n2, "\nWarning: Critical angle is not defined, since n1 <= n2!"
    return math.degrees(math.asin(n2/n1))


def Brewster(n1, n2):
    """Calculate Brewster angle in degrees."""
    return math.degrees(math.atan(n2/n1))


def main(args):
    print("\nstart time:", datetime.now())

    # --------------------------------------------------------------------------
    # physical parameters characterizing light source and interface characteristics
    # (must be adjusted - eihter here or via command line interface (CLI))
    # --------------------------------------------------------------------------
    interface = args.interface
    s_pol = args.s_pol
    ref_medium = args.ref_medium

    n1 = args.n1
    n2 = args.n2

    kw_0 = args.kw_0
    kr_w = args.kr_w

    kr_c = args.kr_c

    # angle of incidence
    chi_deg = args.chi_deg
    #chi_deg = 1.0*Critical(n1, n2)
    #chi_deg = 0.95*Brewster(n1, n2)
    
    # --------------------------------------------------------------------------
    # specific Meep parameters (may need to be adjusted - either here or via CLI)
    # --------------------------------------------------------------------------
    # TODO: add short comments for every parameter
    sx = 5
    sy = 5
    pml_thickness = 0.25
    freq = 12
    runtime = 10
    pixel = 10
    source_shift = -2.15
    
    # --------------------------------------------------------------------------
    # derived Meep parameters (do not change)
    # --------------------------------------------------------------------------
    k_vac = 2 * math.pi * freq
    k1 = n1 * k_vac
    n_ref = (1  if ref_medium == 0 else
             n1 if ref_medium == 1 else
             n2 if ref_medium == 2 else math.nan)
    rw = kr_w / (n_ref * k_vac)  # TODO: rw --> r_w
    w_0 = kw_0 / (n_ref * k_vac)
    r_c = kr_c / (n_ref * k_vac)
    shift = source_shift + rw
    
    # --------------------------------------------------------------------------
    # placement of the dielectric interface within the computational cell
    # --------------------------------------------------------------------------
    # helper functions
    def alpha(chi_deg):
        """Angle of inclined plane with y-axis in radians."""
        return math.pi/2 - math.radians(chi_deg)

    def Delta_x(alpha):
        """Inclined plane offset to the center of the cell."""
        sin_alpha = math.sin(alpha)
        cos_alpha = math.cos(alpha)
        return (sx/2) * (((math.sqrt(2) - cos_alpha) - sin_alpha) / sin_alpha)

    cell = mp.Vector3(sx, sy, 0)  # geometry-lattice
    
    if interface == "planar":
        default_material = mp.Medium(index=n1)
        geometry = [mp.Block(mp.Vector3(mp.inf, sx*math.sqrt(2), mp.inf),
                         center=mp.Vector3(sx/2 + Delta_x(alpha(chi_deg)), -sy/2),
                         e1=mp.Vector3(1/math.tan(alpha(chi_deg)), 1, 0),
                         e2=mp.Vector3(-1, 1/math.tan(alpha(chi_deg)), 0),
                         e3=mp.Vector3(0, 0, 1),
                         material=mp.Medium(index=n2))]
    elif interface == "concave":
        pass
    elif interface == "convex":
        pass
    
    # --------------------------------------------------------------------------
    # add absorbing boundary conditions and discretize structure
    # --------------------------------------------------------------------------
    pml_layers = [mp.PML(pml_thickness)]
    resolution = pixel * (n1 if n1 > n2 else n2) * freq
    # set Courant factor (mandatory if either n1 or n2 is smaller than 1)
    Courant = (n1 if n1 < n2 else n2) / 2
    
    # --------------------------------------------------------------------------
    # beam profile distribution (field amplitude) at the waist of the beam
    # --------------------------------------------------------------------------
    def Gauss(r, W_y=w_0):
        """Gauss profile."""
        return math.exp(-(r.y / W_y)**2)
    
    # --------------------------------------------------------------------------
    # spectrum amplitude distribution
    # --------------------------------------------------------------------------
    def f_Gauss(k_y, W_y=w_0):
        """Gaussian spectrum amplitude."""
        return math.exp(-(k_y*W_y/2)**2)    
    
    # --------------------------------------------------------------------------
    # plane wave decomposition
    # (purpose: calculate field amplitude at light source position if not
    #           coinciding with beam waist)
    # --------------------------------------------------------------------------
    def integrand(k_y, f, x, y):
        """..."""
        return f(k_y) * sp.exp(1.0j*(x*math.sqrt(k1**2 - k_y**2) + k_y*y))

    def psi(r, f, x):
        """..."""
        result, real_tol, imag_tol = complex_quad(lambda k_y:
                                                  integrand(k_y, f, x, r.y),
                                                  -k1, k1)
        return result
    
    # --------------------------------------------------------------------------
    # display values of physical variables
    # --------------------------------------------------------------------------
    print("\n")
    print("Specified variables and derived values:")
    print("n1:", n1)
    print("n2:", n2)
    print("chi:  ", chi_deg, " [degree]")
    print("incl.:", 90 - chi_deg, " [degree]")
    print("kw_0: ", kw_0)
    if interface != "planar":
        print("kr_c: ", kr_c)
    print("kr_w: ", kr_w)
    print("k_vac:", k_vac)
    print("polarisation:", "s" if s_pol else "p")
    print("interface:", interface)
    print("\n")
    
    # --------------------------------------------------------------------------
    # specify current source, output functions and run simulation
    # --------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-interface',
                        type=interfaceType,
                        default='planar',
                        help=('specify type of interface (concave, convex, '
                              'planar) (default: %(default)s)'))

    parser.add_argument('-n1',
                        type=float,
                        default=1.54,
                        help=('index of refraction of the incident medium '
                              '(default: %(default)s)'))

    parser.add_argument('-n2',
                        type=float,
                        default=1.00,
                        help=('index of refraction of the refracted medium '
                              '(default: %(default)s)'))

    parser.add_argument('-s_pol',
                        type=bool,
                        default=True,
                        help=('True for s-spol, False for p-pol '
                              '(default: %(default)s)'))

    parser.add_argument('-ref_medium',
                        type=int,
                        default=0,
                        help=('reference medium: 0 - free space, '
                              '1 - incident medium, '
                              '2 - refracted medium (default: %(default)s)'))

    parser.add_argument('-kw_0',
                        type=float,
                        default=8,
                        help='beam width (>5 is good) (default: %(default)s)')

    parser.add_argument('-kr_w',
                        type=float,
                        default=60,
                        help=('beam waist distance to interface (30 to 50 is 
                              'good if source position coincides with beam '
                              'waist) (default: %(default)s)'))

    parser.add_argument('-kr_c',
                        type=float,
                        default=150,
                        help=('radius of curvature (if interface is either '
                              'concave of convex) (default: %(default)s)'))

    parser.add_argument('-chi_deg',
                        type=float, default=45,
                        help='incidence angle in degrees (default: %(default)s)')

    args = parser.parse_args()
    main(args)
