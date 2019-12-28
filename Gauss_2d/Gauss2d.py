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

    print(interface)


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
