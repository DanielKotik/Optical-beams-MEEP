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
import numpy as np
import sys

from datetime import datetime
from scipy.integrate import quad

print("Meep version:", mp.__version__)


def interfaceType(string):
    """Provide interface argument type."""
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
        return quad(lambda x: np.real(func(x)), a, b, **kwargs)

    def imag_integral():
        return quad(lambda x: np.imag(func(x)), a, b, **kwargs)

    result = real_integral()[0] + 1j * imag_integral()[0]
    real_tol = real_integral()[1]
    imag_tol = imag_integral()[1]

    return result, real_tol, imag_tol


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
    # physical parameters characterizing light source and interface characteris-
    # tics (must be adjusted - eihter here or via command line interface (CLI))
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

    test_output = args.test_output

    # --------------------------------------------------------------------------
    # specific Meep parameters (may need to be adjusted)
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
    r_w = kr_w / (n_ref * k_vac)
    w_0 = kw_0 / (n_ref * k_vac)
    r_c = kr_c / (n_ref * k_vac)
    shift = source_shift + r_w
    chi_rad = math.radians(chi_deg)

    # --------------------------------------------------------------------------
    # placement of the dielectric interface within the computational cell
    # --------------------------------------------------------------------------
    # helper functions
    def alpha(chi_rad):
        """Angle of inclined plane with y-axis in radians."""
        return math.pi/2 - chi_rad

    def Delta_x(alpha):
        """Inclined plane offset to the center of the cell."""
        sin_alpha = math.sin(alpha)
        cos_alpha = math.cos(alpha)
        return (sx/2) * (((math.sqrt(2) - cos_alpha) - sin_alpha) / sin_alpha)

    cell = mp.Vector3(sx, sy, 0)  # geometry-lattice

    if interface == "planar":
        default_material = mp.Medium(index=n1)
        # located at lower right edge for 45 degree
        geometry = [mp.Block(size=mp.Vector3(mp.inf, sx*math.sqrt(2), mp.inf),
                             center=mp.Vector3(+sx/2 + Delta_x(alpha(chi_rad)),
                                               -sy/2),
                             e1=mp.Vector3(1/math.tan(alpha(chi_rad)), 1, 0),
                             e2=mp.Vector3(-1, 1/math.tan(alpha(chi_rad)), 0),
                             e3=mp.Vector3(0, 0, 1),
                             material=mp.Medium(index=n2))]
    elif interface == "concave":
        default_material = mp.Medium(index=n2)
        # move center to the right in order to ensure that the point of impact
        # is always centrally placed
        geometry = [mp.Cylinder(center=mp.Vector3(-r_c*math.cos(chi_rad),
                                                  +r_c*math.sin(chi_rad)),
                                height=mp.inf,
                                radius=r_c,
                                material=mp.Medium(index=n1))]
    elif interface == "convex":
        default_material = mp.Medium(index=n1)
        # move center to the right in order to ensure that the point of impact
        # is always centrally placed
        geometry = [mp.Cylinder(center=mp.Vector3(+r_c*math.cos(chi_rad),
                                                  -r_c*math.sin(chi_rad)),
                                height=mp.inf,
                                radius=r_c,
                                material=mp.Medium(index=n2))]

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

    if test_output:
        print("Gauss spectrum:", f_Gauss(0.2, w_0))
        sys.exit()

    # --------------------------------------------------------------------------
    # plane wave decomposition
    # (purpose: calculate field amplitude at light source position if not
    #           coinciding with beam waist)
    # --------------------------------------------------------------------------
    def psi(r, f, x):
        """Field amplitude function."""
        try:
            getattr(psi, "called")
        except AttributeError:
            psi.called = True
            print("Calculating inital field configuration. "
                  "This will take some time...")

        def phi(k_y, x, y):
            """Phase function."""
            return x*math.sqrt(k1**2 - k_y**2) + k_y*y

        try:
            (result,
             real_tol,
             imag_tol) = complex_quad(lambda k_y:
                                      f(k_y) * np.exp(1j*phi(k_y, x, r.y)),
                                      -k1, k1)
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return result

    # --------------------------------------------------------------------------
    # display values of physical variables
    # --------------------------------------------------------------------------
    print()
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
    print()

    # --------------------------------------------------------------------------
    # specify current source, output functions and run simulation
    # --------------------------------------------------------------------------
    force_complex_fields = False          # default: False
    eps_averaging = True                  # default: True
    filename_prefix = None

    sources = [mp.Source(src=mp.ContinuousSource(frequency=freq, width=0.5),
                         component=mp.Ez if s_pol else mp.Ey,
                         size=mp.Vector3(0, 2, 0),
                         center=mp.Vector3(source_shift, 0, 0),
                         #amp_func=lambda r: Gauss(r, w_0)
                         amp_func=lambda r: psi(r, lambda k_y: f_Gauss(k_y, w_0), shift)
                         )
               ]

    sim = mp.Simulation(cell_size=cell,
                        boundary_layers=pml_layers,
                        default_material=default_material,
                        Courant=Courant,
                        geometry=geometry,
                        sources=sources,
                        resolution=resolution,
                        force_complex_fields=force_complex_fields,
                        eps_averaging=eps_averaging,
                        filename_prefix=filename_prefix
                        )

    sim.use_output_directory(interface)  # put output files in a separate folder

    def eSquared(r, ex, ey, ez):
        """Calculate |E|^2.

        With |.| denoting the complex modulus if 'force_complex_fields?'
        is set to true, otherwise |.| gives the Euclidean norm.
        """
        return mp.Vector3(ex, ey, ez).norm()**2

    def output_efield2(sim):
        """Output E-field intensity."""
        name = "e2_s" if s_pol else "e2_p"
        func = eSquared
        cs = [mp.Ex, mp.Ey, mp.Ez]
        return sim.output_field_function(name, cs, func, real_only=True)

    sim.run(mp.at_beginning(mp.output_epsilon),
            mp.at_end(mp.output_efield_z if s_pol else mp.output_efield_y),
            mp.at_end(output_efield2),
            until=runtime)

    print("\nend time:", datetime.now())


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
                        help=('beam waist distance to interface (30 to 50 is '
                              'good if source position coincides with beam '
                              'waist) (default: %(default)s)'))

    parser.add_argument('-kr_c',
                        type=float,
                        default=150,
                        help=('radius of curvature (if interface is either '
                              'concave of convex) (default: %(default)s)'))

    parser.add_argument('-chi_deg',
                        type=float,
                        default=45,
                        help='incidence angle in degrees (default: %(default)s)')

    parser.add_argument('-test_output',
                        action='store_true',
                        default=False,
                        help='switch to enable test print statements')

    args = parser.parse_args()
    main(args)
