#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file:    LaguerreGauss3d.py
brief:   Python configuration input file for the FDTD solver Meep simulating the
         scattering of a polarised Laguerre-Gaussian beam at a planar dielectric
         interface (3d)
author:  Daniel Kotik
version: 1.4.0
release date: 12.02.2020
creation date: 10.01.2020


invocation:

 a) launch the serial version of meep with specified polarisation (p)

        python3 LaguerreGauss3d.py -e_z 0 -e_y 1

 b) launch the parallel version of meep using 8 cores

        mpirun -quiet -np 8 python3 LaguerreGauss3d.py


coordinate system in meep (defines center of computational cell):

                        --|-----> x
                          |
                          |
                          v y


visualisation:

 - slice within the plane of incidence (x-y plane)
          h5topng -S2 -0 -z 0 -c hot [HDF5FILE]

 - slice transversal to the incident propagation axis (INDEX specifies slice index)
          h5topng -S2 -x [INDEX] -c hot [HDF5FILE]

 - full 3D simulation (creating a VTK file to be opened e.g., with MayaVi or ParaView)
          h5tovtk [HDF5FILE]

As input HDF5FILE choose between, for example, 'e_real2_p-___.h5' and
'e_imag2_p-___.h5' (these are proportional to the electric field energy
density) or the sum of both 'e2.h5' (which is proportional to the complex modulus
of the complex electric field) obtained by
          h5math -e "d1 + d2" e2.h5 e_real2_p-___.h5 e_imag2_p-___.h5

"""
import argparse
import cmath
import math
import meep as mp
import sys

from datetime import datetime
from scipy.integrate import dblquad


print("Meep version:", mp.__version__)


def complex_dblquad(func, a, b, gfun, hfun, **kwargs):
    """Integrate real and imaginary part of the given function."""
    def real_func(x, y):
        return func(x, y).real

    def imag_func(x, y):
        return func(x, y).imag

    real, real_tol = dblquad(real_func, a, b, gfun, hfun, **kwargs)
    imag, imag_tol = dblquad(imag_func, a, b, gfun, hfun, **kwargs)

    return real + 1j*imag, real_tol, imag_tol


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
    e_z = args.e_z
    e_y = args.e_y

    m_charge = args.m_charge
    ref_medium = args.ref_medium

    n1 = args.n1
    n2 = args.n2

    kw_0 = args.kw_0
    kr_w = args.kr_w

    # angle of incidence
    chi_deg = args.chi_deg
    #chi_deg = 1.0*Critical(n1, n2)
    #chi_deg = 0.95*Brewster(n1, n2)

    test_output = args.test_output

    # --------------------------------------------------------------------------
    # specific Meep parameters (may need to be adjusted)
    # --------------------------------------------------------------------------
    sx = 5   # size of cell including PML in x-direction
    sy = 5   # size of cell including PML in y-direction
    sz = 4   # size of cell including PML in z-direction
    pml_thickness = 0.25   # thickness of PML layer
    freq = 5       # vacuum frequency of source (default 5)
    runtime = 10   # runs simulation for 10 times freq periods

    # number of pixels per wavelength in the denser medium (at least 10,
    # 20 to 30 is a good choice)
    pixel = 10

    # source position with respect to the center (point of impact) in Meep
    # units (-2.15 good); if equal -r_w, then source position coincides with
    # waist position
    source_shift = -2.15

    # --------------------------------------------------------------------------
    # derived (Meep) parameters (do not change)
    # --------------------------------------------------------------------------
    k_vac = 2 * math.pi * freq
    k1 = n1 * k_vac
    n_ref = (1  if ref_medium == 0 else
             n1 if ref_medium == 1 else
             n2 if ref_medium == 2 else math.nan)
    r_w = kr_w / (n_ref * k_vac)
    w_0 = kw_0 / (n_ref * k_vac)
    shift = source_shift + r_w
    chi_rad = math.radians(chi_deg)
    s_pol = True if (e_z == 1 and e_y == 0) else False
    p_pol = True if (e_z == 0 and e_y == 1) else False

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

    cell = mp.Vector3(sx, sy, sz)  # geometry-lattice

    default_material = mp.Medium(index=n1)
    # located at lower right edge for 45 degree tilt
    geometry = [mp.Block(size=mp.Vector3(mp.inf, sx*math.sqrt(2), mp.inf),
                         center=mp.Vector3(+sx/2 + Delta_x(alpha(chi_rad)),
                                           -sy/2),
                         e1=mp.Vector3(1/math.tan(alpha(chi_rad)), 1, 0),
                         e2=mp.Vector3(-1, 1/math.tan(alpha(chi_rad)), 0),
                         e3=mp.Vector3(0, 0, 1),
                         material=mp.Medium(index=n2))]

    # --------------------------------------------------------------------------
    # add absorbing boundary conditions and discretize structure
    # --------------------------------------------------------------------------
    pml_layers = [mp.PML(pml_thickness)]
    resolution = pixel * (n1 if n1 > n2 else n2) * freq
    # set Courant factor (mandatory if either n1 or n2 is smaller than 1)
    Courant = (n1 if n1 < n2 else n2) / 3

    # --------------------------------------------------------------------------
    # 2d-beam profile distribution (field amplitude) at the waist of the beam
    # --------------------------------------------------------------------------
    def Gauss(r, W_y=w_0):
        """Gauss profile."""
        return math.exp(-((r.y**2 + r.z**2) / W_y**2))

    # --------------------------------------------------------------------------
    # some test outputs
    # --------------------------------------------------------------------------
    if test_output:
        print("Gauss 2d beam profile:", Gauss(mp.Vector3(0, 0.5, 0.2), w_0))
        print()

    # --------------------------------------------------------------------------
    # spectrum amplitude distribution(s)
    # --------------------------------------------------------------------------

    # cartesian coordinates (not recommmended) -------------------------
    def phi(k_y, k_z):
        """Azimuthal angle.

        Part of coordinate transformation from k-space to (theta, phi)-space.
        """
        return math.atan2(k_y, -k_z)

    def theta(k_y, k_z, k):
        """Polar angle.

        Part of coordinate transformation from k-space to (theta, phi)-space.
        """
        return math.acos(cmath.sqrt(k**2 - k_y**2 - k_z**2).real / k)

    def f_Gauss_cartesian(k_y, k_z, W_y=w_0):
        """2d-Gaussian spectrum amplitude.

        Impementation for Cartesian coordinates.
        """
        return math.exp(-W_y**2 * (k_y**2 + k_z**2)/4)

    def f_Laguerre_Gauss_cartesian(k_y, k_z, W_y=w_0, m=m_charge, k=k1):
        """Laguerre-Gaussian spectrum amplitude.

        Impementation for Cartesian coordinates.
        """
        return f_Gauss_cartesian(k_y, k_z, W_y) * \
            cmath.exp(1j*m*phi(k_y, k_z)) * theta(k_y, k_z, k)**abs(m)

    # spherical coordinates --------------------------------------------
    def f_Gauss_spherical(sin_theta, theta, phi, W_y=w_0, k=k1):
        """2d-Gaussian spectrum amplitude.

        Impementation for spherical coordinates.
        """
        return math.exp(-(k*W_y*sin_theta/2)**2)

    def f_Laguerre_Gauss_spherical(sin_theta, theta, phi, W_y=w_0, m=m_charge, k=k1):
        """Laguerre-Gaussian spectrum amplitude.

        Impementation for spherical coordinates.
        """
        return f_Gauss_spherical(sin_theta, theta, phi, W_y, k) * theta**abs(m) * \
            cmath.exp(1j*m*phi)

    # --------------------------------------------------------------------------
    # some test outputs
    # --------------------------------------------------------------------------
    if test_output:
        k_y, k_z = 1.0, 5.2
        theta_ = theta(k_y, k_z, k1)
        phi_ = phi(k_y, k_z)

        print("Gauss spectrum (cartesian):",
              f_Gauss_cartesian(k_y, k_z, w_0))
        print("Gauss spectrum (spherical):",
              f_Gauss_spherical(math.sin(theta_), theta_, phi_, w_0, k1))
        print()
        print("L-G spectrum   (cartesian):",
              f_Laguerre_Gauss_cartesian(k_y, k_z, w_0, m_charge, k1))
        print("L-G spectrum   (spherical):",
              f_Laguerre_Gauss_spherical(math.sin(theta_), theta_, phi_, w_0,
                                         m_charge, k1))
        print()

    # --------------------------------------------------------------------------
    # plane wave decomposition
    # (purpose: calculate field amplitude at light source position if not
    #           coinciding with beam waist)
    # --------------------------------------------------------------------------
    def psi_cartesian(r, f, x, k=k1):
        """Field amplitude function.

        Integration in Cartesian coordinates.
        """
        try:
            getattr(psi_cartesian, "called")
        except AttributeError:
            psi_cartesian.called = True
            print("Calculating inital field configuration. "
                  "This will take some time...")

        def phase(k_y, k_z, x, y, z):
            """Phase function."""
            return x*cmath.sqrt(k**2 - k_y**2 - k_z**2).real + y*k_y + z*k_z

        try:
            (result,
             real_tol,
             imag_tol) = complex_dblquad(lambda k_y, k_z:
                                         f(k_y, k_z) * \
                                         cmath.exp(1j*phase(k_y, k_z, x, r.y, r.z)),
                                         -k, k, -k, k)
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return result

    def psi_spherical(r, f, x, k=k1):
        """Field amplitude function.

        Integration in spherical coordinates.
        """
        try:
            getattr(psi_spherical, "called")
        except AttributeError:
            psi_spherical.called = True
            print("Calculating inital field configuration. "
                  "This will take some time...")

        def phase(theta, phi, x, y, z):
            """Phase function."""
            sin_theta, sin_phi = math.sin(theta), math.sin(phi)
            cos_theta, cos_phi = math.cos(theta), math.cos(phi)

            return k*(sin_theta*(y*sin_phi - z*cos_phi) + cos_theta*x)

        try:
            (result,
             real_tol,
             imag_tol) = complex_dblquad(lambda theta, phi:
                                         math.sin(theta) * math.cos(theta) * \
                                         f(math.sin(theta), theta, phi) * \
                                         cmath.exp(1j*phase(theta, phi, x, r.y, r.z)),
                                         0, 2*math.pi, 0, math.pi/2)
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return k**2 * result

    # --------------------------------------------------------------------------
    # some test outputs (uncomment if needed)
    # --------------------------------------------------------------------------
    if test_output:
        k_y, k_z = 1.0, 5.2
        x, y, z = -2.15, 0.3, 0.5
        r = mp.Vector3(0, y, z)

        print("phi:", phi(k_y, k_z))
        print()
        print("psi            (cartesian):",
              psi_cartesian(r, f_Laguerre_Gauss_cartesian, x, k1))
        print("psi            (spherical):",
              psi_spherical(r, f_Laguerre_Gauss_spherical, x, k1))
        print("psi       (origin, simple):", Gauss(r))
        sys.exit()

    # --------------------------------------------------------------------------
    # display values of physical variables
    # --------------------------------------------------------------------------
    print()
    print("Expected output file size:",
          round(8*(sx*sy*sz*resolution**3)/(1024**2)), "MiB")
    print()
    print("Specified variables and derived values:")
    print("n1:", n1)
    print("n2:", n2)
    print("chi:  ", chi_deg, " [degree]")
    # interface inclination with respect to the x-axis
    print("incl.:", 90 - chi_deg, " [degree]")
    print("kw_0: ", kw_0)
    print("kr_w: ", kr_w)
    print("k_vac:", k_vac)
    print("vortex charge:", m_charge)
    print("Jones vector components: "
          "(e_z=", e_z, ", e_y=", e_y, ")", sep="", end="")
    print(" --->", ("s-" if s_pol else
                    "p-" if p_pol else
                    "mixed-") + "polarisation")
    print("degree of linear   polarisation at pi/4:",
          2*(-e_z.conjugate()*e_y).real)
    print("degree of circular polarisation:", 2*(-e_z.conjugate()*e_y).imag)

    # --------------------------------------------------------------------------
    # exploiting symmetries to reduce computational effort
    # (only possible for beams without intrinsic orbital angular momentum, i.e.
    #  no vortex charge)
    # --------------------------------------------------------------------------

    # The plane of incidence (x-y-plane) is a mirror plane which is characterised
    # to be orthogonal to the z-axis (symmetry of the geometric structure).
    # Symmetry of the sources must be ensured simultaneously, which is only
    # possible for certain cases. If I am not mistaken this can only be achieved
    # for vortex free beams with pure s- or p-polarisation, i.e. where either
    # the Ez or Ey component is specified.
    symmetries = []
    if m_charge == 0:
        if s_pol:
            symmetries.append(mp.Mirror(mp.Z, phase=-1))
        if p_pol:
            symmetries.append(mp.Mirror(mp.Y, phase=+1))

    # --------------------------------------------------------------------------
    # specify current source, output functions and run simulation
    # --------------------------------------------------------------------------
    force_complex_fields = True           # default: True
    eps_averaging = True                  # default: True

    sources = []

    if e_z != 0:
        source_Ez = mp.Source(src=mp.ContinuousSource(frequency=freq, width=0.5),
                              component=mp.Ez,
                              amplitude=e_z,
                              size=mp.Vector3(0, 3, 3),
                              center=mp.Vector3(source_shift, 0, 0),
                              #amp_func=lambda r: Gauss(r, w_0)
                              #amp_func=lambda r: psi_cartesian(r, f_Laguerre_Gauss_cartesian, shift)
                              amp_func=lambda r: psi_spherical(r, f_Gauss_spherical, shift) if m_charge == 0 else
                                                 psi_spherical(r, f_Laguerre_Gauss_spherical, shift)
                              )
        sources.append(source_Ez)

    if e_y != 0:
        source_Ey = mp.Source(src=mp.ContinuousSource(frequency=freq, width=0.5),
                              component=mp.Ey,
                              amplitude=e_y,
                              size=mp.Vector3(0, 3, 3),
                              center=mp.Vector3(source_shift, 0, 0),
                              #amp_func=lambda r: Gauss(r, w_0)
                              #amp_func=lambda r: psi_cartesian(r, f_Laguerre_Gauss_cartesian, shift)
                              amp_func=lambda r: psi_spherical(r, f_Gauss_spherical, shift) if m_charge == 0 else
                                                 psi_spherical(r, f_Laguerre_Gauss_spherical, shift)
                              )
        sources.append(source_Ey)

    sim = mp.Simulation(cell_size=cell,
                        boundary_layers=pml_layers,
                        symmetries=symmetries,
                        default_material=default_material,
                        Courant=Courant,
                        geometry=geometry,
                        sources=sources,
                        resolution=resolution,
                        force_complex_fields=force_complex_fields,
                        eps_averaging=eps_averaging
                        )

    sim.use_output_directory()   # put output files in a separate folder

    def efield_real_squared(r, ex, ey, ez):
        """Calculate |Re E|^2.

        With |.| denoting the complex modulus if 'force_complex_fields?'
        is set to true, otherwise |.| gives the Euclidean norm.
        """
        return ex.real**2 + ey.real**2 + ez.real**2

    def efield_imag_squared(r, ex, ey, ez):
        """Calculate |Im E|^2.

        With |.| denoting the complex modulus if 'force_complex_fields?'
        is set to true, otherwise |.| gives the Euclidean norm.
        """
        return ex.imag**2 + ey.imag**2 + ez.imag**2

    def output_efield_real_squared(sim):
        """Output E-field (real part) intensity."""
        name = "e_real2_s" if s_pol else "e_real2_p" if p_pol else "e_real2_mixed"
        func = efield_real_squared
        cs = [mp.Ex, mp.Ey, mp.Ez]
        return sim.output_field_function(name, cs, func, real_only=True)

    def output_efield_imag_squared(sim):
        """Output E-field (imag part) intensity."""
        name = "e_imag2_s" if s_pol else "e_imag2_p" if p_pol else "e_imag2_mixed"
        func = efield_imag_squared
        cs = [mp.Ex, mp.Ey, mp.Ez]
        return sim.output_field_function(name, cs, func, real_only=True)

    run_args = [#mp.at_beginning(mp.output_epsilon),    # output of dielectric function
                #mp.at_end(mp.output_efield_x),         # output of E_x component
                #mp.at_end(mp.output_efield_z),         # output of E_y component
                #mp.at_end(mp.output_efield_y),         # output of E_z component
                mp.at_end(output_efield_real_squared)  # output of electric field intensity
                ]

    if force_complex_fields:
        run_args.append(mp.at_end(output_efield_imag_squared))

    sim.run(*run_args, until=runtime)

    print("\nend time:", datetime.now())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-e_z',
                        type=eval,
                        default=1,
                        help='z-component of Jones vector (default: %(default)s)')

    parser.add_argument('-e_y',
                        type=eval,
                        default=0,
                        help='y-component of Jones vector (default: %(default)s)')

    parser.add_argument('-m_charge',
                        type=int,
                        default=2,
                        help=('vortex charge (azimuthal quantum number)'
                              ' (default: %(default)s)'))

    parser.add_argument('-n1',
                        type=float,
                        default=1.00,
                        help=('index of refraction of the incident medium '
                              '(default: %(default)s)'))

    parser.add_argument('-n2',
                        type=float,
                        default=1.54,
                        help=('index of refraction of the refracted medium '
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
                        default=0,
                        help=('beam waist distance to interface (30 to 50 is '
                              'good if source position coincides with beam '
                              'waist) (default: %(default)s)'))

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
