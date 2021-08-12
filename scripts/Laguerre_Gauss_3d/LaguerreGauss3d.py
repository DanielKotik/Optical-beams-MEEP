#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file:    LaguerreGauss3d.py
brief:   Python configuration input file for the FDTD solver Meep simulating the
         scattering of a polarised Laguerre-Gaussian beam at a planar dielectric
         interface (3d)
author:  Daniel Kotik
version: 1.5-beta
release date: xx.xx.2020
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
import math
import meep as mp
import optbeam as op

from datetime import datetime

print("Meep version:", mp.__version__)


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
    #chi_deg = 1.0*op.critical(n1, n2)
    #chi_deg = 0.95*op.brewster(n1, n2)

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

    params = dict(W_y=w_0, m=m_charge, k=k1)

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

    # specify optical beam
    LGbeam = op.LaguerreGauss3d(x=shift, params=params)

    if e_z != 0:
        source_Ez = mp.Source(src=mp.ContinuousSource(frequency=freq, width=0.5),
                              component=mp.Ez,
                              amplitude=e_z,
                              size=mp.Vector3(0, 3, 3),
                              center=mp.Vector3(source_shift, 0, 0),
                              amp_func=LGbeam.profile
                              )
        sources.append(source_Ez)

    if e_y != 0:
        source_Ey = mp.Source(src=mp.ContinuousSource(frequency=freq, width=0.5),
                              component=mp.Ey,
                              amplitude=e_y,
                              size=mp.Vector3(0, 3, 3),
                              center=mp.Vector3(source_shift, 0, 0),
                              amp_func=LGbeam.profile
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

    args = parser.parse_args()
    main(args)
