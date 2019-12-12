#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file:    Airy2d.py
brief:   ...
author:  Daniel Kotik
version: X.X.X
date:    30.11.2019
"""
import math
import meep as mp
import sys
from datetime import datetime

# TODO: parse arguments with argparse

def Critical(n1, n2):
    """Calculate critical angle in degrees."""
    assert n1 > n2, "\nWarning: Critical angle is not defined, since n1 <= n2!"
    return math.degrees(math.asin(n2/n1))


def Brewster(n1, n2):
    """Calculate Brewster angle in degrees."""
    return math.degrees(math.atan(n2/n1))


print("\nstart time:", datetime.now())

# -----------------------------------------------------------------------------
# physical parameters characterizing light source and interface characteristics
# (must be adjusted - either here or via command line interface (CLI))
# -----------------------------------------------------------------------------
s_pol = True
ref_medium = 0

n1 = 1.0
n2 = 0.65
kw_0 = 12
kr_w = 0

# Airy beam parameters
M = 0
W = 4

# angle of incidence
try:
    chi_deg = 45
    #chi_deg = 1.0*Critical(n1, n2)
    #chi_deg = 0.95*Brewster(n1, n2)
except Exception as e:
    print(e)
    sys.exit(1)


# -----------------------------------------------------------------------------
# specific Meep parameters (may need to be adjusted - either here or via CLI)
# -----------------------------------------------------------------------------
sx = 10
sy = 10
pml_thickness = 0.25
freq = 5
runtime = 90
pixel = 6
source_shift = 0

# -----------------------------------------------------------------------------
# derived Meep parameters (do not change)
# -----------------------------------------------------------------------------
k_vac = 2 * math.pi * freq
k1 = n1 * k_vac
n_ref = 1.0
rw = kr_w / (n_ref * k_vac)
w_0 = kw_0 / (n_ref * k_vac)
shift = source_shift + rw

# -----------------------------------------------------------------------------
# placement of the dielectric interface within the computational cell
# -----------------------------------------------------------------------------
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
default_material = mp.Medium(index=n1)
geometry = [mp.Block(mp.Vector3(mp.inf, 1, mp.inf),
                     center=mp.Vector3(),
                     material=mp.Medium(index=n2))]


# -----------------------------------------------------------------------------
# add absorbing boundary conditions and discretize structure
# -----------------------------------------------------------------------------
pml_layers = [mp.PML(pml_thickness)]
resolution = pixel * n1 * freq
# set Courant factor (mandatory if either n1 or n2 is smaller than 1)
Courant = (n1 if n1 < n2 else n2) / 2


# -----------------------------------------------------------------------------
# beam profile distribution (field amplitude) at the waist of the beam
# -----------------------------------------------------------------------------
def Gauss(r, W_y=w_0):
    """Gauss profile."""
    return math.exp(-(r.y / W_y)**2)


# -----------------------------------------------------------------------------
# spectrum amplitude distribution
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# plane wave decomposition
# (purpose: calculate field amplitude at light source position if not 
#           coinciding with beam waist)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# display values of physical variables
# -----------------------------------------------------------------------------
print("\n")
print("Specified variables and derived values:")
print("chi:  ", chi_deg, " [degree]")
print("incl.:", 90 - chi_deg, " [degree]")
print("kw_0: ", kw_0)
print("kr_w: ", kr_w)
print("k_vac:", k_vac)
print("polarisation:", "s" if s_pol else "p")
print("\n")

# -----------------------------------------------------------------------------
# specify current source, output functions and run simulation
# -----------------------------------------------------------------------------
sources = [mp.Source(src=mp.ContinuousSource(frequency=freq, width=0.5),
                     component=mp.Ez if s_pol else mp.Ey,
                     size=mp.Vector3(0, 9, 0),
                     center=mp.Vector3(source_shift, 0, 0),
                     amp_func=lambda r: Gauss(r, w_0))
           ]

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    default_material=default_material,
                    Courant=Courant,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)

sim.run(mp.at_beginning(mp.output_epsilon),
        mp.at_end(mp.output_efield_z if s_pol else mp.output_efield_y),
        until=runtime)


print("\nend time:", datetime.now())
