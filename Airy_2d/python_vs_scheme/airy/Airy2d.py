#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file:    Airy2d.py
brief:   ...
author:  Daniel Kotik
version: X.X.X
date:    30.11.2019
"""
import meep as mp
import math

from datetime import datetime


def ciritial_deg(n1, n2):
    """Calculates the critical angle in degrees."""
    raise NotImplementedError
    
def brewster_deg(n1, n2):
    """calculates the Brewster angle in degrees."""
    raise NotImplementedError

    
print("\nstart time:", datetime.now())

#------------------------------------------------------------------------------
# physical parameters characterizing light source and interface characteristics 
# (must be adjusted - either here or via command line)
#------------------------------------------------------------------------------
polarisation = "s"
ref_medium = 0

n1 = 1.0
n2 = 0.65
kw_0 = 12
kr_w = 0

M = 0
W = 4




#------------------------------------------------------------------------------
# specific Meep paramters (may need to be adjusted - either here or via command line)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# derived Meep parameters (do not change)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# placement of the dielectric interface within the computational cell
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# add absorbing boundary conditions and discretize structure
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# beam profile distribution (field amplitude) at the waist of the beam
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# spectrum amplitude distribution
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# plane wave decomposition 
# (purpose: calculate field amplitude at light source position if not coinciding with beam waist)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# display values of physical variables
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# specify current source, output functions and run simulation
#------------------------------------------------------------------------------


print("\nend time:", datetime.now())
"""
cell = mp.Vector3(16, 8, 0)

geometry = [mp.Block(mp.Vector3(mp.inf, 1, mp.inf),
                     center=mp.Vector3(),
                     material=mp.Medium(epsilon=12))]

sources = [mp.Source(mp.ContinuousSource(frequency=0.15),
                     component=mp.Ez,
                     center=mp.Vector3(-7, 0))]

pml_layers = [mp.PML(1.0)]

resolution = 10

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)

sim.run(until=200)
"""
