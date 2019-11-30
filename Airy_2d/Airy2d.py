#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file:   Airy2d.py
brief:  ...
author: Daniel Kotik
date:   30.11.2019
"""
import meep as mp
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
