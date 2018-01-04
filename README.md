![planar](Gauss_2d/img/planar_cropped_rotated_resized.png)
# Optical-beams-MEEP
Simulation of reflection and refraction of polarized optical beams at plane and curved dielectric interfaces.
*   Collection of Scheme configuration files for Gaussian beams (2d), Laguerre-Gaussian (vortex) beams (3d) and Airy
    beams (2d)
*   Based on [Meep](https://github.com/stevengj/meep) as underlying FDTD simulation software package

![field_pattern_overlay](https://cloud.githubusercontent.com/assets/28047702/26213015/b876ff6e-3bf7-11e7-8da4-9f2dffd5d470.png)
## Invocation
A Scheme configuration file (extension ``.ctl``) may be launched with the parallel version of Meep by executing the command 

``mpirun -np X meep-mpi planar_new.ctl``

 with ``X`` indicating the number of cores.
The generated HDF5 files can be processed by various visualisation toolkits. Meep comes bundled with the h5utils programs ...

``h5topng -Zc dkbluered -a gray -A eps-000000000.h5 ez-000003696.h5``

The Meep Scheme tutorial provides further assistance.
For a detailed explanation of the configuration files and the physical background, please see my dissertation thesis. Available soon.
