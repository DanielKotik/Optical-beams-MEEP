# Optical-beams-MEEP
Simulation of reflection and refraction of polarized opticial beams at plane and curved dielectric interfaces.
*   Collection of Scheme configuration files for Gaussian beams (2d), Laguerre-Gaussian (vortex) beams (3d) and Airy beams (2d)
*   Based on [Meep](https://github.com/stevengj/meep) as underlying FDTD simulation software package

![field_pattern_overlay](https://cloud.githubusercontent.com/assets/28047702/26213015/b876ff6e-3bf7-11e7-8da4-9f2dffd5d470.png)
## Invocation
A Scheme configuration file (extension ``.ctl``) may be launched with the parallel version of Meep by executing the follogwing command

``mpirun -np X meep-mpi planar_new.ctl`` with ``X`` indicating the number of cores.

For a detailed explanation of the script files, please see my dissertation thesis. Available soon.
