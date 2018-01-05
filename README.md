![planar](Gauss_2d/img/planar_cropped_rotated_resized.png) 
# Optical-beams-MEEP
FDTD simulation of reflection and refraction of polarized optical beams at plane and curved dielectric interfaces.
*   Collection of Scheme configuration files for Gaussian beams (2d) impinging upon planar, concave and convex
    dielectric interfaces
*   Based on [Meep](https://github.com/stevengj/meep) as underlying FDTD simulation software package
*   The focus of the beam can be placed anywhere along the propagation direction - independently of the location of the 
    source current distribution

## Invocation
A Scheme configuration file (extension ``.ctl``) may be launched with the parallel version of Meep by executing the command 

``mpirun -np X meep-mpi planar_new.ctl``

 with ``X`` indicating the number of cores.
The generated HDF5 files can be processed by various visualisation toolkits. [Meep](https://github.com/stevengj/meep) 
comes bundled with the [h5utils](https://github.com/stevengj/h5utils) programs ...

``h5topng -S2 -Zc dkbluered -a gray -A eps-000000000.h5 ez-000003696.h5`` (real part of the field pattern)

``h5topng -S2 -c hot -d e2_s.r e2_s-000003696.h5`` (intensity distribution)

The [Meep Scheme tutorial](https://meep.readthedocs.io/en/latest/Scheme_Tutorials/Basics/) provides further useful 
information and assistance.
For a more detailed explanation of our configuration files and the physical background, please see my dissertation 
thesis. Available soon.

## Currently supported _beam - interface - polarisation_ configurations
- [x] Gaussian beams (2d), planar, s- and p-polarisation
- [ ] Gaussian beams (2d), concave, s- and p-polarisation
- [ ] Gaussian beams (2d), convex, s- and p-polarisation
- [ ] Laguerre-Gaussian (vortex) beams (3d)
- [ ] Airy beams (2d)
