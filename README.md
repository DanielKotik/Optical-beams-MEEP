![concave](Gauss_2d/img/concave_intensity_cropped_rotated_resized.png) 
![planar](Gauss_2d/img/planar_intensity_cropped_rotated_resized.png)
![convex](Gauss_2d/img/convex_intensity_cropped_rotated_resized.png)

![snap](Laguerre_Gauss_3d/img/vortex_beam_m_2_longitudinal_resized.png)
![snap](Laguerre_Gauss_3d/img/vortex_beam_m_2_transverse_resized.png)
![snap](Laguerre_Gauss_3d/img/vortex_beam_m_2_3d_half_resized.png)
# Optical-beams-MEEP
FDTD simulation of reflection and refraction of polarised optical beams at plane and curved dielectric interfaces.
*   Collection of Scheme configuration files for Gaussian beams (2d) impinging upon planar, concave and convex
    dielectric interfaces
*   Based on [Meep](https://github.com/stevengj/meep) as underlying FDTD simulation software package
*   The focus of the beam can be placed anywhere along the propagation direction - independently of the location of the 
    source current distribution

## Invocation
A Scheme configuration file (extension ``.ctl``) may be launched with the serial or parallel version of Meep and with parameters specified via command line arguments, e.g. by executing the command 

``mpirun -np X meep-mpi interface='"concave"' Gauss2d.ctl`` (notice the combined single and double quotes)

 with ``X`` indicating the number of cores.
The generated HDF5 files can be processed by different visualisation toolkits. [Meep](https://github.com/stevengj/meep) 
comes bundled with the [h5utils](https://github.com/stevengj/h5utils) programs. Utilising these tools, visualisation is easily performed by issuing the commands (examples)

``h5topng -S2 -Zc dkbluered -a gray -A eps-000000000.h5 ez-000003696.h5`` (real part of the field pattern, optical 
denser material is shaded in grey)

``h5topng -S2 -c hot -a yarg -A eps-000000000.h5 e2_s-000003696.h5`` (intensity distribution, optical 
denser material is shaded in grey)

The [Meep Scheme tutorial](https://meep.readthedocs.io/en/latest/Scheme_Tutorials/Basics/) provides further useful 
information and assistance.
For a more detailed explanation of our configuration files and the physical background, please see my dissertation 
thesis. Available soon.

## Supported _beam - interface - polarisation_ configurations
-   [x] Gaussian beams (2d), planar, s- and p-polarisation
-   [x] Gaussian beams (2d), concave, s- and p-polarisation
-   [x] Gaussian beams (2d), convex, s- and p-polarisation
-   [x] Laguerre-Gaussian (vortex) beams (3d), planar, arbitrary complex polarisation
-   [ ] Airy beams (2d) (free space propagation)
