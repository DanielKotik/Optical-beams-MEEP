![concave](Gauss_2d/img/concave_intensity_cropped_rotated_resized.png) 
![planar](Gauss_2d/img/planar_intensity_cropped_rotated_resized.png)
![convex](Gauss_2d/img/convex_intensity_cropped_rotated_resized.png)

![snap](Laguerre_Gauss_3d/img/vortex_beam_m_2_longitudinal_resized.png)
![snap](Laguerre_Gauss_3d/img/vortex_beam_m_2_transverse_resized.png)
![snap](Laguerre_Gauss_3d/img/vortex_beam_m_2_3d_half_resized.png)
# Optical-beams-MEEP
FDTD simulation of reflection and refraction of polarised optical (vortex) beams at plane and curved dielectric interfaces based on [Meep](https://github.com/stevengj/meep) as underlying FDTD simulation software package.  

The provided files (and features) are:
*   Scheme configuration file for Gaussian beams (2d) impinging upon planar, concave and convex dielectric interfaces
*   Scheme configuration file for Laguerre-Gaussian (vortex) beams (3d) impinging upon a planar dielectric interface
*   Python scripts for enhanced visualisation and analysis of the generated HDF5 output files
*   The focus of the beams can be placed anywhere along the propagation direction - independently of the location of the source current distribution

## Invocation
A Scheme configuration file (extension ``.ctl``) may be launched with the serial or parallel version of Meep and with parameters specified via command line arguments, e.g. by executing the commands:

``mpirun -quiet -np X meep-mpi interface='"concave"' Gauss2d.ctl`` (notice the combined single and double quotes)  
or  
``mpirun -quiet -np X meep-mpi e_z=0 e_y=1 freq=5 LaguerreGauss3d.ctl``

with ``X`` indicating the number of cores. All possible Meep parameters that can be set from the command line are 
defined in expressions beginning with ``(define-param ...`` in the respective configuration files.

## Visualisation
The generated HDF5 files can be processed by different visualisation tools. To get a quick impression of the data 
[Meep](https://github.com/stevengj/meep) comes bundled with the [h5utils](https://github.com/stevengj/h5utils) 
programs. Utilising these tools, visualisation is fairly easy performed by issuing for example the following commands  :

_for 2d simulations_  
``h5topng -S2 -Zc dkbluered -a gray -A eps-000000000.h5 ez-000003696.h5`` (real part of the field pattern, optical 
denser material is shaded in grey)

``h5topng -S2 -c hot -a yarg -A eps-000000000.h5 e2_s-000003696.h5`` (intensity distribution, optical 
denser material is shaded in grey)

The [Meep Scheme tutorial](https://meep.readthedocs.io/en/latest/Scheme_Tutorials/Basics/) provides further useful 
information and assistance.
For a more detailed explanation of our configuration files and the physical background, please see my dissertation thesis. Coming soon.

## Supported _beam - interface - polarisation_ configurations
-   [x] Gaussian beams (2d), planar, s- and p-polarisation
-   [x] Gaussian beams (2d), concave, s- and p-polarisation
-   [x] Gaussian beams (2d), convex, s- and p-polarisation
-   [x] Laguerre-Gaussian (vortex) beams (3d), planar, arbitrary complex polarisation
-   [ ] Airy beams (2d) (free space propagation)
