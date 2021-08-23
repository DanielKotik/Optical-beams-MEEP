"""optbeam: exact and fast calculation of optical beam profiles.

Copyright (c) 2021 Daniel Kotik, Jörg B. Götte.
All rights reserved.

optbeam is a Cython-based package for fast calculation of exact field
configurations of various optical beams in two and three dimensions.
Intended for use together with `PyMeep`_.

.. _PyMeep: https://meep.readthedocs.io/en/latest/

Classes
-------
* :class:`._2d.beams.Gauss2d` : 2d Gaussian beam
* :class:`._2d.beams.IncAiry2d` : 2d incomplete Airy beam
* :class:`._3d.beams.LaguerreGauss3d` : 3d vortex beam

"""

__version__ = "2.1.2"

from math import (degrees as _degrees,
                  asin as _asin,
                  atan as _atan)

from ._2d.beams import Gauss2d, IncAiry2d
from ._3d.beams import LaguerreGauss3d


def critical(n1, n2):
    """Calculate critical angle.

    Parameters
    ----------
    n1 : float
        Index of refraction of the incident medium.
    n2 : float
        Index of refraction of the refractive medium.

    Returns
    -------
    float
        Angle in degrees.

    """
    assert n1 > n2, "\nWarning: Critical angle is not defined, since n1 <= n2!"
    return _degrees(_asin(n2/n1))


def brewster(n1, n2):
    """Calculate Brewster angle.

    Parameters
    ----------
    n1 : float
        Index of refraction of the incident medium.
    n2 : float
        Index of refraction of the refractive medium.

    Returns
    -------
    float
        Angle in degrees.

    """
    return _degrees(_atan(n2/n1))
