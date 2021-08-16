
__version__ = "2.0.3"

from math import (degrees as _degrees,
                  asin as _asin,
                  atan as _atan)

from ._2d.beams import Gauss2d, IncAiry2d
from ._3d.beams import LaguerreGauss3d


def critical(n1, n2):
    """Calculate critical angle in degrees."""
    assert n1 > n2, "\nWarning: Critical angle is not defined, since n1 <= n2!"
    return _degrees(_asin(n2/n1))


def brewster(n1, n2):
    """Calculate Brewster angle in degrees."""
    return _degrees(_atan(n2/n1))
