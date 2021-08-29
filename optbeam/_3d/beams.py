"""Contains classes for beams in 3d."""
import math

try:
    import cython
    cython_imported = True
except ModuleNotFoundError:
    cython_imported = False

# TODO: check if conditional statement is needed at all
if not cython_imported or not cython.compiled:
    from math import (exp as _exp,
                      acos as _acos,
                      atan2 as _atan2)
    from cmath import (exp as _cexp,
                       sqrt as _csqrt)
    from builtins import abs as _abs
    from .helpers import Beam3dSpherical, Beam3dCartesian


if cython_imported and not cython.compiled:
    print("\nPlease consider compiling `%s.py` via Cython: "
          "`$ cythonize -3 -i %s.py`\n" % (__name__, __name__))


# -----------------------------------------------------------------------------
# specific beam implementations based on Beam3dSpherical or Beam3dCartesian
# base classes depending on the choice of a suitable coordinate system
# -----------------------------------------------------------------------------
class LaguerreGauss3dCartesian(Beam3dCartesian):
    """Cartesian implementaion of a 3d Laguerre-Gauss beam."""

    def __init__(self, x, params, called=False):
        """Laguerre-Gauss beam specifc parameters."""
        super().__init__(x, params, called)
        self._W_y = params['W_y']
        self._m = params['m']
        self._norm = 4 * math.pi / (self._W_y**2)

    def profile(self, r):
        """..."""
        if self.x == 0 and self._m == 0:
            return self._norm * _exp(-((r.y**2 + r.z**2) / self._W_y**2))
        else:
            return super().profile(r)

    def spectrum(self, k_y, k_z):
        """Spectrum amplitude function, f."""
        if self._m == 0:
            return self._f_Gauss_cartesian(k_y, k_z, self._W_y)
        else:
            return self._f_Laguerre_Gauss_cartesian(k_y, k_z,
                                                    self._W_y, self._k, self._m)

    def _phi(self, k_y, k_z):
        """Azimuthal angle.

        Part of coordinate transformation from k-space to (theta, phi)-space.
        """
        return _atan2(k_y, -k_z)

    def _theta(self, k_y, k_z, k):
        """Polar angle.

        Part of coordinate transformation from k-space to (theta, phi)-space.
        """
        return _acos(_csqrt(k**2 - k_y**2 - k_z**2).real / k)

    def _f_Gauss_cartesian(self, k_y, k_z, W_y):
        """2d-Gaussian spectrum amplitude.

        Impementation for Cartesian coordinates.
        """
        return _exp(-W_y**2 * (k_y**2 + k_z**2)/4)

    def _f_Laguerre_Gauss_cartesian(self, k_y, k_z, W_y, k, m):
        """Laguerre-Gaussian spectrum amplitude.

        Impementation for Cartesian coordinates.
        """
        return self._f_Gauss_cartesian(k_y, k_z, W_y) * \
            _cexp(1j*m*self._phi(k_y, k_z)) * self._theta(k_y, k_z, k)**_abs(m)


class LaguerreGauss3d(Beam3dSpherical):
    r"""3d Laguerre-Gauss beam.

    The implementation is based on [1]_.

    Parameters
    ----------
    x : type
        Description of parameter `x`.
    params : dict
        Description of parameter `params`.
    called : bool
        Description of parameter `called` (the default is False).


    References
    ----------
    .. [1] Bliokh, K Y, Aiello, A: Goos-Hänchen and Imbert-Fedorov beam shifts:
           an overview, Journal of Optics 15(1), 014001, 2013.

    Examples
    --------
    >>> LGbeam = LaguerreGauss3d(x=x, params=params)
    >>> LGbeam.profile(r)
    >>> LGbeam.spectrum(sin_theta, theta, phi)

    """

    def __init__(self, x, params, called=False):
        """Laguerre-Gauss beam specifc parameters."""
        super().__init__(x, params, called)
        self._W_y = params['W_y']
        self._m = params['m']
        self._norm = 4 * math.pi / (self._W_y**2)

    def profile(self, r):
        r"""Beam profile, psi.

        Parameters
        ----------
        r : type
            Description of parameter `r`.

        Returns
        -------
        complex
            Description of returned object.

        Notes
        -----
        .. math::

           \psi_\text{LG}(x,y,z) = k^2 \int_0^{\pi/2}\mathrm{d}\theta\int_0^{2\pi}\mathrm{d}\phi\, \sin{\theta}\cos{\theta}\,f_\text{LG}(\theta, \phi) \exp\bigl[\mathrm{i}k \bigl(x \cos{\theta} + y\sin{\theta}\sin{\phi} - z \sin{\theta}\cos{\phi}\bigr)\bigr]

        """
        if self.x == 0 and self._m == 0:
            return self._norm * _exp(-((r.y**2 + r.z**2) / self._W_y**2))
        else:
            return super().profile(r)

    def spectrum(self, sin_theta, theta, phi):
        r"""Spectrum amplitude function, f.

        Parameters
        ----------
        sin_theta : type
            Description of parameter `sin_theta`.
        theta : type
            Description of parameter `theta`.
        phi : type
            Description of parameter `phi`.

        Returns
        -------
        complex
            Description of returned object.

        Notes
        -----
        The implementation is based on [1]_.

        The spectrum amplitude is of the form

        .. math::

           f_\text{LG}(\theta,\phi) = \theta^{|m|} \
                                      \exp\biggl[-\Bigl(\frac{kw_0}{2} \
                                                        \sin\theta\Bigr)^2 \
                                          \biggr] \exp(\mathrm{i} m \phi)

        """
        if self._m == 0:
            return self._f_Gauss_spherical(sin_theta, self._W_y, self._k)
        else:
            return self._f_Laguerre_Gauss_spherical(sin_theta, theta, phi,
                                                    self._W_y, self._k, self._m)

    def _f_Gauss_spherical(self, sin_theta, W_y, k):
        """2d-Gaussian spectrum amplitude.

        Implementation for spherical coordinates.
        """
        return _exp(-(k*W_y*sin_theta/2)**2)

    def _f_Laguerre_Gauss_spherical(self, sin_theta, theta, phi, W_y, k, m):
        """Laguerre-Gaussian spectrum amplitude.

        Implementation for spherical coordinates.
        """
        return self._f_Gauss_spherical(sin_theta, W_y, k) * theta**_abs(m) * \
            _cexp(1j*m*phi)
