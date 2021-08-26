"""Contains classes for beams in 3d."""
import math
import sys

try:
    import cython
    cython_imported = True
except ModuleNotFoundError:
    cython_imported = False

# TODO: check if conditional statement is needed at all
if not cython_imported or not cython.compiled:
    from math import (sin as _sin,
                      cos as _cos,
                      exp as _exp,
                      acos as _acos,
                      atan2 as _atan2)
    from cmath import (exp as _cexp,
                       sqrt as _csqrt)
    from builtins import abs as _abs
    from .helpers import Beam3d, _complex_dblquad


if cython_imported and not cython.compiled:
    print("\nPlease consider compiling `%s.py` via Cython: "
          "`$ cythonize -3 -i %s.py`\n" % (__name__, __name__))


class Beam3dSpherical(Beam3d):
    """Reference implementaion of a 3d beam in spherical coordinates."""

    def profile(self, r):
        """Beam profile function, psi."""
        self._ry = r.y
        self._rz = r.z

        if not self.called:
            print("Calculating inital field configuration. "
                  "This will take some time...")
            self.called = True

        try:
            (result,
             real_tol,
             imag_tol) = _complex_dblquad(self if cython_imported and cython.compiled
                                          else self._integrand,
                                          0, 2*math.pi, 0, math.pi/2)
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return self._k**2 * result

    def spectrum(self, sin_theta, theta, phi):
        """Spectrum amplitude function, f."""
        raise NotImplementedError

    def _phase(self, sin_theta, cos_theta, phi, x, y, z):
        """Phase function."""
        sin_phi = _sin(phi)
        cos_phi = _cos(phi)

        return self._k*(sin_theta*(y*sin_phi - z*cos_phi) + cos_theta*x)

    def _integrand(self, theta, phi):
        """Integrand function."""
        sin_theta = _sin(theta)
        cos_theta = _cos(theta)

        return sin_theta * cos_theta * self.spectrum(sin_theta, theta, phi) * \
            _cexp(1j*self._phase(sin_theta, cos_theta,
                                 phi, self.x, self._ry, self._rz))


class Beam3dCartesian(Beam3d):
    """Reference implementaion of a 3d beam in Cartesian coordinates."""

    def profile(self, r):
        """Beam profile function, psi."""
        self._ry = r.y
        self._rz = r.z

        if not self.called:
            print("Calculating inital field configuration. "
                  "This will take some time...")
            self.called = True

        try:
            (result,
             real_tol,
             imag_tol) = _complex_dblquad(self if cython_imported and cython.compiled
                                          else self._integrand,
                                          -self._k, self._k, -self._k, self._k)
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return result

    def spectrum(self, k_y, k_z):
        """Spectrum amplitude function, f."""
        raise NotImplementedError

    def _phase(self, k_y, k_z, x, y, z):
        """Phase function."""
        return x*_csqrt(self._k**2 - k_y**2 - k_z**2).real + y*k_y + z*k_z

    def _integrand(self, k_y, k_z):
        """Integrand function."""
        return self.spectrum(k_y, k_z) * \
            _cexp(1j*self._phase(k_y, k_z, self.x, self._ry, self._rz))


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
    params : type
        Description of parameter `params`.
    called : type
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
        type
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
        type
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
