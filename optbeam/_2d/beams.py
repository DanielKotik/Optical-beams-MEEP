"""Contains classes for beams in 2d."""

import math
import sys

try:
    import cython
    cython_imported = True
except ModuleNotFoundError:
    cython_imported = False

if not cython_imported or not cython.compiled:
    from math import (exp as _exp,
                      sqrt as _sqrt)
    from cmath import exp as _cexp
    from .helpers import Beam2d, _complex_quad

if cython_imported and not cython.compiled:
    print("\nPlease consider compiling `%s.py` via Cython: "
          "`$ cythonize -3 -i %s.py`\n" % (__name__, __name__))


class Beam2dCartesian(Beam2d):
    """..."""

    def profile(self, r):
        """Field amplitude function.

        Plane wave decomposition: calculate field amplitude at light source
        position if not coinciding with beam waist.
        """
        self._ry = r.y
        self._rz = r.z

        if not self.called:
            print("Calculating inital field configuration. "
                  "This will take some time...")
            self.called = True

        try:
            (result,
             real_tol,
             imag_tol) = _complex_quad(self if cython_imported and cython.compiled
                                       else self._integrand,
                                       -self._k, self._k)
        except Exception as e:
            print(type(e).__name__ + ":", e)
            sys.exit()

        return result

    def spectrum(self, k_y):
        """Spectrum amplitude distribution function, f."""
        raise NotImplementedError

    def _phase(self, k_y, x, y):
        """Phase function."""
        return x*_sqrt(self._k**2 - k_y**2) + k_y*y

    def _integrand(self, k_y):
        """Integrand function."""
        return self.spectrum(k_y) * _cexp(1j*self._phase(k_y, self.x, self._ry))


# -----------------------------------------------------------------------------
# specific beam implementations based on Beam2dCartesian base class
# -----------------------------------------------------------------------------
class Gauss2d(Beam2dCartesian):
    """2d Gauss beam."""

    def __init__(self, x, params, called=False):
        """..."""
        self._W_y = params['W_y']
        self._norm = 2 * _sqrt(math.pi) / self._W_y
        super().__init__(x, params, called)

    def profile(self, r):
        """..."""
        # beam profile distribution (field amplitude) at the waist of the beam
        if self.x == 0:
            return self._norm * _exp(-(r.y / self._W_y)**2)
        else:
            return super().profile(r)

    def spectrum(self, k_y):
        """Spectrum amplitude function, f."""
        return self._f_Gauss(k_y, self._W_y)

    def _f_Gauss(self, k_y, W_y):
        """Gaussian spectrum amplitude."""
        return _exp(-(k_y*W_y/2)**2)


class IncAiry2d(Beam2dCartesian):
    """2d incomplete Airy beam.

    The implementation is based on [2]_.

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
    .. [2] Ring, J D, Howls, C J, Dennis, M R: Incomplete Airy beams:
           finite energy from a sharp spectral cutoff , Optics Letters 38(10),
           OSA, 1639–1641, May 2013.

    """

    def __init__(self, x, params, called=False):
        """..."""
        self._W_y = params['W_y']
        self._M = params['M']
        self._W = params['W']
        super().__init__(x, params, called)

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
        The beam profile at waist is defined by the incomplete Airy function,
        see [2]_

        .. math::

           \psi^\text{Airy}_{M,W}(x, z=0) = \int_{M-W}^{M+W}\mathrm{d}\xi\, \
           \mathrm{exp}\left[\mathrm{i}\left(\frac{1}{3} \xi^3 + \
           \xi \frac{x}{w_0}\right)\right]

        """
        if self.x == 0:
            # adjust integration boundaries
            self._a = self._M-self._W
            self._b = self._M+self._W

        return super().profile(r)

    def spectrum(self, k_y):
        r"""Spectrum amplitude function, f.

        Parameters
        ----------
        k_y : type
            Description of parameter `k_y`.

        Returns
        -------
        f_Airy: complex
            Spectrum amplitude function

        Notes
        -----
        .. math::

           \begin{align*}
           f^\text{Airy}_{M,W}(k_x) &= \frac{1}{2 \pi} \int_{-\infty}^\infty\mathrm{d}x\, \psi^\text{Airy}_{M,W}(x, z=0) \exp(-\mathrm{i}k_x x)\\
                                    &= \int_{M-W}^{M+W}\mathrm{d}\xi\, \exp\Bigl(\mathrm{i}\frac{1}{3}\xi^3 \Bigr) \delta\Bigl(\frac{\xi}{w_0} -k_x\Bigr)\\
                                    &=\begin{cases}
                                       w_0 \exp\Bigl[\mathrm{i} \frac{1}{3} \bigl(w_0k_x \bigr)^3\Bigr] & M-W < w_0k_x <M+W \\
                                       0                                                                  & \text{otherwise}
                                      \end{cases}
           \end{align*}

        """
        return self._f_Airy(k_y, self._W_y, self._M, self._W)

    def _f_Airy(self, k_y, W_y, M, W):
        """Airy spectrum amplitude."""
        return W_y * _cexp(1.0j*(-1/3)*(k_y*W_y)**3) \
            * self._heaviside(W_y * k_y - (M - W)) \
            * self._heaviside((M + W) - W_y * k_y)

    def _heaviside(self, x):
        """Theta (Heaviside step) function."""
        return 0 if x < 0 else 1

    def _integrand(self, k_y):
        """..."""
        if self.x == 0:
            xi = k_y
            return _cexp(1.0j*(-(xi**3)/3 + (xi * self._ry/self._W_y)))
        else:
            # next line needs _integrand declared cpdef without nogil attribute,
            # and will execute slower than repeating the super class integration
            # routine function here
            #return super(IncAiry2d, self)._integrand(k_y)
            return self.spectrum(k_y) * _cexp(1j*self._phase(k_y, self.x, self._ry))
