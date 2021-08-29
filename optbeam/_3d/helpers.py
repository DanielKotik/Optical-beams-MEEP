
import math
import sys

try:
    import cython
    cython_imported = True
except ModuleNotFoundError:
    cython_imported = False

if cython_imported:
    if cython.compiled:
        from scipy import LowLevelCallable
    else:
        print("\nPlease consider compiling `%s.py` via Cython: "
              "`$ cythonize -3 -i %s.py`\n" % (__name__, __name__))

if not cython_imported or not cython.compiled:
    from math import (sin as _sin,
                      cos as _cos)
    from cmath import (exp as _cexp,
                       sqrt as _csqrt)

from scipy.integrate import dblquad
from types import MappingProxyType


def _real_2d_func(x, y, func):
    """Return real part of a 2d function."""
    return func(x, y).real


def _imag_2d_func(x, y, func):
    """Return imag part of a 2d function."""
    return func(x, y).imag


def _imag_2d_func_c(n, arr, func_ptr):
    """Return imag part of a 2d function.

    Cython implementation.
    """
    # pure python formulation of:
    # return (<Beam3d>func_ptr)(arr[0], arr[1]).imag
    return cython.cast(Beam3d, func_ptr)._integrand(arr[0], arr[1]).imag


def _real_2d_func_c(n, arr, func_ptr):
    """Return real part of a 2d function.

    Cython implementation.
    """
    # pure python formulation of:
    # return (<Beam3d>func_ptr)(arr[0], arr[1]).real
    return cython.cast(Beam3d, func_ptr)._integrand(arr[0], arr[1]).real


def _complex_dblquad(func, a, b, gfun, hfun, kwargs={}):
    """Integrate real and imaginary part of the given function."""
    if cython_imported and cython.compiled:
        # pure python formulation of: cdef void *f_ptr = <void*>func
        f_ptr = cython.declare(cython.p_void, cython.cast(cython.p_void, func))

        func_capsule = PyCapsule_New(f_ptr, cython.NULL, cython.NULL)

        current_module = sys.modules[__name__]

        ll_real_2d_func_c = LowLevelCallable.from_cython(current_module,
                                                         '_real_2d_func_c',
                                                         func_capsule)
        ll_imag_2d_func_c = LowLevelCallable.from_cython(current_module,
                                                         '_imag_2d_func_c',
                                                         func_capsule)
        real, real_tol = dblquad(ll_real_2d_func_c, a, b, gfun, hfun, **kwargs)
        imag, imag_tol = dblquad(ll_imag_2d_func_c, a, b, gfun, hfun, **kwargs)
    else:
        real, real_tol = dblquad(
            _real_2d_func, a, b, gfun, hfun, (func,), **kwargs)
        imag, imag_tol = dblquad(
            _imag_2d_func, a, b, gfun, hfun, (func,), **kwargs)

    return real + 1j*imag, real_tol, imag_tol


class Beam3d:
    """Abstract base class."""

    def __init__(self, x, params, called=False):
        """..."""
        self.x = x   # TODO: rename x to x_shift
        self._k = params['k']
        self._params = MappingProxyType(params)  # read-only view of a dict
        self.called = called

    @property
    def params(self):
        """Beam specific parameters.

        This is a read-only property.
        """
        return self._params

    def _integrand(self, x, y):
        """Integrand function over two coordinates x and y."""
        raise NotImplementedError


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
