
import cython
import sys

if cython.compiled:
    from scipy import LowLevelCallable

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
    if cython.compiled:
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
