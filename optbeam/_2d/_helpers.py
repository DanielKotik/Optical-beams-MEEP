
import cython
import sys

if cython.compiled:
    from scipy import LowLevelCallable

from scipy.integrate import quad

from .beams import Beam2dCartesian


def _real_1d_func(x, func):
    """Return real part of a 1d function."""
    return func(x).real


def _imag_1d_func(x, func):
    """Return imag part of a 1d function."""
    return func(x).imag


def _imag_1d_func_c(n, arr, func_ptr):
    """Return imag part of a 1d function.

    Cython implementation.
    """
    # pure python formulation of:
    # return (<Beam2dCartesian>func_ptr)(arr[0]).imag
    return cython.cast(Beam2dCartesian, func_ptr)._integrand(arr[0]).imag


def _real_1d_func_c(n, arr, func_ptr):
    """Return real part of a 1d function.

    Cython implementation.
    """
    # pure python formulation of:
    # return (<Beam2dCartesian>func_ptr)(arr[0]).real
    return cython.cast(Beam2dCartesian, func_ptr)._integrand(arr[0]).real


def _complex_quad(func, a, b, kwargs={}):
    """Integrate real and imaginary part of the given function."""
    if cython.compiled:
        # pure python formulation of: cdef void *f_ptr = <void*>func
        f_ptr = cython.declare(cython.p_void, cython.cast(cython.p_void, func))

        func_capsule = PyCapsule_New(f_ptr, cython.NULL, cython.NULL)

        current_module = sys.modules[__name__]

        ll_real_1d_func_c = LowLevelCallable.from_cython(current_module,
                                                         '_real_1d_func_c',
                                                         func_capsule)
        ll_imag_1d_func_c = LowLevelCallable.from_cython(current_module,
                                                         '_imag_1d_func_c',
                                                         func_capsule)
        real, real_tol = quad(ll_real_1d_func_c, a, b, **kwargs)
        imag, imag_tol = quad(ll_imag_1d_func_c, a, b, **kwargs)
    else:
        real, real_tol = quad(_real_1d_func, a, b, (func,), **kwargs)
        imag, imag_tol = quad(_imag_1d_func, a, b, (func,), **kwargs)

    return real + 1j*imag, real_tol, imag_tol
