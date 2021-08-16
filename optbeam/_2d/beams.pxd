from .helpers cimport Beam2d, _complex_quad

# -----------------------------------------------------------------------------
# declare C functions as "cpdef" to export them to the module
# -----------------------------------------------------------------------------
cdef extern from "math.h":
    cpdef double _exp "exp" (double x) nogil
    cpdef double _sqrt "sqrt" (double x) nogil

cdef extern from "complex.h":
    cpdef double complex _cexp "cexp" (double complex z) nogil

# -----------------------------------------------------------------------------
# function declarations
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# class declarations
# -----------------------------------------------------------------------------
cdef class Beam2dCartesian(Beam2d):
    cdef:
        double _ry, _rz

    cdef double complex spectrum(self, double k_y) nogil
    cdef double _phase(self, double k_y, double x, double y) nogil
    cdef double complex _integrand(self, double k_y) nogil

cdef class Gauss2d(Beam2dCartesian):
    cdef:
        double _W_y
        double _norm

    cdef double _f_Gauss(self, double k_y, double W_y) nogil
    cdef double complex spectrum(self, double k_y) nogil

cdef class IncAiry2d(Beam2dCartesian):
    cdef:
        double _W_y
        double _M
        double _W
        double xi

    cdef double _heaviside(self, double x) nogil
    cdef double complex _f_Airy(self, double k_y, double W_y, double M, double W) nogil
    cdef double complex spectrum(self, double k_y) nogil
