from .helpers cimport Beam3d, _complex_dblquad

# -----------------------------------------------------------------------------
# declare C functions as "cpdef" to export them to the module
# -----------------------------------------------------------------------------
# TODO: check why we use cpdef and not cdef
cdef extern from "stdlib.h":
    cpdef int _abs "abs" (int n) nogil

cdef extern from "math.h":
    cpdef double _sin "sin" (double x) nogil
    cpdef double _cos "cos" (double x) nogil
    cpdef double _exp "exp" (double x) nogil
    cpdef double _acos "acos" (double x) nogil
    cpdef double _atan2 "atan2" (double y, double x) nogil

cdef extern from "complex.h":
    cpdef double complex _cexp "cexp" (double complex z) nogil
    cpdef double complex _csqrt "csqrt" (double complex z) nogil

# -----------------------------------------------------------------------------
# function declarations
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# class declarations
# -----------------------------------------------------------------------------
cdef class Beam3dSpherical(Beam3d):
    cdef:
        double _ry, _rz

    cdef double _phase(self, double sin_theta, double cos_theta, double phi,
                       double x, double y, double z) nogil
    cdef double complex spectrum(self, double sin_theta, double theta, double phi) nogil
    cdef double complex _integrand(self, double theta, double phi) nogil

cdef class Beam3dCartesian(Beam3d):
    cdef:
        double _ry, _rz

    cdef double _phase(self, double k_y, double k_z, double x, double y, double z) nogil
    cdef double complex spectrum(self, double k_y, double k_z) nogil
    cdef double complex _integrand(self, double k_y, double k_z) nogil

cdef class LaguerreGauss3d(Beam3dSpherical):
    cdef:
      int _m
      double _W_y
      double _norm

    cdef double complex _f_Gauss_spherical(self, double sin_theta, double _W_y, double k) nogil
    cdef double complex _f_Laguerre_Gauss_spherical(self, double sin_theta, double theta, double phi,
                                                    double _W_y, double k, int _m) nogil
    cdef double complex spectrum(self, double sin_theta, double theta, double phi) nogil

cdef class LaguerreGauss3dCartesian(Beam3dCartesian):
    cdef:
      int _m
      double _W_y
      double _norm

    cdef double _phi(self, double k_y, double k_z) nogil
    cdef double _theta(self, double k_y, double k_z, double k) nogil
    cdef double complex _f_Gauss_cartesian(self, double k_y, double k_z, double _W_y) nogil
    cdef double complex _f_Laguerre_Gauss_cartesian(self, double k_y, double k_z,
                                                    double _W_y, double k, int _m) nogil
    cdef double complex spectrum(self, double k_y, double k_z) nogil
