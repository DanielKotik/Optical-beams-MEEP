cimport cython
from cpython.pycapsule cimport PyCapsule_New
from cpython cimport bool

# -----------------------------------------------------------------------------
# function declarations
# -----------------------------------------------------------------------------
cdef double _imag_2d_func_c(int n, double *arr, void *func_ptr)
cdef double _real_2d_func_c(int n, double *arr, void *func_ptr)

@cython.locals(real=cython.double, imag=cython.double, real_tol=cython.double,
               imag_tol=cython.double)
cpdef (double complex, double, double) _complex_dblquad(Beam3d func,
                                                       double a, double b,
                                                       double gfun, double hfun,
                                                       dict kwargs=*)

# -----------------------------------------------------------------------------
# class declarations
# -----------------------------------------------------------------------------
cdef class Beam3d:
    cdef:
        dict __dict__
        double x, _k
        public bool called

    cdef double complex _integrand(self, double x, double y) nogil
