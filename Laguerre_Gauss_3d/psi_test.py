"""
Check Meep's integration routine.
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from scipy.integrate import dblquad

def complex_dblquad(func, a, b, gfun, hfun, **kwargs):
    """Integrates the real and imaginary part of the given function."""
    def real_func(x, y):
        return sp.real(func(x, y))
    def imag_func(x, y):
        return sp.imag(func(x, y))
    
    def real_integral():
        return dblquad(real_func, a, b, gfun, hfun, **kwargs)
    def imag_integral():
        return dblquad(imag_func, a, b, gfun, hfun, **kwargs)
    
    return (real_integral()[0] + 1j * imag_integral()[0], real_integral()[1:], imag_integral()[1:])


## test paramters
freq  = 4.0
k_vac = 2 * np.pi * freq
w_0   = 8.0 / (2 * np.pi * 4)

def f_Gauss(W_y, k_y, k_z): 
    """
    """
    return sp.exp(-(W_y ** 2.0) * (k_y ** 2.0 + k_z ** 2.0) / 4.0)

def integrand(x, y, z, k, k_y, k_z):
    """
    """
    ## first variant (taking the real part of the suqare root)
    return f_Gauss(w_0, k_y, k_z) * sp.exp(1.0j * (x * (sp.sqrt(k**2 - k_y**2 - k_z**2) + y * k_y + z * k_z).real))
    
    ## second variant (leave square root as is, but perform second integration with non-constant bounds)
    #return f_Gauss(w_0, k_y, k_z) * sp.exp(1.0j * (x * (sp.sqrt(k**2 - k_y**2 - k_z**2) + y * k_y + z * k_z)))

def psi(x, y, z, k):
    """
    """
    integrand_ = lambda k_y, k_z: integrand(x, y, z, k, k_y, k_z)
    
    ## constant integration bounds (appropriate for first variant)
    return complex_dblquad(integrand_, -k, k, lambda x: -k, lambda x: k)[0]
    
    ## non-constant integration bounds (appropriate for second variant)
    #return complex_dblquad(integrand_, -k, k, lambda x: -np.sqrt(k**2 - x**2), lambda x: np.sqrt(k**2 - x**2))[0]
    
    
x_shift = -2.0
print "integrand: ", integrand(x_shift, 0.3, 0, k_vac, 4, 0)
print "psi:       ", psi(x_shift, 0, 0, k_vac)

K_y = np.linspace(-k_vac, k_vac, 200)

#INTEGRAND = integrand(x_shift, 0.0, 0.0, k_vac, K_y, 0.0)
#plt.plot(K_y, INTEGRAND.real)
#plt.plot(K_y, INTEGRAND.imag)
#plt.show()