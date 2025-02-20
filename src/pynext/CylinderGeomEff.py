import numpy as np
import math
from scipy.integrate import dblquad

def barrel_detection_efficiency(R, L):
    """
    Compute the fraction (eta) of isotropically emitted photons 
    that first intersect the barrel (cylindrical side at r=R),
    rather than an end cap at z=±L.
    
    Based on the 2D integral:
      eta = (1 / (R^2 * L)) ∫_{z=-L to L} ∫_{r=0 to R} 
                r * f_side(r,z) dr dz,
    where 
      f_side(r,z) = 1/2 * [1 / sqrt(1 + ((R-r)/(L - z))^2)]
                   + 1/2 * [1 / sqrt(1 + ((R-r)/(L + z))^2)].
    """
    
    def f_side(r, z):
        """
        Fraction of directions (out of 4π) at point (r,z)
        that hit barrel first, split into upward/downward hemispheres.
        """
        # Upward fraction (cos(theta) > 0)
        # Avoid division by zero if L == z
        denom_up = (L - z)
        if abs(denom_up) < 1e-14:
            f_up = 0.0
        else:
            arg_up = (R - r) / denom_up
            f_up = 0.5 / math.sqrt(1.0 + arg_up*arg_up)
        
        # Downward fraction (cos(theta) < 0)
        # Avoid division by zero if L == -z
        denom_down = (L + z)
        if abs(denom_down) < 1e-14:
            f_down = 0.0
        else:
            arg_down = (R - r) / denom_down
            f_down = 0.5 / math.sqrt(1.0 + arg_down*arg_down)
        
        return f_up + f_down

    def integrand(z, r):
        # dblquad expects the integrand as integrand(z, r)
        # where r is the outer variable, z is the inner variable.
        # We'll just call f_side(r,z) and multiply by r.
        return r * f_side(r, z)

    # Perform the double integral:
    # r in [0, R], z in [-L, L].
    # dblquad syntax: dblquad(func, r_min, r_max, z_min, z_max)
    #  BUT note the function signature is func(z, r).
    # We integrate w.r.t. z first (inner integral), then r (outer integral).
    res, err = dblquad(
        integrand,
        0,   # r_min
        R,   # r_max
        lambda r: -L,  # z_min
        lambda r:  L   # z_max
    )

    # Multiply by the factor 1 / (R^2 * L)
    eta = res / (R**2 * L)
    return eta

