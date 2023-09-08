
import numpy  as  np 
# import mpmath as mp
# import sympy  as smp
from scipy.optimize import fsolve


# ------------------------------------------------- ------------------- 
# Define relevant constants for use in calculations
# ------------------------------------------------- -------------------

sig_sb = 5.6703e-5 # Stefan-Boltzmann const  [erg/cm^2/K^4/s]
m_sun = 1.9884e33 # Sun mass [g] 
temp_sun = 5.777e3 # effective temperature of the sun [K] 
l_s = 3.828e33 # solar intensity [erg/s] defined by IAU 
r_sol = 6.96e10 # solar radius [cm] 
G = 6.67430e-8 # Universal gravitational constant changing.... 
mu = 2.353 # Average molecular weight 
k_b = 1.380649e-16 # Boltzmann constant defined
m_p = 1.672621898e-24 # proton mass [g]
x_rim = 2.0  # Constant in range 2 - 6

K, T_rim, M, L = [k_b, 1500.0, 0.5 * m_sun, 0.7 * l_s]

# ------------------------------------------------------------------------------
# Mathematical Calculations now that constants are defined
# ------------------------------------------------------------------------------

h_const = K * T_rim / (mu * m_p * G * M)
r_const = np.sqrt(0.25 * L / (np.pi  * (temp_sun ** 4) * sig_sb))


def r_rim(r_i):
    """Set equation for R_rim."""
    return (r_const * np.sqrt(1 + x_rim * np.sqrt(h_const * r_i))) - r_i


# Create initial guess
r_0 = 0.01 * r_sol

# Numerically solve for R_rim
eff_r_rim = fsolve(func=r_rim, x0=r_0, factor=100)

# Calculate h_rim using eff_r_rim, and subsequent H_rim
h_rim = np.sqrt(h_const * (eff_r_rim ** 3))
H_rim = x_rim * h_rim