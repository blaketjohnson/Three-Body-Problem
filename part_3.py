import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp
import sys

# Import equations from parts 1 and 2
sys.path.append('/Users/blakejohnson/Documents/Thesis/Three Body Problem')
from part_2 import *

"""
Solving for the Jacobi Constant
I defined a new function using the same equations and inputs from the three body equations of motion.
However this time I set the function to return the Jacobi Constant
"""
def jacobi_constant(state, mu):
    x, y, z, vx, vy, vz = state
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    U = 0.5 * (x**2 + y**2 + z**2) + (1 - mu) / r1 + mu / r2

    # murray p. 68 (3.28)
    C = 2*U - (vx**2 + vy**2 + vz**2)
    return C

# Calculate Jacobi constant for initial conditions
C = jacobi_constant(state0, mu)
print(f"Jacobi constant: {C}")
