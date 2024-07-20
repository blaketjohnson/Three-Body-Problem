import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp
import sys

# Import constants from your specific path
sys.path.append('/Users/blakejohnson/Documents/Thesis/Three Body Problem')
from constants import *

"""
Part 1: Tritons Orbit using the two body problem
1. Using these values to determine Tritons Orbit
3. Use the Two Body Problem to find the Equations of Motion and plot Tritons orbit around Neptune

"""


""" 
1.1 Calculate the values for Tritons Elliptical orbit around Neptune

Using Orbital Mechanics for Engineering Students by Curtis, I calculated Tritons semi-major axis, semi-minor axis, angular momentum,
tangental velocity at perigee, specific orbital energy, and orbital period.

These values will help me plot Tritons Orbit in the next section

"""


# Triton's elliptical orbit

# (Curtis p. 63 2.21)
mu_2d = G * (m1 + m2)  # Gravitational parameter in m^3/s^2

# curtis p. 82 (2.73)
rp2_km = a2_km * (1 - e2)  # Radius at perigee in km
rp2_m = rp2_km * 1000  # Radius at perigee in meters

# curtis p. 82 (2.76)
b2_km = a2_km * math.sqrt(1 - e2**2)  # Semi-minor axis in km


#h2_p = math.sqrt(mu_2d * rp2_m * (1 + e2))  # Angular momentum of Triton in m^2/s

# curtis p. 81 (2.71)
h2_p = math.sqrt(mu_2d * a2_m * (1 - e2**2))  # Angular momentum of Triton in m^2/s


# curtis p 69 (2.31)
v2_p = h2_p / rp2_m  # Velocity of Triton at perigee in m/s

# curtis p.83 (2.80)
E2_sp = -mu_2d / (2 * a2_m)  # Specific orbital energy in J/kg

#curtis p.84 (2.83)
T2_orbit = (2 * math.pi / math.sqrt(mu_2d)) * (a2_m ** (3 / 2))  # Orbital period in seconds

# Output the results
print(f"Radius at perigee: {rp2_km:.2f} km")
print(f"Radius at perigee: {rp2_m:.2f} m")
print(f"Semi-minor axis: {b2_km:.2f} km")
print(f"Angular momentum: {h2_p:.2e} m^2/s")
print(f"Velocity at perigee: {v2_p:.2f} m/s")
print(f"Specific orbital energy: {E2_sp:.2e} J/kg")
print(f"Orbital period: {T2_orbit:.2f} s")


"""  

1.2 Calculate the Two Body Equations of Motion and plot Tritons orbit

"""

"""
Two Body Equation of Motion:
We know that the time deriviative of a particles position is its tangental velocity.
Additionally we know that the time derivative of a particles velocity is its acceleration.
Using pythagorean theorm we can say that position is a function of x,y,z with: r**2 = x**2 + y**2 + z**2
Additionally the equation for acceleeration can be found on curtis p. 23 (2.5)
We can use these values to create a variable that will find the velocity and acceleration of Triton

"""

"""
Defining the equations of motion:

Lets define the function dfdt since the function is returning the time derivatives
I set the inputs to be of a variable f and t
f represents the cartesian positions and velocities of Triton
t represents time

Then I set the equations of motion for the two body problem:
 I use the pythagorean theorem to define r
 I said that the time derivatives for x, y, and z are the cartesian velocities
 I set the variables for acceleration equal to equation 2.5 from murray INCLUDING the unit vector

The function will return the velocities and accelerations in Cartesian Coordinates over the given time.

"""


# Define the two-body problem differential equations
def dfdt(f, t):
    x, y, z, vx, vy, vz = f
    r = np.sqrt(x**2 + y**2 + z**2)
    dx = vx
    dy = vy
    dz = vz
    
    #murray p.23 (2.5)
    ddx = (-mu_2d / r**3) * x
    ddy = (-mu_2d / r**3) * y
    ddz = (-mu_2d / r**3) * z
    return [dx, dy, dz, ddx, ddy, ddz]

"""

Now that the function is defined and the equations are set up, I created some initial conditions.
I used the values when Triton is at perigeee as the starting values.
time at the start is 0s.
Since Triton is at perigee the radius is entirely in the x-direction. I set the position as a vector < a,0,0 > (this is why it is an array)
Additionally the velocity of Triton is entirley in the y direction when it is at perigee, so the veloctiy is <0,v,0>

"""

# Initial conditions
t0 = 0  # Set initial time to 0
r0 = np.array([rp2_km * 1000, 0, 0])  # Initial position in meters
v0 = np.array([0, v2_p, 0])  # Initial velocity in m/s


"""
Next I insert the intial conditions into the variable f so that I can run it in the function. I list f0 to represent intial conditions.
The first variables in my position and velocity array's above represent the variable in the x-axis. These are the first numbers in the array
which python represents as [0]. Similarly variables in the y-axis are [1] and z-axis are [2]. 

Below I am stating that f0 is a combination of my inital conditions, r0[0] = initial position along the x-axis, r0[1] = initial condition in the y-axis,etc.
"""

# Initial conditions for odeint: [x, y, z, vx, vy, vz]
f0 = [r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]]


"""

Determine the time span for the orbit:
Let tf be the time we are measuring and I want multiple orbits to ensure a fully complete plot.
I used Tritons orbital period in seconds T2_s and multiplied by 10 to run 10 full orbits
dt is the time step in seconds which I set to 100 so that there is consistant readings without overstressing python.
tspan defines the integration. I am saying that I want to evaluate the orbit from 0-tf and I want to take a measurement every 100s.
"""

# Time span for the integration
tf = T2_s   # Integrate for multiple orbits
dt = 100  # Time step for integration
tspan = np.arange(t0, tf, dt)

"""
Integrate the Equations of motion:
I used the odeint function to integrate the equations in the dfdt function above.
I am using intial conditions f0, and I am integrating the values for a range of tspan.

"""

# Integrate the equations of motion
sol = odeint(dfdt, f0, tspan)


# Convert position from meters to kilometers for plotting
sol_km = sol[:, :3] / 1000


"""
Plot the results:

Finally I plotted the results. I set the position in the x-axis as sol_km[0] and y-position as sol_km[1].
I set Triton and its orbit to be a brown color and Neptune to be a Blue Color. I want this to be consistent during the entire project
I added axis and units and a legend in the bottom.
"""
# Plotting the trajectory
plt.figure()
plt.plot(sol_km[:, 0], sol_km[:, 1], color='saddlebrown', label="Triton's Orbit")
plt.scatter(0, 0, color='blue', s=50, label='Neptune')  # Neptune at the origin
plt.scatter(sol_km[0, 0], sol_km[0, 1], color='saddlebrown', s=30, label='Triton at Perigee')  # Triton at initial position

# Plotting settings
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.legend(loc='upper right')
plt.title("Triton's Orbit around Neptune")
plt.grid()
plt.axis()
plt.show()




