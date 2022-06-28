""" 1D Oscillators system resolution
"""
import logging
import random as r
from math import ceil
from enum import Enum

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from cycler import cycler
from scipy.fft import fftshift, rfft, rfftfreq

class LimitCondition(Enum):
    """Enum for the limit conditions.
    """
    FREE = 0
    STATIC = 1
    PERIODIC = 2

def randomize_positions(index_start: int, index_end: int, sigma: float):
    """Gaussian randomisation of the positions of the particles.

    Args:
        index_start (int): index after which position will be randomized.
        index_end (int): index before which position will be randomized.
        sigma (float): standard deviation of the gaussian distribution.
    """
    for i in range(index_start, index_end):
        positions[i] += r.gauss(0,sigma)

logging.basicConfig() # Configure logging


##---------------------------------------Parameters---------------------------------------
logger = logging.getLogger(__name__) # Create a logger for our main file
LIMIT = LimitCondition.FREE
NB_MASS = 20  # Number of masses
D0 = 1 # Initial distance between masses (just for plotting)
dt = 1e-3 # Time step
NB_PAS = int(5e4) # Number of time steps
K = .5 # Spring constant
logger.setLevel(logging.DEBUG) # Set your logger level (DEBUG, INFO, WARNING, ERROR, CRITICAL)

FILE = "results.dat" # Name of thefile to save the results

##---------------------------------------Initialization---------------------------------------
dt2 = dt**2 # dt^2
t = 0 # Time
masses = np.ones(NB_MASS)

# Loop to set mass of 2
for i in range(0, NB_MASS, 2):
    masses[i] = 1
inverse_masses = 1/masses # Inverse of masses, multiplication is faster than division

times = np.zeros(NB_PAS)
positions = np.zeros(NB_MASS)
velocities = np.zeros(NB_MASS)
randomize_positions(2, 3, 5e-1)
logger.debug(f"Initial particles positions: {positions}\n")

logger.debug(f"Initial particles masses: {masses}")

sum_forces = np.zeros(NB_MASS)
sum_forces_1 = np.zeros(NB_MASS)

result = np.zeros((NB_PAS, NB_MASS)) # Array to store the postions at each time step

## Boundary conditions
# m for minus (particles on the left), p for plus (particles on the right)
if LIMIT == LimitCondition.PERIODIC:
    positions_p = np.concatenate((positions[1:], [positions[0]]), axis = None)
    positions_m = np.concatenate(([positions[-1]], positions[:-1]), axis = None)
elif LIMIT == LimitCondition.STATIC:
    ## Static boundary conditions, (like a wall)
    # We set right and left boundaries to 0
    positions_p = np.concatenate((positions[1:], [0]), axis = None)
    positions_m = np.concatenate(([0], positions[:-1]), axis = None)
elif LIMIT == LimitCondition.FREE:
    ## Free boundary conditions
    # We set right and left boundaries to virtual particles that follow the same motion as the first and last particles
    positions_p = np.concatenate((positions[1:], [positions[-1]]), axis = None)
    positions_m = np.concatenate(([positions[0]], positions[:-1]), axis = None)

sum_forces[:] = K*(positions_p[:] - 2 * positions[:] + positions_m[:]) * inverse_masses
logger.info(f"Starting simulation with {NB_PAS} time steps")

# Main loop
for i in range(NB_PAS):
    if i % ceil(NB_PAS/10) == 0:
        # Simulation progress
        logger.info(f"{i} time steps done. {(i/NB_PAS)*100:.2f}%")

    t += dt
    times[i] = t

    ## Velocity-Verlet algorithm
    positions[:] += dt * velocities[:] + .5 * dt2 * sum_forces[:]
    result[i] = positions

    # Boundary conditions
    if LIMIT == LimitCondition.PERIODIC:
        positions_p = np.concatenate((positions[1:], [positions[0]]), axis = None) # Positions of the particles on the right side
        positions_m = np.concatenate(([positions[-1]], positions[:-1]), axis = None) # Positions of the particles on the left side
    elif LIMIT == LimitCondition.STATIC:
        positions_p = np.concatenate((positions[1:], [0]), axis = None)
        positions_m = np.concatenate(([0], positions[:-1]), axis = None)
    elif LIMIT == LimitCondition.FREE:
        ## Free boundary conditions
        # We set right and left boundaries to virtual particles that follow the same motion as the first and last particles
        positions_p = np.concatenate((positions[1:], [positions[-1]]), axis = None)
        positions_m = np.concatenate(([positions[0]], positions[:-1]), axis = None)

    sum_forces_1[:] = K*(positions_p[:] - 2 * positions[:] + positions_m[:]) * inverse_masses
    velocities[:] += .5 * dt * (sum_forces[:] + sum_forces_1[:])
    sum_forces[:] = sum_forces_1[:]

logger.info(f"Simulation done. {(i/NB_PAS)*100:.2f}%")

for i in range(NB_MASS):
    result[:, i] += D0 * i # Add the initial distance between masses

with open(FILE, "w") as f:
    for i in range(NB_PAS):
        f.write(f"{times[i]:.6f}, ")
        for j in range(NB_MASS):
            f.write(f"{result[i, j]:.8f}, ")
        f.write("\n")

mpl.rcParams['axes.prop_cycle'] = cycler(color=['r', 'g', 'b', 'k'])
mpl.rcParams['figure.figsize'] = (7, 7)
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['axes.labelsize'] = 35
mpl.rcParams['legend.fontsize'] = 15
mpl.rcParams['lines.linewidth'] = 2

# Plotting positions
fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(times, result, '-', label = 'Position')
ax.set_xlabel('Time')
ax.set_ylabel('Position')
plt.tight_layout()

# Plotting frequencies
fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(fftshift(rfftfreq(NB_PAS, d=dt)), rfft(result[:,0]), '-', label = 'Position')
ax.set_xlim(240, 260)
ax.set_xlabel('Frequence')
ax.set_ylabel('Amplitude')
plt.tight_layout()
plt.show()
