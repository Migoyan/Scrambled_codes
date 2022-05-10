from re import T
import numpy as np
import matplotlib.pyplot as plt

NB_MASS = 4
D0 = 2e-1
dt = 1e-3
dt2 = dt**2
t = 0
NB_PAS = int(5e5)
K = 1
L = 1

times = np.zeros(NB_PAS)
positions = np.linspace(0, NB_MASS * D0, NB_MASS)
velocities = np.zeros(NB_MASS)
sum_forces = np.zeros(NB_MASS)
sum_forces_1 = np.zeros(NB_MASS)

result = np.zeros((NB_PAS, NB_MASS))

positions_p = np.concatenate((positions[1:], [(NB_MASS+1) * D0]), axis = None)
positions_m = np.concatenate(([-D0], positions[:-1]), axis = None)

sum_forces[:] = K*(positions_p[:] - 2 * positions[:] + positions_m[:])
for i in range(NB_PAS):
    t += dt
    times[i] = t

    # Velocity-Verlet algorithm
    positions[:] += dt * velocities[:] + .5 * dt2 * sum_forces[:]
    result[i] = positions

    positions_p = np.concatenate((positions[1:], [(NB_MASS+1) * D0]), axis = None)
    positions_m = np.concatenate(([-D0], positions[:-1]), axis = None)
    
    sum_forces_1[:] = K*(positions_p[:] - 2 * positions[:] + positions_m[:])
    velocities[:] += .5 * dt * (sum_forces[:] + sum_forces_1[:])
    sum_forces[:] = sum_forces_1[:]

fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(times, result, '-', label = 'Position')
plt.show()
