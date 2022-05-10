from re import T
import random as r
import numpy as np
import matplotlib.pyplot as plt

def randomize_positions(i_0, i_1):
    for i in range(i_0, i_1):
        positions[i] += r.gauss(0,2.5e-1)

NB_MASS = 100
D0 = 5e-1
dt = 1e-3
dt2 = dt**2
t = 0
NB_PAS = int(5e4)
K = 1
L = 1

times = np.zeros(NB_PAS)
positions = np.linspace(0, NB_MASS * D0, NB_MASS, endpoint=False)
velocities = np.zeros(NB_MASS)
randomize_positions(int(NB_MASS/2), int(NB_MASS/2)+1)
print(f"Positions des particules: {positions}")

sum_forces = np.zeros(NB_MASS)
sum_forces_1 = np.zeros(NB_MASS)

result = np.zeros((NB_PAS, NB_MASS))

positions_p = np.concatenate((positions[1:], [positions[0] + NB_MASS * D0]), axis = None)
positions_m = np.concatenate(([positions[-1] - NB_MASS * D0], positions[:-1]), axis = None)

sum_forces[:] = K*(positions_p[:] - 2 * positions[:] + positions_m[:])
for i in range(NB_PAS):
    t += dt
    times[i] = t

    # Velocity-Verlet algorithm
    positions[:] += dt * velocities[:] + .5 * dt2 * sum_forces[:]
    result[i] = positions

    positions_p = np.concatenate((positions[1:], [positions[0] + NB_MASS * D0]), axis = None)
    positions_m = np.concatenate(([positions[-1] - NB_MASS * D0], positions[:-1]), axis = None)
    
    sum_forces_1[:] = K*(positions_p[:] - 2 * positions[:] + positions_m[:])
    velocities[:] += .5 * dt * (sum_forces[:] + sum_forces_1[:])
    sum_forces[:] = sum_forces_1[:]

fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(times, result, '-', label = 'Position')
plt.show()
