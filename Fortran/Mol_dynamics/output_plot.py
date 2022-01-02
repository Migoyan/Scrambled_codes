# %% Imports
import numpy as np
import matplotlib.pyplot as plt

# %%
files = [
    'base_result/output_seq_nopara.log',
    'base_result/output_seq_E_4th.log',
    'base_result/output_seq_4th.log'
]
keys = ['no_para', 'E_4th', '4th']
colors = ['r', 'g', 'b']

dict_time = {}
dict_energies = {}

for file, key in zip(files, keys):
    time = []
    energies = []
    with open(file, 'r') as f:
        for i in range(6):
            f.readline()
        for line in f:
            line = line.split()
            time += [float(line[0])]
            energies += [line[1:]]
    dict_time[key] = np.array(time)
    dict_energies[key] = np.array(energies, dtype=float).T
    
print(energies)

# %% Potenial energies
plt.figure('E_pot')
plt.title('Potential energies')
plt.xlabel('Steps')
plt.ylabel('Energy')
for key, color in zip(keys, colors):
    plt.plot(dict_time[key], dict_energies[key][0], color, label = 'Epot ' + key)
plt.legend()

# %% Kinetic energies
plt.figure('E_kin')
plt.title('Kinetic energies')
plt.xlabel('Steps')
plt.ylabel('Energy')
for key, color in zip(keys, colors):
    plt.plot(dict_time[key], dict_energies[key][2], color, label = 'Ek ' + key)
plt.legend()

# %% Total energies
plt.figure('E_tot')
plt.title('Total energies')
plt.xlabel('Steps')
plt.ylabel('Energy')
for key, color in zip(keys, colors):
    plt.plot(dict_time[key], dict_energies[key][3], color, label = 'E ' + key)
plt.legend()

#%% Temperature
plt.figure('Temp')
plt.title('Temperature')
plt.xlabel('Steps')
plt.ylabel('Temperature (K)')
for key, color in zip(keys, colors):
    plt.plot(dict_time[key], dict_energies[key][1], color, label = 'T ' + key)
plt.legend()
plt.show()
# %%
