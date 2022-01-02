import matplotlib.pyplot as plt
import numpy as np

x_n = []
x_n_1 = []

with open("./result/rng.dat", 'r') as file:
    line = " "
    i = 0
    while line != "":
        i += 1
        line = file.readline()
        data = line.split(';')
        if data[0] != "":
            x_n += [float(data[0])]
            x_n_1 += [float(data[1])]

x_n = np.array(x_n)
x_n_1 = np.array(x_n_1)

plt.figure()
plt.xlabel('X_n')
plt.ylabel('X_(n + 1)')
plt.title('X_n')
plt.plot(x_n, x_n_1, '.')
plt.show()
