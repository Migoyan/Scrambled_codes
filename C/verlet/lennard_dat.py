import matplotlib.pyplot as plt
import numpy as np

dt = np.array([10., 5., 2., 1., .5, .2, .1, 5e-1, 2e-1, 1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 1e-6, 1e-7, 1e-8])

dE = []

with open("data.dat", "r") as file:
	for i in range(len(dt)):
		line = file.readline()
		data = line.split(';')
		del data[-1]
		data = [float(j) for j in data]
		dE += [float(np.mean(data))]

dE = np.array(dE)

dt = np.log(dt)
dE= np.log(dE)

plt.figure()
plt.title("Variation d'énergie moyenne lié à l'erreur en fonction du pas de temps")
plt.plot(dt, dE, "rx", label = "Variation moyenne d'énergie totale")
plt.xlabel("Ln(dt)")
plt.ylabel("ln(dE)")
plt.legend()
plt.show()