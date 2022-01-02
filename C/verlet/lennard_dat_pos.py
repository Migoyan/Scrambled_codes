import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

nb_p = 0
time = []
x, y = [], []

with open("pos_data.dat", "r") as file:
	line = " "
	while line != "":
		x_t, y_t = [], []
		line = file.readline()
		if line == "":
			break
		data = line.split(';')
		time += [float(data[-1])]
		del data[-1]
		for i in range(len(data)):
			if i%2 == 0:
				x_t += [float(data[i])]
			else:
				y_t += [float(data[i])]
		x += [x_t]
		y += [y_t]

fig = plt.figure()
plt.title("Positions des particules")
plt.xlim(-1, 17)
plt.ylim(-1, 17)
plt.xlabel("x")
plt.ylabel("y")

line, = plt.plot([],[], ".")
def init():
    line.set_data([],[])
    return line,

def animate(n): 
	line.set_data(x[n], y[n])
	return line,
 
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(time), blit=True, interval=20, repeat=False)
anim.save("particle.gif", writer='PillowWriter')

plt.show()