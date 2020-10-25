
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D


ear = list(open("earth_verlet_case5_2.000000_250_5.txt").readlines())[1:]
mer = list(open("mercury_verlet_case5_2.000000_250_5.txt").readlines())[1:]
sun = list(open("sun_verlet_case5_2.000000_250_5.txt").readlines())[1:]
ven = list(open("venus_verlet_case5_2.000000_250_5.txt").readlines())[1:]
mar = list(open("mars_verlet_case5_2.000000_250_5.txt").readlines())[1:]
ura = list(open("uranus_verlet_case5_2.000000_250_5.txt").readlines())[1:]
sat = list(open("saturn_verlet_case5_2.000000_250_5.txt").readlines())[1:]
jup = list(open("jupiter_verlet_case5_2.000000_250_5.txt").readlines())[1:]
nep = list(open("neptune_verlet_case5_2.000000_250_5.txt").readlines())[1:]
plu = list(open("pluto_verlet_case5_2.000000_250_5.txt").readlines())[1:]


def return_values(n):
	return([list(map(float, ear[n].split()))[1:], list(map(float, sun[n].split()))[1:], list(map(float, nep[n].split()))[1:], list(map(float, jup[n].split()))[1:],\
		list(map(float, mar[n].split()))[1:], list(map(float, ven[n].split()))[1:], list(map(float, mer[n].split()))[1:], list(map(float, sat[n].split()))[1:],\
		list(map(float, ura[n].split()))[1:], list(map(float, plu[n].split()))[1:]])


values = return_values(0)

fig, ax = plt.subplots(figsize=(14,13))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim3d([-3, 3]); ax.set_xlabel('X'); ax.set_ylim3d([-3, 3]); ax.set_ylabel('Y'); ax.set_zlim3d([-0.25, 0.25]); ax.set_zlabel('Z')
ax.set_facecolor('xkcd:black')
ax.view_init(elev=+40., azim=+120)
ax.set_axis_off()
ax.axis('off')

scatters = [ax.scatter(values[1][0], values[1][1], values[1][2], color='yellow', s=500, alpha= 0.95), ax.scatter(values[0][0], values[0][1], values[0][2], color='aqua', s=5, alpha= 0.75)\
	, ax.scatter(values[2][0], values[2][1], values[2][2], color='royalblue', s=20, alpha= 0.75) , ax.scatter(values[3][0], values[3][1], values[3][2], color='burlywood', s=50, alpha= 0.75)\
	, ax.scatter(values[4][0], values[4][1], values[4][2], color='orangered', s=3, alpha= 0.75), ax.scatter(values[5][0], values[5][1], values[5][2], color='tan', s=5, alpha= 0.75) \
	, ax.scatter(values[6][0], values[6][1], values[6][2], color='slategray', s=1.5, alpha= 0.75), ax.scatter(values[7][0], values[7][1], values[7][2], color='peachpuff', s=40, alpha= 0.75)\
	, ax.scatter(values[8][0], values[8][1], values[8][2], color='lightblue', s=25, alpha= 0.75), ax.scatter(values[9][0], values[9][1], values[9][2], color='brown', s=1, alpha= 0.75)]


def plot(i):
	values = return_values(i)
	scatters[1]._offsets3d = ([values[0][0]], [values[0][1]], [values[0][2]])
	scatters[0]._offsets3d = ([values[1][0]], [values[1][1]], [values[1][2]])
	scatters[2]._offsets3d = ([values[2][0]], [values[2][1]], [values[2][2]])
	scatters[3]._offsets3d = ([values[3][0]], [values[3][1]], [values[3][2]])
	scatters[4]._offsets3d = ([values[4][0]], [values[4][1]], [values[4][2]])
	scatters[5]._offsets3d = ([values[5][0]], [values[5][1]], [values[5][2]])
	scatters[6]._offsets3d = ([values[6][0]], [values[6][1]], [values[6][2]])
	scatters[7]._offsets3d = ([values[7][0]], [values[7][1]], [values[7][2]])
	scatters[8]._offsets3d = ([values[8][0]], [values[8][1]], [values[8][2]])
	scatters[9]._offsets3d = ([values[9][0]], [values[9][1]], [values[9][2]])
	print(i)
	return scatters


iterations = None
ani = animation.FuncAnimation(fig, plot, iterations, interval=20, blit=False, repeat=True)
# ani.save('Solar.gif', dpi=80, writer='imagemagick')
plt.show()
