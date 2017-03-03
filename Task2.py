import matplotlib.pyplot as plt
import numpy as np

Nt = 1000

for i in range(Nt):
	if i % 100 == 0:
		solution = np.loadtxt("Task2_Time_Step" + str(i) + ".txt")
		x = solution[:,0]
		y = solution[:,1]
		plt.plot(x,y, label='Solution at time step' + str(i))
		plt.hold(True)
	else:
		continue


plt.legend()
plt.title('Plot of Displacement in the Vertical Direction', fontsize=20)
plt.xlabel('Distance from Beginning of Beam [mm]', fontsize=14)
plt.ylabel('Displacement [mm] in Vertical Direction', fontsize=14)
plt.grid(True)
plt.show()
