import matplotlib.pyplot as plt
import numpy as np

# Nt = 10000

# for i in range(Nt):
# 	if i % (Nt/10) == 0:
# 		solution = np.loadtxt("Task2_Time_Step" + str(i) + ".txt")
# 		x = solution[:,0]
# 		y = solution[:,1]
# 		plt.plot(x,y, label='Solution at time step' + str(i))
# 		plt.hold(True)
# 	else:
# 		continue

solution = np.loadtxt("Task2.txt")
x = solution[:,0]
y = solution[:, 1]
plt.plot(x,y)


# plt.title('Displacment of point L/2 over time', fontsize=20)
# plt.xlabel('Time Steps [x/10000 seconds]', fontsize=14)
# plt.ylabel('Displacement in Vertical Direction [mm]', fontsize=14)
# plt.grid(True)
plt.show()
