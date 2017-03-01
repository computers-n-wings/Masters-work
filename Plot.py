import matplotlib.pyplot as plt
import numpy as np

solution = np.loadtxt('displacements.txt')
x = solution[:,0]
y = solution[:,1]

qy = 1.0
L = 10000.0
E = 2.1e+06
I = 1.44e+07
Fy = 1.0e+03

u = np.zeros_like(x)

for i in range(len(x)):
	if x[i] <= L/2:
		u[i] += Fy*(3*L-4*x[i])*x[i]**2/(48*E*I) + (x[i]**2*qy*(L-x[i])**2)/(24*E*I)
	else:
		u[i] = u[-i-1]


plt.plot(x,u, 'r--', label='True Solution')
plt.hold(True)
plt.plot(x,y, 'bs', label='Computed Solution')
plt.legend()
plt.title('Plot of y displacement in the steel beam', fontsize=20)
plt.xlabel('Distance from beginning of beam', fontsize=14)
plt.ylabel('Displacement in y axis', fontsize=14)
plt.grid(True)
plt.show()

v = max(abs(u-y))
print v