import numpy as np
from matplotlib import pyplot as plt

err_w_dt = np.loadtxt('crank_err_dt.txt')
dt = err_w_dt[:,1]
err_dt = err_w_dt[:,2]
deg1_fit = np.polyfit(dt, err_dt, 1)

plt.figure()
plt.plot(dt, err_dt, 'ro')
plt.plot(dt, deg1_fit[0]*dt + deg1_fit[1])
plt.show()

err_w_dx = np.loadtxt('crank_err_dx.txt')
dx = err_w_dx[:,0]
err_dx = err_w_dx[:,2]
deg2_fit = np.polyfit(dx, err_dx, 2)

plt.figure()
plt.plot(dx, err_dx, 'ro')
plt.plot(dx, deg2_fit[0]*dx**2 + deg2_fit[1]*dx + deg2_fit[2] )
plt.show()