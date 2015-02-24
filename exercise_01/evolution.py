import numpy as np
import matplotlib.pyplot as plt

nSpecies = 64

chain = np.random.rand(nSpecies)

mAct = []

T = 4000
t = 0

while(t < T):
	minBarrier = min(chain)
	for i in range(0, nSpecies):
		if chain[i] == minBarrier:
			mAct[len(mAct):] = [i]
			
			chain[i] = np.random.rand()
			
			if i == 0:
				chain[-1] = np.random.rand()
			else:
				chain[i-1] = np.random.rand()
			
			if i == nSpecies - 1:
				chain[0] = np.random.rand()  
			else:
				chain[i+1] = np.random.rand()
			
			break
	t += 1

# Plot of where in the chain the mutations happens vs time
time = range(0, T)
plt.plot(mAct, time, 'k.')
plt.xlim(0, 64)
plt.ylim(0, T)
plt.xlabel('Mutation activity')
plt.ylabel('Time')
plt.show()
