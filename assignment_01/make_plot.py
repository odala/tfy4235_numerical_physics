# Reads and plots data from:
file_name = 'coordinates';

import numpy as np
import matplotlib
matplotlib.use('Agg');
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
import scipy
from scipy.integrate import quad


def read_data():
	fin = open(file_name + '.txt', 'r')
	dim1, dim2 = [int(x) for x in fin.readline().split()] 	# read first line
	dt = [float(x) for x in fin.readline().split()]			# read second line
	dt = dt[0]
	alfa = [float(x) for x in fin.readline().split()]		# read third line
	alfa = alfa[0]
	tau = [float(x) for x in fin.readline().split()]		# read third line
	tau = tau[0]
	dU = [float(x) for x in fin.readline().split()]			# read third line
	dU = dU[0]
	L = [float(x) for x in fin.readline().split()]			# read third line
	L = L[0]
	kT = [float(x) for x in fin.readline().split()]			# read third line
	kT = kT[0]
	data = [float(number) for line in fin for number in line.split()]
	minimum = min(data)
	maximum = max(data)
	data = np.reshape(data, (dim1, dim2))
	fin.close()
	return (data, dim1, dim2, dt, alfa, tau, dU, L, kT, minimum, maximum)
	
def potential(x, alfa):

	x = x - np.floor(x/1.0)*1.0
	if (x >= 0.0 and x < alfa):
		U = x/alfa
	else:
		U = (1.0-x)/(1.0-alfa)
	return U

def plot_trajectory(data, dim1, dim2, alfa, tau, minimum, maximum):
	trajectory = plt.figure()
	define_color_cycle(trajectory)
	
	# Plot of the particles	
	y = np.linspace(0, (dim2 - 1), scale_time*dim2)
	for i in range(0, dim1):
		plt.plot(scale_length*data[i][:], y, linewidth = 0.3)
	plt.ylabel(r'Timestep}', fontsize=15)
	plt.xlabel(r'Position [\mu m]}', fontsize=15)
	plt.axis([minimum*scale_length, maximum*scale_length, 0.0, y[-1]])
	plt.title(r'Trajectories of the particles with $\tau$ = %f'%(tau), fontsize = 20)
	
	# Plot of the potential
	nPoints = 250
	x = np.linspace(minimum, maximum, nPoints)
	pot = np.zeros(nPoints)
	for i in range(0, len(x)):
		pot[i] = scale_potential*potential(x[i], alfa)
	#plt.plot(x*scale_length, 0.2*y[-1]/max(pot) * pot, 'k')		# The potential as a graph

	stuff = np.linspace(0, y[-1], nPoints)
	X,Y = np.meshgrid(x*scale_length, stuff)
	Z = np.zeros((nPoints,nPoints))
	for i in range(0, nPoints):
		Z[i] = pot

	plt.contourf(X, Y, Z, 50, alpha=.75, cmap='binary', vmin=abs(Z).min(), vmax=abs(Z).max() + 0.5*(abs(Z).max() - abs(Z).min()))
	cbar = plt.colorbar()
	cbar.set_label(r'Potential in factors of $k_B T$', fontsize = 15)
	
	# Saving figure
	plt.savefig('trajectory' + '.png');

def define_color_cycle(fig_name):
	# Define color cycle for the figure fig_name
	NUM_COLORS = 10
	cm = plt.get_cmap('nipy_spectral')	#'terrain'
	ax = fig_name.add_subplot(111)
	ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])


def make_boltzmann_histogram(data, minimum, maximum,  dU, kT):
	hist, bins = np.histogram(data, bins=1000, range=(minimum, maximum), normed=True)
	
	# Plot histogram
	width = 0.7 * (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	histogram = plt.figure()
	plt.title(r'Boltzmann Distribution for $dU / k_BT=$ %d'%(dU/kT), fontsize = 20)
	plt.bar(center, hist, align='center', width=width)			#plt.plot(center, hist,'o')
	
	# Plot Boltzmann distribution
	boltzmann_dist = lambda U: np.exp(-U/kT) / (kT*(1-np.exp(-dU/kT)))
	U = np.linspace(minimum, maximum, 1000)
	plt.plot( U, boltzmann_dist(U) )
	print 'The integral of the Boltzmann distribution from 0.0 to 5.0E-20 = ', scipy.integrate.quad(boltzmann_dist, 0.0, 5.0e-20)
	
	# Saving figure
	plt.savefig('boltzmann_distribution' + '.png');
	
	
	
	
############################################################################		
### ------------------------------- MAIN ------------------------------- ###

make_boltzmann = False

# Read in data [ALL IS IN REDUCED UNITS]
data, dim1, dim2, dt, alfa, tau, dU, L, kT, minimum, maximum = read_data()

# Convert from reduced units to real units with these scaling factors
scale_length = L/1.0e-6		# micrometer of real length
scale_time = 1.0			# time steps
scale_potential = dU/kT		# potential in factors of kT

plot_trajectory(data, dim1, dim2, alfa, tau, minimum, maximum)

if make_boltzmann:
	U = np.zeros(len(data[0][:]))
	for i in range(0, len(data[0][:])):				# convert to real potential
		U[i] = dU*potential(data[0][i], alfa)
	
	make_boltzmann_histogram(U, min(U), max(U), dU, kT)
