import numpy as np
import matplotlib
matplotlib.use('Agg');
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
import scipy
from scipy.integrate import quad

# Read some constants from the file "constants.txt"
def read_constants():
	fin = open("constants" + '.txt', 'r')
	N, nTimeSteps = [int(x) for x in fin.readline().split()] 	# read first line
	dt = [float(x) for x in fin.readline().split()]			# read second line
	dt = dt[0]
	alfa = [float(x) for x in fin.readline().split()]		# read third line
	alfa = alfa[0]
	tau = [float(x) for x in fin.readline().split()]		# read third line
	tau = tau[0]
	omega = [float(x) for x in fin.readline().split()]		# read third line
	omega = omega[0]
	dU = [float(x) for x in fin.readline().split()]			# read third line
	dU = dU[0]
	L = [float(x) for x in fin.readline().split()]			# read third line
	L = L[0]
	kT = [float(x) for x in fin.readline().split()]			# read third line
	kT = kT[0]
	fin.close()
	return (N, nTimeSteps, dt, alfa, tau, omega, dU, L, kT)


# The potential
def potential(x, alfa):

	x = x - np.floor(x/1.0)*1.0
	if (x >= 0.0 and x < alfa):
		U = x/alfa
	elif (x >= alfa and x < 1.0):
		U = (1.0-x)/(1.0-alfa)
	return U

# Plot the trajectory of the particle(s) and the potential they are in
def plot_trajectory(alfa, tau):
	data = np.loadtxt('trajectory.txt')
	dim1, dim2 = data.shape
	minimum = data.min()
	maximum = data.max()
	
	trajectory = plt.figure()
	define_color_cycle(trajectory)
	
	# Plot of the particles	
	y = np.linspace(0, (dim2 - 1), dim2)
	for i in range(0, dim1):
		plt.plot(scale_length*data[i][:], scale_time*y, linewidth = 0.3)
	plt.ylabel(r'Time [s]', fontsize=15)
	plt.xlabel(r'Position [\mu m]', fontsize=15)
	plt.axis([minimum*scale_length, maximum*scale_length, 0.0, y[-1]*scale_time])
	plt.title(r'Trajectories of the particles with $\tau$ = %f'%(tau), fontsize = 20)
	
	# Plot of the potential
	nPoints = 250
	x = np.linspace(minimum, maximum, nPoints)
	pot = np.zeros(nPoints)
	for i in range(0, len(x)):
		pot[i] = scale_potential/kT*potential(x[i], alfa)			# Potential in factors of kT
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

# Defines a cosumised color cycle for the plotting
def define_color_cycle(fig_name):
	# Define color cycle for the figure fig_name
	NUM_COLORS = 10
	cm = plt.get_cmap('nipy_spectral')	#'terrain'
	ax = fig_name.add_subplot(111)
	ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])


# Makes a histogram of where in the potential the particles have been and
# compares it to the true Boltzmann distribution from statistical mechanics
def make_boltzmann_histogram(data,  dU, kT):
	# Plot histogram of Boltzmann distribution
	minimum = data.min()
	maximum = data.max()
	hist, bins = np.histogram(data, bins=100, range=(minimum, maximum), normed=True)
	center = (bins[:-1] + bins[1:]) / 2
	histogram = plt.figure()		#plt.bar(center, hist, align='center', width=bins[1] - bins[0])			
	plt.plot(center, hist,'ro')
	plt.title(r'Boltzmann Distribution for $dU / k_BT=$ %d'%(dU/kT), fontsize = 20)
	
	# Plot exact Boltzmann distribution
	boltzmann_dist = lambda U: np.exp(-U/kT) / (kT*(1-np.exp(-dU/kT)))
	U = np.linspace(minimum, maximum, 1000)
	plt.plot( U, boltzmann_dist(U) )
	print 'The integral of the Boltzmann distribution from 0.0 to 5.0E-20 = ', scipy.integrate.quad(boltzmann_dist, 0.0, 5.0e-20)
	
	# Saving figure
	plt.savefig('boltzmann_distribution' + '.png');
	

# Checks if the gaussian Box-Muller algorithm really are a normal distribution
def plot_gaussian():
	normal_dist = lambda x: 1/np.sqrt(2*np.pi)*np.exp(-x**2/2)
	x = np.linspace(-5, 5, 1000)
	rand_gauss = np.loadtxt('check_gaussian.txt')
	hist, bins = np.histogram(rand_gauss, bins=100, normed=True)
	center = (bins[:-1] + bins[1:]) / 2
	gaussian = plt.figure()
	#plt.bar(center, hist, align='center', width=bins[1] - bins[0])			
	plt.plot(center, hist,'ro')
	plt.plot(x, normal_dist(x), 'b')
	plt.title(r'Standard Gaussian Random Number Generator', fontsize = 20)
	
	# Saving figure
	plt.savefig('normal_distribution' + '.png');
	

# Plots the drift velocity versus flashing period	
def plot_drift_velocity():
	temp_arr = np.loadtxt('drift_velocity.txt')
	taus = temp_arr[:, 0]
	drift_velocity = temp_arr[:, 1]
	drift_vel = plt.figure()
	plt.plot(taus, drift_velocity, 'x')
	
	# Saving figure
	plt.savefig('drift_velocity' + '.png');
	
	
	
############################################################################		
### ------------------------------- MAIN ------------------------------- ###
############################################################################

do_check_gaussian = False
do_plot_trajectory = True
do_plot_boltzmann = False
do_plot_drift_velocity = True

# Read in data [ALL IS IN REDUCED UNITS]
N, nTimeSteps, dt, alfa, tau, omega, dU, L, kT = read_constants()

if do_plot_trajectory:
	# Convert from reduced units to real units with these scaling factors
	scale_length = L/1.0e-6			# micrometer of real length
	scale_time = 10*dt/omega		# time steps
	scale_potential = dU			# potential

	plot_trajectory(alfa, tau)

if do_plot_drift_velocity:
	plot_drift_velocity()
	
if do_plot_boltzmann:
	U = np.zeros(len(data[0][:]))
	for i in range(0, len(data[0][:])):				# convert to real potential
		U[i] = dU*potential(data[0][i], alfa)
	
	make_boltzmann_histogram(U, dU, kT)
	
if do_check_gaussian:
	plot_gaussian()
