import numpy as np
import matplotlib
matplotlib.use('Agg');
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm 
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import scipy
from scipy.integrate import quad


# ---------------------------------------------------------------------------
# Read the constants from the file "constants.txt"
# ---------------------------------------------------------------------------
def read_constants(filename):
    fin = open(filename, 'r')
    N, nSteps = [int(x) for x in fin.readline().split()]     # read first line
    delta_t = [float(x) for x in fin.readline().split()]     # read second line
    delta_t = delta_t[0]
    alfa = [float(x) for x in fin.readline().split()]    
    alfa = alfa[0]
    tau = [float(x) for x in fin.readline().split()]
    tau = tau[0]
    zeta = [float(x) for x in fin.readline().split()]
    zeta = zeta[0]
    r = [float(x) for x in fin.readline().split()]
    r = r[0]
    dU = [float(x) for x in fin.readline().split()]    
    dU = dU[0]
    L = [float(x) for x in fin.readline().split()]    
    L = L[0]
    kT = [float(x) for x in fin.readline().split()]        
    kT = kT[0]
    fin.close()
    return (N, nSteps, delta_t, alfa, tau, zeta, r, dU, L, kT)


# ---------------------------------------------------------------------------
# The potential
# ---------------------------------------------------------------------------
def potential(x, alfa):
    x = x/scale_length - np.floor(x/scale_length/1.0)*1.0
    if (x >= 0.0 and x < alfa):
        U = x/alfa
    elif (x >= alfa and x < 1.0):
        U = (1.0-x)/(1.0-alfa)
    return U


# ---------------------------------------------------------------------------
# Plot the trajectory of the particle(s) and the potential they are in
# ---------------------------------------------------------------------------
def plot_trajectory(data, alfa, tau, realisation):
    dim1, dim2 = data.shape
    minimum = -30e-6#data.min()
    maximum = 720e-6#data.max()
    
    trajectory = plt.figure()
    define_color_cycle(trajectory)
    
    # Plotting particles    
    y = np.linspace(0, (dim2 - 1), dim2)
    for i in range(0, dim1):
        plt.plot(1.0e6*data[i][:], scale_time*y, linewidth = 0.3)
    plt.ylabel(r'Time [s]', fontsize=15)
    plt.xlabel(r'Position [\mu m]', fontsize=15)
    plt.axis([minimum*1.0e6, maximum*1.0e6, 0.0, y[-1]*scale_time])
    plt.title(r'Trajectories of the particles with $\tau$ = %f'%(tau), fontsize = 20)
    
    # Plotting potential
    nPoints = 250
    x = np.linspace(minimum, maximum, nPoints)
    pot = np.zeros(nPoints)
    for i in range(0, len(x)):
        pot[i] = scale_potential/kT*potential(x[i], alfa)            # Potential in factors of kT

    stuff = np.linspace(0, y[-1], nPoints)
    X,Y = np.meshgrid(x*1.0e6, stuff)
    Z = np.zeros((nPoints,nPoints))
    for i in range(0, nPoints):
        Z[i] = pot

    plt.contourf(X, Y, Z, 50, alpha=.75, cmap='binary', vmin=abs(Z).min(), vmax=abs(Z).max() + 0.5*(abs(Z).max() - abs(Z).min()))
    cbar = plt.colorbar()
    cbar.set_label(r'Potential in factors of $k_B T$', fontsize = 15)
    
    # Saving figure
    plt.savefig('fig/trajectory' + realisation + '.png');
    

# ---------------------------------------------------------------------------
# Plots the motion of the ensemble of particles and how the particle density
# behaves in time, 3D-plot
# ---------------------------------------------------------------------------
def plot_ensemble_in_time(data, realisation):
    dim1, dim2 = data.shape

    ensemble_motion = plt.figure()
    ax = ensemble_motion.gca(projection='3d')
    
    start = np.floor(data.min()/scale_length)- (1.0-alfa)
    stop = np.ceil(data.max()/scale_length) + alfa
    print '[start, stop]', start, stop
    nbins = np.ceil((stop-start))
    start = scale_length*start
    stop = scale_length*stop
    print 'nbins', nbins

    n = 4
    
    t = np.linspace(10, dim2 - 1, n)
    t = np.floor(t)
    
    verts = []
    verts2 = []
    zs = np.linspace(0.0, 1.0*(n-1), n)
    hist_max = 0.0
    
    D = kT/gama

    vd = 0.0
    for i in range(0, N):
        vd = vd + (data[i, -1] - data[i, 0])/T/N        

    for z in zs:
        hist, x_bins = np.histogram(data[:, t[z]], bins=nbins, range=(start, stop), normed=False)
        if (max(hist) > hist_max):
            hist_max = max(hist)
        x_center = (x_bins[:-1] + x_bins[1:]) / 2
        print 'TIME: ', t[z]*scale_time
        x_teo = x_center - vd*t[z]*scale_time
        hist_teo = N/np.sqrt(4*np.pi*D*t[z]*scale_time)*np.exp(-np.power(x_teo, 2)/(4*D*t[z]*scale_time))
        print 'MAX: ', max(hist_teo)
        x_center = 1.0e6*x_center     # is now in micrometer
        verts.append(zip(x_center, hist))
        verts2.append(zip( x_center,  hist_teo))
        ax.bar(x_center, hist, align='center', zs=t[z]*scale_time, zdir='y', alpha=0.4, width = (x_center[1] - x_center[0]))
        
    
    # Convert verts to poly
    cc = lambda c: colorConverter.to_rgba(c, alpha=0.6)
    #poly = PolyCollection(verts)
    #poly.set_alpha(0.7)
    #ax.add_collection3d(poly, zs=t, zdir = 'y')
    poly2 = PolyCollection(verts2)
    poly2.set_alpha(0.7)
    ax.add_collection3d(poly2, zs=t*scale_time, zdir = 'y')
    
    ax.set_xlabel('Position [$\mu$m/s]')
    ax.set_xlim3d(min(x_center), max(x_center))
    ax.set_ylabel('Time [s]')
    ax.set_ylim3d(min(t)*scale_time, max(t)*scale_time)
    ax.set_zlabel('Number of particles')
    ax.set_zlim3d(0.0, hist_max)
    
    # Saving figure
    plt.savefig('fig/ensemble_motion' + realisation + '.png');
    
    return vd


# ---------------------------------------------------------------------------
# Defines a costumised color cycle for the plotting
# ---------------------------------------------------------------------------
def define_color_cycle(fig_name):
    # Define color cycle for the figure fig_name
    NUM_COLORS = 10
    cm = plt.get_cmap('nipy_spectral')    #'terrain'
    ax = fig_name.add_subplot(111)
    ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])


# ---------------------------------------------------------------------------
# Makes a histogram of where in the potential the particle have been and
# compares it to the true Boltzmann distribution from statistical mechanics
# ---------------------------------------------------------------------------
def make_boltzmann_histogram(data,  dU, kT):
    # Plot histogram of Boltzmann distribution
    minimum = data.min()
    maximum = data.max()
    hist, bins = np.histogram(data, bins=100, range=(minimum, maximum), normed=True)
    center = (bins[:-1] + bins[1:]) / 2
    histogram = plt.figure()        #plt.bar(center, hist, align='center', width=bins[1] - bins[0])            
    plt.plot(center, hist,'ro')
    plt.title(r'Boltzmann Distribution for $dU / k_BT=$ %d'%(dU/kT), fontsize = 20)
    
    # Plot exact Boltzmann distribution
    boltzmann_dist = lambda U: np.exp(-U/kT) / (kT*(1-np.exp(-dU/kT)))
    U = np.linspace(minimum, maximum, 1000)
    plt.plot( U, boltzmann_dist(U) )
    #print 'The integral of the Boltzmann distribution from 0.0 to 5.0E-20 = ', scipy.integrate.quad(boltzmann_dist, 0.0, 5.0e-20)
    
    # Saving figure
    plt.savefig('fig/boltzmann_distribution' + '.png');
    

# ---------------------------------------------------------------------------
# Checks if the gaussian Box-Muller algorithm really are a normal distribution
# ---------------------------------------------------------------------------
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
    plt.savefig('fig/normal_distribution' + '.png');
    

# ---------------------------------------------------------------------------
# Plots the drift velocity versus flashing period   
# --------------------------------------------------------------------------- 
def plot_drift_velocity(filename):
    temp_arr = np.loadtxt(filename)     # drift velocity in m/s
    temp_arr = temp_arr[temp_arr[:,0].argsort()]
    taus = temp_arr[:, 0]
    vd = 1.0e6 * temp_arr[:, 1]
    std_dev = 1.0e6 * temp_arr[:, 2]
    tau_w_max_vel = taus[vd.argmax()]
    drift_vel = plt.figure()
    plt.plot(taus, vd, 'b*')
    plt.plot(taus, vd+std_dev, 'r*')
    plt.plot(taus, vd-std_dev, 'r*')
    plt.errorbar(taus, vd, yerr=std_dev, fmt=None, color='b')
    plt.axis([0.0, max(taus), 0.0, 1.05*(max(vd)+max(std_dev))])
    plt.ylabel(r'Average drift velocity, $\left< v_d \right>$ [$\mathrm{\mu}$m/s]', fontsize=15)
    plt.xlabel(r'Period of flash, $\tau$ [s]', fontsize=15)
    
    # Saving figure
    plt.savefig('fig/' + filename + '.png');
    
    return tau_w_max_vel
    
    
############################################################################        
### ------------------------------- MAIN ------------------------------- ###
############################################################################

do_plot_trajectory = True
do_plot_ensemble_in_time = True
do_plot_drift_velocity = True
do_check_gaussian = False
do_plot_boltzmann = False

# Read in data [dt IS IN REDUCED UNITS!]

std_tau = 0.1
std_N   = 100

realisation = '_N1000_rnm12.0_tau.50'
file_const = 'constants' + realisation + '.txt'
file_trajectory = 'trajectory' + realisation + '.txt'
file_vd = 'drift_velocity' + realisation + '.txt'

N, nSteps, delta_t, alfa, tau, zeta, r, dU, L, kT = read_constants(file_const)
gama = 6.0*np.pi*zeta*r
omega = dU/(gama*L**2)
T = nSteps*delta_t/omega

# Convert from reduced units to real units with these scaling factors
scale_length = L               # micrometer of real length
scale_time = delta_t/omega     # time steps
scale_potential = dU           # potential

data = scale_length*np.loadtxt(file_trajectory)

if do_plot_trajectory:
    plot_trajectory(data, alfa, tau, realisation)
    
if do_plot_ensemble_in_time:
    vd = plot_ensemble_in_time(data, realisation)
    print 'The average drift velocity is: ', vd*1.0e6, 'micrometer/s'

if do_plot_drift_velocity:
    tau_w_max_vel = plot_drift_velocity(file_vd)
    print 'The tau which results in the maximum drift velocity is tau = ', tau_w_max_vel
    
if do_plot_boltzmann:
    U = np.zeros(len(data[0][:]))
    for i in range(0, len(data[0][:])):         # convert to real potential
        U[i] = dU*potential(data[0][i], alfa)
    
    make_boltzmann_histogram(U, dU, kT)
    
if do_check_gaussian:
    plot_gaussian()
