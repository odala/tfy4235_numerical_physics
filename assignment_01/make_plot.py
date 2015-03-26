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
from matplotlib import animation
import scipy
from scipy.integrate import quad

# ---------------------------------------------------------------------------
# Read the constants from the file "constants.txt"
# ---------------------------------------------------------------------------
def read_constants(realisation):
    filename = 'constants' + realisation + '.txt'
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
def potential(x, t):
    x = x/t.scale_length - np.floor(x/t.scale_length/1.0)*1.0
    if (x >= 0.0 and x < t.alfa):
        U = x/t.alfa
    elif (x >= t.alfa and x < 1.0):
        U = (1.0-x)/(1.0-t.alfa)
    return U


# ---------------------------------------------------------------------------
# Plot the trajectory of the particle(s) and the potential they are in
# ---------------------------------------------------------------------------
def plot_trajectory(trajectories):
    trajectory = plt.figure()

    minimum = min(trajectories[i].min for i in range(len(trajectories)))
    maximum = max(trajectories[i].max for i in range(len(trajectories)))    #-100e-6 # 720e-6
    min_T = min(trajectories[i].T for i in range(len(trajectories)))
    
    # Plotting particles   
    colors = cm.nipy_spectral(np.linspace(0.1, 0.9, len(trajectories)))
    colors = ['k', 'b', 'r']
    for j in range(len(trajectories)):
        t = trajectories[j]
        y = np.linspace(0, (t.nSteps - 1), t.nSteps)
        for i in range(0, t.N):
            plt.plot(1.0e6*t.trajectory[i][:], t.scale_time*y, linewidth = 0.5, color=colors[j], label= (r'$\frac{\Delta U}{k_BT} = $'+str(t.dU/t.kT)) if i == 0 else "")
    
    # Plotting potential
    nPoints = 250
    x = np.linspace(minimum, maximum, nPoints)
    pot = np.zeros(nPoints)
    for i in range(0, len(x)):
        pot[i] = potential(x[i], t)            # *t.scale_potential/t.kT Potential in factors of kT

    stuff = np.linspace(0, y[-1], nPoints)
    X,Y = np.meshgrid(x*1.0e6, stuff)
    Z = np.zeros((nPoints,nPoints))
    for i in range(0, nPoints):
        Z[i] = pot

    plt.contourf(X, Y, Z, 50, alpha=.75, cmap='binary', vmin=abs(Z).min(), vmax=abs(Z).max() + 0.5*(abs(Z).max() - abs(Z).min()))
    cbar = plt.colorbar()
    cbar.set_label(r'Potential / $\Delta U$', fontsize = 15)

    # Layout of plot
    plt.ylabel('Time [s]', fontsize=15)
    plt.xlabel('Position [$\mu$m]', fontsize=15)
    plt.axis([minimum*1.0e6, maximum*1.0e6, 0.0, min_T])
    plt.legend(loc='lower right')
    print 'Plotted the trajectories of the particles with tau = ', t.tau
    
    # Saving figure
    plt.savefig('fig/trajectory' + t.realisation + '.png');


# ---------------------------------------------------------------------------
# Plots the motion of the ensemble of particles and how the particle density
# behaves in time, 3D-plot
# ---------------------------------------------------------------------------
def plot_ensemble_in_time(data, realisation):
    ensemble_motion = plt.figure()
    ax = ensemble_motion.gca(projection='3d')
    
    # initialise start, stop and nbins so that each bin is one ratchet
    start = np.floor(data.min()/scale_length)- (1.0-alfa)
    stop = np.ceil(data.max()/scale_length) + alfa
    nbins = np.ceil((stop-start))
    start = scale_length*start
    stop = scale_length*stop

    nt = 4
    ts = np.linspace(10, nSteps - 1, nt)
    ts = np.floor(ts)
    
    verts = []
    hist_max = 0.0

    D = kT/gama

    vd = 0.0
    for i in range(0, N):
        vd = vd + (data[i, -1] - data[i, 0])/T/N        

    for z in range(0,nt):
        t = ts[z]*scale_time
        hist, x_bins = np.histogram(data[:, ts[z]], bins=nbins, range=(start, stop), normed=False)
        if (max(hist) > hist_max):
            hist_max = max(hist)
        x_center = 1.0e6*(x_bins[:-1] + x_bins[1:]) / 2     # is now in micrometer
        ax.bar(x_center, hist, align='center', zs=t, zdir='y', alpha=0.4, width = (x_center[1] - x_center[0]))
        if (np.floor(dU*6.24150934e18+0.5) == 0):
            x_0 = x_center*1e-6 - vd*t
            if (t==0.0):
                n = np.zeros(len(x_center))
                n[np.where(np.abs(x_center)==np.abs(x_center).argmin())] = N
            else:
                n = N/np.sqrt(4*np.pi*D*t)*np.exp(-np.power(x_0, 2)/(4*D*t))*(x_0[1] - x_0[0])
            verts.append(zip( x_center,  n))

    
    # Convert verts to poly
    cc = lambda c: colorConverter.to_rgba(c, alpha=0.6)
    poly = PolyCollection(verts, facecolors=cc('r'))
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=ts*scale_time, zdir = 'y')
    
    ax.set_xlabel('Position [$\mu$m]')
    ax.set_xlim3d(min(x_center), max(x_center))
    ax.set_ylabel('Time [s]')
    ax.set_ylim3d(min(ts)*scale_time, max(ts)*scale_time)
    ax.set_zlabel('Number of particles')
    ax.set_zlim3d(0.0, hist_max)
    
    # Saving figure
    plt.savefig('fig/ensemble_motion' + realisation + '.png');

    return vd

# ---------------------------------------------------------------------------
# Plots the motion of the ensemble of particles and how the particle density
# behaves in time, 3D-plot, TAKES IN ARRAY OF OBJECTS!!! 
# ---------------------------------------------------------------------------
def plot_objects_in_time(o):
    ensemble_motion = plt.figure()
    ax = ensemble_motion.gca(projection='3d')

    minimum, maximum = 0, 0
    T = o[0].T
    for i in range(0, len(o)):
        minimum = o[i].min if o[i].min < minimum else minimum
        maximum = o[i].max if o[i].max > maximum else maximum
        T       = o[i].T if o[i].T < T else T
    
    # initialise start, stop and nbins so that each bin is one ratchet
    start = np.floor(minimum/o[0].scale_length) - (1.0-o[0].alfa)
    stop  = np.ceil(maximum/o[0].scale_length) + o[0].alfa
    nbins = np.ceil((stop-start))
    start = o[0].scale_length*start
    stop  = o[0].scale_length*stop

    nt = 4
    ts = np.linspace(10, T - o[0].delta_t*o[0].scale_time, nt)
    
    hist_max = 0.0

    D = o[0].kT/o[0].gama

    vd = 0.0   

    cs = ['b', 'g', 'y', 'r']  

    for z in range(0,nt):
        for i in range(0, len(o)):
            tStep = np.ceil(ts[z]/o[0].scale_time)
            hist, x_bins = np.histogram(o[i].trajectory[:, tStep], bins=nbins, range=(start, stop), normed=False)
            hist_max = max(hist) if (max(hist) > hist_max) else hist_max
            x_center = 1.0e6*(x_bins[:-1] + x_bins[1:]) / 2     # is now in micrometer
            ax.bar(x_center, hist, align='center', zs=ts[z], zdir='y', alpha=0.4, width=(x_center[1]-x_center[0]), color=cs[i])
    
    ax.set_xlabel('Position [$\mu$m]')
    ax.set_xlim3d(minimum*1e6, maximum*1e6)
    ax.set_ylabel('Time [s]')
    ax.set_ylim3d(min(ts), max(ts))
    ax.set_zlabel('Number of particles')
    ax.set_zlim3d(0.0, hist_max)
    
    # Saving figure
    plt.savefig('fig/objects_in_time' + '.png');

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
def compare_to_boltzmann(data,  dU, kT):
    # Plot histogram of energy distribution
    minimum = data.min()
    maximum = data.max()
    hist, bins = np.histogram(data, bins=100, range=(minimum, maximum), normed=True)
    U = (bins[:-1] + bins[1:]) / 2
    histogram = plt.figure()           
    plt.plot(6.24150934e18*U, hist,'ro')
    
    # Plot Boltzmann distribution
    boltzmann_dist = lambda U: np.exp(-U/kT) / (kT*(1-np.exp(-dU/kT)))
    U = np.linspace(minimum, maximum, 1000)
    b_dist = boltzmann_dist(U)
    plt.plot(6.24150934e18*U, b_dist)
    #print 'The integral of the Boltzmann distribution from 0.0 to 5.0E-20 = ', scipy.integrate.quad(boltzmann_dist, 0.0, 5.0e-20)
    
    plt.xlabel(r'Energy $U$ [eV]', fontsize=15)
    plt.axis([6.24150934e18*minimum, 6.24150934e18*maximum, 0, 1.1*max(b_dist)])
    print 'Generated a plot of the Boltzmann Distribution for dU / k_BT = ', dU/kT

    # Saving figure
    plt.savefig('fig/boltzmann_distribution_' + str(dU/kT) + '.png');
    

# ---------------------------------------------------------------------------
# Standard Gaussian Random Number Generator
# Checks if the gaussian Box-Muller algorithm really are a normal distribution
# ---------------------------------------------------------------------------
def plot_gaussian():
    normal_dist = lambda x: 1/np.sqrt(2*np.pi)*np.exp(-x**2/2)
    x = np.linspace(-5, 5, 1000)
    rand_gauss = np.loadtxt('check_gaussian.txt')
    hist, bins = np.histogram(rand_gauss, bins=100, normed=True)
    center = (bins[:-1] + bins[1:]) / 2
    gaussian = plt.figure()           
    plt.plot(center, hist,'ro')
    plt.plot(x, normal_dist(x), 'b')
    
    # Saving figure
    plt.savefig('fig/normal_distribution' + '.png');


# ---------------------------------------------------------------------------
# Plot with standard deviation
# returns the maximum and minimum value of y with corresponding x
# --------------------------------------------------------------------------- 
def plot_w_std_dev(x, y, std_dev, xlabel, ylabel, savename, errorbar='errorbar'):
    fig = plt.figure()
    plt.plot(x, y, 'bo')
    if errorbar == 'errorbar':
        plt.errorbar(x, y, yerr=std_dev, fmt=None, color='r')
    elif errorbar == 'filled':
        plt.fill_between(x, y+std_dev, y-std_dev, color='r', alpha='0.5')
    elif errorbar == 'line':
        plt.plot([x, x], [y+std_dev, y-std_dev], 'r')
    else:
        print 'Could not plot std_dev. You have to specify if you want filled, line or errorbar.'
    y_max, y_min = max(y), min(y)
    xy_max, xy_min = x[y.argmax()], x[y.argmin()]
    plt.axis([min(x), max(x), 1.05*(min(y)-max(std_dev)), 1.05*(max(y)+max(std_dev))])
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    
    # Saving figure
    plt.savefig(savename);

    return xy_max, y_max, xy_min, y_min

class Trajectory:
    def __init__(self, realisation):
        self.realisation = realisation
        # Read in data and set some constants
        self.N, self.nSteps, self.delta_t, self.alfa, self.tau, self.zeta, self.r, self.dU, self.L, self.kT = read_constants(realisation)
        self.gama = 6.0*np.pi*self.zeta*self.r
        self.D = self.kT/self.gama
        self.omega = self.dU/(self.gama*self.L**2)
        self.T = self.nSteps*self.delta_t/self.omega
        # Convert from reduced units to real units with these scaling factors
        self.scale_length = self.L                  # reduced units -> meter
        self.scale_time = self.delta_t/self.omega   # time steps    -> seconds
        self.scale_potential = self.dU              # reduced units -> joule
        # Read in trajectory and convert to real length
        self.trajectory = self.scale_length*np.loadtxt('trajectory'+realisation+'.txt')
        self.max, self.min = self.trajectory.max(), self.trajectory.min()


############################################################################        
### ------------------------------- MAIN ------------------------------- ###
############################################################################

# Preferences
do_plot_trajectory = True
do_plot_ensemble_in_time = False
do_plot_drift_velocity = False
do_check_gaussian = False
do_plot_boltzmann = False

# The different realisations
t_poff  = Trajectory('_N1000_rnm12.0_tau.00_dU.0000')       # potential off
t_std   = Trajectory('_N1000_rnm12.0_tau.50_dU80.0000')     # std
t_std2  = Trajectory('_N1000_rnm36.0_tau.50_dU80.0000')     # std med r=36nm
t_foff  = Trajectory('_N1000_rnm12.0_tau.00_dU80.0000')     # flashing off

t_1     = Trajectory('_N1_rnm12.0_tau.00_dU.2600')          # boltzmann-test
t_2     = Trajectory('_N1_rnm12.0_tau.00_dU.0260')          # boltzmann-test
t_3     = Trajectory('_N1_rnm12.0_tau.00_dU.0026')          # boltzmann-test

t_d1    = Trajectory('_N10_rnm12.0_tau.00_dU.2600')
t_d2    = Trajectory('_N10_rnm12.0_tau.00_dU.0260')
t_d3    = Trajectory('_N10_rnm12.0_tau.00_dU.0026')

t_b = t_2

if do_plot_trajectory:
    plot_trajectory([t_d1, t_d2, t_d3])
    
if do_plot_ensemble_in_time:
    #vd = plot_ensemble_in_time(trajectory, realisation)
    plot_objects_in_time([t_1, t_2, t_3])
    #print 'The average drift velocity is: ', vd*1.0e6, 'micrometer/s'

if do_plot_drift_velocity:
    data = np.loadtxt('drift_velocity_N1000_rnm12.0_tau.50_dU80.0000'+'.txt')                  # drift velocity in m/s
    taus = data[:, 0]
    vd = 1.0e6 * data[:, 1]
    std_dev = 1.0e6 * data[:, 2]
    tau_w_max_vd, vd_max, tau_w_min_vd, vd_min = plot_w_std_dev(taus, vd, std_dev, r'Period of flash, $\tau$ [s]', r'Average drift velocity, $\left< v_d \right>$ [$\mathrm{\mu}$m/s]', 'fig/test.png', errorbar = 'filled')
    print 'The tau which results in the maximum drift velocity is tau = ', tau_w_max_vd

if do_plot_boltzmann:
    U = np.zeros(len(t_b.trajectory[:]))
    for i in range(0, len(t_b.trajectory[:])):         # convert to real potential
        U[i] = t_b.dU*potential(t_b.trajectory[i], t_b)
    
    compare_to_boltzmann(U, t_b.dU, t_b.kT)
    
if do_check_gaussian:
    plot_gaussian()