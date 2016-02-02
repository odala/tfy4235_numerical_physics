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
# --- Class that contains all the simulation parameters
#     and the trajectories.
# ---------------------------------------------------------------------------
class Particles:
    def __init__(self, input_file, trajectory_file):
        # --- Open file.
        fin = open(input_file, 'r')
        
        # --- Read parameters from file.
        Nparameters = 13
        parameters = np.zeros(Nparameters)
        i = 0
        for line in fin:
            line = line.partition(',')[0]
            line = line.rstrip()
            row = line.split()
            if i > 0 and i <= Nparameters:
                parameters[i-1] = float(row[2])
            i = i + 1

        # --- Close file.
        fin.close()
        
        # --- Assign parameters.         
        self.Nparticles  = int(parameters[0])
        self.Nsteps      = int(parameters[1])
        self.timestep    = parameters[2]
        self.radius      = parameters[3]
        self.x0          = parameters[4]
        self.tau         = parameters[5]
        self.dU          = parameters[6]*1.60217657E-19
        self.L           = parameters[7]
        self.alpha       = parameters[8]
        self.zeta        = parameters[9]
        self.kT          = parameters[10]*1.60217657E-19
        self.fraction_off= parameters[11]
        self.dSteps      = parameters[12]

        # --- Calculate additional parameters.
        self.gamma = 6.0*np.pi*self.zeta*self.radius
        self.D = self.kT/self.gamma
        self.omega = self.dU/(self.gamma*self.L**2)
        self.total_time = self.Nsteps*self.timestep

        # --- Read in trajectory and convert to real length.
        self.trajectory = self.L*np.loadtxt(trajectory_file) 
        self.max, self.min = self.trajectory.max(), self.trajectory.min()
        self.Nmeasurements = len(self.trajectory[0][:])

    def potential(self, x):
        x = x - np.floor(x/self.L)*self.L
        U = np.zeros(len(x))
        for i in range(len(x)):
            if (x[i] >= 0.0 and x[i] < self.alpha*self.L):
                U[i] = x[i]/(self.alpha*self.L)*self.dU
            elif (x[i] >= self.alpha*self.L and x[i] < self.L):
                U[i] = (self.L-x[i])/(self.L*(1.0-self.alpha))*self.dU
        return U


# -----------------------------------------------------------------------
# Plot the trajectory of the particle(s) and the potential they are in
# -----------------------------------------------------------------------
def plot_trajectory(species):

    # --- Calculate number of species.
    Nspecies = len(species)
    print('#Species: ', Nspecies)

    # --- Calculate limits for time and position.
    minimum_x = min(species[i].min for i in range(Nspecies))
    maximum_x = max(species[i].max for i in range(Nspecies))    #-100e-6 # 720e-6
    minimum_total_time = min(species[i].total_time for i in range(Nspecies))

    # --- Plot trajectory.
    trajectory = plt.figure()
    colors = cm.nipy_spectral(np.linspace(0.1, 0.9, Nspecies))      #colors = ['k', 'b', 'r']
    colors = ['b', 'g', 'r', 'y']
    for j in range(Nspecies):
        s = species[j]
        t = s.timestep*np.linspace(0, (s.Nsteps - 1), s.Nsteps/s.dSteps)
        if s.Nparticles == 1:
            plt.plot(1.0e6*s.trajectory[0][:], t, linewidth = 0.5, color=colors[j], label= (r'$\frac{\Delta U}{k_BT} = $'+str(s.dU/s.kT)))
        else:
            for i in range(0, s.Nparticles):
                plt.plot(1.0e6*s.trajectory[i][:], t, linewidth = 0.5, color=colors[j], label= (r'$\frac{\Delta U}{k_BT} = $' + str(s.dU/s.kT)) if i == 0 else "")
        
    # --- Plot potential.
    Npoints = 250; x = np.linspace(minimum_x, maximum_x, Npoints); pot = species[0].potential(x)/species[0].dU
    stuff = np.linspace(0, t[-1], Npoints)
    X,T = np.meshgrid(x*1e6, stuff)
    U = np.zeros((Npoints,Npoints))
    for i in range(0, Npoints):
        U[i] = pot
    plt.contourf(X, T, U, 50, alpha=.75, cmap='binary', vmin=abs(U).min(), vmax=abs(U).max() + 0.5*(abs(U).max() - abs(U).min()))

    # --- Layout of plot.
    plt.xlabel('Position [$\mu$m]', fontsize=20)
    plt.ylabel('Time [s]', fontsize=20) 
    cbar = plt.colorbar(); cbar.set_label(r'Potential / $\Delta U$', fontsize = 20)
    plt.axis([minimum_x*1.0e6, maximum_x*1.0e6, 0.0, 1.0*minimum_total_time])
    plt.legend(loc='lower right')
    
    # --- Saving figure.
    plt.savefig('fig/trajectory.png')
    print 'Plotted the trajectories of the particles with tau = ', s.tau

# ---------------------------------------------------------------------------
# Plots the motion of the ensemble of particles and how the particle density
# behaves in time, 3D-plot
# ---------------------------------------------------------------------------
def plot_ensembles_in_time(species):

    # --- Calculate limits for time and position.
    minimum, maximum = 0, 0
    total_time = species[0].total_time
    for i in range(0, len(species)):
        minimum     = min(species[i].min, minimum)
        maximum     = max(species[i].max, maximum)
        total_time  = min(species[i].total_time, total_time)
    
    # --- Initialise start, stop and Nbins so that 
    #     each bin is one potential well.
    start = np.floor(minimum/species[0].L) - (1.0-species[0].alpha)
    stop  = np.ceil(maximum/species[0].L) + species[0].alpha
    Nbins = np.ceil(stop-start)
    start = start*species[0].L
    stop  = stop*species[0].L

    # --- Initialise time array.
    Nts = 4
    ts = np.linspace(5, species[0].Nsteps/species[0].dSteps - 1, Nts)
    ts = np.floor(ts)

    # --- Plot histograms (3-Dimensional).
    ensemble_motion = plt.figure()
    ax = ensemble_motion.gca(projection='3d')
    colors = ['b', 'g', 'r', 'y'] 
    hist_max = 0.0
    verts = []
    
    # --- Iterate through the different ensembles.
    for j in range(len(species)):

        # --- If the potential is off:
        #     Calculate the average drift velocthity.
        if (species[j].fraction_off == 1.0):
            avg_vd = 0.0
            for p in range(0, species[j].Nparticles):
                avg_vd = avg_vd + (species[j].trajectory[p][-1] - species[j].trajectory[p][0])/species[j].total_time/species[j].Nparticles  

        # --- Plot histograms of the ensemble in time.
        for i in ts:
            t = i*species[j].timestep*species[j].dSteps
            hist, x_bins = np.histogram(species[j].trajectory[:,i], bins=Nbins, range=(start, stop), normed=False)
            hist_max = max(max(hist), hist_max)
            x_center = (x_bins[:-1] + x_bins[1:]) / 2
            ax.bar(1e6*x_center, hist, align='center', zs=t, zdir='y', alpha=0.4, width = 1e6*(x_center[1] - x_center[0]), color=colors[j])

            # --- If the potential is off:
            #     Plot the theoretical number of particles as a function of t and x.
            if (species[j].fraction_off == 1.0):
                x_0 = x_center - avg_vd*t
                if (t==0.0):
                    n = np.zeros(len(x_center))
                    n[np.where(np.abs(x_center)==np.abs(x_center).argmin())] = species[j].Nparticles
                else:
                    n = species[j].Nparticles/np.sqrt(4*np.pi*species[j].D*t)*np.exp(-np.power(x_0, 2)/(4*species[j].D*t))*(x_0[1] - x_0[0])
                verts.append(zip( 1e6*x_center, n))
     
    # Convert verts to poly
    cc = lambda c: colorConverter.to_rgba(c, alpha=0.6)
    poly = PolyCollection(verts, facecolors=cc('r'), alpha=0.7)
    ax.add_collection3d(poly, zs=ts*species[0].timestep*species[0].dSteps, zdir = 'y')

    # --- Plot layout.
    ax.set_xlabel('Position [$\mu$m]')
    ax.set_xlim3d(minimum*1e6, maximum*1e6)
    ax.set_ylabel('Time [s]')
    ax.set_ylim3d(min(ts)*species[0].timestep*species[0].dSteps, max(ts)*species[0].timestep*species[0].dSteps)
    ax.set_zlabel('Number of particles')
    ax.set_zlim3d(0.0, hist_max/2)

    # Saving figure
    plt.savefig('fig/ensemble_motion.png');

# ---------------------------------------------------------------------------
# Makes a histogram of where in the potential the particle have been and
# compares it to the true Boltzmann distribution from statistical mechanics
# ---------------------------------------------------------------------------
def compare_to_boltzmann(data, dU, kT):

    # --- Divide data into 100 bins.
    minimum = data.min()
    maximum = data.max()
    print(minimum, maximum, len(data))
    hist, bins = np.histogram(data, bins=100, range=(minimum, maximum), normed=True)
    bin_values = (bins[:-1] + bins[1:]) / 2

    # --- Plot histogram.
    histogram = plt.figure()           
    plt.plot(bin_values, hist,'ro')
    
    # --- Plot Boltzmann distribution.
    boltzmann_dist = lambda U: np.exp(-U/kT) / (kT*(1-np.exp(-dU/kT)))
    U = np.linspace(minimum, maximum, 1000); b_dist = boltzmann_dist(U)
    plt.plot(U, b_dist)
    print 'The integral of the Boltzmann distribution from 0.0 to 0.5 = ', scipy.integrate.quad(boltzmann_dist, 0.0, 0.5)
    
    plt.xlabel(r'Energy $U$ [eV]', fontsize=20)
    plt.ylabel(r'$P(U)$', fontsize=20)
    plt.axis([minimum, maximum, 0, 1.1*max(b_dist)])
    print 'Generated a plot of the Boltzmann Distribution for dU / k_BT = ', dU/kT

    # --- Saving figure.
    plt.savefig('fig/boltzmann_distribution.png')

# -----------------------------------------------------------------------
# --- Plot y(x) with standard deviation.
#     Input : xvalues, yvalues, xlabel, ylabel, savename, errorbar
#     Output: x and y values at function maximum and functino minimum
# -----------------------------------------------------------------------
def plot_w_std_dev(x, y, std_dev, xlabel, ylabel, savename, errorbar='errorbar'):
    
    # --- Plot y as a function of x.
    fig = plt.figure()
    plt.plot(x, y, 'bo')

    # --- Plot errorbars.
    if errorbar == 'errorbar':
        plt.errorbar(x, y, yerr=std_dev, fmt=None, color='r')
    elif errorbar == 'filled':
        plt.fill_between(x, y+std_dev, y-std_dev, color='#808080', alpha='0.5')
    elif errorbar == 'line':
        plt.plot([x, x], [y+std_dev, y-std_dev], 'r')
    else:
        print 'Could not plot std_dev. You have to specify if you want filled, line or errorbar.'
    
    # --- Find the maximum and minimum of y with corresponding x-values.
    fmax_y, fmin_y = max(y), min(y)
    fmax_x, fmin_x = x[y.argmax()], x[y.argmin()]
    
    # --- Layout of plot.
    plt.axis([min(x), max(x), 1.05*(min(y)-max(std_dev)), 1.05*(max(y)+max(std_dev))])
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    
    # --- Saving figure.
    plt.savefig(savename)

    return fmax_x, fmax_y, fmin_x, fmin_y
    
# -----------------------------------------------------------------------        
# -------------------------------- MAIN ---------------------------------
# -----------------------------------------------------------------------

# TEST BLOCK.
#test_run = Particles('input.txt', 'output.txt')
#plot_trajectory([test_run])

# TASK 7.
# --- Assign parameters from file and plot the trajectory.
'''
run1 = Particles('input_x1e1.txt', 'output_x1e1.txt')
run2 = Particles('input_x1e0.txt', 'output_x1e0.txt')
run3 = Particles('input_x1e-1.txt', 'output_x1e-1.txt')
plot_trajectory([run1, run2, run3])
'''

# --- Compare with Boltzmann distribution.
'''
b_run = Particles('input_boltzmann.txt', 'output_boltzmann.txt')
compare_to_boltzmann(b_run.potential(b_run.trajectory[:]), b_run.dU, b_run.kT)
'''

# TASK 9.
# --- Drift velocity with respect to flashing period.
v_run = Particles('input_vd.txt', 'output_x1e1.txt')
data = np.loadtxt('drift_velocities.txt')
tau = data[:, 0]
avg_vd = data[:, 1]*v_run.L*v_run.omega
std_dev = data[:, 2]*v_run.L*v_run.omega
fmax_x, fmax_y, fmin_x, fmin_y = plot_w_std_dev(tau, 1e6*avg_vd, 1e6*std_dev, r'Period of flash, $\tau$ [s]', r'Average drift velocity, $\left< v_d \right>$ [$\mathrm{\mu}$m/s]', 'fig/drift_velocities.png', errorbar = 'filled')
print 'The tau which results in the maximum drift velocity is tau = ', fmax_x

# TASK 12.
'''
#ensemble1 = Particles('input_ensemble_check.txt', 'output_ensemble_check.txt')
#ensemble2 = Particles('input_ensemble_check2.txt', 'output_ensemble_check2.txt')
#ensemble3 = Particles('input.txt', 'output.txt')
ensemble_rad12 = Particles('input_ensemble_rad12.txt', 'output_ensemble_rad12.txt')
ensemble_rad36 = Particles('input_ensemble_rad36.txt', 'output_ensemble_rad36.txt')
#run3 = Particles('input_x1e-1.txt', 'output_x1e-1.txt')
#plot_ensembles_in_time([ensemble_rad12, ensemble_rad36])
plot_trajectory([ensemble_rad12, ensemble_rad36])
'''