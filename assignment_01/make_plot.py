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
    plt.xlabel('Position [$\mu$m]', fontsize=labelsize)
    plt.ylabel('Time [s]', fontsize=labelsize) 
    plt.tick_params(axis='both', which='major', labelsize=ticksize)
    cbar = plt.colorbar(); cbar.set_label(r'Potential / $\Delta U$', fontsize = labelsize)
    cbar.ax.tick_params(labelsize=ticksize) 
    plt.axis([minimum_x*1.0e6, maximum_x*1.0e6, 0.0, 1.0*minimum_total_time])
    
    # --- Saving figure.
    plt.tight_layout()
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

    # --- Divide data into Nbins bins.
    Nbins = 100
    minimum = data.min(); maximum = data.max()
    print(minimum, maximum, len(data))
    hist, bins = np.histogram(data, bins=Nbins, range=(minimum, maximum), normed=True)
    bin_values = (bins[:-1] + bins[1:]) / 2

    # --- Plot histogram.
    histogram = plt.figure()
    plt.scatter(bin_values/dU, hist*dU, s=50, facecolors='none', edgecolors='r')
    #plt.plot(bin_values/dU, hist*dU,'ro')
    
    # --- Plot Boltzmann distribution.
    boltzmann_dist = lambda U: np.exp(-U/kT) / (kT*(1-np.exp(-dU/kT)))
    U = np.linspace(minimum, maximum, Nbins); b_dist = boltzmann_dist(U)
    plt.plot(U/dU, b_dist*dU)
    print 'The integral of the Boltzmann distribution from 0.0 to 0.5 = ', scipy.integrate.quad(boltzmann_dist, minimum, maximum)
    
    plt.xlabel(r'Relative energy $\frac{U}{\Delta U}$', fontsize=labelsize)
    plt.ylabel(r'$P(U)$', fontsize=labelsize)
    plt.tick_params(axis='both', which='major', labelsize=ticksize)
    plt.axis([minimum/dU, maximum/dU, 0, 1.1*max(b_dist*dU)])
    print 'Generated a plot of the Boltzmann Distribution for dU / k_BT = ', dU/kT

    # --- Saving figure.
    plt.tight_layout()
    plt.savefig('fig/boltzmann_distribution.png')

# -----------------------------------------------------------------------
# --- Plot potential
#     Input : ensemble
#     Output: 
# -----------------------------------------------------------------------
def plot_potential(ensemble):

    # --- Plot potential as a function of x.
    fig = plt.figure()
    x = np.array([-ensemble.L, -(1-ensemble.alpha)*ensemble.L, 0, ensemble.alpha*ensemble.L, ensemble.L])
    x_temp = np.array([-ensemble.L, -0.9*(1-ensemble.alpha)*ensemble.L, 0, ensemble.alpha*ensemble.L, ensemble.L])
    U = ensemble.potential(x)
    my_xticks = ['$-L$',r'$-L+\alpha L$','$0$',r'$\alpha L$','$L$']
    my_yticks = ['$0$', r'$\Delta U$', '$0$', r'$\Delta U$', '$0$']
    plt.xticks(x_temp, my_xticks, fontsize = ticksize)
    plt.yticks(U, my_yticks, fontsize = ticksize)
    plt.plot(x, U, 'k-', linewidth = 1.5)

    # --- 
    plt.plot([x[1], x[1]], [0, max(U)], 'k--', linewidth = 1.5)
    plt.plot([x[3], x[3]], [0, max(U)], 'k--', linewidth = 1.5)

    # --- Layout of plot.
    plt.axis([min(x), max(x), 0, 1.2*max(U)])
    plt.xlabel('Position $x$', fontsize=labelsize)
    plt.ylabel('Potential $U(x)$', fontsize=labelsize) 

    # --- Saving figure.
    plt.tight_layout()
    plt.savefig('fig/potential.png')

# -----------------------------------------------------------------------
# --- Plot y(x) with standard deviation.
#     Input : xvalues, yvalues, xlabel, ylabel, savename, errorbar
#     Output: x and y values at function maximum and functino minimum
# -----------------------------------------------------------------------
def plot_w_std_dev(x, y, std_dev, xlabel, ylabel, savename, colorline, errorbar='errorbar'):
    
    Nplots = len(x)

    # --- Allocate arrays for minima and maxima.
    fmax_y, fmin_y = np.zeros(Nplots), np.zeros(Nplots)
    fmax_x, fmin_x = np.zeros(Nplots), np.zeros(Nplots)
    fmax_std, fmin_std = np.zeros(Nplots), np.zeros(Nplots)

    fig = plt.figure()
    for i in range(Nplots):
        # --- Plot y as a function of x.
        plt.plot(x[i], y[i], colorline[i])

        # --- Plot errorbars.
        if (min(std_dev[i]) == 0 and max(std_dev[i]) == 0):
            print('No std_dev.')
        elif errorbar == 'errorbar':
            plt.errorbar(x[i], y[i], yerr=std_dev[i], fmt=None, color='r')
        elif errorbar == 'filled':
            plt.fill_between(x[i], y[i]+std_dev[i], y[i]-std_dev[i], color='#808080', alpha='0.5')
        elif errorbar == 'line':
            plt.plot([x[i], x[i]], [y[i]+std_dev[i], y[i]-std_dev[i]], 'r')
        else:
            print 'Could not plot std_dev. You have to specify if you want filled, line or errorbar.'
    
        # --- Find the maximum and minimum of y with corresponding x-values.
        fmax_y[i], fmin_y[i] = max(y[i]), min(y[i])
        fmax_x[i], fmin_x[i] = x[i][y[i].argmax()], x[i][y[i].argmin()]
        fmax_std[i], fmin_std[i] = std_dev[i][y[i].argmax()], std_dev[i][y[i].argmin()]

    # --- Layout of plot.
    plt.axis([0.0, 2.0, min(fmin_y), 1.05*max(fmax_y)])
    plt.xlabel(xlabel, fontsize=labelsize)
    plt.ylabel(ylabel, fontsize=labelsize)
    
    # --- Saving figure.
    plt.tight_layout()
    plt.savefig(savename)

    return fmax_x, fmax_y, fmax_std, fmin_x, fmin_y, fmin_std

    
# -----------------------------------------------------------------------        
# -------------------------------- MAIN ---------------------------------
# -----------------------------------------------------------------------

labelsize = 25
ticksize = 15

# TEST BLOCK.
#test_run = Particles('input.txt', 'output.txt')
#plot_trajectory([test_run])

# TASK 7.
# --- Assign parameters from file and plot the trajectory.
#     Compare with Boltzmann distribution.
'''
run1 = Particles('input_dU_kTe+1.txt', 'output_dU_kTe+1.txt')
run2 = Particles('input_dU_kTe-1.txt', 'output_dU_kTe-1.txt')
plot_trajectory([run1])
plot_trajectory([run2])
compare_to_boltzmann(run1.potential(run1.trajectory.flatten()), run1.dU, run1.kT)
compare_to_boltzmann(run2.potential(run2.trajectory.flatten()), run2.dU, run2.kT)
'''

# --- Plot potential.
'''
test = Particles('input_x1e1.txt', 'output_x1e1.txt')
plot_potential(test)
'''

# --- Plot trajectories of different tau values.
'''
high    = Particles('input_tau_0point01.txt', 'output_tau_0point01.txt')
medium  = Particles('input_tau_0point5.txt', 'output_tau_0point5.txt')
low     = Particles('input_tau_40point0.txt', 'output_tau_40point0.txt')
print('Plot for taus: ', high.tau, medium.tau, low.tau)
plot_trajectory([low, medium, high])
'''

# TASK 9.
# --- Drift velocity with respect to flashing period.
'''
data12nm = np.loadtxt('drift_velocities.txt')
tau12nm = data12nm[:, 0]
avg_vd12nm = data12nm[:, 1]
std_dev12nm = data12nm[:, 2]
data36nm = np.loadtxt('drift_velocities.txt')
tau36nm = data36nm[:, 0]
avg_vd36nm = data36nm[:, 1]*0.5
std_dev36nm = data36nm[:, 2]
fmax_x, fmax_y, fmax_std, fmin_x, fmin_y, fmin_std = plot_w_std_dev([tau12nm, tau12nm*3, tau36nm], [1e6*avg_vd12nm, 1e6*avg_vd12nm/3, 1e6*avg_vd36nm], [1e6*std_dev12nm, np.zeros(len(tau12nm)), 1e6*std_dev36nm], r'Period of flash, $\tau$ [s]', r'Average drift velocity, $\left< v_d \right>$ [$\mathrm{\mu}$m/s]', 'fig/drift_velocities.png', ['b.', 'r', 'r.'], errorbar = 'filled')
print 'The tau which results in the maximum drift velocity is tau = ', fmax_x, ' with drift velocity ', fmax_y, ' and standard deviation ', fmax_std
'''

# TASK 12.
# --- Motion of two ensembles in time.

#ensemble1 = Particles('input_ensemble_check.txt', 'output_ensemble_check.txt')
#ensemble2 = Particles('input_ensemble_check2.txt', 'output_ensemble_check2.txt')
#ensemble3 = Particles('input.txt', 'output.txt')
#ensemble_rad12 = Particles('input_ensemble_rad12.txt', 'output_ensemble_rad12.txt')
#ensemble_rad36 = Particles('input_ensemble_rad36.txt', 'output_ensemble_rad36.txt')
#run3 = Particles('input_x1e-1.txt', 'output_x1e-1.txt')
#plot_ensembles_in_time([ensemble_rad12, ensemble_rad36])
#plot_trajectory([ensemble_rad12, ensemble_rad36])