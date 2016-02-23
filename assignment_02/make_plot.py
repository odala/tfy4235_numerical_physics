import numpy as np
import scipy
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg');
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
from mpl_toolkits.mplot3d import axes3d
#from mpl_toolkits.mplot3d.axes3d import Axes3D
#from matplotlib import cm 
#from matplotlib.collections import PolyCollection
#from matplotlib.colors import colorConverter
#from matplotlib import animation

# --- Define global constants.
e   = 1.602176565E-19           # elementary charge (C)
m_e = 9.10938291E-31            # electron mass (kg)
m_p = 1.672621777E-27           # proton mass (kg)

# ---------------------------------------------------------------------------
# --- Class that contains all the simulation parameters
#     and the trajectories.
# ---------------------------------------------------------------------------
class Particle:
    def __init__(self, input_file):
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
        self.dSteps      = int(parameters[2])
        self.timestep    = parameters[3]
        self.sgn         = parameters[4]
        self.x0          = parameters[5]
        self.y0          = parameters[6]
        self.z0          = parameters[7]
        self.u0          = parameters[8]
        self.v0          = parameters[9]
        self.w0          = parameters[10]
        self.B           = parameters[11]
        self.E           = parameters[12]

        self.mass = m_p if (self.sgn > 0) else m_e
        
        # --- Calculate additional parameters.
        self.omega = e * self.B / self.mass
        self.r = np.sqrt(self.u0**2 + self.v0**2) / self.omega
        self.total_time = self.Nsteps*self.timestep

def plot_xy(data, savename):

    fig = plt.figure()

    # --- Plot points in 2D.
    colors = ['r', 'b']
    linestyles = ['dashed', 'solid']

    for i, item in enumerate(data):
        plt.plot(item[0], item[1], alpha = 0.7, c=colors[i], ls = linestyles[i])

    # --- Layout.
    plt.xlabel(r'$x$', fontsize=labelsize)
    plt.ylabel(r'$y$', fontsize=labelsize)
    plt.axis([-1.5, 15, -2.5, 0.5])

    # --- Saving figure.
    plt.tight_layout()
    plt.savefig(savename);

def plot_xyz(data, savename):

    fig = plt.figure()
    ax  = fig.gca(projection='3d')

    # --- Plot points in 3D.
    colors = ['r', 'b']
    linestyles = ['dashed', 'solid']

    for i, item in enumerate(data):
        ax.plot(item[0], item[1], item[2], alpha = 0.7, c=colors[i], ls = linestyles[i])

    # --- Layout.
    ax.set_xlabel(r'$x$', fontsize=labelsize)
    ax.set_ylabel(r'$y$', fontsize=labelsize)
    ax.set_zlabel(r'$z$', fontsize=labelsize)

    # --- Saving figure.
    plt.tight_layout()
    plt.show()
    plt.savefig(savename);

def plot_test_2():
    # generate 3D sample data
    X,Y,Z = axes3d.get_test_data(0.05)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X,Y,Z,alpha=0.5)

    plt.savefig('fig/test2.png');
        
# -----------------------------------------------------------------------        
# -------------------------------- MAIN ---------------------------------
# -----------------------------------------------------------------------

plt.switch_backend('TkAgg')
#print(matplotlib.get_backend())

labelsize = 20
ticksize = 15

particle = Particle('input.txt')

v_perp = np.sqrt(particle.u0**2 + particle.v0**2)
v_parall = particle.w0
delta0 = np.pi/2 #np.arctan(particle.u0 / particle.v0)

# --- Analytical in reduced units.
ts_red = np.linspace(0, particle.total_time, 100000)
Axs = - (np.cos(ts_red + delta0) - np.cos(delta0)) + particle.x0/particle.r
Ays = (np.sin(ts_red + delta0) - np.sin(delta0)) + particle.y0/particle.r - particle.E/particle.B*ts_red
Azs = v_parall/v_perp*ts_red + particle.z0/particle.r

# --- Read in trajectory.
data = np.loadtxt('e_output.txt')
positions = data[0:3, :]

plot_xy([ [Axs, Ays, Azs], positions], 'fig/e.png')
print('one')

# --- Read in trajectory.
data = np.loadtxt('mp_output.txt')
positions = data[0:3, :]

plot_xy([ [Axs, Ays, Azs], positions], 'fig/mp.png')
print('two')

# --- Read in trajectory.
data = np.loadtxt('rk4_output.txt')
positions = data[0:3, :]

plot_xy([ [Axs, Ays, Azs], positions], 'fig/rk4.png')
print('three')

# --- Analytical in real units.
ts_real = ts_red / particle.omega
Axs = - particle.r*(np.cos(particle.omega*ts_real + delta0) - np.cos(delta0)) + particle.x0
Ays = particle.r*(np.sin(particle.omega*ts_real + delta0) - np.sin(delta0)) + particle.y0 - particle.E/particle.B*ts_real
Azs = v_parall*ts_real + particle.z0
plot_xyz([ [Axs/particle.r, Ays/particle.r, Azs/particle.r], positions], 'fig/test2.png')