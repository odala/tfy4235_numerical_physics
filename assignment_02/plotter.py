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
        Nparameters = 15
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
        self.Bz           = parameters[11]
        self.Ex           = parameters[12]
        self.Ey           = parameters[13]
        self.Ez           = parameters[14]

        self.mass = m_p if (self.sgn > 0) else m_e
        
        # --- Calculate additional parameters.
        self.omega = e * self.Bz / self.mass
        self.r = np.sqrt((self.u0-self.Ey/self.Bz)**2 + self.v0**2) / self.omega
        self.total_time = self.Nsteps*self.timestep

def plot_xy(data, savename, xlabel = r'$x$', ylabel = r'$y$'):

    fig = plt.figure()

    # --- Plot points in 2D.
    colors  = ['k', 'b', 'r', 'g']
    markers = ['s', 'D', 'o', 'x']

    for i, item in enumerate(data):
        plt.scatter(item[0], item[1], alpha = 0.9, facecolors='none', edgecolors=colors[i], marker=markers[i])

    # --- Layout.
    plt.xlabel(xlabel, fontsize=labelsize)
    plt.ylabel(ylabel, fontsize=labelsize)
    #plt.axis([-1.5, 15, -2.5, 0.5])

    # --- Saving figure.
    plt.tight_layout()
    plt.savefig(savename + '.png');

def plot_xyz(data, savename):

    fig = plt.figure()
    ax  = fig.gca(projection='3d')

    # --- Plot points in 3D.
    colors  = ['k', 'b', 'r', 'g']
    markers = ['s', 'D', 'o', 'x']

    for i, item in enumerate(data):
        ax.scatter(item[0], item[1], item[2], alpha = 0.9, facecolors='none', edgecolors=colors[i], marker = markers[i])

    if False:
        theta = np.linspace(0, 2*np.pi, 100)
        x     = 0.2*np.cos(theta)
        y     = 0.3*np.sin(theta)
        z1    = [1]*100
        z2    = [-1]*100
        ax.scatter(x, y, z1, alpha = 0.7, c='r', marker = 'o' )
        ax.scatter(x, y, z2, alpha = 0.7, c='r', marker = 'o' )
    
    # --- Layout.
    ax.set_xlabel(r'$x$', fontsize=labelsize)
    ax.set_ylabel(r'$y$', fontsize=labelsize)
    ax.set_zlabel(r'$z$', fontsize=labelsize)

    # --- Saving figure.
    plt.tight_layout()
    '''ax.view_init(elev=90., azim=270)
    plt.savefig(savename + '_xy' + '.png') #OK!
    ax.view_init(elev=0., azim=270)
    plt.savefig(savename + '_xz' + '.png') # OK!
    ax.view_init(elev=0., azim=0)
    plt.savefig(savename + '_zy' + '.png')'''
    #plt.savefig(savename + '.png')
    plt.show()

def plot_test_2():
    # generate 3D sample data
    X,Y,Z = axes3d.get_test_data(0.05)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X,Y,Z,alpha=0.5)

    plt.savefig('fig/test2' + '.png');
        
# -----------------------------------------------------------------------        
# -------------------------------- MAIN ---------------------------------
# -----------------------------------------------------------------------

# --- Some configuration for the plotting.
plt.switch_backend('TkAgg')
#print(matplotlib.get_backend())
labelsize = 20
ticksize = 15

# --- Read in experimental parameters.
particle = Particle('input/input.txt')

# --- Read in exact trajectory (only used for task 1).
data_exact       = np.loadtxt('output/exact_output.txt')
ts_exact         = data_exact[ 0 , :]
positions_exact  = data_exact[1:4, :]
velocities_exact = data_exact[4:7, :]

# --- Read in numerical trajectory.
data       = np.loadtxt('output/rk4_output.txt')
ts         = data[ 0 , :]
positions  = data[1:4, :]
velocities = data[4:7, :]

# --- TASK 1 --- #

# --- Read in trajectories.
positions_e   = np.loadtxt('output/e_output.txt')[1:4, :]
positions_mp  = np.loadtxt('output/mp_output.txt')[1:4, :]
positions_rk4 = np.loadtxt('output/rk4_output.txt')[1:4, :]

# --- Plot trajectories.
plot_xy([ positions_e, positions_mp, positions_rk4], 'fig/projection_xy', r'$\hat{x}$', r'$\hat{y}$')
plot_xy([ positions_e[0:3:2], positions_mp[0:3:2], positions_rk4[0:3:2] ], 'fig/projection_xz', r'$\hat{x}$', r'$\hat{z}$')
plot_xyz([ positions_e, positions_mp, positions_rk4 ], 'fig/projection_xyz')

# --- Read in velocites.
velocities_e   = np.loadtxt('output/e_output.txt')[4:7, :]
velocities_mp  = np.loadtxt('output/mp_output.txt')[4:7, :]
velocities_rk4 = np.loadtxt('output/rk4_output.txt')[4:7, :]

# --- Plot relative variation in kinetic energy.
calculate_kinetic_energy = lambda u, v, w: ( np.sqrt(np.power(u, 2) + np.power(v, 2) + np.power(w, 2)) - np.sqrt(u[0]**2 + v[0]**2 + w[0]**2)) / np.sqrt(u[0]**2 + v[0]**2 + w[0]**2)
plt.semilogy(ts, abs(calculate_kinetic_energy(velocities_e[0], velocities_e[1], velocities_e[2])))
plt.semilogy(ts, abs(calculate_kinetic_energy(velocities_mp[0], velocities_mp[1], velocities_mp[2])))
plt.semilogy(ts, abs(calculate_kinetic_energy(velocities_rk4[0], velocities_rk4[1], velocities_rk4[2])))
#plt.plot(ts, calculate_kinetic_energy(velocities_exact[0], velocities_exact[1], velocities_exact[2]))
plt.show()

#plot_xy([ positions[0:2], positions_exact[0:2]], 'fig/rk4_xy')
#plot_xy([ [positions[0], positions[2]], [positions_exact[0], positions_exact[2]] ], 'fig/rk4_xz', r'$x$', r'$z$')
#plot_xy([ [positions[2], positions[1]], [positions_exact[2], positions_exact[1]] ], 'fig/rk4_zy', r'$z$', r'$y$')
#plot_xyz([ positions, positions_exact ], 'fig/3d_projection_red')


# --- TASK 2 --- #