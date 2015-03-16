"""
This example shows how to use a path patch to draw a bunch of
rectangles for an animated histogram
"""
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.animation as animation
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})

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
# Animate the particle density in time
# ---------------------------------------------------------------------------

fig, ax = plt.subplots()
ax.set_xlabel('$\mathrm{Position\,[}\mu\mathrm{m]}$')
ax.set_ylabel('$\mathrm{Number\,of\,particles}$')

realisation = '_N1000_rnm12.0_tau.00_dU.0000'

# Read in data and set some constants
N, nSteps, delta_t, alfa, tau, zeta, r, dU, L, kT = read_constants(realisation)
gama = 6.0*np.pi*zeta*r
D = kT/gama
omega = dU/(gama*L**2)
T = nSteps*delta_t/omega

# Convert from reduced units to real units with these scaling factors
scale_length = L               # reduced units -> meter
scale_time = delta_t/omega     # time steps    -> seconds
scale_potential = dU           # reduced units -> joule

# the data to be used
trajectory = np.loadtxt('trajectory' + realisation + '.txt')
xmin = np.floor(trajectory.min())- (1.0-alfa)
xmax = np.ceil(trajectory.max()) + alfa
nbins = xmax - xmin
xmax = xmax*scale_length*1e6
xmin = xmin*scale_length*1e6
trajectory = trajectory*scale_length*1e6

# histogram our data with numpy used to initialise corners
data = trajectory[:,0]
n, bins = np.histogram(data, bins=nbins, range=(xmin, xmax))

# get the corners of the rectangles for the histogram
left = np.array(bins[:-1])
right = np.array(bins[1:])
bottom = np.zeros(len(left))
top = bottom + n
nrects = len(left)

# here comes the tricky part -- we have to set up the vertex and path
# codes arrays using moveto, lineto and closepoly
# for each rect: 1 for the MOVETO, 3 for the LINETO, 1 for the
# CLOSEPOLY; the vert for the closepoly is ignored but we still need
# it to keep the codes aligned with the vertices
nverts = nrects*(1+3+1)
verts = np.zeros((nverts, 2))
codes = np.ones(nverts, int) * path.Path.LINETO
codes[0::5] = path.Path.MOVETO
codes[4::5] = path.Path.CLOSEPOLY
verts[0::5,0] = left
verts[0::5,1] = bottom
verts[1::5,0] = left
verts[1::5,1] = top
verts[2::5,0] = right
verts[2::5,1] = top
verts[3::5,0] = right
verts[3::5,1] = bottom

barpath = path.Path(verts, codes)
patch = patches.PathPatch(barpath, facecolor='green', edgecolor='yellow', alpha=0.5)
ax.add_patch(patch)

ax.set_xlim(left[0], right[-1])
ax.set_ylim(bottom.min(), top.max())

time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    return time_text

def animate(i):
    data = trajectory[:,i]
    n, bins = np.histogram(data, bins=nbins, range=(xmin, xmax))
    top = bottom + n
    verts[1::5,1] = top
    verts[2::5,1] = top
    time_text.set_text('time = %6.1f s' % (i*scale_time))
    return time_text

# frames is number of frames, interval is the time delay between frames i ms
anim = animation.FuncAnimation(fig, animate, frames = nSteps, interval=250, repeat=False, blit=True, init_func=init)
anim.save('density_animation'+realisation+'.mp4', fps=100, extra_args=['-vcodec', 'libx264'])
