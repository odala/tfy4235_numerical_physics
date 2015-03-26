"""
This script makes it possible to animate a function which changes in time
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation

# ---------------------------------------------------------------------------
# Animate the particle density in time
# ---------------------------------------------------------------------------

# First set up the figure, the axis, and the plot element we want to animate
N = 1000
nSteps = 1100
alfa = 0.2
realisation = '_N1000_rnm12.0_tau.50_dU80.0000'
data = np.loadtxt('trajectory' + realisation + '.txt')
xmin = np.floor(data.min())- (1.0-alfa)
xmax = np.ceil(data.max()) + alfa

fig = plt.figure()
ax = plt.axes(xlim=(xmin, xmax), ylim=(0, N))
line, = ax.plot([], [], lw=2)

# Initialization function: plot the background of each frame.
def init():
    line.set_data([], [])
    return line,

# Animation function. This is called sequentially in range(0, nSteps)
def particle_density(t):
    hist, xbins = np.histogram(data[:, t], bins=(xmax-xmin), range=(xmin, xmax), normed=False)
    x = (xbins[:-1] + xbins[1:]) / 2
    line.set_data(x, hist)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, particle_density, init_func=init,
                               frames=nSteps, interval=500, blit=True)

# Saves the animation as an mp4. This requires ffmpeg or mencoder to be installed.
# The extra_args ensure that the x264 codec is used, so that the video can be
# embedded in html5. For more info, see: http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('particle_density'+realisation+'.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
