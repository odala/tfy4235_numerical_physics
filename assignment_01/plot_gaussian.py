import numpy as np
import matplotlib
matplotlib.use('Agg');
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# ---------------------------------------------------------------
# Standard Gaussian Random Number Generator
# Checks if the gaussian Box-Muller algorithm really generates 
# numbers which are from a normal distribution.
# ---------------------------------------------------------------
def plot_gaussian(): 

    # --- Normal distribution function.
    normal_dist = lambda x: 1/np.sqrt(2*np.pi)*np.exp(-x**2/2)

    # --- Read data from file.
    rand_gauss = np.loadtxt('check_gaussian.txt')
    hist, bins = np.histogram(rand_gauss, bins=100, normed=True)
    center = (bins[:-1] + bins[1:]) / 2

    # --- Plot figure.
    gaussian = plt.figure()         
    plt.plot(center, hist,'ro')
    x = np.linspace(-5, 5, 1000)  
    plt.plot(x, normal_dist(x), 'b')

    # --- Axis labels.
    plt.xlabel(r'$x$', fontsize=15)
    plt.ylabel(r'$P(x)$', fontsize=15)
    
    # Saving figure.
    plt.savefig('fig/normal_distribution' + '.png');

# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

plot_gaussian()