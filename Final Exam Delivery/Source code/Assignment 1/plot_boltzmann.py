import numpy as np
import matplotlib
matplotlib.use('Agg');
import matplotlib.pyplot as plt

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
    
