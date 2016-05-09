# Hopf's equation and shock #
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as spla
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
import time

# --- Lax-Wendroff scheme to solve the Hopf's equation
#     With periodic boundary conditions.
def laxwendroff(u, dt, dx):
    C = dt/dx
    n     = len(u)
    new_u = np.zeros(n)
    # Internal points.
    for i in range(1, n-1):
        new_u[i] = u[i] - 0.25*C*(u[i+1]**2 - u[i-1]**2) + 0.125*C**2*( (u[i+1] + u[i])*(u[i+1]**2 - u[i]**2) - (u[i] + u[i-1])*(u[i]**2 - u[i-1]**2) )
    # Boundaries.
    new_u[0] = u[0] - 0.25*C*(u[1]**2 - u[n-1]**2) + 0.125*C**2*( (u[1] + u[0])*(u[1]**2 - u[0]**2) - (u[0] + u[n-1])*(u[0]**2 - u[n-1]**2) )
    new_u[n-1] = u[n-1] - 0.25*C*(u[0]**2 - u[n-2]**2) + 0.125*C**2*( (u[0] + u[n-1])*(u[i+1]**2 - u[i]**2) - (u[i] + u[i-1])*(u[i]**2 - u[i-1]**2) )
    return new_u


# Main
def main():

    # --- Initial condition.
    #u0 = lambda x: np.exp(-np.power(x-0.5, 2)/(0.05))
    #u0 = lambda x: 1.0*(np.greater_equal(x, 0.4)) + 1.0*(np.less_equal(x, 0.6)) - 1.0
    u0 = lambda x: np.sin(2*np.pi*x)
    #u0 = lambda x: np.sin(np.pi*x)

    # --- Discretize space.
    xstart = 0.0
    xstop  = 1.0
    Nx     = 101
    x      = np.linspace(xstart, xstop, Nx)
    dx     = x[1] - x[0]

    # --- Discretize time.
    tstart = 0.0
    tstop  = 0.23
    Nt     = 201
    t      = np.linspace(tstart, tstop, Nt)
    dt     = t[1] - t[0]

    # --- Initialize u(x,t=0).
    u = u0(x)

    #plt.ion()
    #fig = plt.figure()
    # --- Evolve in time.
    for n in range(1, Nt+1):
        u = laxwendroff(u, dt, dx)
        if n == round(Nt/3):
            u1 = u
        elif n == round(2*Nt/3):
            u2 = u
        '''fig.clear()
        plt.plot(x, u, label='Time: ' + str(n*dt))
        plt.axis([0.0, 1.0, 1.1*min(u0(x)), 1.1*max(u0(x))])
        plt.xlabel('Position $x$', fontsize=20); plt.ylabel('$u(x,t)$', fontsize=20)
        plt.legend(loc='upper right')
        plt.draw()'''


    # --- Plot.
    plt.ioff()
    plt.figure()
    plt.plot(x, u0(x), 'r', label='$t=0$')
    plt.plot(x, u1, 'g', label='$t=$' + str(dt*Nt/3))
    plt.plot(x, u2, 'k', label='$t=$' + str(dt*Nt*2/3))
    plt.plot(x, u, 'b--', label='$t=$' + str(dt*(Nt-1)))
    plt.axis([0.0, 1.0, 1.1*min(u0(x)), 1.1*max(u0(x))])
    plt.xlabel('Position $x$', fontsize=20); plt.ylabel('$u(x,t)$', fontsize=20)    
    plt.legend(loc='upper right')
    plt.show()

#
if __name__ == "__main__":
    main()


