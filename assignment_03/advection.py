# Advection equation #
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as spla
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
import time

# --- Schemes to solve the advection equation
#     With periodic boundary conditions.
def upwind(old_u, C):
    n     = len(old_u)
    new_u = np.zeros(n)
    # Internal points.
    for i in range(1, n):
        new_u[i] = old_u[i] - C*(old_u[i] - old_u[i-1])
    # Left boundary.
    new_u[0] = old_u[0] - C*(old_u[0] - old_u[-1])
    return new_u

def downwind(old_u, C):
    n     = len(old_u)
    new_u = np.zeros(n)
    # Internal points.
    for i in range(0, n-1):
        new_u[i] = old_u[i] - C*(old_u[i+1] - old_u[i])
    # Right boundary.
    new_u[n-1] = old_u[n-1] - C*(old_u[0] - old_u[n-1])
    return new_u

def explicit(old_u, C):
    n     = len(old_u)
    new_u = np.zeros(n)
    # Internal points.
    for i in range(1, n-1):
        new_u[i] = old_u[i] - 0.5*C*(old_u[i+1] - old_u[i-1])
    # Boundaries.
    new_u[0]   = old_u[0] - 0.5*C*(old_u[1] - old_u[n-1])
    new_u[n-1] = old_u[n-1] - 0.5*C*(old_u[0] - old_u[n-2])
    return new_u

def implicit(old_u, C):
    n   = len(old_u)
    ld  = -0.5*C*np.ones(n-1)
    d   = np.ones(n)
    ud  = 0.5*C*np.ones(n-1)
    ruc = -0.5*C*np.ones(1)
    llc = 0.5*C*np.ones(1)
    A   = sparse.diags([ld,d,ud, ruc, llc], [-1, 0, 1, n-1, -(n-1)], format='csc')
    return spla.spsolve(A, old_u)

def laxfriedrichs(old_u, C):
    n     = len(old_u)
    new_u = np.zeros(n)
    # Internal points.
    for i in range(1, n-1):
        new_u[i] = 0.5*(old_u[i+1] + old_u[i-1]) - 0.5*C*(old_u[i+1] - old_u[i-1])
    # Boundaries.
    new_u[0]   = 0.5*(old_u[1] + old_u[n-1]) - 0.5*C*(old_u[1] - old_u[n-1])
    new_u[n-1] = 0.5*(old_u[0] + old_u[n-2]) - 0.5*C*(old_u[0] - old_u[n-2])
    return new_u

def laxwendroff(old_u, C):
    n     = len(old_u)
    new_u = np.zeros(n)
    # Internal points.
    for i in range(1, n-1):
        new_u[i] = old_u[i] - 0.5*C*(old_u[i+1] - old_u[i-1]) + 0.5*C**2*(old_u[i-1] - 2*old_u[i] + old_u[i+1])
    # Boundaries.
    new_u[0]   = old_u[0] - 0.5*C*(old_u[1] - old_u[n-1]) + 0.5*C**2*(old_u[n-1] - 2*old_u[0] + old_u[1])
    new_u[n-1] = old_u[n-1] - 0.5*C*(old_u[0] - old_u[n-2]) + 0.5*C**2*(old_u[n-2] - 2*old_u[n-1] + old_u[0])
    return new_u


# Main
def main():

    # --- Initial condition.
    #u0 = lambda x: np.exp(-np.power(x-0.5, 2)/(0.01))   # Gauss
    #u0 = lambda x: 1.0*(np.greater_equal(x, 0.4)) + 1.0*(np.less_equal(x, 0.6)) - 1.0
    u0 = lambda x: np.sin(2*np.pi*x)

    # --- Velocity.
    c = 1.0

    # --- Decide Courant number.
    C = 0.1

    # --- Discretize space.
    xstart = 0.0
    xstop  = 1.0
    Nx     = 101
    x      = np.linspace(xstart, xstop, Nx)
    dx     = x[1] - x[0]

    # --- Discretize time.
    tstart = 0.0
    tstop  = 1.00
    dt     = C*dx/c
    Nt     = (tstop - tstart)/dt + 1
    t      = np.linspace(tstart, tstop, Nt)

    # --- Initialize u(x,t=0).
    uu = ud = ue = ui = uf = uw = u0(x)

    #plt.ion()
    #fig = plt.figure()
    # --- Evolve in time.
    for i in t:
        uu = upwind(uu, C)
        ud = downwind(ud, C)
        ue = explicit(ue, C)
        ui = implicit(ui, C)
        uf = laxfriedrichs(uf, C)
        uw = laxwendroff(uw, C)
        '''
        fig.clear()
        plt.plot(x, u, label='Time: ' + str(i))
        plt.axis([0.0, 1.0, 0.0, 1.1])
        plt.xlabel('Position $x$', fontsize=20); plt.ylabel('$u(x,t)$', fontsize=20)
        time.sleep(0.01)
        plt.draw()'''

    # --- Calculate global truncation error and save to file.
    '''gte = sum(abs(uw-u0(x)))/len(x)
    fout = open(filename,'ab')
    np.savetxt(fout, np.atleast_2d([dt, dx, gte]), delimiter=',')
    fout.close()'''

    print(dt, dx)
    
    # --- Plot.
    plt.ioff()
    plt.figure()
    plt.plot(x, u0(x), 'r-', label='Initial condition')
    plt.plot(x, uu, 'b--', label='Upwind scheme')
    #plt.plot(x, ud, 'y--', label='Downwind scheme')
    #plt.plot(x, ue, 'g--', label='Explicit centered scheme')
    plt.plot(x, ui, 'm--', label='Implicit centered scheme')
    plt.plot(x, uf, 'c-.', label='Lax-Friedrichs scheme')
    plt.plot(x, uw, 'k--', label='Lax-Wendroff scheme')
    plt.axis([0.0, 1.0, 1.1*min(u0(x)), 1.1*max(u0(x))])
    plt.xlabel('Position $x$', fontsize=20); plt.ylabel('$u(x,t)$', fontsize=20)    
    plt.legend(loc='upper right')
    plt.savefig('square.png')
    plt.show()


#
if __name__ == "__main__":
    main()

    # --- Plot global truncation error.
    '''data = np.loadtxt('uu_dx.out', delimiter=',')
    dt = data[:,0]
    dx = data[:,1]
    uu_dx = data[:,2]
    ud_dx = np.loadtxt('ud_dx.out', delimiter=',')[:,2]
    ue_dx = np.loadtxt('ue_dx.out', delimiter=',')[:,2]
    ui_dx = np.loadtxt('ui_dx.out', delimiter=',')[:,2]
    uf_dx = np.loadtxt('uf_dx.out', delimiter=',')[:,2]
    uw_dx = np.loadtxt('uw_dx.out', delimiter=',')[:,2]
    plt.figure()
    plt.plot(dx, uu_dx, 'b', label='Upwind scheme')
    plt.plot(dx, ud_dx, 'g', label='Downwind scheme')
    plt.plot(dx, ue_dx, 'c', label='Explicit centered scheme')
    plt.plot(dx, ui_dx, 'm', label='Implicit centered scheme')
    plt.plot(dx, uf_dx, 'y', label='Lax-Friedrichs scheme')
    plt.plot(dx, uw_dx, 'k', label='Lax-Wendroff scheme')
    plt.axis([0.0, 0.12, 0.0, 0.5])
    plt.xlabel('Spatial step size $dx$', fontsize=20); plt.ylabel('$GTE$', fontsize=20)    
    plt.legend(loc='upper right')
    plt.show()'''




