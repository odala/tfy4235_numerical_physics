# Diffusion
# from myIterativeMethods import explicit_euler, implicit_euler, ...

import sys  # for error massages
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
import time
import scipy.linalg
from scipy.special import erf

# --- Import my functions.
from tdma import tdma

# --- Function that iterate through one time step
#     using the explicit Euler method. Replace the 
#     time derivative by a forward difference and 
#     the spatial derivative by a central difference.
def explicit_euler(old_y, dt, dx, bc):

    C = dt/dx**2 * get_diffusivity(0)
    new_y = np.zeros(len(old_y))
    jmax   = len(old_y)

    # --- Set boundary conditons.
    if bc == 'reflective':
        new_y[0]      = old_y[0] + 2*C*(old_y[1] - old_y[0])
        new_y[jmax-1] = old_y[jmax-1] + 2*C*(old_y[jmax-2]-old_y[jmax-1])
    elif bc == 'absorbing':
        new_y[0]      = 0
        new_y[jmax-1] = 0

    # --- Calculate all internal spatial points.
    for i in range(1, jmax-1):
        new_y[i] = old_y[i] + C*(old_y[i+1] - 2*old_y[i] + old_y[i-1])

    return new_y

# --- Function that iterate through one time step
#     using the implicit Euler method. Replace the 
#     time derivative by a backward difference and 
#     the spatial derivative by a central difference.
def implicit_euler(old_y, dt, dx, bc):

    C = dt/dx**2 * get_diffusivity(0)

    jmax  = len(old_y)

    # --- Data strutures for the linear system.
    A = np.zeros((jmax, jmax))
    for i in range(1,jmax-1):
        A[i,i+1] = -C
        A[i,  i] = 1 + 2*C
        A[i,i-1] = -C

    # --- Set boundary conditons.
    A[0,0] = A[jmax-1, jmax-1] = 1 + 2*C
    if bc == 'reflective':
        A[0,1] = -2*C
        A[jmax-1, jmax-2] = -2*C
    elif bc == 'absorbing':
        A[0,1] = 0
        A[jmax-1, jmax-2] = 0

    '''a = np.ones(jmax)*(-C)
    b = np.ones(jmax)*(1+2*C)
    c = np.ones(jmax)*(-C)

    # --- Set boundary conditons.
    if bc == 'reflective':
        c[0] = -2*C
        a[jmax-1] = -2*C
    elif bc == 'absorbing':
        c[0] = 0
        a[jmax-1] = 0
    new_y = tdma(a,b,c,old_y)
    '''

    # --- Calculate all internal spatial points.
    new_y = scipy.linalg.solve(A, old_y)

    return new_y

# --- Function that iterate through one time step
#     using the Crank-Nicolson method. Using a 50-50
#     combination of the explicit and the implicit 
#     Euler method.
def crank_nicolson(old_y, dt, dx, bc):

    C = dt/dx**2
    jmax  = len(old_y)

    # --- Data strutures for the linear system.
    A = np.zeros((jmax, jmax))
    b = np.zeros(jmax)
    for j in range(1,jmax-1):
        D0 = get_diffusivity((j-0.5)/(jmax-1))
        D1 = get_diffusivity((j+0.5)/(jmax-1))
        A[j,j-1] = -0.5*C*D0
        A[j,  j] = 1 + 0.5*C*(D0+D1)
        A[j,j+1] = -0.5*C*D1
        b[j] = 0.5*C*D0*old_y[j-1] + ( 1-0.5*C*(D0+D1) )*old_y[j] + 0.5*C*D1*old_y[j+1]

    # --- Set boundary conditons.
    D0 = get_diffusivity((-0.5)/(jmax-1))
    D1 = get_diffusivity((0.5)/(jmax-1))
    A[0,0] = A[jmax-1, jmax-1] = 1+0.5*C*(D0+D1)
    if bc == 'reflective':
        A[0,1] = -0.5*C*(D0+D1)
        b[0]   = 0.5*C*(D0+D1)*old_y[1] + (1-0.5*C*(D0+D1))*old_y[0]
        D0 = get_diffusivity((jmax-1.5)/(jmax-1))
        D1 = get_diffusivity((jmax-0.5)/(jmax-1))
        A[jmax-1, jmax-2] = -0.5*C*(D0+D1)
        b[jmax-1] = (1-0.5*C*(D0+D1))*old_y[jmax-1] + 0.5*C*(D0+D1)*old_y[jmax-2]
    elif bc == 'absorbing':
        A[0,1] = 0
        A[jmax-1, jmax-2] = 0
        b[0] = b[jmax-1] = 0

    # --- Calculate all internal spatial points.
    new_y = scipy.linalg.solve(A, b)

    return new_y

def exact_unbounded(x, t, u0, xstart, x0, xstop):
    return u0/np.sqrt(4*np.pi*get_diffusivity((x-xstart)/(xstop-xstart))*t)*np.exp(-np.power(x-x0, 2)/(4*get_diffusivity((x-xstart)/(xstop-xstart))*t))

def exact_unbounded_step_diffusivity(xs, t, u0, xstart, x0, xstop, dx):
    u = np.zeros(len(xs))
    for x, j in zip(xs, range(0, len(xs))):
        D0 = get_diffusivity((x-0.5*dx-xstart)/(xstop-xstart))
        D1 = get_diffusivity((x+0.5*dx-xstart)/(xstop-xstart))
        A1 = 2/(1 + erf(x0/np.sqrt(4*D1*t)) + np.sqrt(D0/D1)*np.exp((D1-D0)*x0**2/4/D1/D0/t)*(1-erf(x0/np.sqrt(4*D0*t))) )
        if x >= 0:
            u[j] = u0*A1/np.sqrt(4*np.pi*D1*t)*np.exp(-(x-x0)**2/4/D1/t)
        else:
            A0   = A1*np.sqrt(D0/D1)*np.exp((D1-D0)*x0**2/4/D1/D0/t)
            u[j] = u0*A0/np.sqrt(4*np.pi*D0*t)*np.exp(-(x-x0)**2/4/D0/t)
    return u

def exact_bounded(x, t, u0, xstart, x0, xstop, bc):
    L = xstop-xstart
    xRel = x - xstart
    if bc == 'reflective':
        u = 1/L*u0*np.ones(len(x))
        for n in range(1, 100):
                u += u0 * np.exp(- (n*np.pi/L)**2*get_diffusivity(xRel/L)*t) * 2/L * np.cos(n*np.pi*x0/L) * np.cos(n*np.pi*x/L)
    elif bc == 'absorbing':
        u = np.zeros(len(x))
        for n in range(1, 100):
            #print(np.exp(- (n*np.pi/L)**2*get_diffusivity(xRel/L)*t))
            #print(2/L * np.sin(n*np.pi*x0/L) * np.sin(n*np.pi*x/L))
            u += u0 * np.exp(- (n*np.pi/L)**2*get_diffusivity(xRel/L)*t) * 2/L * np.sin(n*np.pi*x0/L) * np.sin(n*np.pi*x/L)
    return u

# --- Takes in a number between 0.0 and 1.0.
def get_diffusivity(percentil):
    # --- Constant diffusivity
    #return 1.0
    # --- One step
    #return (percentil<0.5)*0.05 + (percentil>=0.5)*0.1
    # --- Two steps
    #return (percentil<0.3)*10 + (percentil>=0.3 and percentil<0.7)*1 + (percentil>=0.7)*10
    # --- Gaussian.
    #return np.exp(-np.power(percentil-0.5, 2)/(0.05))
    # --- Sinusoidal.
    #return np.abs(np.sin(5*np.pi*percentil)) + 1.e-6
    # --- Triangular.
    return 2*abs(percentil-0.5) + 1.e-6
###
def main(jmax, nmax, bc, bestC=False):

    # --- Discretization of space.
    xstart = 0.0
    xstop  = 1.0
    jmax   = jmax
    x      = np.linspace(xstart, xstop, jmax)
    dx     = x[1] - x[0]

    # --- Discretization of time.
    tstart = 0.0
    tstop  = 0.15
    nmax   = nmax
    if bestC == True:
        C = 1.0/6
        dt     = C*dx**2/get_diffusivity(0)
        nmax   = int((tstop-tstart)/dt) + 1
    t      = np.linspace(tstart, tstop, nmax)
    dt     = t[1] - t[0]

    # --- Allocate memory for spatial array.
    u0 = 1.0
    u = np.zeros(jmax)

    # --- Set initial condition u(x,0) = delta(x - x0).
    x0 = 0.5
    #for i in range(0, jmax):
    #    if x[i] + dx > x0 and x[i] - dx < x0:
    #        u[i] = u0/dx
    u = exact_unbounded(x, dt, u0, xstart, x0, xstop)
    mass1 = sum(u[:-1] + u[1:]) / 2 * dx
    flux1 = get_diffusivity(0.0)*(u[0] - u[1])/dx + get_diffusivity(1.0)*(u[len(u)-1] - u[len(u)-2])/dx

    #plt.ion()
    #plt.figure()

    # --- Iterate through time.
    for n in range(1, nmax):
        u = crank_nicolson(u, dt, dx, bc)
        if n == round(nmax/3):
            u1 = u
            #print('2-norm: ', np.sqrt(dx*np.power(sum(u1 - exact_bounded(x, round(nmax/3)*dt, u0, xstart, x0, xstop, bc)),2)))
            #print('T. err: ', max(abs(u1 - exact_bounded(x, round(nmax/3)*dt, u0, xstart, x0, xstop, bc))))
        elif n == round(2*nmax/3):
            u2 = u
            #print('2-norm: ', np.sqrt(dx*np.power(sum(u2 - exact_bounded(x, round(2*nmax/3)*dt, u0, xstart, x0, xstop, bc)),2)))
            #print('T. err: ', max(abs(u2 - exact_bounded(x, round(2*nmax/3)*dt, u0, xstart, x0, xstop, bc))))
        
        # --- Plot time evolution.
        '''plt.axis([xstart, xstop, 0.0, 1.5])
        plt.plot(x, u, label=r'experimental bounded '+bc)
        plt.plot(x, exact_bounded(x, n*dt, u0, xstart, x0, xstop, bc), label = 'theoretical bounded '+bc)
        #plt.plot(x, exact_unbounded(x, n*dt, u0, xstart, x0, xstop), label = r'theoretical unbounded')
        plt.plot(x, exact_unbounded_step_diffusivity(x, n*dt, u0, xstart, x0, xstop, dx), label = r'theoretical unbounded')
        plt.title(n*dt)
        plt.legend()
        plt.draw()
        plt.clf()
        #time.sleep(0.01)'''

        # --- Check if rate of mass change is consistent with outward flux.
        '''mass2 = sum(u[:-1] + u[1:]) / 2 * dx
        mass_change_rate = (mass2 - mass1)/dt
        if bc == 'reflective':
            print('Rate of mass change = ', mass_change_rate)
        elif bc == 'absorbing':
            flux2 = get_diffusivity(0.0)*(u[0] - u[1])/dx + get_diffusivity(1.0)*(u[len(u)-1] - u[len(u)-2])/dx
            print('Rate of mass change - flux at boundaries = ', mass_change_rate - (flux1+flux2)/2)
        mass1 = mass2
        flux1 = flux2'''

    plt.ioff()
    # --- Plot at end time.
    plt.figure()
    plt.axis([xstart, xstop, 0.0, 6.0])
    plt.plot(x, u, 'o', mew=1, mec='b', mfc='None', label=r'experimental bounded '+bc)
    plt.plot(x, exact_bounded(x, nmax*dt, u0, xstart, x0, xstop, bc), 'b', label = 'theoretical bounded '+bc)
    plt.plot(x, u1, 'o', mew=1, mec='g', mfc='None')
    plt.plot(x, exact_unbounded_step_diffusivity(x, round(nmax/3)*dt, u0, xstart, x0, xstop, dx), 'g')
    plt.plot(x, u2, 'o', mew=1, mec='y', mfc='None')
    plt.plot(x, exact_unbounded_step_diffusivity(x, round(2*nmax/3)*dt, u0, xstart, x0, xstop, dx), 'y')
    plt.xlabel('Position $x$', fontsize=20); plt.ylabel('$u(x,t)$', fontsize=20)
    #plt.legend()
    plt.savefig('diffusion_step_diffusivity.png')
    plt.show()

    err = np.sqrt(dx*np.power(sum(u - exact_bounded(x, nmax*dt, u0, xstart, x0, xstop, bc)),2))
    print('2-norm: ', err, 'C = ', dt/dx**2*get_diffusivity(0))

    #with open('diffusion_abs_exp_err_dx_best.out','a') as fout:
    #    np.savetxt(fout, np.array([dx, dt, err]).reshape(1,3))

    #print('dt = ', dt, 'dx = ', dx)

if __name__ == "__main__":
    # --- Make the get_diffusivity function 
    get_diffusivity = np.vectorize(get_diffusivity)

    nmax = 501
    jmax = 101
    main(jmax, nmax, 'absorbing')
    #for jmax in np.logspace(1.0, 3.0, 50):
    #    jmax = round(jmax)
    #    print(jmax)
    #    main(jmax, nmax, 'absorbing', bestC=True)

    '''# --- Plot global truncation errors in dx.
    # Crank-Nicolson.
    data = np.loadtxt('diffusion_abs_crank_err_dx.out')
    dx = data[:,0]
    cerr = data[:,2]
    data = np.loadtxt('diffusion_ref_crank_err_dx.out')
    dx2 = data[:,0]
    cerr2 = data[:,2]
    # Implicit Euler.
    data = np.loadtxt('diffusion_abs_imp_err_dx.out')
    idx = data[:,0]
    ierr = data[:,2]
    data = np.loadtxt('diffusion_ref_imp_err_dx.out')
    idx2 = data[:,0]
    ierr2 = data[:,2]
    # Explicit Euler.
    data = np.loadtxt('diffusion_abs_exp_err_dx.out')
    edx = data[:,0]
    eerr = data[:,2]
    data = np.loadtxt('diffusion_ref_exp_err_dx.out')
    edx2 = data[:,0]
    eerr2 = data[:,2]
    data = np.loadtxt('diffusion_abs_exp_err_dx_best.out')
    edx3 = data[:,0]
    edt3 = data[:,1]
    eerr3 = data[:,2]

    plt.figure()
    plt.loglog(dx, cerr, 'bo', label='Crank-Nicolson scheme (absorbing)')
    plt.loglog(dx2, cerr2, 'bv', label='Crank-Nicolson scheme (reflective)')
    plt.loglog(idx, ierr, 'ro', label='Implicit Euler scheme (absorbing)')
    plt.loglog(idx2, ierr2, 'rv', label='Implicit Euler scheme (reflective)')
    plt.loglog(edx, eerr, 'go', label='Explicit Euler scheme (absorbing)')
    plt.loglog(edx2, eerr2, 'gv', label='Explicit Euler scheme (reflective)')
    plt.axis([min(dx), max(dx), min(cerr2), max(eerr)])
    plt.xlabel('Spatial step size $dt$', fontsize=20); plt.ylabel('Global truncation error', fontsize=20)    
    #plt.legend(loc='lower right')
    plt.savefig('diffusion_gte_dx.png')
    plt.show()

    # --- Plot global truncation errors in dt.
    # Crank-Nicolson.
    data = np.loadtxt('diffusion_abs_crank_err_dt.out')
    dt = data[:,1]
    cerr = data[:,2]
    data = np.loadtxt('diffusion_ref_crank_err_dt.out')
    dt2 = data[:,1]
    cerr2 = data[:,2]
    # Implicit Euler.
    data = np.loadtxt('diffusion_abs_imp_err_dt.out')
    idt = data[:,1]
    ierr = data[:,2]
    data = np.loadtxt('diffusion_ref_imp_err_dt.out')
    idt2 = data[:,1]
    ierr2 = data[:,2]
    # Explicit Euler.
    data = np.loadtxt('diffusion_abs_exp_err_dt.out')
    edt = data[:,1]
    eerr = data[:,2]
    data = np.loadtxt('diffusion_ref_exp_err_dt.out')
    edt2 = data[:,1]
    eerr2 = data[:,2]
    
    plt.figure()
    plt.loglog(dt, cerr, 'bo', label='Crank-Nicolson scheme (absorbing)')
    plt.loglog(dt2, cerr2, 'bv', label='Crank-Nicolson scheme (reflective)')
    plt.loglog(idt, ierr, 'ro', label='Implicit Euler scheme (absorbing)')
    plt.loglog(idt2, ierr2, 'rv', label='Implicit Euler scheme (reflective)')
    plt.loglog(edt, eerr, 'go', label='Explicit Euler scheme (absorbing)')
    plt.loglog(edt2, eerr2, 'gv', label='Explicit Euler scheme (reflective)')
    plt.axis([min(dt), max(dt), min(cerr2), max(eerr)])
    plt.xlabel('Temporal step size $dt$', fontsize=20); plt.ylabel('Global truncation error', fontsize=20)    
    #plt.legend(loc='lower right')
    plt.savefig('diffusion_gte_dt.png')
    plt.show()'''


