# Iterative methods module
# from myIterativeMethods import explicit_euler, implicit_euler, ...

import sys  # for error massages
import numpy as np
from matplotlib import pyplot as plt
#from matplotlib import rc
#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)
import time
import scipy.linalg

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
    A[0,0] = A[jmax-1, jmax-1] = 1 + 2*C

    # --- Set boundary conditons.
    if bc == 'reflective':
        A[0,1] = -2*C
        A[jmax-1, jmax-2] = -2*C
    elif bc == 'absorbing':
        A[0,1] = 0
        A[jmax-1, jmax-2] = 0

    a = np.ones(jmax)*(-C)
    b = np.ones(jmax)*(1+2*C)
    c = np.ones(jmax)*(-C)

    # --- Set boundary conditons.
    if bc == 'reflective':
        c[0] = -2*C
        a[jmax-1] = -2*C
    elif bc == 'absorbing':
        c[0] = 0
        a[jmax-1] = 0

    # --- Calculate all internal spatial points.
    t0 = time.time()
    new_y = scipy.linalg.solve(A, old_y)
    #print('Scipy: ', time.time() - t0)
    #
    t0 = time.time()
    new_y = tdma(a,b,c,old_y)
    #print('TDMA:  ', time.time() - t0)

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

def get_diffusivity(percentil):
    return 0.1
    #return (percentil>=0.5)*0.1 + (percentil<0.5)*0.01


### 

# --- Discretization of space.
xstart = 0.0
xstop  = 1.0
jmax   = 11
x      = np.linspace(xstart, xstop, jmax)
dx     = x[1] - x[0]

# --- Discretization of time.
tstart = 0.0
tstop  = 10.0
nmax   = 101
#dt     = C*dx**2/D
#nmax   = int((tstop-tstart)/dt) + 1
print('nmax: ', nmax)
t      = np.linspace(tstart, tstop, nmax)
dt     = t[1] - t[0]

# --- Allocate memory for spatial array.
u0 = 1.0
u = np.zeros(jmax)
print('dt = ', dt, 'dx = ', dx, 'C(0) = ', dt*dx**2*get_diffusivity(0))

# --- Set initial condition u(x,0) = delta(x - x0).
x0 = 0.5
for i in range(0, jmax):
    if x[i] + dx > x0 and x[i] - dx < x0:
        u[i] = u0/dx
mass1 = sum(u[:-1] + u[1:]) / 2 * dx

plt.ion()
plt.show()

bc = 'reflective'

# --- Iterate through time.
for n in range(1, nmax):
    u = crank_nicolson(u, dt, dx, bc)
    plt.axis([xstart, xstop, 0.0, 1.5])
    plt.plot(x, u, label=r'experimental bounded '+bc)
    plt.plot(x, exact_bounded(x, n*dt, u0, xstart, x0, xstop, bc), label = 'theoretical bounded '+bc)
    plt.plot(x, exact_unbounded(x, n*dt, u0, xstart, x0, xstop), label = r'theoretical unbounded')
    plt.legend()
    plt.draw()
    plt.clf()
    time.sleep(0.01)
    mass2 = sum(u[:-1] + u[1:]) / 2 * dx
    mass_change_rate = (mass2 - mass1)/dt
    '''if bc == 'reflective':
        print('Rate of mass change = ', mass_change_rate)
    elif bc == 'absorbing':
        flux = D*(u[0] - u[1])/dx + D*(u[len(u)-1] - u[len(u)-2])/dx
        print('Rate of mass change - flux at boundaries = ', mass_change_rate - flux)

    mass1 = mass2'''

err = max(abs(exact_bounded(x, nmax*dt, u0, xstart, x0, xstop, bc) - u))
print(err)

with open('test.out','a') as f_handle:
    np.savetxt(f_handle, np.array([dx, dt, err]), newline=' ')




