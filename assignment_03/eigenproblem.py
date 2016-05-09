# Spectral Theory and Eigenvalueproblem #
import numpy as np
import scipy.sparse as sparse
import scipy.linalg as la
import scipy.sparse.linalg as spla
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm #colormap
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
import time

def tridiag(n, a=None, b=None, c=None):
    """
    tridiag  tridiagonal matrix (sparse).
         tridiag(a, b, c) is the sparse tridiagonal matrix with
         subdiagonal x, diagonal y, and superdiagonal z.
         x and z must be vectors of dimension one less than y.
         Alternatively tridiag(n, c, d, e), where c, d, and e are all
         scalars, yields the toeplitz tridiagonal matrix of order n
         with subdiagonal elements c, diagonal elements d, and superdiagonal
         elements e.   This matrix has eigenvalues (todd 1977)
                  d + 2*sqrt(c*e)*cos(k*pi/(n+1)), k=1:n.
         tridiag(n) is the same as tridiag(n,-1,2,-1), which is
         a symmetric positive definite m-matrix (the negative of the
         second difference matrix).

         References:
         J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
           Birkhauser, Basel, and Academic Press, New York, 1977, p. 155.
         D.E. Rutherford, Some continuant determinants arising in physics and
           chemistry---II, Proc. Royal Soc. Edin., 63, A (1952), pp. 232-241.
    """
    try:
        # --- First see if they are arrays
        na, = n.shape
        nb, = a.shape
        nc, = b.shape
        if (nb - na - 1) != 0 or  (nb - nc - 1) != 0:
            raise Higham('Dimensions of vector arguments are incorrect.')
        # Now swap to match above
        c = b
        b = a
        a = n

    except AttributeError:
        # --- They are not arrays
        if n < 2:
            raise Higham("n must be 2 or greater")

        if a == None and b == None and c == None:
            a = -1
            b =  2
            c = -1

        a = a * np.ones(n - 1)
        c = c * np.ones(n - 1)
        b = b * np.ones(n)

    except ValueError:
        raise Higham("a, b, c must be all scalars or 1-D vectors")

    # ---- t = diag(a, -1) + diag(b) + diag(c, 1);  % For non-sparse matrix.
    n  = np.max(np.size(b))
    za = np.zeros(1)
    t = sparse.spdiags(np.vstack((np.hstack((a, za)), b,
                                  np.hstack((za, c)))),
                       np.array([-1, 0, 1]), n, n)
    return t

def plot2D(x, y, p):
    fig = plt.figure(figsize = (11,7), dp1 = 100)
    ax = Axes3D(fig)
    X, Y = np.meshgrid(x,y)
    surf = ax.plot_surface(X,Y,p[:],rstride=1, cstride=1, cmap=cm.coolwarm, 
        linewidth=0, antialiased=True)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_zlim(0,1)
    ax.view_init(30,225)

def initialU(x,y):
    return np.exp(-np.power(np.power(x-0.5,2) + np.power(y-0.5, 2), 2) / 0.001)


def main(NI):
    # --- Define number of interior points.
    Nx = Ny = NI

    # --- Define domain.
    Lx = 1
    Ly = 1
    xmin = 0
    xmax = xmin + Lx
    ymin = 0
    ymax = ymin + Ly
    
    # --- Discretize domain.
    x  = np.linspace(xmin, xmax, Nx+2)
    y  = np.linspace(xmin, xmax, Nx+2)
    dx = x[1]-x[0]
    dy = y[1]-y[0]

    # --- Define simplifying constants multiplied by dx^2 / dy^2.
    aP = 4#2/dx**2 + 2/dy**2
    aW = 1#/dx**2
    aE = 1#/dx**2
    aS = 1#/dy**2
    aN = 1#/dy**2

    # --- Initialize all matrices.
    Dj  = tridiag(Nx, -aW, aP, -aE)
    ABj = tridiag(Nx, -aS, 0, -aN)
    ide = sparse.eye(Ny, Ny)
    A   = sparse.kron(ide, Dj) + sparse.kron(ABj, ide)
    #print(np.shape(A))
    #print(A.todense())
    # 2.15827030e+04 2.15827030e+04   2.16122668e+04
    # --- Solve linear system and find Nk eigenvalues and eigenmodes.
    Nk = 100
    val, vec = spla.eigsh(A, k = Nk, sigma=0, which = 'LM')

    # --- Renormalise eigenvalues.
    val = np.divide(val, dx**2)
    
    # --- Plot the 100 smallest eigenvalues. Compare exact with numerical.
    '''# --- Need to set Nk to 100! And have Nx = Ny > 10
    #with open('test.out','a') as f_handle:
    #    np.savetxt(f_handle, [np.ones(100)*Nx, val])
    i = 0
    val_e = np.zeros(400)
    for kx in range(1, 21):
        for ky in range(1,21):
            val_e[i] = np.pi**2 * (kx**2/(xmax-xmin) + ky**2/(ymax-ymin))
            i = i + 1
    val_e = np.sort(val_e)[0:100]
    val_p = np.loadtxt('test.out')
    lfig = plt.figure(figsize = (11,7))
    plt.semilogy(range(1,101), val_p[1,:], 'p', mfc='none', mec='r', label='Numerical ' + str(int(val_p[0,0])) + 'x' + str(int(val_p[0,0])))
    plt.semilogy(range(1,101), val_p[3,:], 's', mfc='none', mec='r', label='Numerical ' + str(int(val_p[2,0])) + 'x' + str(int(val_p[2,0])))
    plt.semilogy(range(1,101), val_p[5,:], '*', mfc='none', mec='r', label='Numerical ' + str(int(val_p[4,0])) + 'x' + str(int(val_p[4,0])))
    plt.semilogy(range(1,101), val_p[7,:], '.r', label='Numerical ' + str(int(val_p[6,0])) + 'x' + str(int(val_p[6,0])))
    plt.semilogy(range(1,101), val_e, '.b', label='Exact')
    plt.xlabel('Eigenvalue index', fontsize=20)
    plt.ylabel('Eigenvalue $\lambda$', fontsize=20) 
    plt.legend(loc='lower right')
    plt.savefig('eigenproblem_eigenvalues' + '.png')'''

    [X, Y] = np.meshgrid(x,y[::-1])
    ks = np.array([ [1,1], [2,1], [1,2], [2,2],[1,3],[3,1], [3,2], [2,3], [1,4], [4,1], [3,3] ])
    norms = np.zeros(Nk)
    gtes = np.zeros(Nk)
    alphas = np.zeros(Nk)
    test = 0
    for i in range(Nk):
        # --- Adjust so that eigenvector is positive in lower left corner (x=0+, y=0+).
        if (vec[0,i] < 0):
            vec[:,i] = -vec[:,i]

        # --- Reshape numerical eigenmodes to Nx by Ny matrix with zero at the boundaries.
        #     x=0,y=0 in U[Ny+1, 0].
        U = np.zeros((Ny+2, Nx+2))
        for j in range(Ny):
            U[Ny-j,1:Nx+1] = vec[j*Nx:Nx+j*Nx,i]
        normalize = dx*dy*np.sum(np.sum(np.square(U[1:Ny+1, 1:Nx+1]), axis=0)) + 0.5*dx*dy*(np.sum(np.square(U[0,:]))+np.sum(np.square(U[Ny+1,:]))+np.sum(np.square(U[:,0]))+np.sum(np.square(U[:,Nx+1])))
        U = np.divide(U, np.sqrt(normalize))
        vec = np.divide(vec, np.sqrt(normalize))

        # --- Calculate exact eigenmodes with 2-norm of error.
        if Nk >10:
            norm = 1.0
            gte  = 1.0
        else:
            Z    = 2/np.sqrt(Lx*Ly)*np.sin(ks[i,0]*np.pi*X/Lx)*np.sin(ks[i,1]*np.pi*Y/Ly)
            norm = np.sqrt(dx*dy*np.sum(np.sum(np.square(U-Z), axis=0)))
            gte  = np.sum(np.sum(np.abs(U-Z), axis=0))/(Nx+2)/(Ny+2)
            Z2 = 2/np.sqrt(Lx*Ly)*np.sin(ks[i,1]*np.pi*X/Lx)*np.sin(ks[i,0]*np.pi*Y/Ly)
            norm2 = np.sqrt(dx*dy*np.sum(np.sum(np.square(U-Z2), axis=0)))
            gte2  = np.sum(np.sum(np.abs(U-Z2), axis=0))/(Nx+2)/(Ny+2)
            if norm2 < norm:
                Z    = Z2
                norm = norm2
                gte  = gte2

            # --- Plot all the eigenmodes.
            fig = plt.figure(figsize = (20,8))
            ax = Axes3D(fig)
            ax.plot_surface(X,Y,U, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
            ax.plot_wireframe(X,Y,Z, rstride=3, cstride=3, color='red')
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('u')
            ax.view_init(30,220)
            ax.set_xlim([0,1]); ax.set_ylim([0,1]); ax.set_zlim([-2.5,2.5])
            plt.savefig('eigenmode_' + str(i) + '.png')

        norms[i] = norm
        gtes[i]  = gte

        # --- Project initial condition on the set of eigenmodes.
        inU = initialU(X, Y)
        alphas[i] = dx*dy*np.sum(np.sum(np.multiply(inU, U), axis=0))
        test = test + alphas[i]*vec[:,i]

    # --- Save global truncation error / 2-norm to file.
    '''fout = open('eigenproblem_norms.out','ab')
    np.savetxt(fout, np.atleast_2d(norms), delimiter=',')
    fout.close()
    fout = open('eigenproblem_gtes.out','ab')
    np.savetxt(fout, np.atleast_2d(gtes), delimiter=',')
    fout.close()'''

    # --- Plot difference between exact solution and u_in = sum_k (alpha_k vec_k).
    test2 = np.zeros((Ny+2, Nx+2))
    for j in range(Ny):
        test2[Ny-j,1:Nx+1] = test[j*Nx:Nx+j*Nx]
    fig = plt.figure(figsize = (20,8))
    ax = Axes3D(fig)
    ax.plot_surface(X,Y,test2 - inU, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
    plt.show()

    # --- Visulise the evolution of the solutions of: 
    #     diffusion, wave and Schrodinger.
    plt.ion()
    fig = plt.figure(figsize = (20,8))
    ax = Axes3D(fig)
    case = 'wave'
    for t in np.linspace(0,1.0,100):
        print('it: ', t)
        func = 0
        if case == 'diffusion':
            for k in range(Nk):
                phi = np.exp(-val[k]*t)
                func = func + alphas[k]*phi*vec[:,k]
        elif case == 'wave':
            for k in range(Nk): 
                phi = np.cos(np.sqrt(val[k])*t)
                func = func + alphas[k]*phi*vec[:,k]
        elif case == 'schrodinger':
            for k in range(Nk):
                phi = np.exp(1.0j*val[k]*t)
                func = func + alphas[k]*phi*vec[:,k]
            func = np.multiply(func, np.conj(func))

        # --- Reshape
        funcU = np.zeros((Ny+2, Nx+2))
        for j in range(Ny):
            funcU[Ny-j,1:Nx+1] = func[j*Nx:Nx+j*Nx]
        ax.clear()
        ax.set_xlim([0,1]); ax.set_ylim([0,1]); ax.set_zlim([-1.0,1.0])
        ax.plot_surface(X,Y,funcU, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
        plt.draw()

if __name__ == "__main__":
    main(51)

    # --- Plot global truncation error for the 10
    #     eigenmodes with the smallest eigenvalues.
    data = np.loadtxt('eigenproblem_gtes.out', delimiter=',')
    err0 = data[:,0]
    err1 = data[:,1]
    err2 = data[:,2]
    err3 = data[:,3]
    err4 = data[:,4]
    err5 = data[:,5]
    err6 = data[:,6]
    err7 = data[:,7]
    err8 = data[:,8]
    err9 = data[:,9]
    Nsteps = np.logspace(1.1, 2.553, num=100)
    Nsteps = np.divide(1.0, (Nsteps+1))

    fig = plt.figure(figsize = (11,7))
    #plt.loglog(Nsteps, err0, 'o', label='Eigenmode $u_{1,1}$')
    plt.loglog(Nsteps, err1, '-', label='Eigenmode $u_{1,2}$ and $u_{2,1}$')
    #plt.semilogx(Nsteps, err2, 'o', label='Eigenmode $u_{2,1}$')
    #plt.loglog(Nsteps, err3, 'o', label='Eigenmode $u_{2,2}$')
    plt.loglog(Nsteps, err4, '-', label='Eigenmode $u_{1,3}$')
    #plt.semilogx(Nsteps, err5, 'o', label='Eigenmode $u_{3,1}$')
    plt.loglog(Nsteps, err6, '-', label='Eigenmode $u_{2,3}$')
    #plt.semilogx(Nsteps, err7, 'o', label='Eigenmode $u_{3,2}$')
    plt.loglog(Nsteps, err8, '-', label='Eigenmode $u_{1,4}$')
    #plt.semilogx(Nsteps, err9, 'o', label='Eigenmode $u_{4,1}$')
    plt.xlabel('Step size $h$', fontsize=20)
    plt.ylabel('Global truncation error', fontsize=20) 
    plt.legend(loc='lower right')
    #plt.savefig('eigenvectors_err' + '.png')
    plt.show()

    '''# --- 11, 21, 31, 41, 51, 61, 71, 81, 91, 101, 501, 1001
    err = np.matrix([ [  4.65259398e-16,   3.75571477e-01,   3.75571477e-01,
         2.87054010e-15,   2.63389841e-01,   2.63389841e-01,
         5.37047846e-01,   5.37047846e-01,   2.91993712e-01,
         2.91993712e-01], [  9.58782958e-16,   5.88879009e-01,   5.88879009e-01,
         9.54669608e-16,   5.44200281e-02,   5.44200281e-02,
         4.36340249e-01,   4.36340249e-01,   6.23258778e-01,
         6.23258778e-01], [  6.68660661e-16,   6.66354170e-01,   6.66354170e-01,
         1.22486183e-15,   7.03211254e-01,   7.03211254e-01,
         4.99290865e-01,   4.99290865e-01,   6.22143306e-01,
         6.22143306e-01], [  1.64905472e-15,   1.26357868e-02,   1.26357868e-02,
         1.84631769e-15,   7.40773526e-01,   7.40773526e-01,
         6.38482539e-01,   6.38482539e-01,   7.06401551e-01,
         7.06401551e-01], [  3.67517890e-15,   6.63746938e-01,   6.63746938e-01,
         5.64725264e-15,   2.84230379e-01,   2.84230379e-01,
         4.96188951e-01,   4.96188951e-01,   5.38128085e-01,
         5.38128085e-01], [  2.31514810e-15,   2.71136583e-01,   2.71136583e-01,
         3.64884498e-15,   3.30171223e-02,   3.30171223e-02,
         7.55092275e-01,   7.55092275e-01,   2.79353038e-01,
         2.79353038e-01] ])
    
    err = np.matrix([  
        [  4.83478571e-16,   5.62991244e-01,   5.62991244e-01,
         1.70521125e-15,   4.98034175e-01,   1.01750121e+00,
         5.77556994e-01,   5.77556994e-01,   1.19130218e+00,
         2.93579699e-01,   1.05818730e-14],
         [  1.07442328e-15,   5.88879009e-01,   5.88879009e-01,
         9.07368870e-16,   5.44200281e-02,   5.44200281e-02,
         4.13103762e-01,   4.13103762e-01,   9.02208519e-01,
         6.24187196e-01,   6.71191333e-15],
         [  6.65907087e-16,   1.36536077e-01,   1.36536077e-01,
         1.27281063e-15,   6.89810953e-01,   6.89810953e-01,
         4.86765241e-01,   4.86765241e-01,   6.33094060e-01,
         6.33094060e-01,   6.22069753e-15],
         [  1.55162325e-15,   7.83378611e-02,   7.83378611e-02,
         2.68038430e-15,   2.75305122e-01,   2.75305122e-01,
         6.26413840e-01,   6.26413840e-01,   7.14721701e-01,
         8.15443378e-01,   6.66313855e-15],
         [  3.67147617e-15,   6.51259801e-01,   6.51259801e-01,
         5.11559479e-15,   1.85471115e-01,   1.85471115e-01,
         2.42900602e-01,   2.42900602e-01,   1.69749358e-01,
         1.69749358e-01,   6.43077401e-15],
         [  2.22653014e-15,   6.91241255e-01,   6.91241255e-01,
         3.76229514e-15,   3.30171223e-02,   3.30171223e-02,
         6.06447699e-01,   6.06447699e-01,   3.10613292e-01,
         3.10613292e-01,   4.95302509e-15],
         [  1.92852449e-15,   5.10139123e-01,   5.10139123e-01,
         2.05422878e-15,   5.68073186e-01,   5.68073186e-01,
         7.18166879e-01,   7.18166879e-01,   2.82679396e-01,
         2.82679396e-01,   1.14232780e-14],
         [  2.75878546e-15,   5.61593899e-01,   5.61593899e-01,
         3.21197461e-15,   6.93555416e-01,   6.93555416e-01,
         4.84794743e-02,   4.84794743e-02,   3.08753768e-01,
         3.08753768e-01,   1.47668421e-14],
         [  4.13419568e-15,   6.75037642e-01,   6.75037642e-01,
         4.50984556e-15,   6.75414229e-01,   6.75414229e-01,
         3.86923735e-01,   3.86923735e-01,   6.75485771e-01,
         6.75485771e-01,   8.98676029e-15],
         [  1.65189009e-14,   6.26027603e-01,   6.26027603e-01,
         1.78866715e-14,   7.55945877e-01,   7.55945877e-01,
         7.58913930e-01,   7.58913930e-01,   4.70523366e-01,
         4.70523366e-01,   1.40743251e-14],
         [  1.09861148e-13,   5.99797574e-01,   5.99797574e-01,
         8.85295729e-14,   5.58687299e-01,   5.58687299e-01,
         5.83107557e-01,   5.83107557e-01,   5.06328946e-01,
         5.06328946e-01,   1.28557844e-13],
         [  2.15743797e-13,   2.44302002e-01,   2.44302002e-01,
         1.96526310e-13,   2.00313610e-01,   2.00313610e-01,
         3.97228546e-01,   3.97228546e-01,   1.27354955e-01,
         1.27354955e-01,   1.23013073e-13]
        ])'''

    '''fig = plt.figure(figsize = (11,7))
    plt.plot(range(1,11), err[0,:], '.b--', label='Numerical 11x11')
    plt.plot(range(1,11), err[1,:], '.r--', label='Numerical 21x21')
    plt.plot(range(1,11), err[2,:], '.g--', label='Numerical 31x31')
    plt.plot(range(1,11), err[3,:], '.k--', label='Numerical 41x41')
    plt.plot(range(1,11), err[4,:], '.c--', label='Numerical 51x51')
    plt.xlabel('Eigenvalue index', fontsize=20)
    plt.ylabel('Eigenvalue $\lambda$', fontsize=20) 
    plt.legend(loc='lower right')
    plt.savefig('eigenvectors_err' + '.png')'''

    '''fig = plt.figure(figsize = (11,7))
    print(np.shape(err))
    print(np.shape(range(1,6)))
    hs = [1.0/12, 1.0/22, 1.0/32, 1.0/42, 1.0/52, 1.0/62]
    print(hs)
    plt.plot(hs, err[:,0], '.--', label='Eigenvalue (1,1)')
    plt.plot(hs, err[:,1], '.--', label='Eigenvalue (2,1)')
    plt.plot(hs, err[:,2], '.--', label='Eigenvalue (1,2)')
    plt.plot(hs, err[:,3], '.--', label='Eigenvalue (2,2)')
    plt.plot(hs, err[:,4], '.--', label='Eigenvalue (3,1)')
    plt.plot(hs, err[:,5], '.--', label='Eigenvalue (1,3)')
    plt.plot(hs, err[:,5], '.--', label='Eigenvalue (3,2)')
    plt.plot(hs, err[:,5], '.--', label='Eigenvalue (2,3)')
    plt.plot(hs, err[:,5], '.--', label='Eigenvalue (4,1)')
    plt.plot(hs, err[:,5], '.--', label='Eigenvalue (1,4)')
    plt.xlabel('Eigenvalue index', fontsize=20)
    plt.ylabel('Eigenvalue $\lambda$', fontsize=20) 
    plt.legend(loc='lower right')
    plt.savefig('eigenvectors_err_2' + '.png')
    '''
    


