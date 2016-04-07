# TDMA module
import numpy as np

# --- Solution of a linear system of algebraic 
#     equations with a tri-diagonal matrix of coefficients
#     using the Thomas-algorithm (No pivoting).
#     Equation no. i :
#       a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i), i = 1,..,n
#     here n is the number of equations.
#     
#     === Input ===
#     Vectors may be in either row - or column-form
#
#     a[0:n-1] .. Lower diagonal. a[0] is not used
#     b[0:n-1] .. Main diagonal. 
#     c[0:n-1] .. Upper diagonal. c[n-1] is not used
#     d[0:n-1] .. Right hand side of the system.
#
#     === Output ===
#     x[0:n-1] .. The solution vector
def tdma(a,b,c,d):

    # --- Get length of arrays.
    n = len(b)

    # --- Allocate array for solution.
    x = np.zeros(b.shape)
    
    # --- Elimination
    c[0] = c[0] / b[0]
    d[0] = d[0] / b[0]
    for k in range(1,n):
        c[k] = c[k] / (b[k] - a[k]*c[k-1])
        d[k] = ( d[k] - a[k]*d[k-1] ) / ( b[k] - a[k]*c[k-1] )
        #q = a[k]/b[k-1]
        #b[k] = b[k] - c[k-1]*q
        #d[k] = d[k] - d[k-1]*q

    
    # --- Backsubstitution
    #q = d[n-1]/b[n-1]
    #x[n-1] = q
    x[n-1] = d[n-1]
    for k in range(n-2, 0, -1):
        #q = (d[k] - c[k]*q)/b[k]
        #x[k] = q
        x[k] = d[k] - c[k]*x[k+1]

    return x