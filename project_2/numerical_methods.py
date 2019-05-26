''' lam,x = jacobi(a,tol = 1.0e-9).
    Solution of std. eigenvalue problem [a]{x} = lam{x}
    by Jacobi's method. Returns eigenvalues in vector {lam}
    and the eigenvectors as columns of matrix [x].
'''
from __future__ import division
import numpy as np
from numpy import array,identity,diagonal
from math import sqrt

def jacobi(a,tol = 1.0e-8): # Jacobi method

    def maxElem(a): # Find largest off-diag. element a[k,l]
        n = len(a)
        aMax = 0.0
        for i in range(n-1):
            for j in range(i+1,n):
                if abs(a[i,j]) >= aMax:
                    aMax = abs(a[i,j])
                    k = i; l = j
        return aMax,k,l

    def rotate(a,p,k,l): # Rotate to make a[k,l] = 0
        n = len(a)
        aDiff = a[l,l] - a[k,k]

        if abs(a[k,l]) < abs(aDiff)*1.0e-36:
            t = a[k,l]/aDiff

        else:
            phi = aDiff/(2.0*a[k,l])
            t = 1.0/(abs(phi) + sqrt(phi**2 + 1.0))
            if phi < 0.0:
                t = -t

        c = 1.0/sqrt(t**2 + 1.0); s = t*c
        tau = s/(1.0 + c)
        temp = a[k,l]
        a[k,l] = 0.0
        a[k,k] = a[k,k] - t*temp
        a[l,l] = a[l,l] + t*temp
        for i in range(k):      # Case of i < k
            temp = a[i,k]
            a[i,k] = temp - s*(a[i,l] + tau*temp)
            a[i,l] = a[i,l] + s*(temp - tau*a[i,l])
        for i in range(k+1,l):  # Case of k < i < l
            temp = a[k,i]
            a[k,i] = temp - s*(a[i,l] + tau*a[k,i])
            a[i,l] = a[i,l] + s*(temp - tau*a[i,l])
        for i in range(l+1,n):  # Case of i > l
            temp = a[k,i]
            a[k,i] = temp - s*(a[l,i] + tau*temp)
            a[l,i] = a[l,i] + s*(temp - tau*a[l,i])
        for i in range(n):      # Update transformation matrix
            temp = p[i,k]
            p[i,k] = temp - s*(p[i,l] + tau*p[i,k])
            p[i,l] = p[i,l] + s*(temp - tau*p[i,l])

    n = len(a)
    maxRot = 5*(n**2)       # Set limit on number of rotations
    p = identity(n)*1.0     # Initialize transformation matrix
    for i in range(maxRot): # Jacobi rotation loop
        aMax,k,l = maxElem(a)
        if aMax < tol:
            return diagonal(a), p
        rotate(a,p,k,l)
    print 'Jacobi method did not converge'



def buckling_beam_tridiagonal(dim):
    """Function for initializing a tridiagonal (dim x dim) matrix
    for the buckling beam problem"""


    r_max = 1
    step = r_max / (dim+1)


    A = np.zeros((dim, dim))
    a = -1.0 / (step * step)
    d = 2.0 / (step * step)
    A[0, 0] = d
    A[0, 1] = a
    A[dim-1, dim-2] = a
    A[dim-1, dim-1] = d


    for i in range(1, dim-1):

        A[i, i-1] = a
        A[i, i] = d
        A[i, i+1] = a


    return A


def calc_buckling_exact_eigenvalues(n):
    """analytic expression for eigenvalues"""
    r_max = 1
    step = r_max / (dim+1)

    a = -1 / step**2
    d = 2 / step**2
    print(d)

    lambdas = np.zeros(n)
    for j in range(1, n - 1):
        lambdas[j] = d + 2 * a * np.cos((j * np.pi) / (n + 1))

    return lambdas


dim = 20
A = buckling_beam_tridiagonal(dim)
#print(A)

B, p = jacobi(A)
print(np.sort(B))

exact_eigenvalues = calc_buckling_exact_eigenvalues(dim)
print(exact_eigenvalues)
