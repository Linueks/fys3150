from __future__ import division, print_function
import numpy as np


# finding the largest number above the diagonal
def find_a_max(A, N):
    """
    Returning the value and index of the largest matrix element above the diagonal.
    Requires symmetric matrix
    """
    a_max = 0
    for i in range(N):
        for j, value in enumerate(A[i, i + 1:]):
            if abs(value) >= a_max:
                a_max = abs(value)
                k, l = i, j + i + 1
                #print(value**2)
    return a_max, l, k




def rotate(A, S, k, l, n):
    if(A[k, l] != 0):
        tau = (A[l, l] - A[k, k]) / (2 * A[k, l])
        if(tau >= 0):
            t = 1 / (tau + np.sqrt(1 + tau**2))
        else:
            t = -1 / (-tau + np.sqrt(1 + tau**2))

        c = 1 / np.sqrt(1 + t**2)
        s = c * t
    else:
        c = 1
        s = 0

    a_kk = A[k, k]
    a_ll = A[l, l]

    A[k, k] = c**2 * a_kk - 2 * c * s * A[k, l] + s**2 * a_ll
    A[l, l] = s**2 * a_kk + 2 * s * c * A[k, l] + c**2 * a_ll
    A[k, l] = 0
    A[l, k] = 0

    for i in range(n):
        if(i != k and i != l):
            a_ik = A[i, k]
            a_il = A[i, l]
            A[i, k] = c * a_ik - s * a_il
            A[k, i] = A[i, k]
            A[i, l] = c * a_il + s * a_ik
            A[l, i] = A[i, l]

        r_ik = S[i, k]
        r_il = S[i, l]
        S[i, k] = c*r_ik - s*r_il
        S[i, l] = c*r_il + s*r_ik


    return A, S


def jacobi_method(A, n):
    """Function to iterate the Jacobi rotation"""

    A = A.copy()
    S = np.identity(n)

    epsilon = 1e-8
    max_iterations = n**3
    iterations = 0
    max_val, l, k = find_a_max(A, n)
    while(np.abs(max_val**2) > epsilon and iterations < max_iterations):
        max_val, l, k = find_a_max(A, n)
        A, S = rotate(A, S, k, l, n)
        iterations += 1


    return iterations, np.sort(np.diagonal(A))



def calc_buckling_exact_eigenvalues(n):
    """analytic expression for eigenvalues"""
    r_max = 1
    step = r_max / (dim+1)

    a = -1 / step**2
    d = 2 / step**2

    lambdas = np.zeros(n)
    for j in range(1, n - 1):
        lambdas[j] = d + 2 * a * np.cos((j * np.pi) / (n + 1))

    return lambdas


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


def single_electron_tridiagonal(dim):
    """Function for initializing a tridiagonal (dim x dim) matrix
    for the buckling beam problem"""


    r_max = 5
    step = r_max / (dim + 1)

    A = np.zeros((dim, dim))
    V = np.zeros(dim)

    for i in range(1, dim-1):
        a = -1 / step ** 2
        d = 2 / step ** 2 + ((i + 1) * step)**2
        A[i, i-1] = a
        A[i, i] = d
        A[i, i+1] = a

    A[0, 0] = d
    A[0, 1] = a

    A[dim-1, dim-2] = a
    A[dim-1, dim-1] = d

    return A


def compare_iter_dim():
    """Looking at how the number of similarity transformations varies with dimensionality"""


    dimensions = np.arange(200, 210, 10)

    for dim in dimensions:
        A = buckling_beam_tridiagonal(dim)
        iterations, B = jacobi_method(A, dim)
        print(dim, iterations)


#compare_iter_dim()


dim = 50
#analytic_eigenvals = calc_buckling_exact_eigenvalues(dim)
#print(analytic_eigenvals)
#A = buckling_beam_tridiagonal(dim)
A = single_electron_tridiagonal(dim)
#print(A)
iterations, B = jacobi_method(A, dim)
print(B)
eigenvalues, eigenvectors = np.linalg.eig(A)
permute = eigenvalues.argsort()
eigenvalues = eigenvalues[permute]
print(eigenvalues)
