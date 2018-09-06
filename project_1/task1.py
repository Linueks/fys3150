from __future__ import division, print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
from scipy.linalg import lu_solve, lu_factor
style.use('ggplot')


def source_term(x, n):
    """
    Source term for the right hand side of Av=b
    Using special case that has closed form solution for
    one-dimensional Poisson equation
    """
    h = 1/(n + 1)
    value = h**2 * 100 * np.exp(-10 * x)

    return value


def u(x):
    """
    Closed form solution for 1D Poisson Equation
    given our source term: f(x) = 100e^{-10x}
    """

    value = 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)
    return value


def general_tridiagonal(x, n):
    """
    General algorithm to solve tridiagonal matrix equation Av=b
    Equation is solved by row reducing and back substitution -> *Av'=b'
    """


    main_diag = np.array([2 for i in range(n)])
    row_reduced_main_diag = np.zeros(n)
    row_reduced_main_diag[0] = main_diag[0]


    lower_diag = np.array([-1 for i in range(n-1)])
    upper_diag = np.array([-1 for i in range(n-1)])


    right_side = source_term(x, n)
    row_reduced_right_side = np.zeros(n)
    row_reduced_right_side[0] = right_side[0]


    solution = np.zeros(n)


    #Forward substitution
    for i in range(1, n - 1):
        row_reduced_main_diag[i] = main_diag[i] - (lower_diag[i-1] / row_reduced_main_diag[i-1]) * upper_diag[i-1]
        row_reduced_right_side[i] = right_side[i] - (lower_diag[i-1] / row_reduced_main_diag[i-1]) * row_reduced_right_side[i-1]


    #Backward substitution
    for i in range(n - 2, 0, -1):
        solution[i] = (row_reduced_right_side[i] - upper_diag[i] * solution[i+1]) / row_reduced_main_diag[i]


    return solution


def relative_error_calc(x, n, numerical_solution):


    analytic_solution = u(x)
    error = np.log10(np.abs((numerical_solution[1:] - analytic_solution[1:]) / analytic_solution[1:]))
    max_error = np.max(error)


    return max_error


def lu_decomposition(x, n):
    """Generating matrix and running scipy lu decomposition to solve system of equations"""


    matrix = np.zeros((n, n))
    main_diag = np.array([2 for i in range(n)])
    lower_diag = np.array([-1 for i in range(n-1)])
    upper_diag = np.array([-1 for i in range(n-1)])
    right_side = source_term(x, n)

    a = -1
    b = 2
    c = -1

    for i in xrange(n):

        matrix[i][i] = b

        if i != 0:
            matrix[i][i-1] = a
        if i != n - 1:
            matrix[i][i+1] = c


    solution = (lu_solve(lu_factor(matrix), right_side))

    return solution


if __name__=='__main__':


    if sys.argv[1] == 'general':


        for exp in xrange(1, 4):
            n = 10**exp
            h = 1 / (n + 1)
            x = np.array([i*h for i in range(n)])

            numerical_solution = general_tridiagonal(x, n)
            relative_error = relative_error_calc(x, n, numerical_solution)


            plt.plot(x, numerical_solution, linestyle='dashed', label='Numerical solution, n:%i' % n)


        plt.plot(x, u(x), label='Analytic solution')
        plt.legend()
        plt.show()


    elif sys.argv[1] == 'lu':


        for exp in xrange(1, 4):
            n = 10**exp
            h = 1 / (n + 1)
            x = np.array([i*h for i in range(n)])

            lu_solution = lu_decomposition(x, n)

            plt.plot(x, lu_solution, label='LU decomposition, n:%i' % n)


        plt.plot(x, u(x), label='Analytic solution')
        plt.legend()
        plt.show()
