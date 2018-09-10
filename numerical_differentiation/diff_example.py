"""
Easy example to get back into coding,
program will numerically calculate the
derivative of the exponential function
"""
from __future__ import division, print_function
import numpy as np



def second_derivative(number_of_steps, x, initial_step, step_modifier):
    """
    Function supposed to loop over different values for the step size to
    see how the numerical accuracy follows from the step size
    """


    h = initial_step
    h_step = np.zeros(number_of_steps)
    derivative = np.zeros(number_of_steps)


    for step in xrange(number_of_steps):
        h_step[step] = 


second_derivative(400, 5, 0.1, 0.01)
