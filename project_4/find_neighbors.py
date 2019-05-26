import numpy as np


def initial_lattice(dim):
    state = np.random.choice([1, -1], size(dim, dim))
    return state


def find_neighbors(state, i, j):
    dim = len(state)

    left = (i, j - 1)
    right = (i, (j + 1) % dim)
    top = (i - 1, j)
    bottom = ((i + 1) % dim, j)
