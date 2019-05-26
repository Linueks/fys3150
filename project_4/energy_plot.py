from ising import metropolis, initial_lattice
import matplotlib.pyplot as plt
import numpy as np

dim = 20
state = initial_lattice(dim)
state, energies, spins = metropolis(state, 1)

energy = 0
cumulative_energy = []
for e in energies:
    energy += np.abs(e)
    cumulative_energy.append(energy)

plt.plot(range(len(energies)), energies)
plt.show()
