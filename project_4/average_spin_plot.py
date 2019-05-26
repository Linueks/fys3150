from ising import metropolis, initial_lattice
import matplotlib.pyplot as plt

dim = 20
state = initial_lattice(dim)
temperatures = [0.5, 2.27, 5.0]


for t in temperatures:
    state, energies, spins = metropolis(state, t)
    spins = spins / dim**2
    plt.plot(range(len(spins)), spins, label = f'T = {t}')


plt.legend(loc = 'best')
plt.xlabel('Number of Sweeps')
plt.ylabel('Average Spin')
plt.ylim(-1.2, 1.2)
plt.savefig('spin_plot.png')
plt.show()
