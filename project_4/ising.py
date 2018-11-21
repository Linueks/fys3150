import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import math as math
import argparse
import multiprocessing as mp
import time
plt.style.use('ggplot')
np.random.seed(0)


dim = 40



@jit
def initial_lattice(dim, ground_state=False):
    if ground_state:
        state = np.zeros((dim, dim), dtype=np.int8) + 1
    else:
        state = np.random.choice([1, -1], size=(dim, dim))

    return state


def periodic_boundary(i, limit, add):
    return (i + limit + add) % limit


def calculate_magnetization(state):
    magnetization = np.sum(state)

    return magnetization


def calculate_energy(state):
    dim = len(state)
    energy = 0

    for j in range(dim):
        for i in range(dim):
            energy -= state[i, j] * (state[i, periodic_boundary(j, dim, -1)] + state[periodic_boundary(i, dim, -1), j])

    return energy


def metropolis(temperature):
    """
        Monte Carlo move using Metropolis algorithm
    Args:
        state (np.array): N * N array representing the spins of a collection of
                            particles
        temperature (float): temperature...
        number_of_sweeps (int): Number of Monte Carlo sweeps

    Returns:
        state (np.array): ensamble after n Monte Carlo sweeps
        energies (np.array): energies for the flips
    """

    state = initial_lattice(dim, ground_state=False)
    number_of_sweeps = 50000
    energy = calculate_energy(state)
    magnetization = calculate_magnetization(state)
    average_energy = 0
    average_magnetization = 0
    average_energy_squared = 0
    average_magnetization_squared = 0
    abs_magnetization = 0
    accepted_configs = 0
    energies = []
    magnetizations = []
    configs = []

    w = np.zeros(17,np.float64)
    for de in range(-8,9,4): #include +8
        w[de+8] = math.exp(-de/temperature)

    for cycle in range(1, number_of_sweeps):
        for i in range(dim):
            for j in range(dim):
                rand_x = np.random.randint(dim)
                rand_y = np.random.randint(dim)

                # Periodic Boundary
                change_in_energy = 2 * state[rand_x, rand_y]\
                    * (state[periodic_boundary(rand_x, dim, -1), rand_y] + state[rand_x, periodic_boundary(rand_y, dim, -1)]\
                    + state[periodic_boundary(rand_x, dim, 1), rand_y] + state[rand_x, periodic_boundary(rand_y, dim, 1)])

                if np.random.uniform() <= w[change_in_energy + 8] or change_in_energy <= 0:
                    #Accept!
                    state[rand_x, rand_y] *= -1
                    magnetization += 2 * state[rand_x, rand_y]
                    energy += change_in_energy
                    accepted_configs += 1

        average_energy += energy
        average_energy_squared += energy**2
        average_magnetization += magnetization
        average_magnetization_squared += magnetization**2
        abs_magnetization += np.abs(magnetization)
        #print(abs_magnetization)

        energies.append(average_energy / cycle)
        magnetizations.append(abs_magnetization / cycle)
        configs.append(accepted_configs)


        #if cycle % 100:
        #    print(cycle)

    average_energy /= number_of_sweeps
    average_energy_squared /= number_of_sweeps
    average_magnetization /= number_of_sweeps
    average_magnetization_squared /= number_of_sweeps
    abs_magnetization /= number_of_sweeps

    energy_variance = (average_energy_squared - average_energy**2)# / (dim**2 * temperature**2)
    magnetization_variance = (average_magnetization_squared - abs_magnetization**2)# / (dim**2 * temperature)

    heat_capacity = energy_variance / temperature**2
    susceptibility = magnetization_variance / temperature

    if number_of_sweeps % 50 == 0:
        print(number_of_sweeps)

    #return (energies, magnetizations, configs)
    return (average_energy, abs_magnetization, heat_capacity, susceptibility, accepted_configs)


def plot_values_vs_temp(energy, specific_heat, magnetization, susceptibility):
    figure = plt.figure()

    top_right = figure.add_subplot(2, 2, 1)
    plt.plot(temperatures, energy, c='royalblue')
    plt.xlabel('Temperature')
    plt.ylabel('Energy')

    top_left = figure.add_subplot(2, 2, 2)
    plt.plot(temperatures, magnetization, c='indianred')
    plt.xlabel("Temperature")
    plt.ylabel("Magnetization")

    bottom_right =  figure.add_subplot(2, 2, 3)
    plt.plot(temperatures, specific_heat, c='sandybrown')
    plt.xlabel("Temperature")
    plt.ylabel("Specific Heat ")

    bottom_left =  figure.add_subplot(2, 2, 4)
    plt.plot(temperatures, susceptibility, c='forestgreen')
    plt.xlabel("Temperature")
    plt.ylabel("Susceptibility")

    plt.show()


def plot_values_vs_cycles(energy_ordered, abs_magnetization_ordered):#, energy_random, abs_magnetization_random):
    figure = plt.figure()
    top_right = figure.add_subplot(2, 1, 1)
    plt.title('Ordered Spins')
    plt.plot(range(len(energy_ordered)), energy_ordered, color='royalblue')
    #plt.plot(range(len(energy_random)), energy_random, 'o', markersize=3)\
    plt.ylabel('Energy')

    top_left = figure.add_subplot(2, 1, 2)
    plt.plot(range(len(abs_magnetization_ordered)), abs_magnetization_ordered, color='indianred')
    #plt.plot(range(len(abs_magnetization_random)), abs_magnetization_random, 'o', markersize=3)
    plt.xlabel("Monte Carlo Cycles")
    plt.ylabel("Magnetization")


    plt.show()


def parallel_sim(cycles):
    start_time = time.time()
    p = mp.Pool()
    results = p.map(metropolis, cycles)
    p.close()
    p.join()
    print(time.time() - start_time)

    energies = np.transpose(results)[0]
    magnetization = np.transpose(results)[1]
    heat_capacity = np.transpose(results)[2]
    susceptibility = np.transpose(results)[3]
    accepted_configs = np.transpose(results)[4]

    np.save('results_temperatures40_30000.npy', results)


    #plot_values_vs_cycles(energies, magnetization, energies, magnetization)
    plot_values_vs_temp(energies, heat_capacity, magnetization, susceptibility)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('lattice_dimension',
                        type=int,
                        help='dimension of spin lattice')
    args = parser.parse_args()
    dim = args.lattice_dimension
    cycles = np.arange(1, 500)
    temperatures = np.arange(2.0, 2.3, 0.05)
    parallel_sim(temperatures)

    #energies, magnetization, heat_capacity, susceptibility, accepted_configs = metropolis(1)

    #plot_values_vs_temp(energies, heat_capacity, magnetization, susceptibility)


    #plot_values_vs_cycles(np.array(energies), np.array(magnetization))
    #plt.plot(range(len(accepted_configs)), accepted_configs)


    #for cycle in range(len(cycles)):
    #    energy[cycle], magnet[cycle], heat_capacity[cycle], susceptibility[cycle] = metropolis(state, 2.4, number_of_sweeps=cycles[cycle])
    #    energy_ordered[cycle], magnet_ordered[cycle], heat_capacity_ordered[cycle], susceptibility_ordered[cycle] = metropolis(state_ordered, 2.4, number_of_sweeps=cycles[cycle])
    #    if cycles[cycle] % 100 == 0:
    #        print(f'done with {cycles[cycle]}')

    #plot_values_vs_cycles(energy, magnet, energy_ordered, magnet_ordered)

    #temperatures = np.linspace(a, args.final_temp, args.temp_divisions)
    #for temp in range(args.temp_divisions):
    #    energy[temp], magnetization_abs[temp], magnetization[temp], heat_capacity[temp], susceptibility[temp] = metropolis(state, temperatures[temp], number_of_sweeps=10000)


    #plot_values_vs_temp(energy, heat_capacity, magnetization, susceptibility, magnetization_abs)
