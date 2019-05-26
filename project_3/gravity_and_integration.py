from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const


def gravitational_acceleration(position_vector):


    solar_mass = 2 * 10**30     # [kg]
    earth_mass = 6 * 10**24 / solar_mass    # [M_S, solar masses]
    gravitational_force = 4 * np.pi**2 * earth_mass / (np.linalg.norm(position_vector))**2     # [AU**3 * M_S / (yr**2 * m**2), newtons]
    gravitational_accel = gravitational_force / earth_mass
    r_unit = -position_vector / np.linalg.norm(position_vector)


    return gravitational_accel * r_unit


def forward_euler(position, velocity, acceleration, n):


    for i in range(n - 1):
        acceleration[i, :] = gravitational_acceleration(position[i, :])
        velocity[i + 1, :] = velocity[i, :] + acceleration[i, :] * time_step
        position[i + 1, :] = position[i, :] + velocity[i, :] * time_step

    return position, velocity, acceleration


def velocity_verlet(position, velocity, acceleration, n):


    acceleration[0, :] = gravitational_acceleration(position[0, :])

    for i in range(n - 1):

        position[i + 1, :] = position[i, :] + velocity[i, :] * time_step + 0.5 * acceleration[i, :] * time_step**2
        acceleration[i + 1, :] = gravitational_acceleration(position[i + 1, :])
        velocity[i + 1, :] = velocity[i, :] + 0.5 * (acceleration[i + 1, :] + acceleration[i, :]) * time_step

    return position, velocity, acceleration


number_of_points = 10000
total_time = 1      # [years]
time_step = total_time / number_of_points


time = np.linspace(0, total_time, number_of_points)
acceleration = np.zeros((number_of_points, 2))
velocity = np.zeros((number_of_points, 2))
position = np.zeros((number_of_points, 2))

position[0, :] = 1, 0
velocity[0, :] = 0, 2 * np.pi

pos, vel, accel = forward_euler(position, velocity, acceleration, number_of_points)
verlet_pos, verlet_vel, verlet_accel = velocity_verlet(position, velocity, acceleration, number_of_points)


plt.plot(verlet_pos[:, 0], verlet_pos[:, 1])
plt.show()
