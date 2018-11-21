from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt


class SolarSystem():

    def __init__(self, total_time, time_step):
        self.solver = Solver(total_time, time_step)
        self.planet = Planet(self.solver, 6e24, (1, 0), (0, 2 * np.pi))
        self.sun = Planet(self.solver, 2e30, (0, 0), (0, 0))
        self.planet.accelerations[0, :] = self.planet.gravitational_acceleration([self.sun], 0)

    def run(self):
        self.solver.velocity_verlet(self.planet, self.sun)
        plt.plot(self.planet.positions[:, 0], self.planet.positions[:, 1])
        plt.show()


class Planet():

    def __init__(self, solver, mass, init_pos, init_vel):
        self.positions = np.zeros((solver.number_of_points, 2))
        self.velocities = np.zeros((solver.number_of_points, 2))
        self.accelerations = np.zeros((solver.number_of_points, 2))
        self.positions[0, :] = init_pos
        self.velocities[0, :] = init_vel


    def gravitational_acceleration(self, planets, t):

        for planet in planets:
            if planet != self:
                position_vector = planet.positions[t, :] - self.positions[t, :]
                gravitational_accel = 4 * np.pi**2  / np.linalg.norm(position_vector)**2
                r_unit = position_vector / np.linalg.norm(position_vector)

        return gravitational_accel * r_unit


class Solver():

    def __init__(self, total_time, time_step):

        self.total_time = total_time
        self.time_step = time_step
        self.number_of_points = int(total_time / time_step)


    def forward_euler(self, planet, bodies):
        for t in range(self.number_of_points - 1):
            planet.accelerations[t, :] = planet.gravitational_acceleration([sun], t)
            planet.velocities[t + 1, :] = planet.velocities[t, :] + planet.accelerations[t, :] * self.time_step
            planet.positions[t + 1, :] = planet.positions[t, :] + planet.velocities[t, :] * self.time_step


    def velocity_verlet(self, planet, sun):
        for t in range(self.number_of_points - 1):
            planet.positions[t + 1, :] = planet.positions[t, :] + planet.velocities[t, :] * self.time_step + 0.5 * planet.accelerations[t, :] * self.time_step**2
            planet.accelerations[t + 1, :] = planet.gravitational_acceleration([sun], t + 1)
            planet.velocities[t + 1, :] = planet.velocities[t, :] + 0.5 * (planet.accelerations[t + 1, :] + planet.accelerations[t, :]) * self.time_step




solarsystem = SolarSystem(1, 1e-3)
solarsystem.run()
