import numpy as np
import matplotlib.pyplot as plt

results = np.load('results_temperatures40_30000.npy')

energy = np.transpose(results)[0]
magnet = np.transpose(results)[1]
heat = np.transpose(results)[2]
suscept = np.transpose(results)[3]
accepted_configs = np.transpose(results)[4]
temps = np.arange(2.0, 2.3, 0.05)


figure = plt.figure()
top_right = figure.add_subplot(2, 2, 1)
plt.plot(temps, energy, c='royalblue')
plt.xlabel('Temperature')
plt.ylabel('Energy')

top_left = figure.add_subplot(2, 2, 2)
plt.plot(temps, magnet, c='indianred')
plt.xlabel("Temperature")
plt.ylabel("Magnetization")

bottom_right =  figure.add_subplot(2, 2, 3)
plt.plot(temps, heat, c='sandybrown')
plt.xlabel("Temperature")
plt.ylabel("Specific Heat ")

bottom_left =  figure.add_subplot(2, 2, 4)
plt.plot(temps, suscept, c='forestgreen')
plt.xlabel("Temperature")
plt.ylabel("Susceptibility")

plt.show()
