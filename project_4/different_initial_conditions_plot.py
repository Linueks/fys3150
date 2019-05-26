import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')


results1 = np.load('1-ordered.npy')
results2 = np.load('2.4-ordered.npy')
results3 = np.load('1-unordered.npy')
results4 = np.load('2.4-unordered.npy')

energies1 = results1[0]
energies2 = results2[0]
energies3 = results3[0]
energies4 = results4[0]

magnet1 = results1[1]
magnet2 = results2[1]
magnet3 = results3[1]
magnet4 = results4[1]


plt.plot(range(len(magnet1)), magnet1, 'o', c='royalblue', markersize=1)
plt.plot(range(len(magnet2)), magnet2, 'o', c='indianred', markersize=1)
plt.plot(range(len(magnet3)), magnet3, 'o', c='sandybrown', markersize=1)
plt.plot(range(len(magnet4)), magnet4, 'o', c='forestgreen', markersize=1)
plt.legend(['$T=2.4$, disordered', '$T=2.4$, ordered', '$T=1$, disordered', '$T=1$, ordered'], loc='upper right')
plt.xlabel('Monte Carlo Cycles')
plt.ylabel('Mean Magnetization')
plt.show()




"""
plt.plot(range(len(energies1)), energies1, 'o', c='royalblue', markersize=1)
plt.plot(range(len(energies2)), energies2, 'o', c='indianred', markersize=1)
plt.plot(range(len(energies3)), energies3, 'o', c='sandybrown', markersize=1)
plt.plot(range(len(energies4)), energies4, 'o', c='forestgreen', markersize=1)
plt.legend(['$T=2.4$, disordered', '$T=2.4$, ordered', '$T=1$, disordered', '$T=1$, ordered'])
plt.xlabel('Monte Carlo Cycles')
plt.ylabel('Mean Energy')
plt.show()
"""
