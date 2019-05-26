import numpy as np
import matplotlib.pyplot as plt


L_40 = np.load('L40.npy')
L_60 = np.load('L60.npy')
L_80 = np.load('L80.npy')
L_100 = np.load('L100.npy')
temperatures = np.linspace(2.0, 2.3, 60)
energy, magnetization, heat_capacity, susceptibility = L_40
energy2, magnetization2, heat_capacity2, susceptibility2 = L_60
energy3, magnetization3, heat_capacity3, susceptibility3 = L_80
energy4, magnetization4, heat_capacity4, susceptibility4 = L_100

results1 = np.load('1-ordered.npy')
results2 = np.load('2.4-ordered.npy')
results3 = np.load('1-unordered.npy')
results4 = np.load('2.4-unordered.npy')

#print(np.shape(results1))




plt.subplot(2, 2, 1)
plt.plot(temperatures, energy / 40**2,
         temperatures, energy2 / 60**2,
         temperatures, energy3 / 80**2,
         temperatures, energy4 / 100**2)
plt.ylabel('Mean Energy')
plt.xlabel('Temperature')
plt.legend(['$L=40$', '$L=60$', '$L=80$', '$L=100$'], loc='upper right')
plt.subplot(2, 2, 2)
plt.plot(temperatures, magnetization / 40**2,
         temperatures, magnetization2 / 60**2,
         temperatures, magnetization3 / 80**2,
         temperatures, magnetization4 / 100**2)
plt.legend(['$L=40$', '$L=60$', '$L=80$', '$L=100$'], loc='upper right')
plt.ylabel('Mean Magnetization')
plt.xlabel('Temperature')
plt.subplot(2, 2, 3)
plt.plot(temperatures, heat_capacity / 40**2,
         temperatures, heat_capacity2 / 60**2,
         temperatures, heat_capacity3 / 80**2,
         temperatures, heat_capacity4 / 100**2)
plt.legend(['$L=40$', '$L=60$', '$L=80$', '$L=100$'], loc='upper right')
plt.ylabel('Heat Capacity')
plt.xlabel('Temperature')
plt.subplot(2, 2, 4)
plt.plot(temperatures, susceptibility / 40**2,
         temperatures, susceptibility2 / 60**2,
         temperatures, susceptibility3 / 80**2,
         temperatures, susceptibility4 / 100**2)
plt.legend(['$L=40$', '$L=60$', '$L=80$', '$L=100$'], loc='upper right')
plt.ylabel('Susceptibility')
plt.xlabel('Temperature')
plt.show()


#plt.plot(temperatures, energy4 / 100**2,
#         temperatures, magnetization4 / 100**2,
#         temperatures, heat_capacity4 / 100**2,
#         temperatures, susceptibility4 / 100**2)
#plt.show()

"""
plt.subplot(2, 2, 1)
plt.title('$L=40$')
plt.plot(temperatures, energy,
            temperatures, magnetization,
            temperatures, heat_capacity,
            temperatures, susceptibility)
plt.legend(['energy', 'magnetization', 'heat capacity', 'susceptibility'])
plt.subplot(2, 2, 2)
plt.title('$L=60$')
plt.plot(temperatures, energy2,
            temperatures, magnetization2,
            temperatures, heat_capacity2,
            temperatures, susceptibility2)
plt.legend(['energy', 'magnetization', 'heat capacity', 'susceptibility'])\
plt.subplot(2, 2, 3)
plt.title('$L=80$')
plt.plot(temperatures, energy3,
            temperatures, magnetization3,
            temperatures, heat_capacity3,
            temperatures, susceptibility3)
plt.legend(['energy', 'magnetization', 'heat capacity', 'susceptibility'])
plt.subplot(2, 2, 4)
plt.title('$L=100$')
plt.plot(temperatures, energy4,
            temperatures, magnetization4,
            temperatures, heat_capacity4,
            temperatures, susceptibility4)
plt.legend(['energy', 'magnetization', 'heat capacity', 'susceptibility'])
plt.show()
"""


#print(np.shape(results))
#print(results)
