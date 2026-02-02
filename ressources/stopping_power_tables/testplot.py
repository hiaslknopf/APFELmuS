import matplotlib.pyplot as plt
import numpy as np

def read_nist_astar(filename):
    energies = []
    stopping_powers = []
    ranges = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines[5:]:  # Skip header lines
            parts = line.strip().split(',')
            energies.append(float(parts[0]))
            ranges.append(float(parts[1]))
            stopping_powers.append(float(parts[2]))
    return np.array(energies), np.array(stopping_powers), np.array(ranges)

FILE1 = "NIST_helium_silicon.csv"
FILE2 = "NIST_helium_water.csv"

energies_si, stopping_powers_si, ranges_si = read_nist_astar(FILE1)
energies_water, stopping_powers_water, ranges_water = read_nist_astar(FILE2)

plt.figure(figsize=(6, 4))
plt.plot(energies_si, stopping_powers_si, label='Helium in Silicon', color='blue')
plt.plot(energies_water, stopping_powers_water, label='Helium in Water', color='orange')
plt.legend()
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Energy (keV)')
plt.ylabel('Stopping Power (keV/um * cm3/g)')
plt.tight_layout()
plt.grid()
plt.show()