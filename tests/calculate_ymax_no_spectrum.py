import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from lib import Calibrate

""" Just get the ymax value for a given chord length, material and particle type. """

#----------------------------------------------------

#The calibration functions can be used without any measurement
mean, max, _ =Calibrate.get_chord_length('None', 'sphere', 10, True)
print(f'mean chord length: {mean} um')
print(f'max chord length: {max} um\n')

ymax, Lmax, _ = Calibrate.get_stopping_power('None', mean, 'ICRU', 0.01, 'proton', 'silicon', plot=True)

print(f'ymax: {ymax} keV/um')
#----------------------------------------------------