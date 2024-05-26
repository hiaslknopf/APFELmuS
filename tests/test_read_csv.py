import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry

""" Test reading in an already analyzed csv file"""


file = 'exp_148.2MeV_proton_si_10um'
path = f'tests/testdata/{file}.csv'

campaign1 = MicroDosimetry()
campaign1.read_file(path)

print(campaign1.measurements[file].data, '\n')
print(campaign1.measurements[file].particle, '\n')
print(campaign1.measurements[file].original_data, '\n')
print(campaign1.measurements[file].y_axis) #For some reason this does not show up in the df header right away