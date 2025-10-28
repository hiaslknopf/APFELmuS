import sys

# Whereever you stored the repository project folder
projekt_dir = 'C:/Users/knopf/Desktop/APFELmuS'

# Linux version
# projekt_dir = '/home/knopf/Desktop/APFELmuS'

if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Calibrate

""" This is a template demonstrating how to find calibration factors using the edge technique.
    The used MCA spectrum has to contain a clear edge """


file = ''

#----------------------------------------------------

#Read in
campaign1 = MicroDosimetry()
campaign1.read_file(f'analysis/{file}.MCA')

#Calibrate
Calibrate.get_chord_length(campaign1.measurements[file], 'slab', 10)
edge_pos_mV, _ = Calibrate.get_edge_pos(campaign1.measurements[file], marker_point='hTC', check_plot=True)
ymax, Lmax, _ = Calibrate.get_stopping_power(campaign1.measurements[file], 'mean', 'SRIM', 0.01, plot=True)

print(f'Edge position: {edge_pos_mV} mV')
print(f'Maximum energy loss: y = {ymax} keV and L = {Lmax} um')
print(f'Scaling factor: {ymax/edge_pos_mV} keVum-1mV-1')

#----------------------------------------------------