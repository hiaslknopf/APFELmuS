import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Calibrate, Output

""" Test the edge calibration with a real measurement. """

edge_file = 'exp_carbon_edge'
example_file = 'exp_carbon_data_1'

#----------------------------------------------------

#Read in
campaign1 = MicroDosimetry()
campaign1.read_file(f'tests/testdata/{edge_file}.MCA')
campaign1.read_file(f'tests/testdata/{example_file}.MCA')

#Calibrate
Calibrate.get_chord_length(campaign1.measurements[edge_file], 'slab', 10)
Calibrate.get_chord_length(campaign1.measurements[example_file], 'slab', 10)

edge_pos_mV, _ = Calibrate.get_edge_pos(campaign1.measurements[edge_file], marker_point='hTC', fit_bounds=[1500,2500], check_plot=True)
ymax, Lmax, _ = Calibrate.get_stopping_power(campaign1.measurements[edge_file], 'mean', 'SRIM', 0.01, plot=True)

#Scale spectrum
Calibrate.scale_energy_axis(campaign1.measurements[example_file], edge_pos_mV, ymax, 'lineal', 'mean')
Output.plot_single(campaign1.measurements[example_file])
#----------------------------------------------------