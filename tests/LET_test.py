import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Calibrate, Spectrum, Output, Advanced

""" Test the LET analysis (after Kellerer) """

plot = True
file = 'exp_carbon_data_1'

#----------------------------------------------------

#Set up a new measurement campaign
campaign1 = MicroDosimetry()
campaign1.read_folder('tests/testdata')

if plot:
    Output.plot_single(campaign1.measurements[file], 'MCA')

#----------------------------------------------------

#Scale the energy axis
chord_len, _, chord_len_dist = Calibrate.get_chord_length(campaign1.measurements[file], 'slab', 10)
Calibrate.scale_energy_axis_with_factor(campaign1.measurements[file], 3.71, 'energy')
if plot:
    Output.plot_single(campaign1.measurements[file])

Spectrum.cutoff(campaign1.measurements[file], channels=100)
if plot:
    Output.plot_single(campaign1.measurements[file])

#-----------------------------------------------------

#And have a look at the LET distribution
Advanced.get_LET_distribution(campaign1.measurements[file], chord_len_dist)
Output.plot_single(campaign1.measurements[file], 'LET distribution (Kellerer)')