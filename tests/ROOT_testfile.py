import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Calibrate, Spectrum, Output

""" Evaluate simulation data (ROOT format) """

file = 'exp_148.2MeV_proton_si_10um'

#----------------------------------------------------

#Set up a new measurement campaign
campaign1 = MicroDosimetry()
campaign1.read_file(f'tests/testdata/{file}.root')

#Add metainformation to the ROOT data
campaign1.attach_info(campaign1.measurements[file], info_dict = {'particle': 'proton', 'detector': 'silicon', 'num_channels': 1024})
Calibrate.get_chord_length(campaign1.measurements[file], 'slab', 10, plot=False)

#print(campaign1.measurements['si_proton_148.2MeV_11.299cm_0.4um'].data)
Output.plot_single(campaign1.measurements[file], name='test')
Output.csv_output(campaign1.measurements[file], output_path='tests/testdata', name=file)

#----------------------------------------------------
Calibrate.lineal_energy_axis(campaign1.measurements[file])
Output.plot_single(campaign1.measurements[file], name='root')
Spectrum.logarithmic_binning(campaign1.measurements[file], 80)
Spectrum.probability_function(campaign1.measurements[file], 'F')
Spectrum.probability_density(campaign1.measurements[file])

y_F = Output.get_mean_value(campaign1.measurements[file])
Output.plot_single(campaign1.measurements[file], mean=y_F)

#----------------------------------------------------