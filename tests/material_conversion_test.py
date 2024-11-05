import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Output, Advanced, Calibrate

""" Convert between materials using stopping power tables """

#----------------------------------------------------

campaign1 = MicroDosimetry()
campaign1.read_file('tests/testdata/spectra/exp_carbon_data_1.MCA')

Calibrate.get_chord_length(campaign1.measurements['exp_carbon_data_1'], 'slab', 10)
Calibrate.scale_energy_axis_with_factor(campaign1.measurements['exp_carbon_data_1'], 3.71)

Output.plot_single(campaign1.measurements['exp_carbon_data_1'], 'Spectrum before conversion')

Advanced.material_conversion(campaign1.measurements['exp_carbon_data_1'], 'silicon', 'diamond', database='SRIM')

Output.plot_single(campaign1.measurements['exp_carbon_data_1'], 'Spectrum after conversion')

#----------------------------------------------------