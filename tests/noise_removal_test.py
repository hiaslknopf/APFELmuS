import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Output, Advanced

""" Remove noise from a measurement by subtracting a background measurement """

#TODO: Better data, that actually fits together

#----------------------------------------------------

campaign1 = MicroDosimetry()
campaign1.read_file('tests/testdata/exp_carbon_data_1.MCA')
campaign1.read_file('tests/testdata/noise_background.MCA')

Output.plot_campaign(campaign1, 'Raw Spectra')

Advanced.remove_noise(campaign1.measurements['exp_carbon_data_1'], campaign1.measurements['noise_background'])

Output.plot_single(campaign1.measurements['exp_carbon_data_1'], 'Noise Removed')

#----------------------------------------------------