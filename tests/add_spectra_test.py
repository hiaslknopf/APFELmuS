import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Output, Spectrum, Calibrate

""" Demonstration of the summing method for obtaining more distinct edges. """

#----------------------------------------------------

campaign1 = MicroDosimetry()
campaign1.read_file('tests/testdata/spectra/exp_carbon_data_1.MCA')
campaign1.read_file('tests/testdata/spectra/exp_carbon_data_2.MCA')

#----------------------------------------------------

for file in campaign1.measurements:
    Calibrate.scale_energy_axis_with_factor(campaign1.measurements[file], 0.8)
    Spectrum.probability_function(campaign1.measurements[file], 'D')
    Spectrum.probability_density(campaign1.measurements[file])
    Spectrum.weighted_probability_density(campaign1.measurements[file])

#----------------------------------------------------

Output.plot_campaign(campaign1, 'All spec')

#Summing the spectra
sum_list = [campaign1.measurements['exp_carbon_data_1'],
            campaign1.measurements['exp_carbon_data_2']]

sum_spec = Spectrum.add_spectra(sum_list)
campaign1.add_measurement('summed_spectrum', sum_spec)

Output.plot_campaign(campaign1, 'All spec')
Output.plot_single(campaign1.measurements['summed_spectrum'], 'Summed spectrum')