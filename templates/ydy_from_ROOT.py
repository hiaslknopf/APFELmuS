import sys

projekt_dir = 'C:/Users/knopf/Desktop/microdosimetry'
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Calibrate, Spectrum, Output

""" calculate ydy spectra from the ROOT files """

#----------------------------------------------------

campaign1 = MicroDosimetry()
campaign1.read_folder('root_u_data')

#----------------------------------------------------

for file in campaign1.measurements:

    campaign1.attach_info(campaign1.measurements[file], info_dict = {'particle': 'proton', 'detector': 'silicon', 'num_channels': 1024})
    Calibrate.get_chord_length(campaign1.measurements[file], 'slab', 10, plot=False)

    Calibrate.lineal_energy_axis(campaign1.measurements[file])
    Spectrum.probability_function(campaign1.measurements[file], 'D')
    Spectrum.probability_density(campaign1.measurements[file])
    Spectrum.weighted_probability_density(campaign1.measurements[file])
    Spectrum.logarithmic_binning(campaign1.measurements[file], 60)	
    Spectrum.normalize_spectrum(campaign1.measurements[file])