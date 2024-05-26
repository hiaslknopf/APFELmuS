import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Output, FileTranslator
import sys
""" Translation of MAESTRO data to .MCA format and plotting of the data """

#----------------------------------------------------

file = 'exp_carbon_data_1'

#You have to attach a linearization to your measured spectra -> .MCA file, used by the framework
#Either manually, with a loop or use the GUI
input_MAESTRO = f'tests/testdata/original_data/{file}.Spe'
input_linearization = 'tests/testdata/original_data/silicon_linearization_med_gain.csv'
output_path = 'tests/testdata'
info_dict = {'DETECTOR': 'silicon', 'GAIN': 'None', 'PARTICLE': 'carbon'}

FileTranslator.translate_MAESTRO_file(input_MAESTRO, input_linearization, output_path, file, info_dict)

#----------------------------------------------------

#Set up a new measurement campaign
campaign1 = MicroDosimetry()
campaign1.read_file(f'tests/testdata/{file}.MCA')

#----------------------------------------------------
Output.plot_single(campaign1.measurements[file], file, output_path)