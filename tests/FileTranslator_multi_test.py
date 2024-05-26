import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Output, FileTranslator
import sys
""" Translation of several MAESTRO data to .MCA format """

#----------------------------------------------------

files = ['sandra_low', 'sandra_mid', 'sandra_high']
gains = ['LOW', 'MID', 'HIGH']

#You have to attach a linearization to your measured spectra -> .MCA file, used by the framework
#Either manually, with a loop or use the GUI

for i, file in enumerate(files):
    input_MAESTRO = f'tests/testdata/original_data/{file}.Spe'
    input_linearization = 'tests/testdata/original_data/exp_linearization_3gains.csv'
    output_path = 'tests/testdata'
    info_dict = {'DETECTOR': 'TEPC', 'GAIN': gains[i], 'PARTICLE': 'proton'}

    FileTranslator.translate_MAESTRO_file(input_MAESTRO, input_linearization, output_path, file, info_dict)

#----------------------------------------------------

#Set up a new measurement campaign
campaign1 = MicroDosimetry()
campaign1.read_file(f'tests/testdata/sandra_low.MCA')

#----------------------------------------------------
Output.plot_single(campaign1.measurements['sandra_low'])