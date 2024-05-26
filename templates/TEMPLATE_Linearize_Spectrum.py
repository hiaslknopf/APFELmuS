import sys

projekt_dir = 'C:/Users/knopf/Desktop/microdosimetry'
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from lib import FileTranslator
from glob import glob

""" This is a template demonstrating how to attach one or more linearization measurements to MAESTRO .Spe files and save this combination as a .MCA file for further analysis. """

#----------------------------------------------------

input_folder = ''
input_calibration_file = ''

input_folder_2 = ''

# It is also advisable to add some meta information about the measurement to the file
meta_info_dict = {'DETECTOR': 'silicon', 'GAIN': 'None', 'PARTICLE': 'carbon'}

output_file_path = ''

#----------------------------------------------------

#for input in glob(input_folder + '/*.Spe'):
#    output_file_name = input.split('/')[-1]
#    FileTranslator.translate_MAESTRO_file(input, input_calibration_file, output_file_path, output_file_name, meta_info_dict)

for input in glob(input_folder_2 + '/*.Spe'):
    output_file_name = input.split('/')[-1]
    FileTranslator.translate_MAESTRO_file(input, input_calibration_file, output_file_path, output_file_name, meta_info_dict)