__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

from MicroDosimetry import Measurement
from lib.Spectrum import probability_function, probability_density, weighted_probability_density
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import rfft, rfftfreq, irfft
from sklearn.preprocessing import normalize

""" Module for more experimental, advanced analysis methods not part of the standard workflow.
    Take this with a grain of salt, there's some of half-baked stuff in here. """

def remove_noise(measurement, background_measurement):
    """ Remove noise from a measurement by subtracting a background measurement (scaled by the relative live time) """
    raise NotImplementedError

def material_conversion(measurement, from_material, to_material, density=1.0, tables='SRIM'):
    """ Convert a spectrum from one material to another using a binwise comparison of stopping power tables """
    raise NotImplementedError

def get_LET_distribution(measurement, chord_length_dist, LET_spec='pdf'):
    """ Calculate an LET distribution from a calibrated f(E) pdf
    Following an algorithm according to Kellerer 1972, Phys.Med.Biol., 17, NO.2, 232-240
    
    Args:
        measurement: The spectrum to be analyzed
        chord_length_dict: Chord length distribution information
        LET_spec: 'pdf' for t(LET) or 'weighted' for LET*t(LET)
    """

    testplot = False

    #TODO: Check if chord_length is already attached

    if measurement.x_axis != 'ENERGY' or measurement.y_axis not in ['COUNTS', 'F(E)', 'f(E)', 'Ef(E)']:
        raise ValueError('This calculation only works for energy calibrated spectra')
    if measurement.mean_chord_length == 0.0:
        raise ValueError('The chord length function has to be called first')

    #Transform pdf data into right format Ef(E)
    if measurement.y_axis == 'COUNTS':
        probability_function(measurement, 'F')
        probability_density(measurement)
        weighted_probability_density(measurement)
    elif measurement.y_axis == 'F(E)':
        probability_density(measurement)
        weighted_probability_density(measurement)
    elif measurement.y_axis == 'f(E)':
        weighted_probability_density(measurement)
    
    #Get pdf data
    E = measurement.data[measurement.x_axis].to_numpy()
    phi = measurement.data[measurement.y_axis].to_numpy()

    #Normalize distributions for Fourier Transform
    gamma = np.multiply(chord_length_dist['pdf'], chord_length_dist['chord_len'])
    norm_gamma = normalize(gamma[:,np.newaxis], axis=0).ravel()
    norm_phi = normalize(phi[:,np.newaxis], axis=0).ravel()

    if testplot:
        fig, ax = plt.subplots()
        ax.plot(chord_length_dist['chord_len'], norm_gamma, label='chord length dist g(l)')
        ax.set_xlabel('Chord length [µm]')
        ax.legend()
        plt.show()

        fig, ax = plt.subplots()
        ax.plot(E, norm_phi, label='normalize phi = Ef(E)')
        ax.set_xlabel('Energy [keV]')
        ax.legend()
        plt.show()

    #Fourier transform
    phi_ft = rfft(norm_phi)
    gamma_ft = rfft(norm_gamma)
    rfftx = rfftfreq(measurement.num_channels, 1 / 2*measurement.num_channels)

    if testplot:
        fig, ax = plt.subplots()
        ax.plot(rfftx, np.abs(phi_ft), label='phi*')
        ax.plot(rfftx, np.abs(gamma_ft), label='gamma*')
        ax.set_xlabel('Discrete step [a.u.]')
        ax.legend()
        plt.show()

    #Calculate tau
    tau_ft = np.divide(phi_ft, gamma_ft, out=np.zeros_like(phi_ft), where=gamma_ft!=0)

    if testplot:
        fig, ax = plt.subplots()
        ax.plot(rfftx, np.abs(tau_ft), label='tau*')
        ax.set_xlabel('Discrete step [a.u.]')
        ax.legend()
        plt.show()

    #Backtrafo
    tau = irfft(tau_ft)
    LET_pdf = np.divide(tau, chord_length_dist['chord_len'], out=np.zeros_like(tau), where=chord_length_dist['chord_len']!=0)

    #Get rid of the negative values
    LET_pdf[LET_pdf < 0] = 0
    tau[tau < 0] = 0

    if testplot:
        fig, ax = plt.subplots()
        ax.plot(LET_pdf, label='LET_pdf = t(LET)')
        ax.plot(tau, label='tau = LET*t(LET)')
        ax.set_xlabel('Discrete step [a.u.]')
        ax.legend()
        plt.show()

    measurement.data.drop(columns=measurement.x_axis, inplace=True)
    measurement.data.drop(columns=measurement.y_axis, inplace=True)

    #TODO: Something is not right about the x-axis
    #TODO: Kellerer mentions this: Intervall of periodiciallity has to be chosen accordingly, that all spectra fit
    measurement.data[measurement.x_axis] = np.divide(E, chord_length_dist['chord_len'], out=np.zeros_like(E), where=chord_length_dist['chord_len']!=0)
    measurement.x_axis = 'LET'

    if LET_spec == 'pdf':
        measurement.data[measurement.y_axis] = LET_pdf
        measurement.y_axis = 't(LET)'
    elif LET_spec == 'weighted':
        measurement.data[measurement.y_axis] = tau
        measurement.y_axis = 'LETt(LET)'
    
    print(f'LET distribution calculated from energy calibrated spectrum: {measurement.name}')