__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

from MicroDosimetry import Measurement
from lib.Spectrum import probability_function, probability_density, weighted_probability_density
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.fft import rfft, rfftfreq, irfft
from sklearn.preprocessing import normalize

import os

ressources_path = 'ressources/stopping_power_tables'


""" Module for more experimental, advanced analysis ideas not part of the standard workflow.
    Take this with a grain of salt, there's some of half-baked stuff in here. """

def remove_noise(measurement, background_measurement):
    """ Remove noise from a measurement by subtracting a background measurement (scaled by the relative live time)
        Should be used directly on the linear spectrum f(Channel)"""
    
    if measurement.y_axis != 'COUNTS' or background_measurement.y_axis != 'COUNTS':
        raise ValueError('Make sure to use the original linear spectrum f(Channel)')
    
    if measurement.live_time == 0.0 or background_measurement.live_time == 0.0:
        raise ValueError('There is no live time information available')
    
    if measurement.real_time == 0.0 or background_measurement.real_time == 0.0:
        raise ValueError('There is no real time information available')

    scaling_factor = measurement.real_time / background_measurement.real_time

    measurement.data[measurement.y_axis] = measurement.data[measurement.y_axis] - scaling_factor * background_measurement.data[background_measurement.y_axis]
    measurement.data.loc[measurement.data[measurement.y_axis] < 0, measurement.y_axis] = 0

    print(f'Noise removed from spectrum: {measurement.name}')
    print(f'Total area removed: {np.sum(scaling_factor * background_measurement.data[background_measurement.y_axis])}')

def material_conversion(measurement, from_material, to_material, database='SRIM'):
    """ Convert a spectrum from one material to another using a binwise comparison of stopping power tables 
    
    Inspired by Magrin et al. 2018
    """

    densities = {'silicon':2.32, 'diamond': 3.52, 'sic': 3.21, 'water': 1} #g/cm3

    if measurement.x_axis not in ['ENERGY', 'LINEAL ENERGY'] or measurement.y_axis != 'COUNTS':
        raise ValueError('This calculation only works for energy calibrated frequency spectra: f(E), f(y)')
    
    if from_material == to_material:
        raise ValueError('The materials for the conversion are the same')
    
    x = measurement.data[measurement.x_axis].to_numpy()
    y = measurement.data[measurement.y_axis].to_numpy()

    # Convert lineal energy to energy
    if measurement.x_axis == 'LINEAL ENERGY':

        if measurement.mean_chord_length == 0.0:
            raise ValueError('The chord length has to be set first')
        x = x / measurement.mean_chord_length

    if not os.path.exists(f"{ressources_path}/{database}_{measurement.particle}_{from_material}.csv"):
            raise ValueError(f"Stopping power table for {measurement.particle} in {from_material} not found in {ressources_path}")
    if not os.path.exists(f"{ressources_path}/{database}_{measurement.particle}_{to_material}.csv"):
            raise ValueError(f"Stopping power table for {measurement.particle} in {to_material} not found in {ressources_path}")

    table_from = pd.read_csv(f"{ressources_path}/{database}_{measurement.particle}_{from_material}.csv",
                            sep=',', skiprows=5, names=['energy', 'range', 'dE/dx'])
    table_to = pd.read_csv(f"{ressources_path}/{database}_{measurement.particle}_{to_material}.csv",
                            sep=',', skiprows=5, names=['energy', 'range', 'dE/dx'])
    
    density_scaling_factor = densities[to_material] / densities[from_material]

    # Linear interpolation
    precision_lower = 0.01
    precision_upper = 1

    table_from = pd.concat([pd.DataFrame([[0, 0, 0]], columns=['energy', 'range', 'dE/dx']), table_from], ignore_index=True)
    table_to = pd.concat([pd.DataFrame([[0, 0, 0]], columns=['energy', 'range', 'dE/dx']), table_to], ignore_index=True)

    new_energy_values_from = np.arange(0.0, 100, precision_lower)
    new_energy_values_from = np.append(new_energy_values_from, np.arange(100, 1000, precision_upper))
    new_table_from = pd.DataFrame(new_energy_values_from, columns=['energy'])
    merged_table_from = pd.merge(new_table_from, table_from, on='energy', how='outer', sort=True)
    interpolated_table_from = merged_table_from.interpolate(method='linear')
    interpolated_table_from.sort_values(by='energy', inplace=True)

    new_energy_values_to = np.arange(0.0, 100, precision_lower)
    new_energy_values_to = np.append(new_energy_values_to, np.arange(100, 1000, precision_upper))
    new_table_to = pd.DataFrame(new_energy_values_to, columns=['energy'])
    merged_table_to = pd.merge(new_table_to, table_to, on='energy', how='outer', sort=True)
    interpolated_table_to = merged_table_to.interpolate(method='linear')
    interpolated_table_to.sort_values(by='energy', inplace=True)

    table_from = interpolated_table_from
    table_to = interpolated_table_to
        
    # Binwise comparison:
    # 1) Find the closest dE/dx in table_from
    # 2) Take the value for the corresponding energy in table_to
    y_new = np.zeros_like(y)

    max_idx_from = np.argmax(table_from['dE/dx'])

    for i in range(len(x)):
        #TODO: Unsure if this is set up correctly

        # Find the closest energy loss value in the table (2 closest values = both sides of the curve)
        idx_from = np.argpartition(np.abs(table_from['dE/dx'] - x[i]), 2)[:2]

        if idx_from[0] > max_idx_from:
            idx_from = idx_from[1]
        else:
            idx_from = idx_from[0]

        # Take energy loss value of the other material for the same energy
        idx_to = np.argmin(np.abs(table_to['energy'] - table_from['energy'][idx_from]))

        """ if i > 1000:
            plt.plot(table_from['energy'], table_from['dE/dx'], label=f'{from_material}', marker='x')
            plt.plot(table_to['energy'], table_to['dE/dx'], label=f'{to_material}', marker='x')
            plt.title(f'Event {i} with energy loss {x[i]}')
            plt.axhline(y=x[i], color='r', linestyle='--', label=f'Spectrum energy loss {x[i]}')
            plt.axvline(x=table_from['energy'][idx_from], color='g', linestyle='--', label=f'idx_from-{idx_from}')
            plt.axhline(y=table_from['dE/dx'][idx_from], color='g', linestyle='--', label=f'dE/dx_from-{table_from["dE/dx"][idx_from]}')
            plt.axvline(x=table_to['energy'][idx_to], color='b', linestyle='--', label=f'idx_to-{idx_to}')
            plt.axhline(y=table_to['dE/dx'][idx_to], color='b', linestyle='--', label=f'dE/dx_to-{table_to["dE/dx"][idx_to]}')
            plt.xlim(0, 5000)
            plt.legend()
            plt.show() """
            
        # Scale the counts accordingly
        y_new[i] = y[i] * table_to['dE/dx'][idx_to] / table_from['dE/dx'][idx_from] * density_scaling_factor
    
    measurement.data.drop(columns=measurement.x_axis, inplace=True)
    measurement.data.drop(columns=measurement.y_axis, inplace=True)

    measurement.data['ENERGY'] = x
    measurement.x_axis = 'ENERGY'
    measurement.data[measurement.y_axis] = y_new

    print(f'Spectrum {measurement.name} converted from {from_material} to {to_material}')


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