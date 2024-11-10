__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

import numpy as np
import pandas as pd
from bisect import *
import math
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

""" Collection of functions concerning the calibration, geometry and the energy axis of spectra """

def _fermi_func(h, A, B, C):
      return A / (1+np.exp(B * (h-C)))

def _get_nearest_range(table, array, value):
    #Find closest range bigger than a given value in stopping power table

    lower = array[bisect_left(array, value) - 1]
    closest_ind = table.index[table['range'] == lower].tolist()[0]

    if value <= array[0]:
        closest_ind = 0
    
    closest_range = table.iloc[closest_ind]['range']
    after_closest_range = table.iloc[closest_ind+1]['range']
    closest_dEdx = table.iloc[closest_ind]['dE/dx']
    after_closest_dEdx = table.iloc[closest_ind+1]['dE/dx']
    closest_energy = table.iloc[closest_ind]['energy']
    after_closest_energy = table.iloc[closest_ind+1]['energy']

    return {'closest_range': closest_range, 'after_closest_range': after_closest_range,
            'closest_dEdx': closest_dEdx, 'after_closest_dEdx': after_closest_dEdx,
            'closest_energy': closest_energy, 'after_closest_energy': after_closest_energy}

def _plot_edge_calibration(dict, ax=None):

    h = dict['h']
    weighted_pdf = dict['weighted_pdf']
    popt = dict['popt']
    fit_bounds = dict['fit_bounds']
    hFlex = dict['hFlex']
    hDD = dict['hDD']
    hTC = dict['hTC']

    channel_axis = dict['channel_axis']

    label_F = f"Sigmoid fit"
    label_hFlex = f"hFlex"
    label_hDD = f"hDD"
    label_hTC = f"hTC"
    label_bounds =f"fit bounds"

    print('hFlex:', hFlex, 'hDD:', hDD, 'hTC:', hTC)

    if ax is None:
        fig, ax = plt.subplots()

    ax.set_xlabel('mV')
    ax.set_ylabel('hd(h)')

    # Add a second x-axis showing the channels
    # This does not start at 0 because of the offset
    ax2 = ax.twiny()
    #ax2.set_xticks(channel_axis[::100])
    #ax2.set_xticklabels(channel_axis[::100])
    ax2.set_xlabel('channel')
    #ax2.set_xscale('log')

    # Add zero to beginning of channel axis and delete last element
    channel_axis = np.insert(channel_axis, 0, 0)
    channel_axis = np.delete(channel_axis, -1)
    
    ax2.set_xlim(min(channel_axis), max(channel_axis)+1, auto=False)
    ax2.grid()

    # Save the output as csv
    #pd.DataFrame({'channel': channel_axis, 'mV': h, 'pdf': weighted_pdf}).to_csv('cal.csv', index=False)

    ax.step(h, weighted_pdf, color='red', label='Spectrum')
    ax.plot(h, _fermi_func(h, *popt), color='blue', linestyle='--', label=label_F, linewidth=2)

    #To check the correspondence between mV and channel
    # Unfortunately this is not quite working out in the plot bc of rounding errors
    #ax2.step(channel_axis, weighted_pdf, color='green', label='Spectrum')

    ax.axvline(hFlex, color='darkviolet', linestyle='--', label=label_hFlex)
    ax.axvline(hDD, color='deeppink', linestyle='--', label=label_hDD)
    ax.axvline(hTC, color='magenta', linestyle='--', label=label_hTC)

    ax.fill_between(h[fit_bounds[0]:fit_bounds[1]], [-1e3]*len(h[fit_bounds[0]:fit_bounds[1]]), [1e3]*len(h[fit_bounds[0]:fit_bounds[1]]), color='orange', alpha=0.1, label=label_bounds)
    ax.axvline(h[fit_bounds[0]], color='orange', linestyle='--')
    ax.axvline(h[fit_bounds[1]], color='orange', linestyle='--')

    ax.axhline(0, color='black', linestyle='--')

    #ax.set_xscale('log')
    ax.set_title('Edge Calibration')
    ax.legend(framealpha=1)

    ax.set_xlim(min(h), max(h), auto=False)
    ax.set_ylim(0, max(weighted_pdf)*1.1)

    plt.show()

def _plot_stopping_power(dict, ax=None):

    table = dict['table']
    interp_range = dict['interp_range']
    interp_dEdx = dict['interp_dEdx']
    range_lower = dict['range_lower']
    range_upper = dict['range_upper']
    database = dict['database']
    particle = dict['particle']
    material = dict['material']

    if ax is None:
        fig, ax = plt.subplots()

    # Save the output as csv
    #pd.DataFrame({'range': interp_range, 'dEdx': interp_dEdx}).to_csv('stopping_power.csv', index=False)	

    ax.plot(interp_range, interp_dEdx, color='blue', label='Linear Interpolation', marker='.')
    ax.plot(table['range'], table['dE/dx'], color='red', label=f'{database} Data', marker='x', linestyle='')

    ax.axvline(range_upper, color='orange', linestyle='--')
    ax.axvline(range_lower, color='orange', linestyle='--', label='Particle Track\nImparting Highest Energy')
    ax.fill_betweenx(np.linspace(-1,5000), range_lower, range_upper, color='orange', alpha=0.1)

    if min(table['range']) < 1e-2:
        ax.set_xlim(1e-2, max(table['range']) * 1.1)
    else:
        ax.set_xlim(min(table['range']), max(table['range']) * 1.1)
    ax.set_ylim(0, max(table['dE/dx']) * 1.1)
    ax.set_xscale('log')	
    ax.set_title(f'{database} Stopping Power for {particle} in {material}')
    ax.set_xlabel('Residual Range [µm]')
    ax.set_ylabel('Electronic Stopping Power [keV/µm]')
    ax.legend(loc='upper right')
    ax.grid()
    plt.show()

def _plot_chord_length_dist(dict, ax=None):
     
    chord_length_dist = dict['chord_length_dist']
    shape = dict['shape']
    dimension = dict['dimension']
    
    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(chord_length_dist['chord_len'], chord_length_dist['pdf'], label='chord length distribution', marker='x')
    ax.set_ylabel('frequency [a.u.]')
    ax.set_xlabel('chord length (µm)')
    ax.set_title(f'Chord length distribution for {shape} of {dimension} µm')
    ax.set_xlim(0, dimension)
    ax.set_ylim(0,np.max(chord_length_dist['pdf']) * 1.1)
    plt.show()

def get_chord_length(measurement, shape:str, dimension:float, plot:bool=False):
    """ Get chord lengths for a given geometry (shape, dimensions)
    This can be used without a measurement by giving 'None'
    Or attached to the measurement object

    Currently only implemented for unidirectional beam and normal incidence.
    Other chord length distributions will be added in the future.

    Args:
        measurement: Spectrum to be analyzed
        shape: Available options: 'slab', 'sphere'
        dimensions: Dimensions of the detector in µm (just float bc. shape=slab, sphere or cylinder(h=d))
        plot (optional): Show a plot of the chord length distribution

    Returns:
        mean: Mean chord length in µm
        max: Maximum chord length in µm
        chord_length_distribution: Distribution of chord lengths in µm (for the given number of channels in the measurement)
    """

    SAMPLING = 1000

    #TODO: Implement other shapes (sphere, cylinder) -> See Kellerer or Bradley PhD Thesis for formulas
    #unidirectional beam and normal incidence
    chord_length_dist = {}

    if measurement == 'None':
        chord_length_dist['chord_len'] = np.arange(0, dimension, dimension/SAMPLING)    
    else:
        chord_length_dist['chord_len'] = np.arange(0, dimension, dimension/measurement.num_channels)

    if shape == 'slab':
        chord_length_dist['pdf'] = np.array([dimension] * measurement.num_channels)
        mean = dimension
        max = dimension
    elif shape == 'sphere':
         chord_length_dist['pdf'] = np.array([(2*i)/(dimension**2) for i in chord_length_dist['chord_len']])
         mean = 0.666 * dimension
         max = dimension
    else:
        raise ValueError(f"Shape {shape} is not a valid option (slab, sphere)")
    
    #If a measurement object is given, attach the chord lengths
    if measurement != 'None':       
        measurement.mean_chord_length = mean
        measurement.max_chord_length = max        

    chord_length_dict = {'mean': mean, 'max': max, 'chord_length_dist': chord_length_dist,
                        'shape': shape, 'dimension': dimension}

    if plot:
        _plot_chord_length_dist(chord_length_dict)
    
    return mean, max, chord_length_dist

def get_stopping_power(measurement, chord_length:str='mean', database:str='SRIM', precision:float=0.05, particle:str='None', material:str='None', plot:bool=False):
        """ Get maximum lineal energy or LET with stopping power tables.
        Algorithm for unidirectional beam and normal incidence.

        This can be used without a measurement by giving 'None'
        Or attached to the measurement object. In this case, the particle and material are taken from the measurement object.

        Args:
            measurement: Spectrum to be analzed - Can be used with a measurement object by giving 'None'
            chord_length: Chord length of the given arrangement in µm
            database: Available options 'ICRU', 'SRIM'
            precision (optional): Precision of the interpolation in keV/µm
            particle (optional): Particle to be used for the calibration (proton, carbon, helium)
            material (optional): Material to be used for the calibration (silicon, diamond, sic, water)
            plot (optional): Show a plot of the interpolation

        Returns:
            ymax: Maximum lineal energy in keV/µm
            Lmax: Maximum LET in keV/µm
            interpol_dict: Dictionary with all relevant information about the interpolation
        """
        # Include whereever this script is located
        ressources_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/ressources/stopping_power_tables/'

        densities = {'silicon':2.32, 'diamond': 3.52, 'sic': 3.21, 'water': 1} #g/cm3

        if not database in ['ICRU', 'SRIM', 'NIST']:
            raise ValueError(f"Database {database} is not a valid option (ICRU, SRIM, NIST)")

        #Option with measurement object
        if measurement != 'None':

            if 'None' in measurement.detector:
                raise ValueError('No detector attached to measurement object. Do this with campaign.attach_info()') 
            if 'None' in measurement.particle:
                raise ValueError('No particle attached to measurement object. Do this with campaign.attach_info()')

            particle = measurement.particle
            material = measurement.detector

            print(f"Using particle {particle} and material {material} from measurement object")

            if chord_length == 'mean':
                chord_length = measurement.mean_chord_length
            elif chord_length == 'max':
                chord_length = measurement.max_chord_length
            else:
                raise ValueError(f"Chord length {chord_length} is not a valid option (mean, max)")            

        #Option without measurement object
        else:
            if particle == 'None' or material == 'None':
                raise ValueError('No particle or material is given')

            chord_length = float(chord_length)
            
        if particle == 'proton':
            nucleons = 1
        elif particle == 'carbon':
            nucleons = 12
        elif particle == 'helium':
            nucleons = 4
        else:
            raise ValueError(f"Calibration not implemented for {particle}. Current options: proton, carbon, helium")
        
        if not os.path.exists(f"{ressources_path}/{database}_{particle}_{material}.csv"):
            raise ValueError(f"Stopping power table for {particle} in {material} not found in {ressources_path}")
        
        ################################
        # PREPARE DATA
        ################################
        table = pd.read_csv(f"{ressources_path}/{database}_{particle}_{material}.csv",
                            sep=',', skiprows=5, names=['energy', 'range', 'dE/dx'])
        table['dE/dx'] = table['dE/dx'].multiply(densities[material])

        interp_range = np.arange(0, 100 + precision, precision) #From 0 to 100 µm
        interp_dEdx = []

        range_values = table['range'].to_numpy()

        for ran in interp_range:
            values = _get_nearest_range(table, range_values, ran)

            if ran < range_values[0]:
                interp_dEdx.append(math.nan)
            else:
                interp_dEdx.append((ran - values['closest_range']) / (values['after_closest_range'] - values['closest_range']) * (values['after_closest_dEdx']-values['closest_dEdx']) + values['closest_dEdx'])
            
        interpol_table = pd.DataFrame({'interp_range': interp_range, 'interp_dEdx': interp_dEdx})

        mass_thickness = chord_length * densities[material]
        num_steps = int(mass_thickness / precision)

        ################################
        # START OF CALCULATIONS
        ################################

        #Sum up dEdx values -> Integrate dE/dx numerically (Stammfunktion)
        sum_dEdx = []
        for i in range(len(interp_range)):
             sum_dEdx.append(np.nansum(interpol_table['interp_dEdx'][i:i+num_steps+1]))
        interpol_table['sum_dEdx'] = sum_dEdx

        #Find residual ranges -> Also R' (für das das Integral maximal ist) und R'+t
        range_lower = precision * (interpol_table['sum_dEdx'].idxmax())
        range_upper = range_lower + mass_thickness

        #print(f"Residual range: {range_lower} µm")
        #print(f"Residual range + chord length: {range_upper} µm")

        #Find energy at R'
        values = _get_nearest_range(table, range_values, range_lower)
        exit_energy = (range_lower - values['closest_range']) / (values['after_closest_range'] - values['closest_range']) * (values['after_closest_energy']-values['closest_energy']) + values['closest_energy']
        if database == 'ICRU':
            exit_energy *= nucleons

        #print(values)
        #print(f"Energy at R': {exit_energy} keV")

        #Find energy at R'+t
        values = _get_nearest_range(table, range_values, range_upper)
        entry_energy = (range_upper - values['closest_range']) / (values['after_closest_range'] - values['closest_range']) * (values['after_closest_energy']-values['closest_energy']) + values['closest_energy']
        if database == 'ICRU':
            entry_energy *= nucleons

        #print(values)
        #print(f"Energy at R'+t: {entry_energy} keV")

        #Calculate ymax, Lmax
        Lmax_tables = interpol_table['sum_dEdx'].max() / num_steps / densities[material]
        Lmax = (entry_energy - exit_energy) / (range_upper - range_lower)
        ymax = (entry_energy - exit_energy) / chord_length

        #print(range_lower, range_upper, entry_energy, exit_energy, Lmax_tables, Lmax, ymax)

        interpol_dict = {'ymax': ymax, 'Lmax': Lmax, 'interpol_table': interpol_table,
                         'range_lower': range_lower, 'range_upper': range_upper,
                         'table': table, 'interp_range': interp_range, 'interp_dEdx': interp_dEdx,
                         'database': database, 'particle': particle, 'material': material}

        if plot:
            _plot_stopping_power(interpol_dict)
        
        return ymax, Lmax, interpol_dict

def get_edge_pos(measurement, marker_point:str, fit_bounds:list=None, check_plot:bool=False):
        """ Find edge position for a given spectrum. This is done via a Fermi function following:        
        
        [1] D. Moro, et al., Radiat Prot Dosimetry, vol. 166, no. 1-4, pp. 233-237, Sep. 2015, doi: 10.1093/rpd/ncv153.
        [2] V. Conte, et al., AIP Conference Proceedings, vol. 1530, no. 1, pp. 171-178, Jul. 2013, doi: 10.1063/1.4812920.

        Args:
            measurement: Spectrum to be analzed
            marker_point: Which curve metric is to be used as the edge (hFlex, hDD, hTC)
            fit_bounds (optional): Channel bounds for the fit function (Default: [100,2500])
            check_plot (optional): Show a plot of the fit

        Returns:
            edge_position_mV: The selected metric for the particle edge
            cal_dict: Dictionary with all relevant information about the calibration
         """
        
        if measurement.x_axis != 'CHANNEL' or measurement.y_axis != 'COUNTS':
               raise ValueError('This manipulation only works for not calibrated MCA spectra [channel, counts]')

        #Get mV per channel from MCA file
        h = measurement.data['mV'].to_numpy()
        counts = measurement.data[measurement.y_axis].to_numpy()

        #Calculate D(h)
        weighted_num_p_above_cutoff = float(np.sum(np.multiply(counts, h)))
        D_temp = [0.0]
        for i in range(1,len(measurement.data)):
            mul = np.multiply(counts[:i], h[:i])
            D_temp.append(np.sum(mul))
        D_norm = np.divide(D_temp, weighted_num_p_above_cutoff, out=np.zeros_like(D_temp), where=weighted_num_p_above_cutoff!=0)

        #Calculate d(h)
        pdf = [0.0]
        for i in range(1,len(measurement.data)):
            delta_y = h[i] - h[i-1]
            #Just a temporary sanity check
            if h[i] < h[i-1]:
                raise ValueError('Something is wrong with the data: y[i] < y[i-1]')
            else:
                delta_F = D_norm[i] - D_norm[i-1]
                pdf.append(np.divide(delta_F,delta_y, out=np.zeros_like(delta_F), where=delta_y!=0))

        #Calculate hd(h) for fitting
        for i in range(1,measurement.num_channels):
            weighted_pdf = np.multiply(h, pdf)

        #Fitting
        if fit_bounds is None:
            fit_bounds = [100,3500] #Default channel bounds for fitting

        h_fit = h[fit_bounds[0]:fit_bounds[1]]
        weighted_pdf_fit = weighted_pdf[fit_bounds[0]:fit_bounds[1]]
        #For a smoother fit -> Start on a plateau
        weighted_pdf_fit[:10] = weighted_pdf_fit[11]

        # Initial guesses (regarding for fit bounds)
        A_guess = 0.5
        B_guess = 0.1
        C_guess = h_fit[0] + (h_fit[-1] - h_fit[0]) / 2

        print(f"Initial guesses: A={A_guess}, B={B_guess}, C={C_guess}")

        bounds = [[0.01, 0.01, -25], [10, 1000, 4096]] #A,B,C min ; A,B,C max
        popt, pcov = curve_fit(_fermi_func, h_fit, weighted_pdf_fit,
                               p0=[A_guess, B_guess, C_guess], bounds=(bounds[0], bounds[1]))
        
        print(f"Optimal parameters: A={popt[0]}, B={popt[1]}, C={popt[2]}")

        #calculate marker points
        hFlex = popt[2]
        hDD = (math.log(2+math.sqrt(3))) / popt[1] + popt[2]
        hTC = 2/popt[1] + popt[2]

        cal_dict = {'hFlex': hFlex, 'hDD': hDD, 'hTC': hTC,
                    'h': h, 'weighted_pdf': weighted_pdf, 'popt': popt, 'fit_bounds': fit_bounds,
                    'channel_axis': measurement.data['CHANNEL'].to_numpy()}
            
        if check_plot:
            _plot_edge_calibration(cal_dict)
        
        #Return edge position
        if marker_point == 'hFlex':
             edge_pos_mV = popt[2]
        elif marker_point == 'hDD':
             edge_pos_mV = (math.log(2+math.sqrt(3))) / popt[1] + popt[2]
        elif marker_point == 'hTC':
             edge_pos_mV = 2/popt[1] + popt[2]
        else:
             raise ValueError(f"Marker point {marker_point} is not a valid option (hFlex, hDD, hTC)")

        return edge_pos_mV, cal_dict

def scale_energy_axis(measurement, edge_pos_mV:float, max_energy_transfer:float, energy_axis:str='lineal', chord_length:str='mean'):
        """ Scales x-axis from CHANNEL to keVµm-1
        The edge position has to be specified in mV. This can be done via the get_edge_pos method.
        A chord length has to be attached before scaling the energy axis. This can be done via the get_chord_length method.

        Args:
            measurement: Spectrum to be analzed
            edge_pos_mV: Position of the Particle edge (Returned from Calibrate.get_edge_pos)
            max_energy_transfer: Maximum energy transfer (ymax, Lmax returned from Calibrate.get_stopping_power)
            energy_axis: 'lineal' (keV/um) or 'energy' (keV)
            chord_length: Choose, if mean or max chord length is to be considered
        """

        if measurement.mean_chord_length == 0.0 or measurement.max_chord_length == 0.0:
            raise ValueError('No chord length attached to measurement object. Do this with Calibrate.get_chord_length()')

        if chord_length == 'mean':
             chord = measurement.mean_chord_length
        elif chord_length == 'max':
                chord = measurement.max_chord_length
        else:
            raise ValueError(f"Chord length {chord_length} is not a valid option (mean, max)")

        if measurement.x_axis != 'CHANNEL' or measurement.y_axis != 'COUNTS':
               raise ValueError('This manipulation only works for not calibrated MCA spectra [channel, counts]. For simulation data use: lineal_energy_axis()')
        
        if energy_axis == 'lineal':
            scaling_coeff = max_energy_transfer / edge_pos_mV

            measurement.data[measurement.x_axis] = measurement.data['mV'].multiply(scaling_coeff)
            print('changed axis lineal')
            measurement.data.drop(columns='mV', inplace=True)
            measurement.x_axis = 'LINEAL ENERGY'
                
        elif energy_axis == 'energy':
            scaling_coeff = max_energy_transfer / edge_pos_mV * chord

            measurement.data[measurement.x_axis] = measurement.data['mV'].multiply(scaling_coeff)
            print('changed axis linear')
            measurement.data.drop(columns='mV', inplace=True)
            measurement.x_axis = 'ENERGY'
        else:
            raise ValueError(f"Axis {energy_axis} is not a valid option")
        
        #print(f'Energy axis scaled to {measurement.x_axis} axis')

def scale_energy_axis_with_factor(measurement, factor:float, energy_axis:str='lineal'):
        """ Scales x-axis with a given conversion factor
        Factor has to be given in terms of mV -> keV/um or keV

        Args:
            measurement: Spectrum to be analzed
            factor: Scaling factor for the energy axis
            energy_axis: 'lineal' (keV/um) or 'energy' (keV)
        """

        if measurement.x_axis != 'CHANNEL' or measurement.y_axis != 'COUNTS':
               raise ValueError('This manipulation only works for not calibrated MCA spectra [channel, counts]. For simulation data use: lineal_energy_axis()')

        if energy_axis == 'lineal':
            measurement.data[measurement.x_axis] = measurement.data['mV'].multiply(factor)
            print('changed axis lineal')
            measurement.data.drop(columns='mV', inplace=True)
            measurement.x_axis = 'LINEAL ENERGY'
                
        elif energy_axis == 'energy':
            measurement.data[measurement.x_axis] = measurement.data['mV'].multiply(factor)
            print('changed axis linear')
            measurement.data.drop(columns='mV', inplace=True)
            measurement.x_axis = 'ENERGY'
        else:
            raise ValueError(f"Axis {energy_axis} is not a valid option")

def lineal_energy_axis(measurement, chord_length:str='mean'):
    """ Scale an energy axis to a lineal energy axis by giving the mean_chord_length

    Args:
        measurement: Spectrum to be analzed
        chord_length: Choose, if mean or max chord length is to be considered
    """

    if measurement.mean_chord_length == 0.0 or measurement.max_chord_length == 0.0:
            raise ValueError('No chord length attached to measurement object. Do this with Calibrate.get_chord_length()')

    if chord_length == 'mean':
        calc_chord_length = measurement.mean_chord_length
    elif chord_length == 'max':
        calc_chord_length = measurement.max_chord_length

    if measurement.x_axis != 'ENERGY':
               raise ValueError(f'This manipulation only works for energy x-axis in keV - not {measurement.x_axis}')
    else:
        measurement.data[measurement.x_axis] = measurement.data[measurement.x_axis].divide(calc_chord_length)
        print('changed axis lineal')
        measurement.x_axis = 'LINEAL ENERGY'
    
    print(f'Energy axis scaled to {measurement.x_axis}')
