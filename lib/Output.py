__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

import matplotlib.pyplot as plt
import math
import numpy as np
import os

""" This module contains functions for plotting and saving the results of the microdosimetric spectra. """

labels_dict = {
    'x_labels':{
        'CHANNEL': 'MCA channel',
        'LINEAL ENERGY': 'y [keVµm$^{-1}$]',
        'ENERGY': 'E [keV]',
        'LET': 'LET [keVµm$^{-1}$]'
        },
    'y_labels': {
        'COUNTS': 'Counts [a.u.]',
        'F(y)': 'F(y) [keV$^{-1}$µm]',
        'F(E)': 'F(E) [keV$^{-1}$]',
        'D(y)': 'D(y) [a.u.]',
        'D(E)': 'D(E) [a.u.]',
        'f(y)': 'f(y)  [a.u.]',
        'f(E)': 'f(E)  [a.u.]',
        'd(y)': 'd(y)  [a.u.]',
        'd(E)': 'd(E)  [a.u.]',
        'yf(y)': 'yf(y)  [a.u.]',
        'Ef(E)': 'Ef(E)  [a.u.]',
        'yd(y)': 'yd(y)  [a.u.]',
        'Ed(E)': 'Ed(E)  [a.u.]',
        't(LET)': 't(LET) [a.u.]',
        'LETt(LET)': 'LETt(LET) [a.u.]'
        },
    'mean_labels': {
        'f(y)': 'y_F',
        'd(y)': 'y_D',
        'Ef(E)': 'EF',
        'Ed(E)': 'ED',
        'yf(y)': 'yF',
        'yd(y)': 'yD',
        't(LET)': 'Mean LET'
    }
}

titles_dict = {
    'COUNTS': 'Pulse height spectrum',
    'F(y)': 'Probability F(y)',
    'F(E)': 'Probability F(E)',
    'D(y)': 'Probability D(y)',
    'D(E)': 'Probability D(E)',
    'f(y)': 'Probability density f(y)',
    'f(E)': 'Probability density f(E)',
    'd(y)': 'Probability density d(y)',
    'd(E)': 'Probability density d(E)',
    'yf(y)': 'Weighted probability density yf(y)',
    'Ef(E)': 'Weighted probability density Ef(E)',
    'yd(y)': 'Weighted probability density yd(y)',
    'Ed(E)': 'Weighted probability density Ed(E)',
    't(LET)': 'LET distribution t(LET)',
    'LETt(LET)': 'Weighted LET distribution LETt(LET)'
}

def __on_pick(event, graphs, fig):
    """ Function for interactive plot control

    Args:
        event: Event that triggers the function
        graphs: Dictionary with line objects
        fig: Figure object to be updated
    """

    legend = event.artist
    isVisible = legend.get_visible()

    graphs[legend].set_visible(not isVisible)
    legend.set_visible(not isVisible)

    fig.canvas.draw()

def get_mean_value(measurement):
    """ Calculate mean value of a given spectrum via integration

    Warning: This is not yet implemented for all possible combinations and not tested!
    
    Args:
        measurement: The spectrum to be analyzed
    Returns:
        mean_F: Frequency mean value
        mean_D: Dose mean value
    """

    x = measurement.data[measurement.x_axis].to_numpy()
    y = measurement.data[measurement.y_axis].to_numpy()

    if np.trapz(x*y, x=x, dx=1.0) == 1.0:
        raise ValueError('The spectrum cannot be normalized when taking an average. The integral of the spectrum is 1.0. Please check the data.')

    if measurement.y_axis == 'f(y)':
        mean_F = np.trapz(x*y, x=x, dx=1.0)
        mean_D = np.trapz(x*x*y, x=x, dx=1.0) / mean_F
    
    elif measurement.y_axis == 'd(y)':
        mean_F = 0.0
        mean_D = np.trapz(x*y, x=x, dx=1.0)
    
    elif measurement.y_axis == 'yf(y)':
        mean_F = np.trapz(y, x=x, dx=1.0)
        mean_D = np.trapz(x*y, x=x, dx=1.0)
    
    elif measurement.y_axis == 'yd(y)':
        mean_F = 0.0
        mean_D = np.trapz(y, x=x, dx=1.0)
        
    if measurement.x_axis == 'ENERGY':
        unit = 'keV'
    elif measurement.x_axis == 'LINEAL ENERGY':
        unit = 'keVµm$^{-1}$'
    elif measurement.x_axis == 'CHANNEL':
        unit = 'channels'
    elif measurement.x_axis == 'LET':
        unit = 'keVµm$^{-1}$'

    if mean_F == 0.0:
        print('Frequency mean value cannot be calculated from the given spectrum')
    else:
        print(f'Frequency mean value of {measurement.y_axis} = {mean_F: .3f} {unit}')
    if mean_D == 0.0:
        print('Dose mean value cannot be calculated from the given spectrum')
    else:   
        print(f'Dose mean value of {measurement.y_axis} = {mean_D: .3f} {unit}')

    return mean_F, mean_D

def plot_single(measurement, name:str=False, output_path:str=False,
                mean:float=False, step:bool=False, scale:str='log', xlim:list=[0.1, 1000],
                show_plot:bool=True, interactive:bool=True, gui=False, fig=None, ax=None):
    """ View and save a single spectrum plot

    Args:
        measurement: The spectrum to be analyzed
        name: Title and filename for the plot (Default: Spectrum name)
        output_path: For the figure to be saved (optional)
        mean: Include a mean value in the plot (optional: float)
        step: Use step plot instead of line plot (optional, default=False)
        scale: Scaling of x-axis (optional: 'lin', default: 'log')
        show_plot: Immediately show test plot (optional, default=True)
        interactive: Enable interactive plot control (Turn on and off lines) (optional, default=True)
        gui: Use existing figure and axis objects generated by GUIs (only for internal use)
    """

    lines = []

    if not name:
        name = measurement.name

    type_of_plot = titles_dict[measurement.y_axis]

    if not gui:
        fig, ax = plt.subplots()

    x = measurement.data[measurement.x_axis]
    y = measurement.data[measurement.y_axis]

    ax.axhline(0, color='black', lw=0.5, linestyle='--')

    if measurement.gain == 'SIMULATION':
        data_label = 'Geant4 Simulation'
    else:
        data_label = 'Measurement'

    if step:
        line, = ax.step(x, y, color='red', where='mid', label=data_label)
    else:
        line, = ax.plot(x, y, color='red', label=data_label)
    
    lines.append(line)

    if measurement.x_axis == 'CHANNEL' or scale == 'lin':
        if xlim != [1, 1000]:
            ax.set_xlim(xlim[0],xlim[1])
        else:
            ax.set_xlim(1,measurement.num_channels)
    elif scale == 'log':
        ax.set_xscale('log')
        ax.set_xlim(xlim[0],xlim[1])
    
    if math.isnan(max(y)) or math.isinf(max(y)):
        pass
    else:
        ax.set_ylim(0.0,1.1*max(y))
    
    if mean:
        mean = ax.axvline(mean, color='blue', linestyle='--', label=labels_dict['mean_labels'][measurement._y_axis])
        lines.append(mean)

    ax.set_xlabel(labels_dict['x_labels'][measurement._x_axis])
    ax.set_ylabel(labels_dict['y_labels'][measurement._y_axis])
    ax.set_title(f"{name}")

    legend = plt.legend(loc='upper left')
    plt.tight_layout()

    if interactive:
        graphs = {}
        for i, line in enumerate(range(len(lines))):
            line = legend.get_lines()[i]
            line.set_picker(True)
            line.set_pickradius(10)

            graphs[line] = lines[i]

    if show_plot:
        if interactive:
            plt.connect('pick_event', lambda event: __on_pick(event, graphs, fig))
        plt.show()

    if output_path:
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        fig.savefig(f"{output_path}/{name}.png", dpi=600, bbox_inches='tight', orientation='portrait')
        print(f'"{name}".png saved to "{output_path}"')

def plot_campaign(campaign, name:str=False, output_path:str=False, xlim:list=[1, 1000],
                  step:bool=False, scale:str='log', skip_noise:int = 50, show_plot:bool=True, interactive:bool=True):
    """ View and save multiple spectra in one plot
    
    Args:
        campaign: Set of measurements to be plotted (MicroDosimetry object)
        name: Title and filename for the plot (Default: "Campaign Testplot")
        output_path: For the figure to be saved (optional)
        step: Use step plot instead of line plot (optional, default=False)
        scale: Scaling of x-axis (optional: 'lin', default: 'log')
        skip_noise: Skip initial tail of the spectrum for finding the optimal ymax for the plot (optional, default=50)
        show_plot: Immediately show test plot (optional, default=True)
        interactive: Enable interactive plot control (Turn on and off lines) (optional, default=True)
    """

    meas_names = list(campaign.measurements.keys())
    type_of_plot = titles_dict[campaign.measurements[meas_names[0]].y_axis]

    # All datasets have to have the same x and y axis
    for measurement in campaign.measurements.values():
        if measurement.x_axis != campaign.measurements[meas_names[0]].x_axis or measurement.y_axis != campaign.measurements[meas_names[0]].y_axis:
            raise ValueError(f'All measurements in a campaign have to have the same x and y axis. {measurement.name} does not match ({campaign.measurements[meas_names[0]].x_axis}({campaign.measurements[meas_names[0]].y_axis})))')
        if measurement.num_channels != campaign.measurements[meas_names[0]].num_channels:
            raise ValueError(f'All measurements in a campaign have to have the same number of channels. {measurement.name} does not match ({campaign.measurements[meas_names[0]].num_channels})')

    if not name:
        name = 'Campaign Testplot'

    fig, ax = plt.subplots()
    ymax = 0

    # Set up the colormap
    col = plt.colormaps['gist_rainbow'](np.linspace(0, 1, len(meas_names)))

    lines = []
    for i, meas in enumerate(meas_names):
        measurement = campaign.measurements[meas]

        x = measurement.data[measurement.x_axis]
        y = measurement.data[measurement.y_axis]

        if step:
            line, = ax.step(x, y, color=col[i], where='mid', label=meas_names[i], alpha=0.5)
        else:
            line, = ax.plot(x, y, color=col[i], label=meas_names[i], alpha=0.5)
        
        lines.append(line)

        ymaxnew = max(y[skip_noise:]) # To skip initial tail
        if ymaxnew > ymax:
            ymax = ymaxnew

    if measurement.x_axis == 'CHANNEL' or scale == 'lin':
        ax.set_xlim(1,measurement.num_channels)
    elif scale == 'log':
        ax.set_xscale('log')
        ax.set_xlim(xlim)
        
    if math.isnan(ymax) or math.isinf(ymax):
        pass
    else:
        ax.set_ylim(0.0,1.1*ymax)

    ax.set_xlabel(labels_dict['x_labels'][campaign.measurements[meas_names[0]]._x_axis])
    ax.set_ylabel(labels_dict['y_labels'][campaign.measurements[meas_names[0]]._y_axis])
    ax.set_title(f"{type_of_plot}: {name}")

    legend = plt.legend(loc='upper right')
    plt.tight_layout()

    if interactive:
        graphs = {}
        for i, line in enumerate(range(len(lines))):
            line = legend.get_lines()[i]
            line.set_picker(True)
            line.set_pickradius(10)

            graphs[line] = lines[i]

    if show_plot:
        if interactive:
            plt.connect('pick_event', lambda event: __on_pick(event, graphs, fig))
        plt.show()

    if output_path:
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        fig.savefig(f"{output_path}/{name}.png", dpi=600, bbox_inches='tight', orientation='portrait')
        print(f'"{name}".png saved to "{output_path}"')

def csv_output(measurement, output_path, name=False):
    """ Save a single spectrum as csv file

    Args:
        measurement: The spectrum to be saved
        output_path: Where to save the file
        name: Filename (optional, default: Spectrum name)
    """

    if 'None' in measurement.detector:
        print('No detector information available. Please add detector information to the measurement object before saving the data.')

    if not name:
        name = measurement.name
    
    if not os.path.exists(output_path):
            os.makedirs(output_path)
    
    with open(f'{output_path}/{name}.csv', 'w') as f:
        f.write(f"APFELmuS output: {titles_dict[measurement.y_axis]}\n{measurement.particle} on {measurement.detector}\n{measurement.date}\n\n")
        
        xlabel = labels_dict['x_labels'][measurement._x_axis]
        ylabel = labels_dict['y_labels'][measurement._y_axis]

        # Replace LaTeX formatting with plain text
        if '$^{-1}$' in xlabel or 'µm' in xlabel:
            xlabel = xlabel.replace('$^{-1}$', '-1')
            xlabel = xlabel.replace('µm', 'um')
        if '$^{-1}$' in ylabel or 'µm' in ylabel:
            ylabel = ylabel.replace('$^{-1}$', '-1')
            ylabel = ylabel.replace('µm', 'um')

        if 'mV' in measurement._data.columns:
            f.write(f"{xlabel};{ylabel};mV\n")
            print('mV column found in data. Saving as third column in csv file - But theres no point in this, its just an .MCA file.')
        else:
            f.write(f"{xlabel};{ylabel}\n")

    measurement._data.to_csv(f'{output_path}/{name}.csv', index=False, sep=';', mode='a', header=False)

    # Metainfo in footer
    with open(f'{output_path}/{name}.csv', 'a') as f:

        f.write('\n--------------------------------------------------------------------------------\n')

        for i, j in zip(['NAME', 'DATE', 'PARTICLE', 'DETECTOR', 'GAIN',
                        'MEAN_CHORD', 'MAX_CHORD', 'LIVE_TIME', 'REAL_TIME'],
                        [measurement.name, measurement.date, measurement.particle, measurement.detector, measurement.gain,
                        measurement.mean_chord_length, measurement.max_chord_length, measurement.live_time, measurement.real_time]):
            f.write(f"{i}\t{j}\n")

        f.close()

    print(f'"{name}".csv saved to "{output_path}"')