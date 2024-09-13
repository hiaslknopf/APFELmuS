class AnalysisManager():
    """ Manager class to handle the GUI settings, communication and analysis with APFELmuS """

    def __init__(self, campaign_name: str, output_path: str):
        """ One object is created per analysis campaign.
            All information for the spectrum analysis is stored in a temporary json file which can be saved and loaded. """

        self.campaign_name: str = campaign_name
        self.output_path: str = output_path

        self.ana_dict: dict = {}

    def create_ana_dict(self):
        """ Create and save an empty analysis dictionary """

        self.ana_dict = {
            "campaign_name": self.campaign_name,
            "output_path": self.output_path,

            "data": {                           # Currently loaded data for this session
                "data_file_paths": [],          # List of paths to the data files -> To reload the data if necessary
                "microdosimetry_object": None,
                "data_mapping_dict": {},        # Geladenen Daten durchnummerieren, dass in analysis_state gscheid referenziert werden kann (filenames sind so lang)

                # In case any means or other parameters are calculated, they are stored in these lists
                "data_maxpos": [],              # List of maxpos means data
                "data_yF": [],                  # List of yF means data
                "data_yD": [],                  # List of yD means data
                "data_y_star": [],              # List of y* means data"
                "data_RBE": [],                 # List of RBE means data

                "data_param_label": [],         # List of labels for the parametrized analysis
                "data_param_list": [],          # List of values for the parametrized analysis
            },
            
            "analysis_state": {                 # Current state of the analysis (what has been done to each measurement_dict entry)
                                                # Damit kannst zu jeden Zeitpunkt diese Excel Tabelle in jeder GUI erstellen

                "include": [],                  # List of bools, if the data should be included in the analysis (if it is currently ticked in the analysis table)
                "linearized": [],               # List of bools, if the data has already been linearized
                "calibrated": [],               # List of bools, if the data has already been calibrated
            },

            "linearization": {                  # Linearization settings
                "linearization_filepaths": [],  # List of paths to the linearization files (if there are any)
                                                # After linearization you have to save the data as a .MCA file and then load it again into the campaign
            },

            "calibration": {                    # Calibration settings
                                                # Depending on where the calibration data comes from, there are entries in different lists and different functions will ru

                "calibration_methods": [],      # List of calibration method for each spectrum
                "calibration_factors": [],      # List of calibration factors
                "calibration_files": [],        # List of calibration files
                "calibration_edges": [],        # List of calibration edges
                "calibration_chordlengths": [], # List of calibration chordlengths
            },

            "gui_state": {                      # Current state of the GUI (what is displayed, what is selected)
            },

            "plot_settings": {                  # Plot settings for the analysis (xlim, ylim, labels, etc.)
                
                "spectra": {                    # Settings for the spectra plot
                },

                "parametrized": {               # Settings for the parametrized plots
                },
            
            }
        }

        raise NotImplementedError
    
    def load_ana_dict(self, ana_dict_path: str):
        """ Load an existing analysis dictionary """

        # Read from the file and create the ana_dict

        raise NotImplementedError
    
    def save_ana_dict(self, save_path: str):
        """ Save the current analysis dictionary """

        # Save the ana_dict permanently to a file

        raise NotImplementedError
    
    def refresh_ana_dict(self):
        """ Refresh the analysis dictionary from the current GUI changes"""

        # Compare current version with the new version from the GUI

        # Update the analysis dictionary

        raise NotImplementedError
    
    def error_checker(self):
        """ Check if the current analysis dict is valid """

        # Prüf Änderungen auf Widersprüche und Vollständigkeit -> Gscheide Fehlermeldungen

        raise NotImplementedError
    
    def run_spectrum_analysis(self):
        """ Run the analysis with the current analysis dictionary """

        # Fang komplett von vorne an
        # Sollt das Performancmässig oasch rennen berechne immer nur die changes

        raise NotImplementedError
    
    