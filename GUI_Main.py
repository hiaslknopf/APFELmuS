from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QFileDialog, QCheckBox, QComboBox
from PyQt5.QtWidgets import QDialog
from PyQt5.QtGui import QIcon

from guis.GUI_AnalysisManager import AnalysisManager

from guis.GUI_Linearization import Linearization
from guis.GUI_Calibration import Calibration
from guis.GUI_Calibration_popup import CalibrationPopup
from guis.GUI_Linearization_popup import LinearizationPopup

from guis.GUI_SaveCSVPopup import SaveCSVPopup
from guis.GUI_SavePlotsPopup import SavePlotsPopup

class GUI_Main(QWidget):
    """ Main GUI window for the APFELmuS analysis.

    Data is to be analyzed in campaigns. Each campaign is displayed in a new tab in a notebook.
    The state of the analysis and the GUI is managed by the AnalysisManager class. It can be accessed and stored for later use.
    A json file is temporarily stored in the data directory so the state of the analysis can be accessed by several GUIs.

    There are additional GUI windows for performing a calibration or linearization or the conversion of a recorded spectrum to a APFELmuS .MCA file.
    Analysis settings will be chosen in popup windows.
    
    """

    def __init__(self):
        super().__init__()

        self.init_main_window()
        self.init_main_buttons()
        self.init_campaign_notebook()

    def init_main_window(self):
        """ Set up the main window widgets """

        self.setWindowTitle("APFELmuS - Analysis GUI")

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

    def init_main_buttons(self):
        """ Main control buttons on the left side """
        pass

    def init_campaign_notebook(self):
        """ Notebook to display the different analysis campaigns """

        pass
        
    def init_new_notebook_tab(self):
        """ Create a new notebook tab for a new analysis campaign.
            This contains:
                analysis_main_buttons
                analysis_state_table
                analysis_plot_settings
                analysis_plot_notebook
                analysis_info_textbox
                analysis_output_buttons """
        
        raise NotImplementedError
    
    def init_analysis_main_buttons(self):
        """ Main control buttons for the analysis campaign """
        
        raise NotImplementedError
    
    def init_analysis_state_table(self):
        """ Table to display the current state of the analysis """
        
        raise NotImplementedError
    
    def init_analysis_plot_settings(self):
        """ Settings for the plot window """
        
        raise NotImplementedError
    
    def init_analysis_plot_notebook(self):
        """ Create the plot window in a new tab """
        
        raise NotImplementedError
    
    def init_analysis_info_textbox(self):
        """ Textbox to display information """
        
        raise NotImplementedError
    
    def init_analysis_output_buttons(self):
        """ Buttons to save and load data """
        
        raise NotImplementedError