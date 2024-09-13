from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QFileDialog, QCheckBox, QComboBox
from PyQt5.QtWidgets import QDialog

class CalibrationPopup(QDialog):
    """ Popup window to set the calibration parameters. """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.init_calibration_popup()
        self.init_cal_analysis_state_table()
    
    def init_calibration_popup(self):
        """ Set up the calibration popup window widgets """

    def init_cal_analysis_state_table(self):
        """ Table to display the current state of the analysis """
        pass