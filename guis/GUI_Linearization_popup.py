from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QFileDialog, QCheckBox, QComboBox
from PyQt5.QtWidgets import QDialog


class LinearizationPopup(QDialog):
    """ Popup window to set the linearization parameters. """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.init_linearization_popup()
        self.init_lin_analysis_state_table()
    
    def init_linearization_popup(self):
        """ Set up the linearization popup window widgets """

    def init_lin_analysis_state_table(self):
        """ Table to display the current state of the analysis """
        pass