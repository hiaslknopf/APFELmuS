from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QFileDialog, QCheckBox, QComboBox
from PyQt5.QtWidgets import QDialog

class Calibration(QDialog):
    """ Popup GUI to create a new calibration from a recorded spectrum (edge calibration) 
    Can be opened either from the main GUI or from the calibration popup """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.init_calibration_GUI()
    
    def init_calibration_GUI(self):
        """ Set up the calibration window widgets """