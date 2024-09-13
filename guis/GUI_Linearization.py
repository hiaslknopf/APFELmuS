from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QFileDialog, QCheckBox, QComboBox
from PyQt5.QtWidgets import QDialog

class Linearization(QDialog):
    """ Popup GUI to create a new linearization from a recorded spectrum """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.init_linearization_GUI()
    
    def init_linearization_GUI(self):
        """ Set up the linearization window widgets """