import os 
import re
import numpy as np
import vtk 
from DataReader import DataReader
import math 

class ReaderCsv(DataReader):
    def __init__(self, cof, ap, parameter=None, value=None, muw=None, vwall=None, fraction=None, phi = None, I = None):
        super().__init__(cof, ap, parameter, value, muw, vwall, fraction, phi)


    def get_data(self, global_path):
        """
        Read the data from the csv files
        """

        data = np.genfromtxt(self.directory + "time_series/data.csv", delimiter=',', skip_header=1)
        self.time = data[:, 0]
        self.Fx = data[:, 1]
        self.Fy = data[:, 2]
        self.Fz = data[:, 3]
        self.tke = data[:, 4] #translational kinetic energy
        self.rke = data[:, 5] #rotational kinetic energy  

