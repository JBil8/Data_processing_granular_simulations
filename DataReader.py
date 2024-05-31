import os 
import re
import numpy as np
import vtk 
from abc import ABC, abstractmethod

class DataReader:
    def __init__(self, cof, ap, parameter=None, value=None, muw=None, vwall=None, fraction=None, phi = None, I = None):
        """
        Specify the parameters that vary in the parametric study
        cof: coefficient of friction
        ap: aspect ratio
        I: inertial number
        phi: volume fraction
        """
        self.cof = cof
        self.ap = ap
        self.parameter = parameter
        self.value = value
        self.muw = muw
        self.vwall = vwall
        self.fraction = fraction
        self.phi = phi
        self.n_sim = None
        self.step = None
        self.directory = None
        self.file_list = None

    def get_number_of_time_steps(self):
        """
        Deduce number of time steps from the number of files and the step between them
        """
        digits = [int(''.join(re.findall(r'\d+', filename))) for filename in self.file_list]
        self.n_sim = len(self.file_list)
        self.step = int((max(digits) - min(digits)) / (self.n_sim - 1))

    def prepare_data_obstruction(self, global_path, reverse_flag=False):
        if reverse_flag:
            self.directory = global_path + f'reverse_alpha_{self.ap}_mup_{self.cof}_muw_{self.muw}_v_{self.vwall}_frac_{self.fraction}/'
        else:
            self.directory = global_path + f'alpha_{self.ap}_mup_{self.cof}_muw_{self.muw}_v_{self.vwall}_phi_{self.phi}_frac_{self.fraction}/'
        self.file_list = os.listdir(self.directory)


    def prepare_data(self, global_path):
        if self.parameter == "I":
            self.directory = global_path + f'I_{self.value}_ap_{self.ap}_cof_{self.cof}/'
        elif self.parameter == "phi":
            self.directory = global_path + f'phi_{self.value}_ap_{self.ap}_cof_{self.cof}_v_1/'
        else:
            raise ValueError('The parameter must be either I or phi')    
        self.file_list = os.listdir(self.directory)

    @abstractmethod
    def filter_relevant_files(self, prefix='shear_ellipsoids_'):
        pass

    def sort_files_by_time(self):
        """
        Sort the files by time
        """
        digits = [int(''.join(re.findall(r'\d+', filename))) for filename in self.file_list]
        time_idxs = np.argsort(digits)
        self.file_list = list(np.asarray(self.file_list)[time_idxs])

    def get_data(self):
        """
        Read the data from the csv files
        """
        data = np.genfromtxt(self.directory + "time_series/data.csv", delimiter=',', skip_header=1)
        self.time = data[:, 0]
        Fx = data[:, 1]
        Fy = data[:, 2]
        Fz = data[:, 3]
        tke = data[:, 4] #translational kinetic energy
        rke = data[:, 5] #rotational kinetic energy  

        return self.time, Fx, Fy, Fz, tke, rke   
    
    def set_reverse_flag(self):
        self.reverse_flag = True