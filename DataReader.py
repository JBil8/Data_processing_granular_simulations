import os 
import re
import numpy as np
import vtk 

class DataReader:
    def __init__(self, cof, ap, I):
        self.cof = cof
        self.ap = ap
        self.I = I
        self.n_wall_atoms = None
        self.n_sim = None
        self.step = None
        self.directory = None
        self.file_list = None
        self.n_central_atoms = None

    def set_number_wall_atoms(self, n_wall_atoms):
        self.n_wall_atoms = n_wall_atoms

    def get_number_of_time_steps(self):
        digits = [int(''.join(re.findall(r'\d+', filename))) for filename in self.file_list]
        self.n_sim = len(self.file_list)
        self.step = int((max(digits) - min(digits)) / (self.n_sim - 1))

    def prepare_data(self, global_path):
        self.directory = global_path + f'ap_{self.ap}_cof_{self.cof}_I_{self.I}/'
        self.file_list = os.listdir(self.directory)

    def filter_relevant_files(self, prefix='shear_ellipsoids_'):
        self.file_list = [filename for filename in self.file_list if filename[len(prefix) + 1].isdigit()]
    
    def sort_files_by_time(self):
        digits = [int(''.join(re.findall(r'\d+', filename))) for filename in self.file_list]
        time_idxs = np.argsort(digits)
        self.file_list = list(np.asarray(self.file_list)[time_idxs])

    def get_number_of_atoms(self):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.directory + self.file_list[0])
        reader.Update()
        polydata = reader.GetOutput()
        self.n_all_atoms = polydata.GetNumberOfPoints()

    def get_number_central_atoms(self):
        self.n_central_atoms = self.n_all_atoms - self.n_wall_atoms

    def get_velocities(self):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.directory + self.file_list[0])
        reader.Update()
        polydata = reader.GetOutput()
        ids = np.array(polydata.GetPointData().GetArray(0))
        sorted_idxs = np.argsort(ids)
        self.v0 = np.array(polydata.GetPointData().GetArray(3))[sorted_idxs, :][self.n_wall_atoms:, :]

    def read_data(self, global_path, prefix):
        self.prepare_data(global_path)
        self.filter_relevant_files(prefix)
        self.sort_files_by_time()
        self.get_number_of_time_steps()
        self.get_number_of_atoms()
        self.get_number_central_atoms()
        self.get_velocities()