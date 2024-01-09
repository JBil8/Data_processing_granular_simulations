import os 
import re
import numpy as np
import vtk 
from DataReader import DataReader

class ReaderVtk(DataReader):
    def __init__(self, cof, ap, parameter, value):
        super().__init__(cof, ap, parameter, value)
        self.n_wall_atoms = None
        self.n_central_atoms = None

    #to call only after the data have been read
    def set_number_wall_atoms(self):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.directory + self.file_list[0])
        reader.Update()
        polydata = reader.GetOutput()
        coor = np.array(polydata.GetPoints().GetData())
        _, counts = np.unique(coor[:,1], axis=0, return_counts=True)
        self.n_wall_atoms = sum(counts[counts>4])

    def get_number_of_atoms(self):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.directory + self.file_list[0])
        reader.Update()
        polydata = reader.GetOutput()
        self.n_all_atoms = polydata.GetNumberOfPoints()

    def get_number_central_atoms(self):
        self.n_central_atoms = self.n_all_atoms - self.n_wall_atoms

    def get_initial_velocities(self):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.directory + self.file_list[0])
        reader.Update()
        polydata = reader.GetOutput()
        ids = np.array(polydata.GetPointData().GetArray(0))
        sorted_idxs = np.argsort(ids)
        self.v0 = np.array(polydata.GetPointData().GetArray(3))[sorted_idxs, :][self.n_wall_atoms:, :]
        self.y0 = np.array(polydata.GetPoints().GetData())[sorted_idxs, :][self.n_wall_atoms:, 1]
        
    def filter_relevant_files(self, prefix='shear_ellipsoids_'):
        self.file_list = [filename for filename in self.file_list if filename[len(prefix) + 1].isdigit() and filename[-1] == 'k']

    def read_data(self, global_path, prefix):
        self.prepare_data(global_path)
        self.filter_relevant_files(prefix)
        self.sort_files_by_time()
        self.get_number_of_time_steps()
        self.get_number_of_atoms()