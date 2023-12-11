import multiprocessing
import numpy as np
from DataProcessor import DataProcessor

class ProcessorDump(DataProcessor):
    def __init__(self, data):
        super().__init__(data)
        self.n_sim = self.data_reader.n_sim
        self.directory = self.data_reader.directory
        self.file_list = self.data_reader.file_list
        self.force = None
        self.point = None
        self.force_tangential = None
        self.shear = None
        self.area = None
        self.delta = None

    def process_data(self, num_processes):
        self.num_processes = num_processes
        with multiprocessing.Pool(self.num_processes) as pool:
            print("Started multiprocessing")
            results = pool.map(self.process_single_step,
                               [step for step in range(self.n_sim)])
        # Extract the averages from the results
        averages = np.array(results)
        return averages
    
    def process_single_step(self, step):
        """Processing on the data for one single step
        Calling the other methods for the processing"""

        data = np.loadtxt(self.directory+self.file_list[step], skiprows=9)
        self.point = data[:, 3:6] #contact_point
        self.force = data[:, 9:12] #contact_force
        self.force_tangential = data[:, 12:15] #contact_tangential_force
        self.shear = data[:, 15:18] #shear components of the contact force
        self.area = data[:, 19] #interpenetration area
        self.delta = data[:, 20] #interpenetration distance

        self.compute_space_averages()
        return np.concatenate((self.force_space_average,
                               self.force_tangential_space_average,
                               self.shear_space_average,
                               self.area_space_average.reshape(1,),
                               self.delta_space_average.reshape(1,)))

    def compute_space_averages(self):
        self.force_space_average = np.mean(self.force, axis=0)
        self.force_tangential_space_average = np.mean(self.force_tangential, axis=0)
        self.shear_space_average = np.mean(self.shear, axis=0)
        self.area_space_average = np.mean(self.area)
        self.delta_space_average = np.mean(self.delta)


