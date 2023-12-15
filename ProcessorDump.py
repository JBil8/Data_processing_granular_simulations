import multiprocessing
import numpy as np
from DataProcessor import DataProcessor

class ProcessorDump(DataProcessor):
    def __init__(self, data, n_wall_atoms, n_central_atoms):
        super().__init__(data)
        self.n_sim = self.data_reader.n_sim
        self.directory = self.data_reader.directory
        self.file_list = self.data_reader.file_list
        self.n_wall_atoms = n_wall_atoms
        self.n_central_atoms = n_central_atoms
        self.force = None
        self.point = None
        self.force_tangential = None
        self.shear = None
        self.area = None
        self.delta = None
        #self.n_wall_atoms = self.data_reader.n_wall_atoms

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
        #excluding contacts between wall particles
        first_col_check = (data[:,6] > self.n_wall_atoms)
        second_col_check = (data[:,7] > self.n_wall_atoms)
        not_wall_contacts = np.logical_and(first_col_check, second_col_check)
        data = data[not_wall_contacts]  
        self.point = data[:, 19:22] #contact_point
        self.force = data[:, 9:12] #contact_force
        self.force_tangential = data[:, 12:15] #contact_tangential_force
        self.shear = data[:, 15:18] #shear components of the contact force
        self.area = data[:, 18] #interpenetration area
        self.delta = data[:, 22] #interpenetration distance

        self.compute_space_averages()
        return np.concatenate((self.force_space_average,
                               self.force_tangential_space_average,
                               self.shear_space_average,
                               self.area_space_average.reshape(1,),
                               self.delta_space_average.reshape(1,), 
                               self.contact_number.reshape(1,)))
    
    def compute_force_distribution(self):
        """Compute the force distribution"""
        self.force_distribution_x = np.histogram(np.abs(self.force[:, 0]), bins=100)
        self.force_distribution_y = np.histogram(np.abs(self.force[:, 1]), bins=100)
        self.force_distribution_z = np.histogram(np.abs(self.force[:, 2]), bins=100)
        self.force_distribution = np.histogram(np.linalg.norm(self.force, axis=1), bins=100)

    def compute_force_tangential_distribution(self):
        """Compute the tangential force distribution"""
        self.force_tangential_distribution_x = np.histogram(self.force_tangential[:, 0], bins=100)
        self.force_tangential_distribution_y = np.histogram(self.force_tangential[:, 1], bins=100)
        self.force_tangential_distribution_z = np.histogram(self.force_tangential[:, 2], bins=100)
        self.force_tangential_distribution = np.histogram(np.linalg.norm(self.force_tangential, axis=0), bins=100)

    def compute_space_averages(self):
        self.force_space_average = np.mean(self.force, axis=0)
        self.force_tangential_space_average = np.mean(self.force_tangential, axis=0)
        self.shear_space_average = np.mean(self.shear, axis=0)
        self.area_space_average = np.mean(self.area)
        self.delta_space_average = np.mean(np.abs(self.delta))
        self.contact_number = len(self.area)/self.n_central_atoms*2
        
    def plot_force_chain(self, step):
        """Plot the force chain"""
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        data = np.loadtxt(self.directory+self.file_list[step], skiprows=9)
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection='3d')
        print(self.n_wall_atoms)
        for i in range(len(data)):
            if data[i, 6] > self.n_wall_atoms and data[i, 7] > self.n_wall_atoms and data[i, 8] != 1:
                ax.plot([data[i, 0], data[i,3]], [data[i, 2], data[i,5]], [data[i, 1], data[i,4]], linestyle='-', color='k')
        # Set equal axis scaling
        ax.set_box_aspect([np.ptp(coord) for coord in [ax.get_xlim(), ax.get_ylim(), ax.get_zlim()]])
    
        ax.set_xlabel('X ')
        ax.set_ylabel('Z ')
        ax.set_zlabel('Y ')
        ax.set_title('Force chains, ap='+str(self.data_reader.ap)+', cof='+str(self.data_reader.cof)+', phi='+str(self.data_reader.parameter)+', step='+str(step))

        plt.savefig('output_plots/force_chain'+str(step)+'.png')
        plt.show()
        plt.close()