import sys
import numpy as np
import vtk
import matplotlib.pyplot as plt
#from scipy.spatial.transform import Rotation as R
import os
import re
import multiprocessing
import time as tm
from DataProcessor import DataProcessor 
from DataReader import DataReader
from DataPlotter import DataPlotter

if __name__ == "__main__":

    # Define the directory where your files are located
    global_path = '/home/jacopo/Documents/PhD_research/Liggghts_simulations/cluster_simulations/simulations_simple_shear/parametric_studies/'

    # Parameter that vary in the parametric study
    #cof_list = ['0.0', '0.4','1.0', '10.0']
    #ap_list = ['1.0', '1.5','2', '2.5', '3.0']
    I_list = ['0.0001', '0.000501', '0.0025', '0.013', '0.063', '0.32']
    cof_list = ['0.4']
    ap_list = ['3.0']
    
    n_param = len(cof_list)*len(ap_list)*len(I_list)
    num_processes = 6

    # data_read = DataReader(cof_list[1], ap_list[1], I_list[4])
    # data_read.set_number_wall_atoms(1062)
    # data_read.read_data(global_path, 'shear_ellipsoids_')
    # to_process = DataProcessor(data_read)
    # seq = to_process.process_data(num_processes)
    # plotting = DataPlotter(seq)
    # fig = plotting.plot_data()
    # fig.savefig('simple_shear_ap' + ap + '_cof_' + cof + '_I_' + I + '.png')
    plt.ioff()
    for ap in ap_list:
        for cof in cof_list:
            for I in I_list:
                data_read = DataReader(cof, ap, I)
                data_read.set_number_wall_atoms(1062)
                data_read.read_data(global_path, 'shear_ellipsoids_')
                to_process = DataProcessor(data_read)
                seq = to_process.process_data(num_processes)
                plotting = DataPlotter(seq)
                fig = plotting.plot_data()
                fig.suptitle('ap = ' + ap + ', cof = ' + cof + ', I = ' + I)
                fig.savefig('simple_shear_ap' + ap + '_cof_' + cof + '_I_' + I + '.png')
    plt.ion()
    #eulerian = processed.eulerian_velocity(10)
    #trajectory = processed.compute_single_particle_trajectory(10, 2)
    #averages = processed.process_data()
    #print(eulerian)
    #print(trajectory)
    # param_averages = np.zeros((300, 10, len(I_list), len(ap_list)))
    # for id2, ap in enumerate(ap_list):
    #     for id1, I in enumerate(I_list):
    #         processor = DataProcessor(cof_list[2], ap, I, n_sim, num_processes)
    #         averages = processor.process_data()
    #         param_averages[:n_sim, :, id1, id2] = np.array(averages)