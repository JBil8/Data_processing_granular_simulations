import sys
import numpy as np
import vtk
import argparse
import matplotlib.pyplot as plt
#from scipy.spatial.transform import Rotation as R
import os
import re
import multiprocessing
import time as tm
from ProcessorVtk import ProcessorVtk
from ProcessorDump import ProcessorDump
from DataPlotter import DataPlotter
from ReaderVtk import ReaderVtk
from ReaderDump import ReaderDump

if __name__ == "__main__":

    
    parser = argparse.ArgumentParser(description='Process granular simulation.')
    parser.add_argument('-c', '--cof', type=float, help='coefficient of friction')
    parser.add_argument('-a', '--ap', type=float, help='aspect ratio')
    parser.add_argument('-t', '--type', type=str, help='simulation type: either I or phi')
    parser.add_argument('-v', '--value', type=float, help='packing fraction or Inertial number depensing on the type of simulation')
    args = parser.parse_args()

    #parsing command line arguments
    cof = args.cof
    ap = args.ap
    param = args.value
    simulation_type = args.type

    # Define the directory where your files are located
    #global_path = "/scratch/bilotto/simulations_volume_fraction/parametric_studies/" #phi_0.6_ap_1.5_cof_0.4_v_5.6"
    # Parameter that vary in the parametric study
    #cof_list = ['0.0', '0.4','1.0', '10.0']
    #ap_list = ['1.0', '1.5', '2.0', '2.5', '3.0']
    #I_list = ['0.0001', '0.000501', '0.0025', '0.013', '0.063', '0.32']
    #phi_list = ['0.5', '0.6', '0.7', '0.8', '0.9']

    num_processes = 8

    if simulation_type == "I":
        global_path = '/home/jacopo/Documents/PhD_research/Liggghts_simulations/cluster_simulations/simulations_simple_shear/parametric_studies/'
    elif simulation_type == "phi":
        global_path = '/home/jacopo/Documents/PhD_research/Liggghts_simulations/cluster_simulations/parametric_studies/'
    else:
        raise ValueError("simulation_type must be either I or phi") 
    plt.ioff()
    #initialize the vtk reader
    data_read = ReaderVtk(cof, ap, simulation_type, param)
    data_read.read_data(global_path, 'shear_ellipsoids_')
    data_read.set_number_wall_atoms()
    data_read.get_number_central_atoms()
    data_read.get_initial_velocities()
    #intialize the dump reader
    data_dump = ReaderDump(cof, ap, simulation_type, param)
    data_dump.read_data(global_path, 'shear_ellipsoids_contact_data_')
    to_process_vtk = ProcessorVtk(data_read)
    to_process_dump = ProcessorDump(data_dump, data_read.n_wall_atoms, data_read.n_central_atoms)
    
    with multiprocessing.Pool(num_processes) as pool:
        print("Started multiprocessing")
        results_vtk = pool.map(to_process_vtk.process_single_step,
                            [step for step in range(to_process_vtk.n_sim)])
    #     results_dump = pool.map(to_process_dump.process_single_step,
    #                         [step for step in range(to_process_dump.n_sim)])
        
        #Extract the averages from the results
        averages_vtk = {}
        
        for key in results_vtk[0].keys():
            if key == 'trackedGrainsOrientation' or key == 'trackedGrainsPosition':
            # Stack the arrays along a new axis to preserve the structure
                averages_vtk[key] = np.stack([result[key] for result in results_vtk], axis=0)
            else:
                averages_vtk[key] = np.array([result[key] for result in results_vtk])
  
    #     averages_dump = {}
    #     for key in results_dump[0].keys():
    #         if key == "trackedGrainsContactData":
    #             # list of multidimensional arrays for contact data
    #             averages_dump[key] = [result[key] for result in results_dump]
    #         else:
    #             averages_dump[key] = np.array([result[key] for result in results_dump])


    plotter = DataPlotter(ap, cof, simulation_type ,param)
    # plotter.plot_data(averages_vtk, averages_dump)
    # plotter.plot_eulerian_velocities(averages_vtk)

    # #force distrubution
    
    # with multiprocessing.Pool(num_processes) as pool:
    #     print("Started multiprocessing")
    #     results = pool.map(to_process_dump.force_single_step,
    #                         [step for step in range(to_process_dump.n_sim)])
    # # stack all the forces and force tangential to compute the distribution
    # force_normal_stack = np.concatenate([result['force_normal'] for result in results])
    # force_tangential_stack = np.concatenate([result['force_tangential'] for result in results])

    # force_normal_distribution = np.histogram(force_normal_stack, bins=100, density=True)
    # force_tangential_distribution = np.histogram(force_tangential_stack, bins=100, density=True)

    # plotter.plot_force_distribution(force_normal_distribution, force_tangential_distribution)

    #analysis of the tracked grains contact data

    #plot the ellipsoids in 3d
    
    plotter.plot_ellipsoids(0, averages_vtk)

    plt.ion()

