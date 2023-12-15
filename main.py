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
    parser.add_argument('-t', '--type', type=str, help='simulation type')
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
        results_dump = pool.map(to_process_dump.process_single_step,
                            [step for step in range(to_process_dump.n_sim)])
    # Extract the averages from the results
        averages_vtk = np.array(results_vtk)
        averages_dump = np.array(results_dump)
    
    #to_process_dump.process_single_step(200)
    #to_process_dump.compute_force_distribution()
    #to_process_dump.plot_force_chain(200)

    # plt.figure()
    # plt.subplot(2,2,1)
    # plt.plot(to_process_dump.force_distribution_x[1][:-1], to_process_dump.force_distribution_x[0])
    # plt.xlabel('Fx')
    # plt.subplot(2,2,2)
    # plt.plot(to_process_dump.force_distribution_y[1][:-1], to_process_dump.force_distribution_y[0])
    # plt.xlabel('Fy')
    # plt.subplot(2,2,3)
    # plt.plot(to_process_dump.force_distribution_z[1][:-1], to_process_dump.force_distribution_z[0])
    # plt.xlabel('Fz')
    # plt.subplot(2,2,4)
    # plt.plot(to_process_dump.force_distribution[1][:-1], to_process_dump.force_distribution[0])
    # plt.xlabel('F')
    # plt.show()
    #seq_dump = to_process_dump.process_data(num_processes)
    plotter = DataPlotter(averages_vtk, averages_dump, ap, cof, simulation_type ,param)
    plotter.plot_data()
    
    # to implement: eulerian velocities and contact network


    plt.ion()