import sys
import numpy as np
import pickle
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
from DataExporter import DataExporter

if __name__ == "__main__":

    
    parser = argparse.ArgumentParser(description='Process granular simulation.')
    parser.add_argument('-c', '--cof', type=float, help='coefficient of friction')
    parser.add_argument('-a', '--ap', type=float, help='aspect ratio')
    parser.add_argument('-t', '--type', type=str, help='simulation type: either I or phi')
    parser.add_argument('-v', '--value', type=float, help='packing fraction or Inertial number depensing on the type of simulation')
    parser.add_argument('-p', '--postprocess', action='store_false', help='whether to postprocess the data or simply import the pkl file, default is True', default=True)

    args = parser.parse_args()

    #parsing command line arguments
    cof = args.cof
    ap = args.ap
    param = args.value
    simulation_type = args.type
    full_postprocess = args.postprocess

    num_processes = 8

    if full_postprocess == True:

        if simulation_type == "I":
            #global_path = "/scratch/bilotto/simulations_inertial_number/parametric_studies/"
            global_path = "/home/jacopo/Documents/PhD_research/Liggghts_simulations/cluster_simulations/parametric_studies/"
        elif simulation_type == "phi":
            # global_path = "/scratch/bilotto/simulations_volume_fraction/parametric_studies/"
            global_path = "/home/jacopo/Documents/PhD_research/Liggghts_simulations/cluster_simulations/parametric_studies/"
        else:
            raise ValueError("simulation_type must be either I or phi") 
        #initialize the vtk reader
        data_read = ReaderVtk(cof, ap, simulation_type, param)
        data_read.read_data(global_path, 'shear_ellipsoids_')
        data_read.set_number_wall_atoms()
        data_read.get_number_central_atoms()
        data_read.get_initial_velocities()
        particles_volume = data_read.get_particles_volume()
        xz_surface = data_read.get_xz_surface()
        #intialize the dump reader
        data_dump = ReaderDump(cof, ap, simulation_type, param)
        data_dump.read_data(global_path, 'shear_ellipsoids_contact_data_')
        
        #plot force chains
        to_process_dump = ProcessorDump(data_dump, data_read.n_wall_atoms, data_read.n_central_atoms)
        to_process_dump.plot_force_chain(800, param)
        #plotter.plot_ellipsoids(0, averages_vtk)

    else:
        #import the data with pickle
        importer = DataExporter(ap, cof, simulation_type ,param)
        averages_vtk, averages_dump = importer.import_with_pickle()

        plotter = DataPlotter(ap, cof, simulation_type ,param)
        #plotter.plot_data(averages_vtk, averages_dump)
        #plotter.plot_eulerian_velocities(averages_vtk)

        print(averages_dump['trackedGrainsContactData'][0][:2, :])
        # plot ellipsois in 3d
        plotter.plot_ellipsoids(500, averages_vtk, averages_dump)