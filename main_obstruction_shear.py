import sys
import numpy as np
import pickle
import vtk
import math
import argparse
import matplotlib.pyplot as plt
#from scipy.spatial.transform import Rotation as R
import os
import re
import multiprocessing
import time as tm
from ProcessorEulerianVtk import ProcessorEulerianVtk   
from ProcessorDump import ProcessorDump
from DataPlotter import DataPlotter
from ReaderVtk import ReaderVtk
from ReaderDump import ReaderDump
from DataExporter import DataExporter
from Math_operations_grid import MathOperationsGrid

def process_single_step_wrapper(args):
    return to_process_vtk.process_single_step(args)

def compute_eulerian_msd(results_vtk, averages, n_sim):
    # Initialize dictionary to store cumulative square deviations
    square_deviations = {key: np.zeros_like(results_vtk[0][key]) for key in averages.keys()}

    # Compute cumulative square deviations
    for step_data in results_vtk:
        for key in averages.keys():
            deviation = step_data[key] - averages[key]
            squared_deviation = deviation ** 2
            square_deviations[key] += squared_deviation

    # Divide cumulative square deviations by n_sim to get mean square deviation
    mean_square_deviations = {key: value / n_sim for key, value in square_deviations.items()}
    return mean_square_deviations

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process granular simulation.')
    parser.add_argument('-c', '--cof', type=float, help='coefficient of friction particle-particle')
    parser.add_argument('-a', '--ap', type=float, help='aspect ratio')
    #parser.add_argument('-t', '--type', type=str, help='simulation type: either I or phi')
    parser.add_argument('-cw', '--cof_w', type=float, help='coefficient of friction particle-wall')
    #parser.add_argument('-val', '--value', type=float, help='packing fraction or Inertial number depensing on the type of simulation')
    parser.add_argument('-vw', '--vwall', type=float, help='wall velocity')
    parser.add_argument('-p', '--postprocess', action='store_false', help='whether to postprocess the data or simply import the pkl file, default is True', default=True)
    parser.add_argument('-f', '--fraction_obstruction', type=int, help='fraction of the obstruction in the channel', default=2)

    args = parser.parse_args()

    #parsing command line arguments
    cof = args.cof
    ap = args.ap
    muw = args.cof_w
    vwall = args.vwall
    full_postprocess = args.postprocess
    fraction_obstruction = args.fraction_obstruction
    global_path = "/home/jacopo/Documents/PhD_research/Liggghts_simulations/test_simulations/shearing_stl/"

    num_processes = 8

    if full_postprocess == True:

        plt.ioff()
        #initialize the vtk reader
        nx_divisions = 40
        data_read = ReaderVtk(cof, ap, muw=muw, vwall=vwall, fraction = fraction_obstruction)
        data_read.read_data(global_path, 'shearing_')
        data_read.no_wall_atoms()
        data_read.get_number_central_atoms()
        data_read.get_initial_velocities()
        particles_volume = data_read.get_particles_volume()
        data_read.get_box_dimensions()
        ny_divisions, dx, dy = data_read.make_2D_grid(nx_divisions = nx_divisions)
        n_sim = data_read.n_sim
        box_volume = data_read.box_x*data_read.box_y*data_read.box_z
        global_phi = particles_volume/box_volume
        #intialize the dump reader

    #     data_dump = ReaderDump(cof, ap, simulation_type, param)
    #     data_dump.read_data(global_path, 'compacted_contact_data_')
        to_process_vtk = ProcessorEulerianVtk(data_read)
        #to_process_vtk.process_single_step(100)
    #     to_process_dump = ProcessorDump(data_dump, data_read.n_wall_atoms, data_read.n_central_atoms)
        with multiprocessing.Pool(num_processes) as pool:
            print("Started multiprocessing")
           
            # results_vtk = pool.map(to_process_vtk.process_single_step,
            #                        [step for step in range(to_process_vtk.n_sim)])
            
            # Apply asynchronously the function to each argument
            results = [pool.apply_async(to_process_vtk.process_single_step,
                                         (step, nx_divisions, ny_divisions, dx, dy))
                        for step in range(to_process_vtk.n_sim)]
            
            # Wait for all processes to finish
            pool.close()
            pool.join()
            #results_vtk = [result.get() for result in results]

                # Collect results
            results_vtk = []
            max_dicts = []
            min_dicts = []
            for result in results:
                grid_dict, max_values, min_values = result.get()
                results_vtk.append(grid_dict)
                max_dicts.append(max_values)
                min_dicts.append(min_values)

    #         results_dump = pool.map(to_process_dump.process_single_step,
    #                             [step for step in range(to_process_dump.n_sim)])
        #save gif from plots genrated by the vtk processor
        #to_process_vtk.save_gif_from_plots()    
    #        
    #Extract the averages from the results
        averages = {'v': np.zeros_like(results_vtk[0]['v']),
                    'F': np.zeros_like(results_vtk[0]['F']),
                    'phi': np.zeros_like(results_vtk[0]['phi']), 
                    'theta_x' : np.zeros_like(results_vtk[0]['theta_x']),
                    'theta_z' : np.zeros_like(results_vtk[0]['theta_z']),
                    'stress': np.zeros_like(results_vtk[0]['stress'])}

            # Compute cumulative sums and counts
        for step_data in results_vtk:
            for key in averages.keys():
                step_data[key] = np.nan_to_num(step_data[key], nan=0.0)
                averages[key] += step_data[key]  # Add values to cumulative sum
                
       # Compute averages
        averages = {key: value / n_sim for key, value in averages.items()}

        # find max and min values of the quantities
        max_values = {key: np.max([step_data[key] for step_data in results_vtk]) for key in averages.keys()}
        min_values = {key: np.min([step_data[key] for step_data in results_vtk]) for key in averages.keys()}

        # compute the mean square deviation
        msd = compute_eulerian_msd(results_vtk, averages, n_sim)

        # compute divergence and curl of velocity field
        math_operations = MathOperationsGrid(nx_divisions, ny_divisions, dx, dy)
        div_v = math_operations.compute_divergence(averages['v'][:, :, 0], averages['v'][:, :, 1])
        pressure = -(averages['stress'][:,:,0]+averages['stress'][:,:,1] + averages['stress'][:,:,2])/3

    plotter = DataPlotter(ap, cof, muw = muw, vwall =vwall, fraction=fraction_obstruction)
    plotter.plot_space_averages_all_cells(averages['v'], quantity='Velocity', component=0, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['F'], quantity='Force', component=0, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['v'], quantity='Velocity', component=1, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['F'], quantity='Force', component=1, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['v'], quantity='Velocity', component=2, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['F'], quantity='Force', component=2, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    
    plotter.plot_space_averages_all_cells(msd['v'], quantity='MSD_VEL', component=0, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(msd['v'], quantity='MSD_VEL', component=1, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(msd['v'], quantity='MSD_VEL', component=2, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(msd['F'], quantity='MSD_F', component=0, nx_divisions=nx_divisions, ny_divisions=ny_divisions)

    plotter.plot_space_averages_all_cells(averages['stress']/box_volume, quantity='Stress $\sigma_{xx}$', component=0, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['stress']/box_volume, quantity='Stress $\sigma_{yy}$', component=1, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['stress']/box_volume, quantity='Stress $\sigma_{zz}$', component=2, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['stress']/box_volume, quantity='Stress $\sigma_{xy}$', component=3, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['stress']/box_volume, quantity='Stress $\sigma_{xz}$', component=4, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['stress']/box_volume, quantity='Stress $\sigma_{yz}$', component=5, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['stress']/box_volume, quantity='Stress $\sigma_{yz}$', component=5, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(pressure/box_volume, quantity='Pressure', component=None, nx_divisions=nx_divisions, ny_divisions=ny_divisions)


    plotter.plot_space_averages_all_cells(div_v, quantity='Div $\\nabla \cdot u$', component=None, nx_divisions=nx_divisions, ny_divisions=ny_divisions)

    plotter.plot_space_averages_all_cells(averages['phi'], quantity='Volume fraction', component=None, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['theta_x'], quantity='Theta_x', component=None, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    plotter.plot_space_averages_all_cells(averages['theta_z'], quantity='Thetat_z', component=None, nx_divisions=nx_divisions, ny_divisions=ny_divisions)
    averages['v_shearing'] = data_read.v_shearing
    #         averages_vtk['phi'] = particles_volume/(xz_surface*averages_vtk['box_height'])  
    #         averages_vtk.pop('box_height')

    #         averages_dump = {}
    #         for key in results_dump[0].keys():
    #             if key == "trackedGrainsContactData":
    #                 # list of multidimensional arrays for contact data
    #                 averages_dump[key] = [result[key] for result in results_dump]
    #             else:
    #                 averages_dump[key] = np.array([result[key] for result in results_dump])

    #     #export the data with pickle
    #     exporter = DataExporter(ap, cof, simulation_type ,param)
    #     exporter.export_with_pickle(averages_vtk, averages_dump)

    #     plotter = DataPlotter(ap, cof, simulation_type ,param)
    #     plotter.plot_data(averages_vtk, averages_dump, particles_volume, xz_surface)
    #     plotter.plot_eulerian_velocities(averages_vtk)

    #     #force distrubution
        
    #     with multiprocessing.Pool(num_processes) as pool:
    #         print("Started multiprocessing")
    #         results = pool.map(to_process_dump.force_single_step,
    #                             [step for step in range(to_process_dump.n_sim)])
    #     # stack all the forces and force tangential to compute the distribution
    #     force_normal_stack = np.concatenate([result['force_normal'] for result in results])
    #     force_tangential_stack = np.concatenate([result['force_tangential'] for result in results])

    #     force_normal_distribution = np.histogram(force_normal_stack, bins=100, density=True)
    #     force_tangential_distribution = np.histogram(force_tangential_stack, bins=100, density=True)

    #     plotter.plot_force_distribution(force_normal_distribution, force_tangential_distribution)

    #     #export force distribution data
    #     exporter.export_force_distribution(force_normal_distribution, force_tangential_distribution)

    #     #analysis of the tracked grains contact data

    #     #plot the ellipsoids in 3d
        
    #     #plotter.plot_ellipsoids(0, averages_vtk)

    #     plt.ion()

    # else:
    #     #import the data with pickle
    #     importer = DataExporter(ap, cof, simulation_type ,param)
    #     averages_vtk, averages_dump = importer.import_with_pickle()

    #     plotter = DataPlotter(ap, cof, simulation_type ,param)
    #     #plotter.plot_data(averages_vtk, averages_dump)
    #     #plotter.plot_eulerian_velocities(averages_vtk)

    #     print(averages_dump['trackedGrainsContactData'][0][:2, :])
    #     # plot ellipsois in 3d
    #     plotter.plot_ellipsoids(500, averages_vtk, averages_dump)