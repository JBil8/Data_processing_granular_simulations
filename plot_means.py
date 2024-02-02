import numpy as np
from matplotlib import pyplot as plt
import pickle as pkl
import argparse
import os
from DataExporter import DataExporter
from matplotlib.ticker import NullFormatter
from matplotlib.font_manager import FontProperties
from matplotlib.gridspec import GridSpec

parser = argparse.ArgumentParser(description='Process granular simulation.')
parser.add_argument('-c', '--cof', type=float, help='coefficient of friction')
parser.add_argument('-t', '--type', type=str, help='simulation type: either I or phi')

args = parser.parse_args()

#parsing command line arguments
cof = args.cof
simulation_type = args.type
aspect_ratios = [1.0, 1.5, 2.0, 2.5, 3.0]
#aspect_ratios = [1.5]

if simulation_type == "I":
    values = [0.1, 0.0398, 0.0158, 0.0063, 0.0025, 0.001]
    #values = [0.1, 0.0398, 0.0158, 0.0063]
    

elif simulation_type == "phi":
    values = [0.5, 0.6, 0.7, 0.8, 0.9]

# loop over the values of the parameter

avg_dict = {'Z': np.zeros((len(aspect_ratios), len(values))),
            'omega_z': np.zeros((len(aspect_ratios), len(values))),
            'theta_x': np.zeros((len(aspect_ratios), len(values))),
            'percent_aligned': np.zeros((len(aspect_ratios), len(values))),
            'S2': np.zeros((len(aspect_ratios), len(values))),
            'box_height': np.zeros((len(aspect_ratios), len(values))),
            'autocorrelation_v': np.zeros((len(aspect_ratios), len(values))),
            'mu_effective': np.zeros((len(aspect_ratios), len(values))),
            'msd': np.zeros((len(aspect_ratios), len(values))), 
            'eulerian_vx': np.zeros((10, len(aspect_ratios), len(values))), 
            'vel_fluct': np.zeros((10, len(aspect_ratios), len(values)))}

for j, value in enumerate(values):
    # loop over the aspect ratios
    for i, ap in enumerate(aspect_ratios):
    # open pickel file from output_data
        with open('output_data/simple_shear_ap' + str(ap) + '_cof_' + str(cof) + '_' + simulation_type + '_' + str(value) + '.pkl', 'rb') as f:
            data_vtk = pkl.load(f)   

        with open('output_data/simple_shear_ap' + str(ap) + '_cof_' + str(cof) + '_' + simulation_type + '_' + str(value) + '_dump.pkl', 'rb') as f:
            data_dump = pkl.load(f)         


        msd = data_vtk['msd']
        indx = np.where(msd == 0.0)[0][0]
        msd = msd[indx:]
        strain = np.arange(0, msd.shape[0])*18/msd.shape[0]

        # Perform the curve fit
        slope = np.polyfit(strain, msd, 1)[0]/2 # divide by 2 to get the diffusion coefficient in 1D

        # average the data over the last 50% of the simulation (after strain = 10)
        avg_dict['Z'][i,  j] = (np.mean(data_dump['Z'][int(data_dump['Z'].shape[0]/2):]))
        avg_dict['omega_z'][i,  j] = (np.mean(data_vtk['omega_z'][int(data_vtk['omega_z'].shape[0]/2):]))
        avg_dict['theta_x'][i,  j] = np.rad2deg(np.mean(data_vtk['theta_x'][int(data_vtk['theta_x'].shape[0]/2):]))
        avg_dict['percent_aligned'][i,  j] = (np.mean(data_vtk['percent_aligned'][int(data_vtk['percent_aligned'].shape[0]/2):]))
        avg_dict['S2'][i,  j] = (np.mean(data_vtk['S2'][int(data_vtk['S2'].shape[0]/2):]))
        avg_dict['box_height'][i,  j] = (np.mean(data_vtk['box_height'][int(data_vtk['box_height'].shape[0]/2):]))
        avg_dict['autocorrelation_v'][i,  j] = (np.mean(data_vtk['autocorrelation_v'][int(data_vtk['autocorrelation_v'].shape[0]/2):]))
        avg_dict['mu_effective'][i,  j] = (np.mean(data_vtk['mu_effective'][int(data_vtk['mu_effective'].shape[0]/2):]))
        avg_dict['msd'][i,  j] = (slope)
        avg_dict['eulerian_vx'][:, i,  j] = (np.mean(data_vtk['eulerian_vx'][int(data_vtk['eulerian_vx'].shape[0]/2):, :], axis=0))
        avg_dict['vel_fluct'][:, i,  j] = (np.mean(
            (data_vtk['eulerian_vx'][int(data_vtk['eulerian_vx'].shape[0]/2):, :]-avg_dict['eulerian_vx'][:, i,  j])**2, axis=0))

# avg_dict['eulerian_vx'].shape == (10, 5, 6)

font_properties = FontProperties(fname="/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", size=12)

# Plot the data for all values of the parameter on the same plot
fig, axes = plt.subplots(3, 3, figsize=(12, 12), sharex='col', gridspec_kw={'hspace': 0.2, 'wspace': 0.5})
cmap = plt.get_cmap('viridis', len(aspect_ratios))

# Flatten the 2D array of subplots into a 1D array
axes = axes.flatten()


#ylabels = ['Z', '$\omega_z$', '\theta_x', 'percent_aligned', 'S_2', 'box_height [m]', 'autocorrelation_v', '\mu_{eff}', 'D_y [m^2]']
ylabels = ['$Z$', '$\omega_z$', '$\\theta_x$', '$\%$ aligned', '$S_2$', '$h$ [m]', '$C_v$', '$\mu_{eff}$', '$D_y$ [m$^2$]']

# Iterate over properties and plot on each subplot
for property_index, property_name in enumerate(['Z', 'omega_z', 'theta_x', 'percent_aligned', 'S2', 'box_height', 'autocorrelation_v', 'mu_effective', 'msd']):
    ax = axes[property_index]  # Get the current subplot

    # Iterate over aspect ratios and plot on the same subplot
    for i, ap in enumerate(aspect_ratios):
        ax.semilogx(values, avg_dict[property_name][i, :], label=f'ap = {ap}', color=cmap(i))

    # Use text function to add LaTeX-formatted y-axis label
    ax.text(-0.3, 0.5, ylabels[property_index], rotation="vertical", va="center", ha="center", transform=ax.transAxes, fontproperties=font_properties)

    # Hide x-axis values for some subplots
    if property_index // 3 != 2:
        ax.xaxis.set_major_formatter(NullFormatter())

# Set xlabel only on the bottom plot of each column
for ax in axes[-3:]:
    if simulation_type == "I":
        ax.set_xlabel('$I$')
    elif simulation_type == "phi":
        ax.set_xlabel('$\phi$')
# Customize the plot (labels, legend, etc.)
axes[0].legend()

fig.suptitle('cof = ' + str(cof))
fig.savefig('output_plots/parametric_plots/simple_shear_ap' + str(ap) + '_cof_' + str(cof) + '_' + simulation_type + '_all_values.png')





# plot average eulerian velocity
y = np.linspace(0, 1, 10)

# Plot the data for all values of the parameter on the same plot
fig, axes = plt.subplots(2, len(values), figsize=(12, 50), sharey='row', gridspec_kw={'hspace': 0.2, 'wspace': 0.5})
cmap = plt.get_cmap('viridis', len(aspect_ratios))

# Flatten the 2D array of subplots into a 1D array
axes = axes.flatten()

for j, value in enumerate(values):
    # Iterate over aspect ratios and plot on the same subplot
    for i, ap in enumerate(aspect_ratios):
        axes[j].plot(avg_dict['eulerian_vx'][:, i, j], y, label=f'ap = {ap}', color=cmap(i))
        axes[j].set_xlabel('$<V_x>$')
        axes[j+len(values)].plot(avg_dict['vel_fluct'][:, i, j], y, label=f'ap = {ap}', color=cmap(i))
        axes[j+len(values)].set_xlabel('$(V-<V_x>)^2$')
        
    axes[j].set_title(f'{simulation_type} = {value}')

axes[0].set_ylabel('$y/H$')
axes[len(values)].set_ylabel('$y/H$')
axes[0].legend()
fig.suptitle('cof = ' + str(cof))
plt.show()
fig.savefig('output_plots/parametric_plots/simple_shear_ap' + str(ap) + '_cof_' + str(cof) + '_' + simulation_type + '_eulerian_statistics.png')