import numpy as np
import pickle as pkl

class DataExporter:
    def __init__(self, ap , cof, parameter, value):
        self.ap = str(ap) 
        self.cof = str(cof)
        self.parameter = parameter
        self.value = str(value)
    
    def export_with_pickle(self, data_vtk, data_dump):
        # export the data with pickle for further analysis with appropriate name
        with open('output_data/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + '.pkl', 'wb') as f:
            pkl.dump(data_vtk, f)
        with open('output_data/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + '_dump.pkl', 'wb') as f:
            pkl.dump(data_dump, f)

    def import_with_pickle(self):
        # import the data with pickle for further analysis with appropriate name
        with open('output_data/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + '.pkl', 'rb') as f:
            data_vtk = pkl.load(f)
        with open('output_data/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + '_dump.pkl', 'rb') as f:
            data_dump = pkl.load(f)
        return data_vtk, data_dump

    def export_data(self, data_vtk, data_dump):
        # export data to csv files for further analysis
        # export the data from data_vtk and data_dump
        strain = np.arange(1, data_vtk['theta_x'].shape[0]+1)*20/data_vtk['theta_x'].shape[0]
        
        #combine the two dictionaries into one with all the keys
        data_vtk.update(data_dump)
        all_data = [strain]
        key_list = []
        #export the data to csv files for all the keys except the last one in the dictionary
        for key in data_vtk.keys() - ['trackedGrainsContactData']:
            if data_vtk[key].ndim == 1:
                key_list.append(key)
                # for 1D arrays, add the data to the list
                all_data.append(data_vtk[key])

        np.savetxt(
            'output_data/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + '.csv',
            np.transpose(all_data),
            delimiter=',',
            header=','.join(['strain'] + key_list)
        )

        # export eulerian velocities as csv with rows for each time step
        np.savetxt(
            'output_data/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + '_eulerian_velocities.csv',
            data_vtk['eulerian_vx'],
            delimiter=','
        )
        

        # export the tracked grains contact data