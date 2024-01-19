import numpy as np
import matplotlib.pyplot as plt

class DataPlotter:
    def __init__(self, ap , cof, parameter, value):
        self.ap = str(ap) 
        self.cof = str(cof)
        self.parameter = parameter
        self.value = str(value)
    
    def plot_data(self, data_vtk, data_dump):
        strain = np.arange(1, data_vtk['theta_x'].shape[0]+1)*20/data_vtk['theta_x'].shape[0]
        
        fig = plt.figure(figsize=(30, 10))
        ax = fig.add_subplot(331)
        plt.plot(strain, data_dump['Z'])
        plt.ylabel('Z')
        ax = fig.add_subplot(332)
        plt.plot(strain, data_vtk['omega_z'])
        plt.ylabel('Omegaz [rad/s]')
        ax = fig.add_subplot(333)
        plt.plot(strain, -data_vtk['F_x'])
        plt.ylabel('Shearing force [N]')
        ax = fig.add_subplot(334)
        plt.plot(strain, np.rad2deg(data_vtk['theta_x'])) #convert to degrees
        plt.ylabel('theta_x [deg]')
        ax = fig.add_subplot(335)
        plt.plot(strain, data_vtk['percent_aligned']) #convert to degrees
        plt.ylabel('aligned particles [%]')
        ax = fig.add_subplot(336)
        plt.plot(strain, data_vtk['S2'])
        plt.ylabel('S2')
        ax = fig.add_subplot(337)
        plt.plot(strain, data_vtk['autocorrelation_v'])
        plt.ylabel('Autocorrelation V')
        plt.xlabel('Strain')

        ax = fig.add_subplot(338)
        plt.plot(strain, data_vtk['mu_effective'])
        plt.ylabel('Effecive \mu')
        # plt.plot(strain, self.data_dump[:,10])
        # plt.ylabel('Average overlap in contacts [m]')

        ax = fig.add_subplot(339)
        plt.plot(strain, data_vtk['msd'])
        plt.ylabel('Mean square displacement [m^2]')

        fig.suptitle('ap = ' + self.ap + ', cof = ' + self.cof + ', ' + self.parameter +  '=' + self.value)
        fig.savefig('output_plots/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + '.png')
        
        plt.clf()
    
    def plot_eulerian_velocities(self, data_vtk):
        #plot eulerian velocities at n_plots time steps in time
        n_plots = 10
        n_time_steps = data_vtk['theta_x'].shape[0]
        time_interval = int(n_time_steps/n_plots)
        y = np.linspace(0, 1, 10)
        fig2 = plt.figure(figsize=(10, 20))
        for i in range(n_plots):
            #ax = fig2.add_subplot(2, int(n_plots/2), i+1)
            color = plt.cm.viridis(i / n_plots)  # color will now be an RGBA tuple
            plt.plot(data_vtk['eulerian_vx'][i*time_interval,:], y, label = f'strain  = {i*time_interval*20/n_time_steps:.2f}', color=color)
        plt.xlabel('Vx/V')
        plt.ylabel('y/H')
        plt.legend()
        fig2.suptitle('ap = ' + self.ap + ', cof = ' + self.cof + ', ' + self.parameter +  '=' + self.value)
        fig2.savefig('output_plots/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + 'eulerian.png')
        plt.clf()
         #plt.show()
        
    def plot_force_distribution(self, force_normal_distribution, force_tangential_distribution):
    #force distrubution    
        plt.figure()
        plt.subplot(2,1,1)
        plt.title('ap = ' + self.ap + ', cof = ' + self.cof + ', ' + self.parameter +  '=' + self.value)
        plt.plot(force_normal_distribution[1][:-1], force_normal_distribution[0])
        plt.xlabel('F_n')
        plt.ylabel('P(F_n)')
        plt.subplot(2,1,2)
        plt.plot(force_tangential_distribution[1][:-1], force_tangential_distribution[0])
        plt.xlabel('F_t')
        plt.ylabel('P(F_t)')
        plt.savefig('output_plots/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + 'force_distribution.png')
        plt.show()

    def plot_ellipsoids(self, step, data_vtk):
        #plot the ellipsoids in 3d
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(data_vtk['trackedGrainsPosition'][step, :, 0], data_vtk['trackedGrainsPosition'][step, :, 2], data_vtk['trackedGrainsPosition'][step, :, 1], s=1)
        ax.set_xlabel('X ')
        ax.set_ylabel('Z ')
        ax.set_zlabel('Y ')
        ax.set_title('ap = ' + self.ap + ', cof = ' + self.cof + ', ' + self.parameter +  '=' + self.value + ', step='+str(step))
        plt.savefig('output_plots/ellipsoids'+str(step)+'.png')
        plt.show()
        plt.close()