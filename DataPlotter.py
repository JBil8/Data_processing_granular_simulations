import numpy as np
import matplotlib.pyplot as plt

class DataPlotter:
    def __init__(self, data_vtk, data_dump, ap , cof, parameter, value):
          self.data_vtk = data_vtk
          self.data_dump = data_dump
          self.ap = str(ap) 
          self.cof = str(cof)
          self.parameter = parameter
          self.value = str(value)

    
    def plot_data(self):
        strain = np.arange(self.data_vtk.shape[0])*20/self.data_vtk.shape[0]
        
        fig = plt.figure(figsize=(30, 10))
        ax = fig.add_subplot(331)
        plt.plot(strain, self.data_dump[:,11])
        plt.ylabel('Z')
        ax = fig.add_subplot(332)
        plt.plot(strain, self.data_vtk[:,7])
        plt.ylabel('Omegaz [rad/s]')
        ax = fig.add_subplot(333)
        plt.plot(strain, -self.data_vtk[:,8])
        plt.ylabel('Shearing force [N]')
        ax = fig.add_subplot(334)
        plt.plot(strain, self.data_vtk[:,9]*180/np.pi) #convert to degrees
        plt.ylabel('theta [deg]')
        ax = fig.add_subplot(335)
        plt.plot(strain, self.data_vtk[:,10]*180/np.pi) #convert to degrees
        plt.ylabel('theta_z [deg]')
        ax = fig.add_subplot(336)
        plt.plot(strain, self.data_vtk[:,11])
        plt.ylabel('S2')
        ax = fig.add_subplot(337)
        plt.plot(strain, self.data_vtk[:,13])
        plt.ylabel('Autocorrelation V')
        plt.xlabel('Strain')

        # ax = fig.add_subplot(337)
        # plt.plot(strain, np.mean(self.data_dump[:,6:8], axis=1))
        # plt.ylabel('Average shear in contacts')

        # ax = fig.add_subplot(338)
        # plt.plot(strain, self.data_dump[:,8])
        # plt.ylabel('Average area in contacts [m^2]')

        ax = fig.add_subplot(338)
        plt.plot(strain, self.data_dump[:,10])
        plt.ylabel('Average overlap in contacts [m]')

        ax = fig.add_subplot(339)
        plt.plot(strain, self.data_vtk[:,14])
        plt.ylabel('Mean square displacement [m^2]')

        fig.suptitle('ap = ' + self.ap + ', cof = ' + self.cof + ', ' + self.parameter +  '=' + self.value)
        fig.savefig('output_plots/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + '.png')
        
        plt.clf()

        #plot eulerian velocities at n_plots time steps in time
        n_plots = 10
        n_time_steps = self.data_vtk.shape[0]
        time_interval = int(n_time_steps/n_plots)
        y = np.linspace(0, 1, 10)
        fig2 = plt.figure(figsize=(10, 20))
        for i in range(n_plots):
            #ax = fig2.add_subplot(2, int(n_plots/2), i+1)
            color = plt.cm.viridis(i / n_plots)  # color will now be an RGBA tuple
            plt.plot(self.data_vtk[i*time_interval,15:], y, label = f'strain  = {i*time_interval*20/n_time_steps:.2f}', color=color)
        plt.xlabel('Vx/V')
        plt.ylabel('y/H')
        plt.legend()
        fig2.suptitle('ap = ' + self.ap + ', cof = ' + self.cof + ', ' + self.parameter +  '=' + self.value)
        fig2.savefig('output_plots/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + 'eulerian.png')
        plt.clf()
         #plt.show()