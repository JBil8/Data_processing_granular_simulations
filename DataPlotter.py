import numpy as np
import matplotlib.pyplot as plt

class DataPlotter:
    def __init__(self, data, ap , cof, parameter, value):
          self.data = data
          self.ap = ap 
          self.cof = cof
          self.parameter = parameter
          self.value = value

    
    def plot_data(self):
        time = np.arange(self.data.shape[0])
        fig = plt.figure(figsize=(30, 10))
        ax = fig.add_subplot(321)

        plt.plot(time, self.data[:,3])
        plt.ylabel('Vx')
        ax = fig.add_subplot(322)
        plt.plot(time, self.data[:,7])
        plt.ylabel('Omegaz')
        ax = fig.add_subplot(323)
        plt.plot(time, self.data[:,8])
        plt.ylabel('Shearing force')
        ax = fig.add_subplot(324)
        plt.plot(time, self.data[:,9])
        plt.ylabel('$\theta$')
        ax = fig.add_subplot(325)
        plt.plot(time, self.data[:,10])
        plt.ylabel('Box height')
        ax = fig.add_subplot(326)
        plt.plot(time, self.data[:,11])
        plt.ylabel('Autocorrelation V')
        plt.xlabel('Time')

        fig.suptitle('ap = ' + self.ap + ', cof = ' + self.cof + ', ' + self.parameter +  '=' + self.value)
        fig.savefig('output_plots/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + '.png')
        
        plt.clf()

        #plot eulerian velocities at n_plots time steps in time
        n_plots = 10
        n_time_steps = self.data.shape[0]
        time_interval = int(n_time_steps/n_plots)
        y = np.linspace(0, 1, 10)
        fig2 = plt.figure(figsize=(30, 10))
        for i in range(n_plots):
            ax = fig2.add_subplot(2, int(n_plots/2), i+1)
            plt.plot(self.data[i*time_interval,12:22], y)
            plt.xlabel('Vx/V')
            plt.ylabel('y/H')
        fig.suptitle('ap = ' + self.ap + ', cof = ' + self.cof + ', ' + self.parameter +  '=' + self.value)
        fig.savefig('output_plots/simple_shear_ap' + self.ap + '_cof_' + self.cof + '_' + self.parameter + '_' + self.value + 'eulerian.png')
        plt.clf()
        #  #plt.show()
        return fig