import numpy as np
import matplotlib.pyplot as plt

class DataPlotter:
    def __init__(self, data):
          self.data = data

    
    def plot_data(self):
        time = np.arange(self.data.shape[0])
        fig = plt.figure(figsize=(30, 10))
        ax = fig.add_subplot(321)
        if self.data.shape[1] == 12:
             pass
            # plt.plot(time, self.data[:,3])
            # plt.ylabel('Vx')
            # ax = fig.add_subplot(322)
            # plt.plot(time, self.data[:,7])
            # plt.ylabel('Omegaz')
            # ax = fig.add_subplot(323)
            # plt.plot(time, self.data[:,8])
            # plt.ylabel('Shearing force')
            # ax = fig.add_subplot(324)
            # plt.plot(time, self.data[:,9])
            # plt.ylabel('$\theta$')
            # ax = fig.add_subplot(325)
            # plt.plot(time, self.data[:,10])
            # plt.ylabel('Box height')
            # ax = fig.add_subplot(326)
            # plt.plot(time, self.data[:,11])
            # plt.ylabel('Autocorrelation V')
            # plt.xlabel('Time')

        #plot eulerian velocities at n_plots time steps in time
        n_plots = 10
        n_time_steps = self.data.shape[0]
        time_interval = int(n_time_steps/n_plots)
        y = np.linspace(0, 1, 10)
        fig = plt.figure(figsize=(30, 10))
        for i in range(n_plots):
            ax = fig.add_subplot(2, int(n_plots/2), i+1)
            plt.plot(self.data[i*time_interval,12:22], y)

        plt.ylabel('Eulerian velocities')        
        
        #  #plt.show()
        return fig