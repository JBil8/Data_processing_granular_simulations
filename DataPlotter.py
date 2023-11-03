import numpy as np
import matplotlib.pyplot as plt

class DataPlotter:
    def __init__(self, data):
          self.data = data
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
        plt.ylabel('S2')
        ax = fig.add_subplot(325)
        plt.plot(time, self.data[:,10])
        plt.ylabel('Box height')
        ax = fig.add_subplot(326)
        plt.plot(time, self.data[:,11])
        plt.ylabel('Autocorrelation V')
        plt.xlabel('Time')

        #plt.show()
        return fig