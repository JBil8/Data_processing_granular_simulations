import numpy as np
import matplotlib.pyplot as plt
import os

class DataPlotter:
    def __init__(self, ap , cof, parameter= None, value= None,  muw=None, vwall=None, fraction=None, phi = None):   
        self.ap = str(ap) 
        self.cof = str(cof)
        self.parameter = parameter
        self.value = str(value)
        self.muw = muw
        self.vwall = vwall
        self.fraction = fraction
        self.phi = phi
    
    def plot_data(self, data_vtk, data_dump, volume, surface):
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
        if self.parameter == "I":
            plt.plot(strain, data_vtk['phi']) 
            plt.ylabel('phi')
        else:
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
            #plot the eulerian velocity at the given time steps skipping the first 2 strains
            if i != 0:
                plt.plot(data_vtk['eulerian_vx'][i*time_interval,:], y, label = f'strain  = {i*time_interval*20/n_time_steps+2:.2f}', color=color)
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

    def plot_ellipsoids(self, step, data_vtk, data_dump):
        #plot the ellipsoids in 3d
        
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection='3d')
        for i in range(len(data_vtk['trackedGrainsShapeX'][0])):
            center = data_vtk['trackedGrainsPosition'][step, i, :]
            radii = [data_vtk['trackedGrainsShapeX'][step, i],
                     data_vtk['trackedGrainsShapeX'][step, i],
                     data_vtk['trackedGrainsShapeZ'][step, i]]
            rotation_matrix = data_vtk['trackedGrainsOrientation'][step, i, :, :]
            self.plot_single_ellipsoid(ax, center, radii, rotation_matrix)
        
        contact_points = data_dump['trackedGrainsContactData'][step][:, 1:4]
        #self.scatter_contact_point_3D(ax, contact_points)

        ax.set_xlabel('X ')
        ax.set_ylabel('Y ')
        ax.set_zlabel('Z ')
        ax.set_title('ap = ' + self.ap + ', cof = ' + self.cof + ', ' + self.parameter +  '=' + self.value + ', step='+str(step))
        plt.savefig('output_plots/ellipsoids'+str(step)+'.png')
        plt.show()
        plt.close()

    def plot_single_ellipsoid(self, ax, center, radii, rotation_matrix, color='b', alpha=0.1):
        """
        Plot an ellipsoid in 3D.

        Parameters:
        - ax: Axes3D object (matplotlib)
        - center: Center of the ellipsoid (tuple or array)
        - radii: Semi-axes lengths (tuple or array)
        - rotation_matrix: Rotation matrix for orientation (2D array)
        - color: Color of the ellipsoid (default: 'b')
        - alpha: Transparency of the ellipsoid (default: 0.1)
        """
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        # Parametric equations of the ellipsoid
        x = center[0] + radii[0] * np.outer(np.cos(u), np.sin(v))
        y = center[1] + radii[1] * np.outer(np.sin(u), np.sin(v))
        z = center[2] + radii[2] * np.outer(np.ones(np.size(u)), np.cos(v))

        # Rotate the ellipsoid based on the rotation matrix
        for i in range(len(x)):
            for j in range(len(x[i])):
                [x[i, j], y[i, j], z[i, j]] = rotation_matrix@[x[i, j], y[i, j], z[i, j]]

        ax.plot_surface(x, y, z, color=color, alpha=alpha, linewidth=0)
        # Set axis equal
        ax.set_box_aspect([np.ptp(coord) for coord in [x, y, z]])

    def plot_space_averages_all_cells(self, value, quantity= 'Velocity', component = None, nx_divisions = 40, ny_divisions=6):
            """
            Function to plot the space averages of the velocities in eulerian cells
            """
    
            if component == 0:
                axis_name = 'x'
            elif component == 1:
                axis_name = 'y'
            elif component == 2:
                axis_name = 'z'
            elif component == 3:
                axis_name = 'xy'
            elif component == 4:
                axis_name = 'yz'    
            elif component == 5:
                axis_name = 'xz'


            from matplotlib import pyplot as plt
            os.makedirs("eulerian_ap_" + str(self.ap) + "_mup_" + str(self.cof) + "_muw_" + str(self.muw) + 
                        "_vw_" + str(self.vwall), "_phi" + str(self.phi) + "_frac_" + str(self.fraction), exist_ok=True)
            os.chdir("eulerian_ap_" + str(self.ap) + "_mup_" + str(self.cof) + "_muw_" + str(self.muw) + 
                        "_vw_" + str(self.vwall), "_phi" + str(self.phi) + "_frac_" + str(self.fraction))

            # Plot x component of the velocity
            plt.figure(figsize=(20, 12))      
            plt.xlabel('Grid cell x')
            plt.ylabel('Grid cell y')
            #plt.show()
            if component is None:
                plt.imshow(value.T, cmap='plasma', origin='lower', extent=[0, nx_divisions, 0, ny_divisions])
                plt.colorbar(label='$'+ quantity[0] + '$')  # Set ticks for colorbar
                plt.title('Eulerian ' + quantity)
                plt.savefig("eulerian_" + "frac_" + str(self.fraction) + "_" + quantity + ".png")
            
            else:
                plt.imshow(value[:,:,component].T, cmap='plasma', origin='lower', extent=[0, nx_divisions, 0, ny_divisions])
                plt.colorbar(label='$'+ quantity[0]+ '_' + axis_name + '$')  # Set ticks for colorbar
                plt.title('Eulerian ' + axis_name+ ' ' + quantity)
                plt.savefig("eulerian_" + "frac_" + str(self.fraction) + "_" + quantity + "_" + axis_name + ".png")
            
                        # save figure to create gif later include the name of the timestep
            plt.close()
            os.chdir("..")

    def plot_time_series_forces(self, time, fx, fy, fz):
        """
        Function to plot the time series of the given quantity
        """
        plt.figure(figsize=(20, 12))
        plt.plot(time, fx, label='Fx')
        plt.plot(time, fy, label='Fy')
        plt.plot(time, fz, label='Fz')
        plt.xlabel('Time [s]')
        plt.ylabel('Force [N]')
        plt.legend()
        plt.title('fraction = ' + str(self.fraction) + ', ap = ' + str(self.ap) + ', cof = ' + str(self.cof) + ', muw = ' + str(self.muw) + ', vwall = ' + str(self.vwall))
        plt.savefig("time_series_force_" + ', ap = ' + str(self.ap) + ', cof = ' + str(self.cof) + ', muw = ' + str(self.muw) + ', vwall = ' + str(self.vwall) + ".png")
        plt.show()

    def plot_time_series_energy(self, time, tke, rke):
        """
        Function to plot the time series of the given quantity
        """
        plt.figure(figsize=(20, 12))
        plt.plot(time, tke, label='TKE')
        plt.plot(time, rke, label='RKE')
        plt.xlabel('Time [s]')
        plt.ylabel('Energy [J]')
        plt.legend()
        plt.title('fraction = ' + str(self.fraction) + ', ap = ' + str(self.ap) + ', cof = ' + str(self.cof) + ', muw = ' + str(self.muw) + ', vwall = ' + str(self.vwall))
        plt.savefig("time_series_energy_" + ', ap = ' + str(self.ap) + ', cof = ' + str(self.cof) + ', muw = ' + str(self.muw) + ', vwall = ' + str(self.vwall) + ".png")
        plt.show()

    # def scatter_contact_point_3D(self, ax, points):
    #     """
    #     Plot a scatter of points in 3D.

    #     Parameters:
    #     - ax: Axes3D object (matplotlib)
    #     - points: 3D points (tuple or array)
    #     """
    #     ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=4)