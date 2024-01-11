import multiprocessing
import numpy as np
import vtk
from DataProcessor import DataProcessor

class ProcessorVtk(DataProcessor):
    def __init__(self, data):
        super().__init__(data)
        self.n_sim = self.data_reader.n_sim
        self.n_central_atoms = self.data_reader.n_central_atoms
        self.v0 = self.data_reader.v0
        self.y0 = self.data_reader.y0
        self.directory = self.data_reader.directory
        self.file_list = self.data_reader.file_list
        self.n_wall_atoms = self.data_reader.n_wall_atoms
        self.n_central_atoms = self.data_reader.n_central_atoms
        self.polydata = None
        self.polydatapoints = None
        self.ids = None

    def process_data(self, num_processes):
        """
        Processing the data in parallel, returning the averages        
        """
        self.num_processes = num_processes
        with multiprocessing.Pool(self.num_processes) as pool:
            print("Started multiprocessing")
            results = pool.map(self.process_single_step,
                               [step for step in range(self.n_sim)])
        # Extract the averages from the results
        averages = np.array(results)
        return averages

    def process_single_step(self, step):
        """Processing on the data for one single step
        Calling the other methods for the processing
        Concatenating the results in a single array"""

        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.directory+self.file_list[step])
        reader.Update()
        self.polydata = reader.GetOutput()
        self.polydatapoints = self.polydata.GetPointData()
        self.get_ids()
        self.get_data()
        self.compute_space_averages()
        self.compute_box_height() 
        self.compute_alignment()
        self.compute_mean_square_displacement()
        eulerian_velocities = self.eulerian_velocity(10)
        return np.concatenate((self.velocities_space_average,
                               self.vx_space_average.reshape(1,),
                               self.omegas_space_average, 
                               self.omegaz_space_average.reshape(1,),
                               self.shearing_force.reshape(1,), 
                               self.alignment_space_average.reshape(1,),
                               self.alignment_out_of_flow.reshape(1,),
                               self.S2_space_average.reshape(1,), 
                               self.box_height.reshape(1,),
                               self.compute_autocorrelation_vel().reshape(1,), 
                               self.mean_square_displacement.reshape(1,),
                               self.effective_friction.reshape(1,),
                               eulerian_velocities),)
    
    def get_ids(self):
        """
        read and sort the identifiers. Particles making the walls have ids with lower values
        """
        ids = np.array(self.polydatapoints.GetArray("id"))
        self.ids = np.argsort(ids) #ids change at every time step so we need to rearrange them

    def get_data(self):
        """
        Store data from vtk into numpy arrays
        """
        self.coor = np.array(self.polydata.GetPoints().GetData())[self.ids][self.n_wall_atoms:,:]
        self.velocities = np.array(self.polydatapoints.GetArray("v"))[self.ids, :][self.n_wall_atoms:, :]
        self.forces_particles = np.array(self.polydatapoints.GetArray("f"))[self.ids, :][self.n_wall_atoms:, :]
        coordinates_bottom_wall = np.array(self.polydata.GetPoints().GetData())[self.ids][int(self.n_wall_atoms/2):self.n_wall_atoms, :]
        if not np.allclose(coordinates_bottom_wall[:,1], np.min(coordinates_bottom_wall[:,1]), atol=1e-5):
            self.forces_walls = np.array(self.polydatapoints.GetArray("f"))[self.ids, :][int(self.n_wall_atoms/2):self.n_wall_atoms, :]
        else:
            self.forces_walls = np.array(self.polydatapoints.GetArray("f"))[self.ids, :][int(self.n_wall_atoms/2+1):self.n_wall_atoms-1, :]
        self.omegas = np.array(self.polydatapoints.GetArray("omega"))[self.ids, :][self.n_wall_atoms:, :]
        self.torques = np.array(self.polydatapoints.GetArray("tq"))[self.ids, :][self.n_wall_atoms:, :]
        self.orientations = np.array(self.polydatapoints.GetArray("TENSOR"))[self.ids, :][self.n_wall_atoms:, :].reshape(self.n_central_atoms,3,3)
        self.stress = np.concatenate((np.array(self.polydatapoints.GetArray("c_stressAtom[1-3]")), np.array(self.polydatapoints.GetArray("c_stressAtom[4-6]"))), axis=1)

    def compute_autocorrelation_vel(self):
        """
        Compute the autocorrelation of the velocities with respect to the initial velocities -> to change to start at steady state
        """
        autocorrelation_vel = np.mean(np.sum((self.v0-np.mean(self.v0))*(self.velocities-np.mean(self.velocities)), axis=1))
        return autocorrelation_vel

    def compute_space_averages(self):
        """
        Compute the space averages of the velocities, angular velocities and forces
        """
        self.velocities_space_average = np.mean(self.velocities, axis=0)
        self.vx_space_average = np.mean(self.velocities[:,0])
        self.omegas_space_average = np.mean(self.omegas, axis=0)
        self.omegaz_space_average = np.mean(self.omegas[:, 2])
        self.shearing_force = np.mean(self.forces_walls[:,0]) #only x component of the force 
        self.vertical_force = np.mean(self.forces_walls[:,1]) #only y component of the force
        self.effective_friction = self.shearing_force/self.vertical_force
    
    def compute_box_height(self):
        """
        Compute the height of the box from the maximum y coordinate of the particles
        """
        self.box_height = np.max(self.coor[:,1])-np.min(self.coor[:,1])
        

    def compute_mean_square_displacement(self):
        """
        Compute the mean square displacement of the particles in the direction perpendicular to the flow to the flow
        """
        #compute the mean square displacement of the particles
        self.mean_square_displacement = np.mean((self.coor[:,1]-self.y0)**2)

    def compute_alignment(self, store_heads = None):
        """
        Compute the alignment of the particles with respect to the flow direction (angle theta in previous papers)
        I am assuming there is no alignment in the z direction (out of plane)
        Then I compute the nematic order parameter S2 along that direction
        """
        starting_vector = np.array([0,0,1])
        flow_angle = np.zeros((self.n_central_atoms, 1))
        out_flow_angle = np.zeros((self.n_central_atoms, 1))
        S2 = np.zeros((self.n_central_atoms, 1))
        if store_heads is not None:
            arrow_heads = np.zeros((self.n_central_atoms, 3)) #for plotting purposes

        for j in range(self.n_central_atoms):
            rot = self.orientations[j]
            if not np.isclose(rot@rot.T, np.diag(np.ones(3))).all():
                raise SystemExit('Error: One of the points did not store a rotation matrix')
            #project the vector on the plane x and y -> just take the first two components
            main_axis_ellipsoid = rot@starting_vector # longest axis is always stored as the third column
            flow_angle[j] = np.arctan2(main_axis_ellipsoid[1], main_axis_ellipsoid[0])
            out_flow_angle[j] = np.arctan2(main_axis_ellipsoid[2], main_axis_ellipsoid[0])
            
        self.alignment_out_of_flow = np.mean(out_flow_angle) #maybe also compute time average
        self.alignment_space_average = np.mean(flow_angle)
        # nematic_vector = np.array([np.cos(self.alignment_space_average), np.sin(self.alignment_space_average), 0]) #direction of average alignment
        # for j in range(self.n_central_atoms):
        #     rot = self.orientations[j]
        #     if store_heads is not None:
        #         arrow_heads[j] = rot@starting_vector 
        #     grain_vector = rot@starting_vector
        #     grain_vectorxy_plane = np.array([grain_vector[0], grain_vector[1], 0])/np.linalg.norm(grain_vector[:2])
        #     cos_theta = np.dot(grain_vector, nematic_vector)
        #     S2[j] = (3*cos_theta**2-1)/2

        #2d alternative as in the paper
        nematic_vector = np.array([np.cos(self.alignment_space_average), np.sin(self.alignment_space_average)])
        Qx = np.zeros(self.n_central_atoms)
        Qy = np.zeros(self.n_central_atoms)
        for j in range(self.n_central_atoms):
            rot = self.orientations[j]
            grain_vector = rot@starting_vector
            grain_vectorxy_plane = np.array([grain_vector[0], grain_vector[1]])/np.linalg.norm(grain_vector[:2])
            cos_theta = np.dot(grain_vectorxy_plane, nematic_vector)
            theta = np.arccos(cos_theta)
            Qx[j] = np.cos(2*theta)
            Qy[j] = np.sin(2*theta)
            
        self.S2_space_average = np.sqrt(np.mean(Qx)**2+np.mean(Qy)**2) #nematic order parameter in 2 D
        

    def eulerian_velocity(self, n_intervals):
        """
        Compute velocities in eulerian coordinates
        """
        # Layers for velocities eulerian
        velocities_per_layer = np.zeros((n_intervals,))
        #finding domain dimensions
        h_max = max(self.coor[:,1]) #maximum height
        h_min = min(self.coor[:,1]) #minimum height
        scaled_factor = (h_max-h_min)/(n_intervals-1)
        scaled_coor = np.floor((self.coor[:,1]-h_min)/scaled_factor)
        for j in range(n_intervals):
            bin_idxs = np.where(scaled_coor==j)[0] #np.where returns a tuple with indexes in the first component
            velocities_per_layer[j] = np.mean(self.velocities[bin_idxs,0])
        return velocities_per_layer

    #compute single particles space trajectory
    def compute_single_particle_trajectory(self, step, n_sampled_particles):
        #find the particle index
        sampledIDxs = np.linspace(0, self.n_central_atoms-1, n_sampled_particles).astype('int')
        trackedGrains = np.zeros((n_sampled_particles, 3))
    
        #find the coordinates of the particles
        trackedGrains = self.coor[sampledIDxs, :]
        
        return trackedGrains
