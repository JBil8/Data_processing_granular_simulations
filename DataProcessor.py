class DataProcessor:
    def __init__(self, data, num_processes):
        self.data = data
        self.num_processes = num_processes
        

    def process_data(self):
        file_list, n_central_atoms, v0, directory = self.data_reader.read_data()
        with multiprocessing.Pool(self.num_processes) as pool:
            averages = pool.starmap(
                self._process_single_data,
                [(i, n_central_atoms, v0, directory) for i in range(self.n_sim)]
            )
        return averages
    
    def _process_single_data(self, i, n_central_atoms, v0, directory):
        # Implement data processing logic here
        pass