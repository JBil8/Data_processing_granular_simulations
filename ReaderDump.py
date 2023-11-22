import os 
import re
import numpy as np
from DataReader import DataReader


class ReaderDump(DataReader):
    def __init__(self, cof, ap, parameter, value):
        super().__init__(cof, ap, parameter, value)

    def read_data(self, global_path, prefix):
        self.prepare_data(global_path)
        self.filter_relevant_files(prefix)
        self.sort_files_by_time()
        self.get_number_of_time_steps()

    def filter_relevant_files(self, prefix='shear_ellipsoids_'):
        self.file_list = [filename for filename in self.file_list if filename[:len(prefix)]==prefix and filename[-1] == 'p']
    