import multiprocessing
import numpy as np
import vtk
import functools
from abc import ABC, abstractmethod

class DataProcessor:
    def __init__(self, data):
        self.data_reader = data
        self.num_processes = None
    @abstractmethod
    def process_data(self, num_processes):
        pass
        
    @abstractmethod
    def process_single_step(self, step):
        """Processing on the data for one single step
        Calling the other methods for the processing"""
        pass
