# parallelized version of parameter
import os
import sys
import time
import copy
import numpy as np
import matplotlib.pyplot as plt

from mpi4py import MPI
from tqdm import tqdm
from colorama import Fore
from tabulate import tabulate
from scipy.optimize import minimize_scalar
from vachoppy.trajectory import *

# for parallelization
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()

class Calculator:
    def __init__(self,
                 data,
                 index,
                 lattice,
                 interval):
        """
        data : vachoppy.inout.DataInfo
        index : index in data.datainfo
        lattice : vachoppy.trajectory.Lattice
        interval : time interval (ps)
        """
        self.data = data
        self.index = index
        self.lattice = lattice
        self.interval = interval
        self.temp, self.label = self.data.datainfo[self.index]
        self.num_path = len(self.lattice.path_names)
        
        # check interval
        potim = self.data.potim[list(self.data.temp).index(self.temp)]
        if (self.interval * 1000) % potim != 0:
            print(f"unvalid interval ({self.temp}K): interval should be a multiple of potim")
            sys.exit(0)
        else:
            self.step_interval = int(self.interval*1000 / potim)
            
        # quantities
        self.counts = None
        self.unknown = None
        self.t_reside = None
        self.encounter_num = None
        self.encounter_msd = None
        self.encounter_path_names = None
        self.encounter_counts = None
        
        # check success
        self.success = True
        
        # get quantities
        self.get_quantities()
        
    def get_quantities(self):
        xdatcar = os.path.join(self.data.prefix1,
                               f"{self.data.prefix2}.{self.temp}K",
                               f"XDATCAR_{self.label}")
        force = os.path.join(self.data.prefix1,
                             f"{self.data.prefix2}.{self.temp}K",
                             f"FORCE_{self.label}")
        
        # instantiate VacancyHopping
        try:
            traj = VacancyHopping(
                xdatcar=xdatcar,
                lattice=self.lattice,
                force=force if self.data.force is not None else None,
                interval=self.step_interval,
                verbose=False
            )
            traj.correct_multivacancy(start=1)
            traj.check_multivacancy()
        except:
            print(f"Error occured during trajectory analysis : {temp}K, {label}\n")
            self.success = False
            return
        
        if traj.multi_vac is True:
            print(f"Multi-vacancy issue occured : {temp}K, {label}\n")
            self.success = False
            return
        
        if self.data.force is not None:
            traj.correct_transition_state()
        
        # instantiate Analyzer
        anal = Analyzer(
            traj=traj,
            lattice=self.lattice,
            verbose=False
        )
        
        self.counts = anal.counts[:self.num_path]
        self.unknown = anal.path[self.num_path:]
        self.t_reside = anal.total_reside_steps * self.interval
        
        # instantiate Encounter
        enc = Encounter(
            analyzer=anal,
            verbose=False
        )
        self.encounter_num = enc.num_enc
        self.encounter_msd = enc.msd
        self.encounter_path_names = enc.path_names
        self.encounter_counts = enc.path_counts
        

