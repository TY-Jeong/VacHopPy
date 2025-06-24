import os
import sys
import time
from tqdm import tqdm
from colorama import Fore
from tabulate import tabulate
from vachoppy.parameter import *

try:
    from mpi4py import MPI
    PARALELL = True
except:
    PARALELL = False

# color map for tqdm
BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m' # Green color
RED = '\033[91m'   # Red color
RESET = '\033[0m'  # Reset to default color


class DataInfo:
    def __init__(self, 
                 prefix1: str ='traj',
                 prefix2: str ='traj',
                 prefix_pos: str = 'pos',
                 prefix_force: str = 'force',
                 prefix_cond: str = 'cond',
                 verbose: bool =True):
        """
        Args:
            prefix1 (str, optional): 
                name of outer directory. Defaults to 'traj'.
            prefix2 (str, optional): 
                prefix for subdirectories. Defaults to 'traj'.
            prefix_pos (str, optional): 
                prefix for pos.npy. Defaults to 'pos'.
            prefix_force (str, optional): 
                prefix for force.npy. Defaults to 'force'.
            prefix_cond (str, optional): 
                prefix for cond.json. Defaults to 'cond'.
            verbose (bool, optional): 
                verbosity flag. Defaults to True.
        """
        
        self.prefix1 = prefix1
        self.prefix2 = prefix2
        self.prefix_pos = prefix_pos
        self.prefix_cond = prefix_cond
        self.prefix_force = prefix_force
        
        # subdirectories
        self.sub_directory = []
        self.get_sub_directory()
        
        # labels
        self.label = []
        self.get_label()
        
        # integrity of input files
        self.temp = []
        self.pos = []
        self.cond = []
        self.force = []
        self.check_integrity()
        
        
        self.datainfo = [
            [temp, label] for i, temp in enumerate(self.temp)
            for label in self.label[i]
        ]
        
        self.label_all = list(
            {element for sublabel in self.label for element in sublabel}
        )
        self.label_all.sort()
        
        if verbose:
            with open('data.txt', 'w') as f:
                original_stdout = sys.stdout
                sys.stdout = f
                try:
                    self.summary()
                finally:
                    sys.stdout = original_stdout
        
    def get_sub_directory(self):
        list_dir = os.listdir(os.path.join(os.getcwd(), self.prefix1))
        self.sub_directory = [d for d in list_dir if self.prefix2 in d]
        self.sub_directory.sort()
        
    def get_label(self):
        for dir in self.sub_directory:
            # gather cond.json file
            label_i = [
                "".join(f.split(".")[:-1]).split("_")[-1] 
                for f in os.listdir(os.path.join(self.prefix1, dir))
                if f.split("_")[0] == self.prefix_cond
                if f.split(".")[-1] == "json"
            ]
            label_i.sort()
            self.label.append(label_i)
            
    def check_integrity(self):
        for i, dir in enumerate(self.sub_directory):
            temp_i = []
            pos_i = []
            cond_i = []
            force_i = []
            for label in self.label[i]:
                pos_file = os.path.join(
                    self.prefix1, dir, f"{self.prefix_pos}_{label}.npy"
                )
                cond_file = os.path.join(
                    self.prefix1, dir, f"{self.prefix_cond}_{label}.json"
                )
                force_file = os.path.join(
                    self.prefix1, dir, f"{self.prefix_force}_{label}.npy"
                )
                
                # check pos.npy, force.npy
                integrity = True
                if not os.path.isfile(pos_file):
                    print(f"Warning: file absent ({pos_file})")
                    integrity = False
                    
                if not os.path.isfile(force_file):
                    print(f"Warning: file absent ({force_file})")
                    integrity = False
                
                if integrity:
                    with open(cond_file, "r") as f:
                        cond = json.load(f)
                    temp_i.append(cond['temperature'])
                    pos_i.append(pos_file)
                    cond_i.append(cond_file)
                    force_i.append(force_file)
                else:
                    self.label[i].remove(label)
                    
            # check temperature
            temp_first = temp_i[0]
            check_temp = all(abs(x - temp_first) < 1e-6 for x in temp_i[1:])
            if check_temp:
                self.temp.append(temp_first)
                self.pos.append(pos_i)
                self.cond.append(cond_i)
                self.force.append(force_i)
            else:
                print(f"Error: temperatures within the same directory must be the same ({dir})")
            
    def summary(self):
        print("# Sumamry of input data:")
        header = ['label'] + [f"{temp}K" for temp in self.temp]
        data = [
            [label] + ['O' if label in self.label[j] else 'X' for j in range(len(self.temp))]
            for label in self.label_all
        ]
        data.append(
            ['Total'] + [len(label) for label in self.label]
        )
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        print("# List of cond.json")
        for i, temp in enumerate(self.temp):
            print(f"{temp} K :")
            for file in self.cond[i]:
                print(f"  {file}")
        print("")
        
        print("# List of pos.npy")
        for i, temp in enumerate(self.temp):
            print(f"{temp} K :")
            for file in self.pos[i]:
                print(f"  {file}")
        print("")
        
        print("# List of force.npy")
        for i, temp in enumerate(self.temp):
            print(f"{temp} K :")
            for file in self.force[i]:
                print(f"  {file}")
        print("")


class VacancyInfo:
    def __init__(self,
                 data,
                 poscar_lattice: str = "POSCAR_LATTICE"):
        """
        Args:
            data (DataInfo): 
                DataInfo object
            poscar_lattice (str): 
                path to POSCAR_LATTICE. Defaults to 'POSCAR_LATTICE'.
        """
        
        self.data = data
        self.poscar_lattice = poscar_lattice
        
        self.lattice_poscar = []
        self.atom_counts_poscar = {}
        self.read_poscar()
        
        self.symbol_vac = None
        self.number_vac = None
        self.read_cond()
        
    def read_poscar(self):
        check_symbol, check_number = False, False
        with open(self.poscar_lattice, 'r') as f:
            for i, line in enumerate(f):
                if i == 2:
                    self.lattice_poscar.append(line.split()[:3])
                    
                if i == 3:
                    self.lattice_poscar.append(line.split()[:3])
                    
                if i == 4:
                    self.lattice_poscar.append(line.split()[:3])
                
                if i == 5:
                    atom_symbol = line.split()
                    check_symbol = True
                    
                if i == 6:
                    atom_number = list(map(int, line.split()))
                    check_number = True
                    
                if check_symbol and check_number:
                    break
                
        self.lattice_poscar = np.array(self.lattice_poscar, dtype=np.float64)
        
        for symbol, number in zip(atom_symbol, atom_number):
            self.atom_counts_poscar[symbol] = number

    def read_cond(self):
        for cond_i in self.data.cond:
            for cond_file in cond_i:
                with open(cond_file, "r") as f:
                    cond = json.load(f)
                    
                # check lattice parameters
                lattice = np.array(cond["lattice"], dtype=np.float64)
                eps = 1e-3
                var = np.linalg.norm(self.lattice_poscar - lattice)
                if var > eps:
                    print(f"Error: unmatched lattice ({cond_file})")
                    sys.exit(0)
                    
                # check symbol
                if set(cond["atom_counts"].keys()) != set(self.atom_counts_poscar.keys()):
                    print(f"Error: unmatched atom species ({cond_file})")
                
                # get vacancy info
                if self.symbol_vac is None:
                    for symbol in cond["atom_counts"].keys():
                        if cond["atom_counts"][symbol] != self.atom_counts_poscar[symbol]:
                            self.symbol_vac = symbol
                            self.number_vac = self.atom_counts_poscar[symbol] - cond["atom_counts"][symbol]
                else:
                    num_vac = self.atom_counts_poscar[symbol] - cond["atom_counts"][symbol]
                    if num_vac != self.number_vac:
                        print(f"Error: unmatched number of vacancies ({cond_file})")
                        sys.exit(0)
                        
                        
def Automation_serial(data, 
                      lattice, 
                      interval,
                      num_vac,
                      tolerance: float = 1e-3,
                      use_incomplete_encounter: bool = True):
    """
    Arguments
    ---------
    data : DataInfo
        DataInfo object
    lattice : Lattice
        Lattice object
    interval : float
        Time interval (ps)
    num_vac : int
        Number of vacancies
    tolerance : flaot, optional
        Tolerance for numerical accuracy
    use_incomplete_encounter : bool, optional
        If true, incomplete encounters are used together in computations.
    """
        
    results = []
    failure = []
    task_size = len(data.datainfo)
    
    for i in tqdm(range(task_size),
                  bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'%(Fore.GREEN, Fore.RESET),
                  ascii=False,
                  desc=f'{RED}{BOLD}Progress{RESET}'):
        
        temp, label = data.datainfo[i]
        try:
            result = Calculator(
                data=data,
                temp=temp,
                label=label,
                lattice=lattice,
                interval=interval,
                num_vac=num_vac,
                tolerance=tolerance,
                use_incomplete_encounter=use_incomplete_encounter
            )
            
        except SystemExit:
            result = Calculator_fail(
                data=data,
                temp=temp,
                label=label
            )
        
        if result.success:
            results.append(result)
            
        else:
            failure.append(
                f"  T={result.temp}K,  Label={result.label} ({result.fail_reason})"
            )
            
    # sort by (temp, label)   
    index = [data.datainfo.index([result.temp, result.label]) for result in results]
    results = [x for _, x in sorted(zip(index, results))]
    
    # print failed calculations
    if len(failure) > 0:
        print(f"Error occured in :")
        for x in failure:
            print(x)
    print('')
    
    return results


def Automation_parallel(data, 
                        lattice, 
                        interval,
                        num_vac,
                        tolerance: float = 1e-3,
                        use_incomplete_encounter: bool = True):
    """
    Arguments
    ---------
    data : DataInfo
        DataInfo object
    lattice : Lattice
        Lattice object
    interval : float
        Time interval (ps)
    num_vac : int
        Number of vacancies
    tolerance : flaot, optional
        Tolerance for numerical accuracy
    use_incomplete_encounter : bool, optional
        If true, incomplete encounters are used together in computations.
    """
    
    time_i = time.time()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    task_size = len(data.datainfo)
    
    if rank==0:
        task_queue = list(range(task_size))
        print(f"Number of AIMD data : {len(task_queue)}")
        
        results, failure = [], []
        completed_task, terminated_worker, active_workers = 0, 0, size - 1

        while completed_task < task_size or terminated_worker < active_workers:
            status = MPI.Status()
            worker_id, task_result = comm.recv(
                source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status
            )
            
            if status.Get_tag() == 4:
                terminated_worker += 1
                continue
                
            if task_result is not None:
                completed_task += 1
                
                if task_result.success:
                    results.append(task_result)
                    state = 'success'
                    
                else:
                    failure.append(
                        f"  T={task_result.temp}K,  Label={task_result.label} ({task_result.fail_reason})"
                    )
                    state = 'fail'
                    
                print(f"Progress: {completed_task}/{task_size} finished ({state}), " +
                      f"T={task_result.temp}K Label={task_result.label}, " + 
                      f"remaining workers = {active_workers - terminated_worker}/{active_workers}")
                
            if task_queue:
                new_task = task_queue.pop()
                comm.send(new_task, dest=worker_id, tag=1)
                
            else:
                comm.send(None, dest=worker_id, tag=0)
                
        while terminated_worker < active_workers:
            worker_id, _ = comm.recv(source=MPI.ANY_SOURCE, tag=4)
            terminated_worker += 1

    else:
        comm.send((rank, None), dest=0, tag=2)
        while True:
            task = comm.recv(source=0, tag=MPI.ANY_TAG)
            
            if task is None:
                comm.send((rank, None), dest=0, tag=4)
                break
            
            try:
                temp, label = data.datainfo[task]
                result = Calculator(
                    data=data,
                    temp=temp,
                    label=label,
                    lattice=lattice,
                    interval=interval,
                    num_vac=num_vac,
                    tolerance=tolerance,
                    use_incomplete_encounter=use_incomplete_encounter
                )
                
            except SystemExit:
                result = Calculator_fail(
                    data=data,
                    temp=temp,
                    label=label
                )
                
            finally:
                comm.send((rank, result), dest=0, tag=3)
            
    if rank==0:
        index = [data.datainfo.index([result.temp, result.label]) for result in results]
        results = [x for _, x in sorted(zip(index, results))]
        
        time_f = time.time()
        
        if failure:
            print(f"\nError occured in :")
            for x in failure:
                print(x)
        print('')
        print(f"Total time taken: {time_f - time_i} s")
        
        return results
    
    
# class DataInfo:
#     def __init__(self, 
#                  prefix1: str ='traj',
#                  prefix2: str ='traj',
#                  prefix_pos: str = 'pos',
#                  prefix_force: str = 'force',
#                  prefix_cond: str = 'cond',
#                  verbose: bool =False):
        
#         self.prefix1 = prefix1
#         self.prefix2 = prefix2
#         self.prefix_pos = prefix_pos
#         self.prefix_force = prefix_force
#         self.prefix_cond = prefix_cond
#         self.verbose = verbose
        
#         # subdirectories
#         self.sub_directory = []
#         self.get_sub_directory()
        
#         self.label = []
#         self.temp = []
#         self.get_label()
#         self.check_files()
#         self.get_temperautre()
        
#         self.datainfo = [
#             [temp, label] for i, temp in enumerate(self.temp)
#             for label in self.label[i]
#         ]
#         if verbose:
#             with open('data.txt', 'w') as f:
#                 original_stdout = sys.stdout
#                 sys.stdout = f
#                 try:
#                     self.summary()
#                 finally:
#                     sys.stdout = original_stdout
    
#     def get_temperautre(self):
#         self.temp = []
#         for i, dir in enumerate(self.sub_directory):
#             cond_file = os.path.join(
#                 self.prefix1, dir, f"cond_{self.label[i][0]}.json"
#             )
#             with open(cond_file, "r") as f:
#                 condition = json.load(f)
#             self.temp.append(condition['temperature'])
            
#     def check_files(self):
#         check = False
#         for i, dir in enumerate(self.sub_directory):
#             for label in self.label[i]:
#                 # check cond.json
#                 cond_file = os.path.join(
#                     self.prefix1, dir, f"cond_{label}.json"
#                 )
#                 pos_file = os.path.join(
#                     self.prefix1, dir, f"pos_{label}.npy"
#                 )
#                 force_file = os.path.join(
#                     self.prefix1, dir, f"force_{label}.npy"
#                 )
#                 if not os.path.isfile(cond_file):
#                     print(f"Error: {cond_file} is not found.")
#                     check = True
#                 if not os.path.isfile(pos_file):
#                     print(f"Error: {pos_file} is not found.")
#                     check = True
#                 if not os.path.isfile(force_file):
#                     print(f"Error: {force_file} is not found.")
#                     check = True
#         if check:
#             sys.exit(0)
        
#     def get_sub_directory(self):
#         list_dir = os.listdir(os.path.join(os.getcwd(), self.prefix1))
#         self.sub_directory = [d for d in list_dir if self.prefix2 in d]
#         self.sub_directory.sort()   
    
#     def get_label(self):
#         for dir in self.sub_directory:
#             list_file = os.listdir(
#                 os.path.join(self.prefix1, dir)
#             )
#             list_file.sort()
            
#             label = []
#             for file in list_file:
#                 extension = file.split(".")[-1]
#                 if extension == "npy":
#                     _file = "".join(file.split(".")[:-1]).split('_')
#                     if _file[0] == self.prefix_pos and len(_file) > 1:
#                         label.append(_file[1])
#             self.label.append(label)
    
#     def summary(self):
#         label_all = list(
#             {element for sublabel in self.label for element in sublabel}
#         )
#         label_all.sort()
        
#         print("# List of data :")
#         header = ['label'] + [f"{temp}K" for temp in self.temp]
#         data = [
#             [label] + ['O' if label in self.label[j] else 'X' for j in range(len(self.temp))]
#             for label in label_all
#         ]
#         data.append(
#             ['Total'] + [len(label) for label in self.label]
#         )
#         print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
#         print('')
             

# class VacancyInfo:
#     def __init__(self,
#                  data,
#                  poscar_lattice: str):
#         """
#         Arguments
#         ---------
#         data : Data
#             Data object
#         poscar_lattice : str
#             Path to POSCAR_LATTICE (POSCAR of pefect crystalline)
#         """
        
#         self.data = data
#         self.poscar_lattice = poscar_lattice
        
#         self.atom_symbol_poscar = None
#         self.atom_number_poscar = None
#         self.read_poscar()
        
#         self.symbol_vac = None
#         self.number_vac = None
#         self.read_xdatcar()
        
#     def read_poscar(self):
#         # perfect crystal structure
#         check_symbol, check_number = False, False
#         with open(self.poscar_lattice, 'r') as f:
#             for i, line in enumerate(f):
#                 if i == 5:
#                     self.atom_symbol_poscar = line.split()
#                     check_symbol = True
                    
#                 if i == 6:
#                     self.atom_number_poscar = list(map(int, line.split()))
#                     check_number = True
                    
#                 if check_symbol and check_number:
#                     break
                
#     def read_xdatcar(self):
#         xdatcar = self.data.xdatcar
        
#         atom_symbol_xdatcar = []
#         atom_number_xdatcar = []
#         _atom_number_xdatcar = None
        
#         for xdatcar_temp in xdatcar:
#             for xdatcar_i in xdatcar_temp:
#                 check_symbol, check_number = False, False
#                 with open(xdatcar_i, 'r') as f:
#                     for i, line in enumerate(f):
#                         if i == 5:
#                             symbol = line.split()
#                             if symbol != self.atom_symbol_poscar:
#                                 print(f"Error: unmatched atom species ({xdatcar_i}).")
#                                 sys.exit(0)
#                             atom_symbol_xdatcar.append(line.split())
#                             check_symbol = True
                            
#                         if i == 6:
#                             number = list(map(int, line.split()))
#                             if _atom_number_xdatcar is None:
#                                 _atom_number_xdatcar = number
                                
#                                 for j, (n1, n2) in enumerate(
#                                     zip(self.atom_number_poscar, number)
#                                     ):
#                                     if n1 != n2:
#                                         self.number_vac = n1 - n2
#                                         self.symbol_vac = self.atom_symbol_poscar[j]
#                                         break
#                             else:
#                                 if number != _atom_number_xdatcar:
#                                     print(f"Error: unmatched atom numbers ({xdatcar_i}).")
#                                     sys.exit(0)
#                                 atom_number_xdatcar.append(number)
                
#                             check_number = True
                            
#                         if check_symbol and check_number:
#                             break   