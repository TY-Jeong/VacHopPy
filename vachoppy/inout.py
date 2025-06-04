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
                 prefix_xdatcar: str = 'XDATCAR',
                 prefix_force: str = 'FORCE',
                 prefix_outcar: str = 'OUTCAR',
                 verbose: bool =False):
        """
        Arguments
        ---------
        prefix1 : str, optional
            Name of outer directory
        prefix2 : str, optional
            Prefix for subdirectories
        prefix_xdatcar : str, optional
            Prefix for XDATCAR files
        prefix_force : str, optional
            Prefix for XDATCAR files
        prefix_outcar : str, optional
            Name of OUTCAR files
        verbose : bool, optional
            Verbosity flag
        """
        
        # directory information
        if os.path.isdir(prefix1):
            self.prefix1 = prefix1
        else:
            print(f'{prefix1} is not found.')
            sys.exit(0)

        self.prefix2 = prefix2
        self.prefix_xdatcar = prefix_xdatcar
        self.prefix_force = prefix_force
        self.prefix_outcar = prefix_outcar
        self.verbose = verbose
            
        # subdirectories
        self.sub_directory = []
        self.get_sub_directory()
        
        # read outcar
        self.temp = []
        self.potim = []
        self.outcar = []
        self.read_outcar()
        
        # labels
        self.label = []
        self.xdatcar = []
        self.read_label()
        
        # force
        self.force = None
        self.read_force()
        
        self.datainfo = [
            [temp, label] for i, temp in enumerate(self.temp)
            for label in self.label[i]
        ]
        
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
        
    def read_outcar(self):
        for dir in self.sub_directory:
            outcar = os.path.join(self.prefix1, dir, self.prefix_outcar)
            
            if not os.path.isfile(outcar):
                print(f"{outcar} is not found.")
                sys.exit(0)
            else:
                self.outcar.append(outcar)
            
            check_temp, check_potim = False, False
            with open(outcar, 'r') as f:
                for line in f:
                    if "TEBEG " in line:
                        token = line.split()[2]
                        temp = float(''.join(c for c in token if c.isdigit() or c == '.'))
                        self.temp.append(temp)
                        check_temp = True
                        
                    if "POTIM " in line:
                        potim = float(line.split()[2])
                        self.potim.append(potim)
                        check_potim = True
                    
                    if check_temp and check_potim:
                        break
                    
    def read_label(self):
        for dir in self.sub_directory:
            list_xdatcar = []
            list_file = os.listdir(
                os.path.join(self.prefix1, dir)
            )
            list_file.sort()
            
            label = []
            for file in list_file:
                _file = file.split('_')
                if _file[0] == self.prefix_xdatcar and len(_file) > 1:
                    label.append(file.split('_')[1])
                    list_xdatcar.append(os.path.join(
                        os.path.join(self.prefix1, dir, file)
                    ))
                    
            self.xdatcar.append(list_xdatcar)
            self.label.append(label)
            
    def read_force(self):
        self.force = []
        for i, dir in enumerate(self.sub_directory):
            list_force = []
            path_dir = os.path.join(self.prefix1, dir)
            for label in self.label[i]:
                if f"FORCE_{label}" not in os.listdir(path_dir):
                    list_force.append(None)
                else:
                    list_force.append(
                        os.path.join(path_dir, f"FORCE_{label}")
                    )
            self.force.append(list_force)
            
    def summary(self):
        label_all = list(
            {element for sublabel in self.label for element in sublabel}
        )
        label_all.sort()
        
        print("# List of temperature (K) :")
        for temp in self.temp:
            print(str(temp), end=' ')
        print('\n')
        
        print("# List of POTIM (fs) : ")
        header = [f"{temp}K" for temp in self.temp]
        data = [self.potim]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        print("# List of XDATCAR :")
        header = ['label'] + [f"{temp}K" for temp in self.temp]
        data = [
            [label] + ['O' if label in self.label[j] else 'X' for j in range(len(self.temp))]
            for label in label_all
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        print("# List of FORCE :")
        header = ['label'] + [f"{temp}K" for temp in self.temp]
        data = [
            [label] +
            [
                'O' if os.path.join(
                    self.prefix1, self.sub_directory[j], f"{self.prefix_force}_{label}"
                ) in self.force[j] else 'X'
                for j in range(len(self.temp))
            ]
            for label in label_all
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
             

class VacancyInfo:
    def __init__(self,
                 data,
                 poscar_lattice: str):
        """
        Arguments
        ---------
        data : Data
            Data object
        poscar_lattice : str
            Path to POSCAR_LATTICE (POSCAR of pefect crystalline)
        """
        
        self.data = data
        self.poscar_lattice = poscar_lattice
        
        self.atom_symbol_poscar = None
        self.atom_number_poscar = None
        self.read_poscar()
        
        self.symbol_vac = None
        self.number_vac = None
        self.read_xdatcar()
        
    def read_poscar(self):
        # perfect crystal structure
        check_symbol, check_number = False, False
        with open(self.poscar_lattice, 'r') as f:
            for i, line in enumerate(f):
                if i == 5:
                    self.atom_symbol_poscar = line.split()
                    check_symbol = True
                    
                if i == 6:
                    self.atom_number_poscar = list(map(int, line.split()))
                    check_number = True
                    
                if check_symbol and check_number:
                    break
                
    def read_xdatcar(self):
        xdatcar = self.data.xdatcar
        
        atom_symbol_xdatcar = []
        atom_number_xdatcar = []
        _atom_number_xdatcar = None
        
        for xdatcar_temp in xdatcar:
            for xdatcar_i in xdatcar_temp:
                check_symbol, check_number = False, False
                with open(xdatcar_i, 'r') as f:
                    for i, line in enumerate(f):
                        if i == 5:
                            symbol = line.split()
                            if symbol != self.atom_symbol_poscar:
                                print(f"Error: unmatched atom species ({xdatcar_i}).")
                                sys.exit(0)
                            atom_symbol_xdatcar.append(line.split())
                            check_symbol = True
                            
                        if i == 6:
                            number = list(map(int, line.split()))
                            if _atom_number_xdatcar is None:
                                _atom_number_xdatcar = number
                                
                                for j, (n1, n2) in enumerate(
                                    zip(self.atom_number_poscar, number)
                                    ):
                                    if n1 != n2:
                                        self.number_vac = n1 - n2
                                        self.symbol_vac = self.atom_symbol_poscar[j]
                                        break
                            else:
                                if number != _atom_number_xdatcar:
                                    print(f"Error: unmatched atom numbers ({xdatcar_i}).")
                                    sys.exit(0)
                                atom_number_xdatcar.append(number)
                
                            check_number = True
                            
                        if check_symbol and check_number:
                            break   
                        
                        
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