import os
import sys
import copy
import numpy as np
from tqdm import tqdm
from colorama import Fore
import matplotlib.pyplot as plt
from itertools import combinations_with_replacement

from vachoppy.inout import *
from vachoppy.trajectory import *
from vachoppy.parameter import *
from vachoppy.fingerprint import *
from vachoppy.utils import *

try:
    from mpi4py import MPI
    PARALELL = True
except:
    PARALELL = False

BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m' # Green color
RED = '\033[91m'   # Red color
RESET = '\033[0m'  # Reset to default color


# -m t option
class MakeAnimation:
    def __init__(self,
                 data,
                 temp,
                 label,
                 interval,
                 poscar_lattice,
                 update_alpha=0.75,
                 show_index=False,
                 dpi=300,
                 rmax=3.25,
                 tol=1e-3,
                 tolerance=1e-3,
                 verbose=True):
        
        print("VacHopPy is running...\n")

        if int(temp) in data.temp:
            self.temp = temp
             
        else:
            print(f"{temp}K is not valid.")
            sys.exit(0)
            
        if label in data.label[list(data.temp).index(self.temp)]:
            self.label = label
            
        else:
            print(f"{label} is not valid.")
            sys.exit(0)
            
        for i, temp_i in enumerate(data.temp):
            if abs(self.temp - temp_i) < 1e-9:
                index_temp = i
                break
        
        index_label = data.label[index_temp].index(self.label)
        
        self.pos = data.pos[index_temp][index_label]
        self.cond = data.cond[index_temp][index_label]
        self.force = data.force[index_temp][index_label]
        
        # self.potim = data.potim[index_temp]
        self.interval = interval
        self.poscar_lattice = poscar_lattice
        self.update_alpha = update_alpha
        self.show_index = show_index
        self.dpi = dpi
        self.rmax = rmax
        self.tol = tol
        self.tolerance = tolerance
        self.verbose = verbose
        
        vac_info = VacancyInfo(
            data=data,
            poscar_lattice=poscar_lattice
        )
        
        self.symbol = vac_info.symbol_vac
        self.num_vac = vac_info.number_vac
        print(f"Vacancy type : {self.symbol}")
        print(f"Number of vacancies : {self.num_vac}")
        
        if self.verbose:
            file_out = "trajectory.txt"
            with open(file_out, "w", encoding="utf-8") as f:
                original_stdout = sys.stdout
                sys.stdout = f
                try:
                    self.lattice = Lattice(
                        poscar_lattice=self.poscar_lattice,
                        symbol=self.symbol,
                        rmax=self.rmax,
                        tol=self.tol,
                        tolerance=self.tolerance,
                        verbose=self.verbose
                    )
                finally:
                    sys.stdout = original_stdout
        else:
            self.lattice = Lattice(
                poscar_lattice=self.poscar_lattice,
                symbol=self.symbol,
                rmax=self.rmax,
                tol=self.tol,
                tolerance=self.tolerance,
                verbose=self.verbose
            )
        
        self.traj = Trajectory(
            interval=interval,
            num_vac=self.num_vac,
            lattice=self.lattice,
            pos_file=self.pos,
            force_file=self.force,
            cond_file=self.cond,
            verbose=self.verbose
        )
            
        self.save_animation()
        
        if self.verbose:
            with open(file_out, "a", encoding="utf-8") as f:
                original_stdout = sys.stdout
                sys.stdout = f
                
                print('')
                try:
                    anal = TrajectoryAnalyzer(
                        lattice=self.lattice,
                        trajectory=self.traj,
                        tolerance=self.tolerance,
                        verbose=self.verbose
                    )
                finally:
                    sys.stdout = original_stdout
                    
            print(f"{file_out} is created.")
        
        print('VacHopPy is done.')
            
    def save_animation(self):
        simulation_time = self.traj.nsw * self.traj.potim / 1000
        t_interval = self.traj.interval
        
        print(f'{GREEN}{BOLD}\nInformation on animation{RESET}')
        print(f"    Total simulation time  : {simulation_time:.3f} ps")
        print(f"    Time interval per step : {t_interval:.3f} ps")
        print(f"    Total number of steps  : {self.traj.num_step} (={simulation_time:.3f}/{t_interval:.3f})")
        print('')
        
        print(f'Enter the initial and final steps (min: 0, max: {self.traj.num_step})')
        step = input(f'{MAGENTA}{BOLD}Answer{RESET} (example: 0 500) : ')
    
        try:
            step = list(map(int, step.split())) 
        except:
            print('The step number must be integer.')
            sys.exit(0)
            
        if step[-1] > self.traj.num_step:
            print(f'    The final step should be less than {self.traj.num_step}')
            print(f'    The final step is set to {self.traj.num_step}')
            step[-1] = self.traj.num_step
            
        step = 'all' if step[-1] == -1 else np.arange(step[0], step[-1])
        
        print('\nEnter the fps value for animation')
        fps = input(f'{MAGENTA}{BOLD}Answer{RESET} (example: 10) : ')
        
        try:
            fps = int(fps)
        except:
            print('The fps must be integer.')
            sys.exit(0)
        
        print('')
        self.traj.animation(
            step=step,
            foldername='snapshot',
            update_alpha=self.update_alpha,
            fps=fps,
            dpi=self.dpi,
            label=self.show_index
        )
        
        print('snapshot directory was created.')
    

# -m p option
class EffectiveHoppingParameter:
    def __init__(self, 
                 data, 
                 interval,
                 poscar_lattice, 
                 parallel,
                 file_out='parameter.txt', 
                 rmax=3.0,
                 tol=1e-3,
                 tolerance=1e-3,
                 use_incomplete_encounter=True,
                 inset_correlatoin_factor=True,
                 verbose=True):
        
        self.data = data
        self.interval = interval
        self.poscar_lattice = poscar_lattice
        self.parallel = parallel
        self.file_out = file_out
        self.rmax = rmax
        self.tol = tol
        self.tolerance = tolerance
        self.use_incomplete_encounter = use_incomplete_encounter
        self.inset_correlatoin_factor = inset_correlatoin_factor
        self.verbose = verbose
        
        vac_info = VacancyInfo(
            data=data,
            poscar_lattice=poscar_lattice
        )
        self.symbol = vac_info.symbol_vac
        self.num_vac = vac_info.number_vac
        
        self.lattice = Lattice(
            poscar_lattice=self.poscar_lattice,
            symbol=self.symbol,
            rmax=self.rmax,
            tol=self.tol,
            tolerance=self.tolerance,
            verbose=False
        )
        
        if self.parallel:
            # get effective hopping parameters
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            
            if rank == 0:
                print('VacHopPy is running...')
                print('')
                print(f"Vacancy type : {self.symbol}")
                print(f"Number of vacancies : {self.num_vac}")
                
                f = open('VACHOPPY_PROGRESS', 'w', buffering=1, encoding='utf-8')
                original_stdout = sys.stdout
                sys.stdout = f
                if size < 2:
                    print("number of cpu node shoud be >= 2.")
                    MPI.COMM_WORLD.Abort(1)

            self.results = Automation_parallel(
                data=self.data, 
                lattice=self.lattice, 
                interval=self.interval,
                num_vac=self.num_vac,
                tolerance=self.tolerance,
                use_incomplete_encounter=self.use_incomplete_encounter
            )
                
            if rank == 0:
                sys.stdout = original_stdout
                f.close()
                # merge unknown paths in parallel calculculations
                self.results[-1].lattice = self.merge_lattice()
                self.get_parameters()
                print('VacHopPy is done.')
                
        else:
            print('VacHopPy is running...')
            print('')
            print(f"Vacancy type : {self.symbol}")
            print(f"Number of vacancies : {self.num_vac}")
            
            self.results = Automation_serial(
                data=self.data, 
                lattice=self.lattice, 
                interval=self.interval,
                num_vac=self.num_vac,
                tolerance=self.tolerance,
                use_incomplete_encounter=self.use_incomplete_encounter
            )
            self.get_parameters()
            print('VacHopPy is done.')
            
    def append_unknown_path(self, path, lattice):
        for path_unknown in lattice.path_unknown:
            check1 = True if path['site_init'] == path_unknown['site_init'] else False
            check2 = True if path['site_final'] == path_unknown['site_final'] else False
            check3 = True if abs(path['distance'] - path_unknown['distance']) < self.tolerance else False
            
            if check1 and check2 and check3:
                return path['name'], path_unknown['name']
        
        name_before = path['name']
        path['name'] = f"unknown{len(lattice.path_unknown) + 1}"
        lattice.path_unknown.append(path)
    
        return name_before, path['name']
            
    def merge_lattice(self):
        lattice_merged = None
        for result in self.results:
    
            if lattice_merged is None:
                lattice_merged = copy.deepcopy(result.lattice)
                continue
            
            if int(np.sum(result.counts_unknown)) == 0:
                # no unknown path detected
                continue
            
            name_unknown_before = []
            name_unknown_after = []
            
            for path in result.path_unknown:
                name_before, name_after = self.append_unknown_path(path, lattice_merged)
                
                if name_before == name_after:
                    continue
                
                name_unknown_after.append(name_before)
                name_unknown_after.append(name_after)
            
            # update hopping history
            for history in result.hopping_history:
                for hop in history:
                    if hop['name'] in name_unknown_before:
                        index = name_unknown_before.index(hop['name'])
                        hop['name'] = name_unknown_after[index]
            
            # update unknown_name
            for i, name in enumerate(result.unknown_name):
                if name in name_unknown_before:
                    index = name_unknown_before.index(name)
                    result.unknown_name[i] = name_unknown_after[index]
        
        return lattice_merged
            
    def get_parameters(self):
        with open(self.file_out, 'w', encoding='utf-8') as f:
            original_stdout = sys.stdout
            sys.stdout = f
            
            try:
                extractor=ParameterExtractor(
                    data=self.data,
                    results=self.results,
                    tolerance=self.tolerance,
                    verbose=self.verbose,
                    figure=self.verbose,
                    file_out=None,
                    inset_correlatoin_factor=self.inset_correlatoin_factor
                )
                
            finally:
                sys.stdout = original_stdout  
                
        print('')
        print(f"{self.file_out} is created.")
        print("D_rand.png is created.")
        print("f_cor.png is created.")
        print("tau.png is created.")
        

# -m pp option
class PostEffectiveHoppingParameter:
    def __init__(self,
                 file_params='parameter.txt',
                 file_neb='neb.csv',
                 file_out='postprocess.txt',
                 verbose=True):
        
        self.file_params = file_params
        self.file_neb = file_neb
        self.file_out = file_out
        self.verbose = verbose
        
        with open(self.file_out, 'w', encoding='utf-8') as f:
            original_stdout = sys.stdout
            sys.stdout = f
            
            try:
                postprocess = PostProcess(
                    file_params=self.file_params,
                    file_neb=self.file_neb,
                    verbose=self.verbose
                )
                
            finally:
                sys.stdout = original_stdout
                
        print(f"{self.file_out} is created.")


# -m f option
class PhaseTransition:
    def __init__(self,
                 vasprun,
                 interval,
                 Rmax,
                 delta,
                 sigma,
                 parallel,
                 poscar_ref='POSCAR_REF',
                 prefix1='snapshots',
                 prefix2='fingerprints'):
        
        self.vasprun = vasprun
        self.interval = interval
        self.Rmax = Rmax
        self.delta = delta
        self.sigma = sigma
        self.parallel = parallel
        self.poscar_ref = poscar_ref
        self.prefix1 = prefix1
        self.prefix2 = prefix2
        
        if self.parallel:
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            
            if rank == 0:
                print('VacHopPy is running...\n')
                f = open('VACHOPPY_PROGRESS', 'w', buffering=1, encoding='utf-8')
                original_stdout = sys.stdout
                sys.stdout = f
                if size < 2:
                    print("number of cpu node shoud be >= 2.")
                    MPI.COMM_WORLD.Abort(1)
                try:
                    snapshots = Snapshots(
                        interval=self.interval,
                        vasprun=self.vasprun,
                        prefix=self.prefix1
                    )
                    
                except SystemExit:
                    print("Error occured duing instantiating fingerprint.Snapshots.")
                    MPI.COMM_WORLD.Abort(1)
                    
                except Exception as e:
                    print(f"Unexpected error: {e}")
                    MPI.COMM_WORLD.Abort(1)
                
            else:
                snapshots = None
            
            snapshots = comm.bcast(snapshots, root=0)
            
            phase_transition_parallel(
                snapshots=snapshots,
                poscar_ref=self.poscar_ref,
                Rmax=self.Rmax,
                delta=self.delta,
                sigma=self.sigma,
                prefix=self.prefix2
            )
            
            if rank==0:
                sys.stdout = original_stdout
                f.close()
                print('VacHopPy is done.')
        
        else:
            print('VacHopPy is running...\n')
            snapshots = Snapshots(
                interval=self.interval,
                vasprun=self.vasprun,
                prefix=self.prefix1
            )
            phase_transition_serial(
                snapshots=snapshots,
                poscar_ref=self.poscar_ref,
                Rmax=self.Rmax,
                delta=self.delta,
                sigma=self.sigma,
                prefix=self.prefix2
            )
            print('VacHopPy is done.')
            
            
# -u fingerprint option
class GetFingerPrint:
    def __init__(self, 
                 poscar,
                 Rmax,
                 delta,
                 sigma,
                 prefix='fingerprint', 
                 disp=False):
        
        self.poscar = poscar
        self.Rmax = Rmax
        self.delta = delta
        self.sigma = sigma
        self.prefix = prefix
        self.disp = disp
        
        if not os.path.isfile(self.poscar):
            print(f'{self.poscar} is not found')
            sys.exit(0)
        
        # read poscar
        self.atom = None
        self.read_poscar()
        
        self.pair = []
        self.pair.extend(combinations_with_replacement(self.atom, 2))
        
        # fingerprint
        self.fingerprint = []
        self.get_fingerprint()
        self.fingerprint = np.array(self.fingerprint)
        
        # concat fingerprints
        self.fingerprint_concat = self.fingerprint.reshape(1,-1).squeeze()
        
        # save fingerprint
        self.save_fingerprint()
        
        
    def read_poscar(self):
        with open(self.poscar, 'r') as f:
            lines = [line for line in f]
        self.atom = lines[5].split()
    
    def get_pair(self):
        pair = input('input A and B (ex. Hf-O / Hf-Hf,Hf-O / all) : ')
        pair = pair.replace(" ", "")
        
        if pair == 'all':
            self.pair.extend(combinations_with_replacement(self.atom, 2))
            
        else:
            pair = pair.split(',')
            
            for p in pair:
                atoms = p.split('-')
                
                if len(atoms) == 2 and all(atom in self.atom for atom in atoms):
                    self.pair.append(tuple(atoms))
                    
                else:
                    print(f'Invalid pair : {p}')
                    sys.exit(0)
                    
    def get_params(self):
        params = input("input Rmax, delta, and sigma (ex. 15, 0.01, 0.3) : ")
        params = list(map(float, params.replace(" ", "").split(',')))
        self.Rmax, self.delta, self.sigma = params[0], params[1], params[2]
        
    def get_fingerprint(self):
        for (A, B) in self.pair:
            finger = FingerPrint(A, B, 
                                 self.poscar,
                                 self.Rmax, self.delta, self.sigma)
            self.fingerprint.append(finger.fingerprint)
    
    def save_fingerprint(self):
        # save figure
        R = np.linspace(0, self.Rmax, len(self.fingerprint[0]))
        for i in range(len(self.pair)):
            x = R + i*self.Rmax*np.ones_like(R)
            plt.plot(x, self.fingerprint[i], label=f'{self.pair[i][0]}-{self.pair[i][1]}')
        
        plt.axhline(0, 0, 1, color='k', linestyle='--', linewidth=1)
        plt.xlabel("r (Å)", fontsize=13)    
        plt.ylabel('Intensity', fontsize=13)
        plt.legend(fontsize=12)
        plt.savefig(f'{self.prefix}.png', dpi=300)
        print(f'{self.prefix}.png is created.')
        if self.disp:
            plt.show()
        plt.close()
        
        R = np.linspace(0, self.Rmax * len(self.pair), len(self.fingerprint_concat))
        with open(f'{self.prefix}.txt', 'w') as f:
            f.write(f'# Rmax, delta, sigma = {self.Rmax}, {self.delta}, {self.sigma}\n')
            f.write('# pair : ')
            for (A, B) in self.pair:
                f.write(f'{A}-{B}, ')
            f.write('\n')
            for x, y in zip(R, self.fingerprint_concat):
                f.write(f'  {x:2.6f}\t{y:2.6f}\n')
        print(f'{self.prefix}.txt is created.')


# -u cosine_distance option
class GetCosineDistance:
    def __init__(self, fp1, fp2):
        if not os.path.isfile(fp1):
            print(f'{fp1} is not found')
            sys.exit(0)
        if not os.path.isfile(fp2):
            print(f'{fp2} is not found')
            sys.exit(0)
        
        self.fp1 = np.loadtxt(fp1, skiprows=2)[:,1]
        self.fp2 = np.loadtxt(fp2, skiprows=2)[:,1]
        
        if self.fp1.shape != self.fp2.shape:
            print('Size of two fingerprints should be the same')
            sys.exit(0)
        
        self.d_cos = CosineDistance(self.fp1, self.fp2)
        print(f'd_cos = {self.d_cos}')
   
