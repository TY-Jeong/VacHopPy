import os
import sys
import csv
import copy
import numpy as np
from tqdm import tqdm
from colorama import Fore
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from itertools import combinations_with_replacement

from vachoppy.inout import *
from vachoppy.einstein import *
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
                 symbol,
                 correlation=True,
                 update_alpha=0.75,
                 show_index=False,
                 dpi=300,
                 rmax=3.0,
                 tol=1e-3,
                 tolerance=1e-3,
                 verbose=True):
        print("VacHopPy is running...")

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
            
        index_temp = list(data.temp).index(self.temp)
        index_label = data.label[index_temp].index(self.label)
        
        self.xdatcar = data.xdatcar[index_temp][index_label]
        self.outcar = data.outcar[index_temp]
        self.force = None if data.force is None else data.force[index_temp][index_label]
        self.potim = data.potim[index_temp]
        
        self.interval = int(interval * 1000 / self.potim) # ps to iteration
        self.poscar_lattice = poscar_lattice
        self.symbol = symbol
        self.correlation = correlation
        self.update_alpha = update_alpha
        self.show_index = show_index
        self.dpi = dpi
        self.rmax = rmax
        self.tol = tol
        self.tolerance = tolerance
        self.verbose = verbose
        
        self.lattice = Lattice(poscar_lattice=self.poscar_lattice,
                               symbol=self.symbol,
                               rmax=self.rmax,
                               tol=self.tol,
                               tolerance=self.tolerance,
                               verbose=self.verbose)
        
        self.traj = Trajectory(xdatcar=self.xdatcar,
                               lattice=self.lattice,
                               force=self.force,
                               interval=self.interval,
                               verbose=self.verbose)
        
        if self.correlation:
            self.do_correction()
            
        self.save_animation()
        
        if self.verbose:
            anal = TrajectoryAnalyzer(
                traj=self.traj,
                lattice=self.lattice,
                tolerance=self.tolerance,
                verbose=self.verbose
            )
        
        print('VacHopPy is done.')
        
    def do_correction(self):
        check_multivac, check_TS = True, True
        _traj = copy.deepcopy(self.traj)
        # multi-vacancy
        try:
            self.traj.correct_multivacancy(start=1)
            self.traj.check_multivacancy()
            if self.traj.multi_vac == False:
                print('Correction for multi-vacancy : success')
            else:
                print('Correction for multi-vacancy : fail')
                check_multivac = False
        except:
            print('Correction for multi-vacancy : fail')
            check_multivac = False
        
        # TS criteria
        try:
            self.traj.correct_transition_state()
            print('Correction for TS criteria : success')
        except:
            print('Correction for TS criteria : fail')
            check_TS = False
        
        if not(check_multivac and check_TS):
            print('The raw trajectory without corrections will be used.')
            self.traj = _traj
            
    def save_animation(self):
        print('\nInformation on animation')
        print(f'  NSW = {self.traj.nsw} ({self.traj.nsw * self.potim / 1000} ps)')
        print(f'  Interval = {self.interval} step ({self.interval * self.potim / 1000} ps)')
        print(f'  Total step = {self.traj.num_step} (={self.traj.nsw}/{self.traj.interval})')
        print('')
        step = input('Enter init and final steps (int; e.g. 0 100 / 0 -1 for all): ')
        try:
            step = list(map(int, step.split()))
        except:
            print('The step number must be integer.')
            sys.exit(0)
            
        if step[-1] > self.traj.num_step:
            print(f'    The final step should be less than {self.traj.num_step}')
            print(f'    The final step is set to {self.traj.num_step}')
            step[-1] = self.traj.num_step
        step = 'all' if step[-1]==-1 else np.arange(step[0], step[-1])
        
        fps = input('Enter fps (int; e.g. 60): ')
        try:
            fps = int(fps)
        except:
            print('The fps must be integer.')
        
        print('')
        self.traj.animation(step=step,
                            potim=self.potim,
                            foldername='snapshot',
                            update_alpha=self.update_alpha,
                            fps=fps,
                            dpi=self.dpi,
                            label=self.show_index)
        
        print('snapshot directory was created.')
    

# -m p option
class EffectiveHoppingParameter:
    def __init__(self, 
                 data, 
                 interval, 
                 poscar_lattice, 
                 symbol,
                 parallel,
                 file_out='parameter.txt', 
                 rmax=3.0,
                 tol=1e-3,
                 tolerance=1e-3,
                 verbose=True):
        
        self.data = data
        self.interval = interval
        self.poscar_lattice = poscar_lattice
        self.parallel = parallel
        self.symbol = symbol
        self.file_out = file_out
        self.rmax = rmax
        self.tol = tol
        self.tolerance = tolerance
        self.verbose = verbose
        
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
                f = open('VACHOPPY_PROGRESS', 'w', buffering=1, encoding='utf-8')
                original_stdout = sys.stdout
                sys.stdout = f
                if size < 2:
                    print("number of cpu node shoud be >= 2.")
                    MPI.COMM_WORLD.Abort(1)

            self.results = VacancyHopping_parallel(
                self.data, self.lattice,self.interval
            )
                
            if rank == 0:
                sys.stdout = original_stdout
                f.close()
                self.get_parameters()
                print('VacHopPy is done.')
                
        else:
            print('VacHopPy is running...')
            self.results = VacancyHopping_serial(
                self.data, self.lattice, self.interval
            )
            self.get_parameters()
            print('VacHopPy is done.')
            
    def get_parameters(self):
        with open(self.file_out, 'w', encoding='utf-8') as f:
            original_stdout = sys.stdout
            sys.stdout = f
            try:
                extractor=ParameterExtractor(
                    results=self.results,
                    data=self.data,
                    lattice=self.lattice,
                    tolerance=self.tolerance,
                    verbose=self.verbose,
                    figure=self.verbose
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
                postprocess = PostProcess(file_params=self.file_params,
                                          file_neb=self.file_neb,
                                          verbose=self.verbose)
            finally:
                sys.stdout = original_stdout
                
        print(f"{self.file_out} is created.")


# -m f option
class PhaseTransition:
    def __init__(self,
                 xdatcar,
                 outcar,
                 interval,
                 Rmax,
                 delta,
                 sigma,
                 parallel,
                 poscar_ref='POSCAR_REF',
                 prefix1='snapshots',
                 prefix2='fingerprints'):
        
        self.xdatcar = xdatcar
        self.outcar = outcar
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
                print('VacHopPy is running...')
                f = open('VACHOPPY_PROGRESS', 'w', buffering=1, encoding='utf-8')
                original_stdout = sys.stdout
                sys.stdout = f
                if size < 2:
                    print("number of cpu node shoud be >= 2.")
                    MPI.COMM_WORLD.Abort(1)
                try:
                    snapshots = Snapshots(
                        xdatcar=self.xdatcar, 
                        outcar=self.outcar, 
                        interval=self.interval,
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
            print('VacHopPy is running...')
            snapshots = Snapshots(
                xdatcar=self.xdatcar, 
                outcar=self.outcar, 
                interval=self.interval,
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
        
        # input parameters
        # self.pair = []
        # self.get_pair()
        
        # params for fingerprint
        # self.Rmax = None
        # self.delta = None
        # self.sigma = None
        # self.get_params()
        
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
   
            
# -m e option
class MSD:
    def __init__(self,
                 data,
                 tmax,
                 skip=0,
                 start=1,
                 symbol='O',
                 x_vac=1):
        '''
        Calculate MSD at temperatures in data.
        Arguments:
            data : instance of dataInfo
            tmax : x-axis of msd plot (ps)
            skip : steps to be skipped (ps)
            start :  (ps)
        '''
        print(f'{CYAN}{BOLD}Diffusion coefficient from Einstein relation.{RESET}')
        self.data = data
        self.tmax = tmax
        self.skip = skip
        self.start = start
        self.symbol = symbol
        self.x_vac = x_vac
        
        self.msd = []
        self.D = []
        self.runEnsembleEinstein()
        self.saveMSD()
        self.D /= x_vac
        
        self.plot_results()
        
    def runEnsembleEinstein(self):
        desc = 'Einstein'
        for i, T in enumerate(
            tqdm(self.data.temp, 
                 bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'% (Fore.GREEN, Fore.GREEN),
                 ascii=False,
                 desc=f'{GREEN}{BOLD}{desc:>9s}')):
            path_dir = os.path.join(self.data.prefix1, f'{self.data.prefix2}.{T}K')
            step_skip = int(self.skip * 1000 / self.data.potim[i])
            step_tmax = int(self.tmax * 1000 / self.data.potim[i])
            step_start = int(self.start * 1000 / self.data.potim[i])
            step_end = int(self.tmax * 1000 / self.data.potim[i])
            
            if (self.data.nsw[i] - step_skip) % step_tmax != 0:
                print(f'The MD time is not divided by the tmax. (T={T} K)')
                time_tot = self.data.nsw[i] * self.data.potim[i] / 1000
                print(f'  total MD time = {time_tot} ps ({self.data.nsw[i]} step)')
                print(f'  skip = {self.skip} ps ({step_skip} step)')
                print(f'  tmax = {self.tmax} ps ({step_tmax} step)')
                print(f'  ({self.data.nsw[i]} - {step_skip}) % {step_tmax} != 0')
                sys.exit(0)
            else:
                segment = int((self.data.nsw[i] - step_skip) / step_tmax)
                
            msd = EnsembleEinstein(symbol=self.symbol,
                                   prefix=path_dir,
                                   labels=self.data.label[i],
                                   segments=segment,
                                   skip=step_skip,
                                   start=step_start,
                                   end=step_end)
            self.msd.append(msd)
    
    def saveMSD(self):
        for msd in self.msd:
            msd.save_msd()
            self.D.append(msd.diffcoeff)
        self.D = np.array(self.D)
        print('msd directory is created.')
        print('msd raw data is written in msd directory.')

    def plot_results(self):
        kb = 8.61733326e-5
        cmap = plt.get_cmap("Set1")
        
        # msd plot
        plt.rcParams['figure.figsize'] = (3.8, 3.8)
        plt.rcParams['font.size'] = 11
        fig, ax = plt.subplots()
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.2)
            
        for i, msd in enumerate(self.msd):
            x, y = msd.timestep, msd.msd
            x_fit = np.linspace(self.start, self.tmax, 1000)
            slop, intercept = np.polyfit(x, y, deg=1)
            ax.plot(x, y, 
                    c=cmap(i), linestyle='-', linewidth=2.5, alpha=0.35)
            ax.plot(x_fit, slop*x_fit+intercept, 
                    c=cmap(i), linestyle='-', linewidth=2.5, label=str(self.data.temp[i]))
        
        ax.axvline(self.start, 0, 1, c='k', linewidth=1, linestyle=':')
        plt.xlim([-1, self.tmax])
        plt.xlabel('t (ps)', fontsize=14)
        plt.ylabel(r'MSD ($Å^2$)', fontsize=14)
        plt.legend(loc='upper left', fontsize=12,
                   title='T (K)', title_fontsize=13, 
                   fancybox=True, framealpha=1, edgecolor='inherit')
        plt.savefig('msd.svg', transparent=True, dpi=600, bbox_inches="tight")
        # plt.show()
        print('msd.svg is created.')
            
        # Arrhenius plot
        plt.rcParams['figure.figsize'] = (3.8, 3.8)
        plt.rcParams['font.size'] = 11
        fig, ax = plt.subplots()
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.2)
            
        slop, intercept = np.polyfit(1/self.data.temp, np.log(self.D), deg=1)
        x_fit = np.linspace(1/self.data.temp[-1], 1/self.data.temp[0], 100)
        plt.plot(x_fit, slop*x_fit+intercept, 'k:')
        
        for i, D in enumerate(self.D):
            plt.scatter(1/self.data.temp[i], np.log(D), 
                        color=cmap(i), marker='s', s=50, label=str(self.data.temp[i]))
        plt.xlabel('1/T (1/K)', fontsize=14)
        plt.ylabel(r'ln $D$', fontsize=14)
        num_data = len(self.D)
        ncol = int(np.ceil(num_data / 5))
        plt.legend(loc='best', fancybox=True, framealpha=1, edgecolor='inherit',
                   ncol=ncol, labelspacing = 0.3, columnspacing=0.5, borderpad=0.2, handlelength=0.6,
                   fontsize=11, title='T (K)', title_fontsize=11)
        
        if num_data >= 3:
            x = np.array([self.data.temp[0], self.data.temp[int(num_data/2)] ,self.data.temp[-1]])
        else:
            x = self.data.temp
            
        x_str = [f'1/{T}' for T in x]
        x = 1/x
        plt.xticks(x, x_str)
        plt.savefig('Arrhenius.svg', transparent=True, dpi=600, bbox_inches="tight")
        # plt.show()
        print('Arrhenius.svg is created.')
        
        # write results
        with open('Einstein.txt', 'w') as f:
            f.write(f'symbol of moving atom = {self.symbol}\n')
            f.write(f'x_vac = {self.x_vac :.5f}\n')
            if self.x_vac == 1:
                f.write('target : atom\n')
            else:
                f.write('target : vacancy\n')
            f.write('\nParameters for diffusion coefficient : \n')
            f.write(f'  D0   = {np.exp(intercept):.5e} m2/s\n')
            f.write(f'  Ea_D = {-kb * slop:.5f} eV\n')
            f.write('\n')
            f.write('Raw data\n')
            f.write('T(K) \tD(m2/s)\n')
            for temp, D in zip(self.data.temp, self.D):
                f.write(f'{temp} \t{D :.5e}\n')
                   
        print('Einstein.txt is created.')
        print('parameters for diffusion coefficient : ')
        print(f'  D0   = {np.exp(intercept) :.5e} m2/s')
        print(f'  Ea_D = {-kb * slop :.5f} eV')
        print('')