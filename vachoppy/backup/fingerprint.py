import os
import sys
import time
import json
import math
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from lxml import etree
from colorama import Fore
from collections import Counter
from pymatgen.io.vasp.outputs import Vasprun
from itertools import combinations_with_replacement

from vachoppy.utils import CosineDistance

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


class FingerPrint:
    def __init__(self,
                 A,
                 B,
                 poscar,
                 Rmax=10,
                 delta=0.08,
                 sigma=0.03,
                 dirac='g'):
        """
        poscar: path of POSCAR file. (in direct format)
        Rmax: threshold radius
               (valid value: 1-30)
        delta: discretization of fingerprint function 
               (valid value: 0.01-0.2)
        sigma: Gaussian broadening of interatomic distances
               (valid value: 0.01-0.1)
        dirac: 's' for square-form dirac function or 
               'g' for Gaussian broaden dirac function.
        """
        self.A = A
        self.B = B
        self.poscar = poscar
        self.Rmax = Rmax
        self.delta = delta
        self.sigma = sigma
        self.dirac = dirac
        self.R = np.arange(0, self.Rmax, self.delta)
        self.fingerprint = np.zeros_like(self.R)
        self.tolerance = 1e-3
        
        # read poscar
        self.lattice = None
        self.atom_species = None
        self.atom_num = None
        self.position = []
        self.read_poscar()
        
        # Search indices of A and B
        self.idx_A = self.search_atom_index(self.A)
        self.idx_B = self.search_atom_index(self.B)
        
        # volume of unit cell
        self.V_unit = np.abs(
            np.dot(np.cross(self.lattice[0], self.lattice[1]), 
                   self.lattice[2]))
        
        # extended coordinaitons of B
        self.B_ext = None
        self.get_extended_coords_B()

        # calculate fingerprint
        self.get_fingerprint()
    
    
    def read_poscar(self):
        with open(self.poscar, 'r') as f:
            lines = [line.strip() for line in f]
        
        scale = float(lines[1])
        self.lattice = [line.split() for line in lines[2:5]]
        self.lattice = np.array(self.lattice, dtype=float) * scale
        
        self.atom_species = lines[5].split()
        self.atom_num = np.array(lines[6].split(), dtype=int)
        
        # check whether poscar is direct type
        if not lines[7][0].lower() == 'd':
            print('only direct type POSCAR is supported.')
            sys.exit(0)
        
        # parse positions
        start = 8
        for atom, num in zip(self.atom_species, self.atom_num):
            coords = [line.split()[:3] for line in lines[start:start+num]]
            coords = np.array(coords, dtype=float)
            self.position.append({'name': atom, 'num': num, 'coords':coords})
            start += num
           
            
    def gaussian_func(self, x):
        return (1/np.sqrt(2*np.pi*self.sigma**2)) *\
               np.exp(-x**2 / (2*self.sigma**2))
    
    
    def square_func(self, x):
        return np.array(np.abs(x) <= self.sigma)/(2*self.sigma)
    
    
    def dirac_func(self, x):
        if self.dirac[0] == 'g':
            return self.gaussian_func(x)
        elif self.dirac[0] == 's':
            return self.square_func(x)
        else:
            raise ValueError(f"{self.dirac} is not defined.")
            
            
    def search_atom_index(self, atom_name):
        for i, atom in enumerate(self.position):
            if atom['name'] == atom_name:
                return i
        raise ValueError(f"Atom {atom_name} not found in POSCAR file.")
    
    
    def get_extended_coords_B(self):
        # supercells within Rmax
        l_lat = np.linalg.norm(self.lattice, axis=1)
        m = np.floor(self.Rmax / l_lat) + 1
        
        # make 3D grid of supercells
        mx, my, mz = [np.arange(-mi, mi+1) for mi in m]
        shifts = np.array(np.meshgrid(mx, my, mz)).T.reshape(-1, 3)

        # save candidated of B atoms
        coords_B = self.position[self.idx_B]['coords']
        self.B_ext = np.vstack([shifts + coord for coord in coords_B])
             
             
    def get_fingerprint(self):
        for coord_A in self.position[self.idx_A]['coords']:
            self.fingerprint += self.get_fingerprint_i(coord_A)
        
        self.fingerprint *= self.V_unit / self.position[self.idx_A]['num']
        self.fingerprint -= 1
        
        
    def get_fingerprint_i(self, coord_A_i):
        # calculate R_ij
        disp = self.B_ext - coord_A_i
        disp_cart = np.dot(disp[:,:], self.lattice)
        R_ij = np.linalg.norm(disp_cart, axis=1)
        
        # When A=B, i=j should be excluded. (else, diverge)
        if self.idx_A == self.idx_B:
            R_ij[R_ij < self.tolerance] = np.inf
            
        # number of B atoms within Rmax
        N_B = np.sum(R_ij <= self.Rmax)
        fingerprint_i = np.zeros_like(self.fingerprint)

        for idx, r in enumerate(self.R):
            dirac_values = self.dirac_func(r - R_ij)
            valid_indices = R_ij <= self.Rmax
            fingerprint_i[idx] = np.sum(
                dirac_values[valid_indices] / R_ij[valid_indices]**2)\
                      / (4 * np.pi * N_B * self.delta)

        return fingerprint_i
    
    
    def plot_fingerprint(self, 
                         disp=True,
                         save=False,
                         label=None,
                         outdir='./',
                         dpi=300,
                         R=None):
        if R is None:
            R = self.R

        if label is None:
            label = f"{self.A}-{self.B}"

        plt.plot(R, self.fingerprint, label=label)
        plt.axhline(0, 0, 1, color='k', linestyle='--', linewidth=1)
        
        plt.xlabel("Distance (Å)", fontsize=13)
        plt.ylabel('Intensity', fontsize=13)
        
        plt.legend(fontsize=12)

        if save:
            if not os.path.isdir(outdir):
                os.makedirs(outdir, exist_ok=True)
            outfig = os.path.join(
                outdir,f"fingerprint_{self.A}-{self.B}.png")
            plt.savefig(outfig, dpi=dpi)
        if disp:
            plt.show()


class Snapshots:
    def __init__(self, 
                 interval: float,
                 vasprun: str = 'vasprun.xml',
                 prefix: str = 'snapshots'):
        """
        Generate step-wise poscar files
        Args:
            interval (float):
                time interval in ps
            vasprun (str, optional):
                vasprun.xml file. Defaults to 'vasprun.xml'.
            prefix (str, optional):
                name of directory where poscar files will be saved. Defaults to 'snapshots'.
        """
        if os.path.isfile(vasprun):
            self.vasprun = vasprun
        else:
            print(f"Error: {vasprun} is not found.")
            sys.exit(0)
            
        if not os.path.isdir(prefix):
            os.makedirs(prefix)
            print(f'{prefix} directory is created.')
        
        self.interval = interval
        self.prefix = prefix
        
        # read vasprun.xml
        self.nsw = None
        self.pos = None
        self.potim = None
        self.digit = None
        self.lattice = None
        self.nsw_cut = None
        self.num_step = None
        self.num_atom = None
        self.atom_counts = None
        self.interval_nsw = None
        self.read_vasprun()
        
        # atom pairs
        self.pair = []
        self.get_atom_pairs()
        
        # save poscars
        self.save_poscars()
    
    def save_poscars(self):
        for i in range(self.num_step):
            poscar = os.path.join(
                self.prefix, f"POSCAR_{format(i, self.digit)}"
            )
            
            with open(poscar, "w") as f:
                f.write(f"step_{i}. generated by vachoppy.\n")
                f.write("1.0\n")
                for lat in self.lattice:
                    f.write("%.6f %.6f %.6f\n"%(lat[0], lat[1], lat[2]))
                f.write(f"{'  '.join(map(str, self.atom_counts.keys()))}\n")
                f.write(f"{'  '.join(map(str, self.atom_counts.values()))}\n")
                f.write("Direct\n")
                for coord in self.pos[i]:
                    f.write("%.6f %.6f %.6f\n"%(coord[0], coord[1], coord[2]))
                    
    def get_atom_pairs(self):
        atom_species = list(self.atom_counts.keys())
        self.pair.extend(combinations_with_replacement(atom_species, 2))
    
    def read_vasprun(self):
        v = Vasprun(self.vasprun,  
                    parse_dos=False, 
                    parse_eigen=False,
                    parse_potcar_file=False)
        
        self.potim = v.incar.get("POTIM")
        
        structure = v.final_structure
        self.lattice = structure.lattice.matrix.tolist()
        self.num_atom = v.final_structure.num_sites
        atom_symbols = [str(site.specie) for site in structure.sites]
        self.atom_counts = dict(Counter(atom_symbols))
        
        # ionic steps
        iterations = v.ionic_steps
        self.nsw = len(iterations)
        
        # convert interation : ps to iteration
        eps = 1e-9
        val = self.interval * 1000 / self.potim
        if not math.isclose(val, round(val), abs_tol=eps):
            print("ERROR: interval must be a multiple of potim.")
            sys.exit(0)
        self.interval_nsw = round(val)
        
        self.num_step = int(self.nsw / self.interval_nsw)
        self.nsw_cut = self.num_step * self.interval_nsw
        self.digit = f'0{int(np.log10(self.num_step)) + 1}'
        
        # atomic positions
        pos = np.zeros((self.nsw_cut, self.num_atom, 3), dtype=np.float64)
        for step_idx, step in enumerate(iterations[:self.nsw_cut]):
            for i in range(self.num_atom):
                pos[step_idx, i, :] = step["structure"].sites[i].frac_coords
        
        # unwrap position
        displacement = np.zeros_like(pos)
        displacement[1:, :] = np.diff(pos, axis=0)
        displacement[displacement>0.5] -= 1.0
        displacement[displacement<-0.5] += 1.0
        displacement = np.cumsum(displacement, axis=0)
        pos = displacement + pos[0]
        
        # step-averaged coordinatoin
        pos = pos.reshape(self.num_step, self.interval_nsw, self.num_atom, 3)
        pos = np.average(pos, axis=1)
        
        # wrao bacj into cell
        pos = pos - np.floor(pos)
        self.pos = pos

                   
def get_fingerprint(poscar, 
                    filename, 
                    atom_pair, 
                    Rmax, 
                    delta, 
                    sigma):
    # get fingerprint
    fingerprint = []
    for (A, B) in atom_pair:
        fingerprint_i = FingerPrint(A, B, poscar, Rmax, delta, sigma)
        fingerprint.append(fingerprint_i.fingerprint)
    fingerprint = np.array(fingerprint).reshape(1, -1).squeeze()
    
    # save fingerprint
    with open(filename, 'w') as f:
        f.write(f'# Rmax, delta, sigma = {Rmax}, {delta}, {sigma}\n')
        f.write('# pair : ')
        for (A, B) in atom_pair:
            f.write(f'{A}-{B}, ')
        f.write('\n')
        
        R = np.linspace(0, Rmax * len(atom_pair), len(fingerprint))
        for x, y in zip(R, fingerprint):
            f.write(f'  {x:2.6f}\t{y:2.6f}\n')

    return fingerprint


def phase_transition_serial(snapshots,
                            poscar_ref,
                            Rmax,
                            delta,
                            sigma,
                            prefix='fingerprints'):
    
    digit = snapshots.digit
    interval = snapshots.interval
    atom_pair = snapshots.pair
    path_poscar = snapshots.prefix
    task_size = snapshots.num_step
    
    if not os.path.isdir(prefix):
            os.makedirs(prefix)
            print(f'{prefix} directory is created.')
    
    fingerprint_ref = get_fingerprint(
        poscar=poscar_ref,
        filename=os.path.join(prefix, "fingerprint_ref.txt"),
        atom_pair=atom_pair,
        Rmax=Rmax,
        delta=delta,
        sigma=sigma
    )
    
    results = []
    for i in tqdm(list(range(task_size)),
                  bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'%(Fore.GREEN, Fore.RESET),
                  ascii=False,
                  desc=f'{RED}{BOLD}Progress{RESET}'):
        fingerprint = get_fingerprint(
                poscar=os.path.join(path_poscar, f"POSCAR_{format(i, digit)}"),
                filename=os.path.join(prefix, f"fingerprint_{format(i, digit)}.txt"),
                atom_pair=atom_pair,
                Rmax=Rmax,
                delta=delta,
                sigma=sigma
            )
        d_cos = CosineDistance(fingerprint_ref, fingerprint)
        results.append([i, d_cos])
    
    results = np.array(results)
    results[:,0] *= interval
    
    # plot cosine distance
    plt.figure(figsize=(12, 4))
    plt.scatter(results[:,0], results[:,1], s=25)
    plt.xlabel("t (ps)", fontsize=13)
    plt.ylabel('Cosine distnace', fontsize=13)
    plt.savefig('cosine_distance.png', dpi=300)
    plt.close()
    print('cosine_distance.png is created.')
    
    # save cosine distance
    with open('cosine_distance.txt', 'w') as f:
        f.write(f'# Rmax, delta, sigma = {Rmax}, {delta}, {sigma}\n')
        f.write('# pair : ')
        for (A, B) in atom_pair:
            f.write(f'{A}-{B}, ')
        f.write('\n')
        for [x, y] in results:
                f.write(f'  {x}\t{y:.6f}\n')
    print('cosine_distance.txt is created.')
        

def phase_transition_parallel(snapshots,
                              poscar_ref,
                              Rmax,
                              delta,
                              sigma,
                              prefix='fingerprints'):
    """
    snapshots : instant of Snapshots class
    """
    # time estimation
    time_i = time.time()
    
    # parameters from Snapshots
    digit = snapshots.digit
    interval = snapshots.interval
    atom_pair = snapshots.pair
    path_poscar = snapshots.prefix
    task_size = snapshots.num_step
    
    # parallelization
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    # master
    if rank == 0:
        if not os.path.isdir(prefix):
            os.makedirs(prefix)
            print(f'{prefix} directory is created.')
            
        task_queue = list(range(task_size))
        print(f"Number of snapshots : {task_size}")
        
        results = []
        completed_task = 0
        terminated_worker = 0
        active_workers = size - 1
        
        while completed_task < task_size or terminated_worker < active_workers:
            status = MPI.Status()
            worker_id, task_result = comm.recv(
                source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status
            )
            
            # tag=4 : response to termination sign
            if status.Get_tag() == 4:
                terminated_worker += 1
                continue
            
            if task_result is not None:
                completed_task += 1
                results.append(task_result)
                print(f"Progress: {completed_task}/{task_size}, " +
                      f"remaining workers = {active_workers - terminated_worker}/{active_workers}")
            
            if task_queue:
                new_task = task_queue.pop()
                # tag=1 : allocating task
                comm.send(new_task, dest=worker_id, tag=1)
            else:
                # tag=0 : termination sign
                comm.send(None, dest=worker_id, tag=0)
        
        # to confirm all workers are terminated
        while terminated_worker < active_workers:
            worker_id, _ = comm.recv(source=MPI.ANY_SOURCE, tag=4)
            terminated_worker += 1
    
    # workers
    else:
        # fingerprint of reference phase
        fingerprint_ref = get_fingerprint(
            poscar=poscar_ref,
            filename=os.path.join(prefix, "fingerprint_ref.txt"),
            atom_pair=atom_pair,
            Rmax=Rmax,
            delta=delta,
            sigma=sigma
        )
            
        # tag=2 : request for very first task
        comm.send((rank, None), dest=0, tag=2)
        
        while True:
            task = comm.recv(source=0, tag=MPI.ANY_TAG)
            
            if task is None:
                # no remaining task
                comm.send((rank, None), dest=0, tag=4)
                break
            
            # get fingerprint
            fingerprint = get_fingerprint(
                poscar=os.path.join(path_poscar, f"POSCAR_{format(task, digit)}"),
                filename=os.path.join(prefix, f"fingerprint_{format(task, digit)}.txt"),
                atom_pair=atom_pair,
                Rmax=Rmax,
                delta=delta,
                sigma=sigma
            )
            
            # get cosine distnace
            d_cos = CosineDistance(fingerprint_ref, fingerprint)
            
            # tag=3 : send result to master
            comm.send((rank, [task, d_cos]), dest=0, tag=3)
            
    if rank == 0:
        # sort results
        results.sort(key=lambda x: x[0])
        results = np.array(results)
        results[:,0] *= interval
        
        # plot cosine distance
        plt.figure(figsize=(12, 4))
        plt.scatter(results[:,0], results[:,1], s=25)
        plt.xlabel("t (ps)", fontsize=13)
        plt.ylabel('Cosine distnace', fontsize=13)
        plt.savefig('cosine_distance.png', dpi=300)
        plt.close()
        print('cosine_distance.png is created.')
        
        # save cosine distance
        with open('cosine_distance.txt', 'w') as f:
            f.write(f'# Rmax, delta, sigma = {Rmax}, {delta}, {sigma}\n')
            f.write('# pair : ')
            for (A, B) in atom_pair:
                f.write(f'{A}-{B}, ')
            f.write('\n')
            for [x, y] in results:
                 f.write(f'  {x}\t{y:.6f}\n')
        print('cosine_distance.txt is created.')
        
        time_f = time.time()
        print(f"total time taken = {time_f - time_i} s")
        
        
        
# class Snapshots:
#     def __init__(self,
#                  xdatcar: str,
#                  outcar: str,
#                  interval: float,
#                  prefix: str ='snapshots'):
#         """
#         Arguements
#         ----------
#         xdatcar : str
#             Path to XDATCAR file
#         outcar : str
#             Path to OUTCAR file
#         interval : float 
#             Time interval (ps)
#         prefix : str, optional
#             Prefix for output files
#         """
#         self._validate_file(xdatcar, "XDATCAR")
#         self._validate_file(outcar, "OUTCAR")

#         self.xdatcar = xdatcar
#         self.outcar = outcar
#         self.prefix = prefix
#         self.interval = interval
        
#         if not os.path.isdir(self.prefix):
#             os.makedirs(self.prefix)
#             print(f'{self.prefix} directory is created.')
            
#         # read outcar
#         self.potim = None
#         self.read_outcar()
        
#         if self.potim is None:
#             print(f"POTIM does not exist in OUTCAR.")
#             sys.exit(0)
            
#         self.interval_nsw = self._compute_interval_nsw(interval, self.potim)
        
#         # read xdatcar
#         self.digit = None
#         self.nsw_cut = None
#         self.num_step = None
#         self.atom_species = None
#         self.lattice_parameter = None
#         self.position = []
#         self.read_xdatcar()
        
#         self.pair = []
#         self.pair.extend(combinations_with_replacement(self.atom_species, 2))
        
#         # save poscar
#         for i in range(self.num_step):
#             filename = os.path.join(
#                 self.prefix, f"POSCAR_{format(i, self.digit)}"
#             )
#             self.save_poscar(i, filename)
            
#         print(f"AIMD snapshots are saved in {self.prefix} directory.")
    
#     def _validate_file(self, path, label="file"):
#         if not os.path.isfile(path):
#             print(f"ERROR: {label} '{path}' is not found.")
#             sys.exit(0)

#     def _compute_interval_nsw(self, interval, potim):
#         eps = 1e-9
#         val = interval * 1000 / potim
#         if math.isclose(val, round(val), abs_tol=eps):
#             return int(round(val))
#         print("ERROR: interval must be a multiple of potim.")
#         sys.exit(0)
                
#     def read_outcar(self):
#         with open(self.outcar, 'r') as f:
#             for line in f:
#                 if 'POTIM' in line:
#                     self.potim = float(line.split()[2])
#                     break
        
#     def read_xdatcar(self):
        
#         with open(self.xdatcar, 'r') as f:
#             lines = np.array([s.strip() for s in f])
        
#         nsw = find_last_direct_line(self.xdatcar)
#         num_step = int(nsw / self.interval_nsw)
#         nsw_cut = num_step * self.interval_nsw
        
#         lattice_parameter = np.array(
#             [s.split() for s in lines[2:5]], dtype=np.float64
#         ) * np.float64(lines[1])
        
#         atom_species = np.array(lines[5].split())
#         num_atoms = np.array(lines[6].split(), dtype=np.int32)
#         num_atoms_tot = np.sum(num_atoms)
        
#         digit = int(np.log10(num_step)) + 1
#         digit = f'0{digit}'
        
#         for i, spec in enumerate(atom_species):           
#             atom = {}
#             atom['species'] = spec
#             atom['num'] = num_atoms[i]
            
#             traj = np.zeros((atom['num'], num_step, 3)) 

#             for j in range(atom['num']):
#                 start = np.sum(num_atoms[:i]) + j + 8
#                 end = 7 + nsw_cut * (num_atoms_tot + 1)
#                 step = num_atoms_tot + 1
#                 coords = [s.split() for s in lines[start:end:step]]
#                 coords = np.array(coords, dtype=float)
                
#                 displacement = np.zeros_like(coords)
#                 displacement[0,:] = 0
#                 displacement[1:,:] = np.diff(coords, axis=0)

#                 # correction for periodic boundary condition
#                 displacement[displacement>0.5] -= 1.0
#                 displacement[displacement<-0.5] += 1.0
#                 displacement = np.cumsum(displacement, axis=0)
#                 coords = coords[0] + displacement

#                 # averaged coordination
#                 coords = coords.reshape(num_step, self.interval_nsw, 3)
#                 coords = np.average(coords, axis=1)

#                 # wrap back into cell
#                 coords = coords - np.floor(coords)
#                 traj[j] = coords
#             atom['traj'] = traj
            
#             self.position += [atom]
        
#         self.digit = digit
#         self.nsw_cut = nsw_cut
#         self.num_step = num_step
#         self.atom_species = atom_species
#         self.lattice_parameter = lattice_parameter
            
#     def save_poscar(self, 
#                     step, 
#                     filename):
        
#         with open(filename, 'w') as f:
#             f.write(f"step_{step}. generated by vachoppy.\n")
#             f.write("1.0\n")

#             for lat in self.lattice_parameter:
#                 f.write("%.6f %.6f %.6f\n"%(lat[0], lat[1], lat[2]))

#             for atom in self.position:
#                 f.write(f"{atom['species']} ")
#             f.write('\n')
            
#             for atom in self.position:
#                 f.write(f"{atom['num']} ")
#             f.write('\n')
            
#             f.write("Direct\n")
            
#             for atom in self.position:
#                 for traj in atom['traj'][:,step,:]:
#                     f.write("%.6f %.6f %.6f\n"%(traj[0], traj[1], traj[2]))
