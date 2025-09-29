import os
import sys
import json
import math
import copy   
import numpy as np
import matplotlib.pyplot as plt

from PIL import Image
from tqdm import tqdm
from colorama import Fore
from tabulate import tabulate

from collections import defaultdict
from itertools import permutations

from pymatgen.core import Structure
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# For Arrow3D
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

# color map for tqdm
BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m' # Green color
RED = '\033[91m'   # Red color
RESET = '\033[0m'  # Reset to default color

class Arrow3D(FancyArrowPatch):
    def __init__(self, 
                 xs, 
                 ys, 
                 zs, 
                 *args, 
                 **kwargs):
        """
        helper class to drqw 3D arrows
        """
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs
        
    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        return np.min(zs)


class Lattice:
    def __init__(self, 
                 poscar_lattice, 
                 symbol,
                 rmax=3.25,
                 tol=1e-3,
                 tolerance=1e-3,
                 verbose=False):
        # read arguments
        self.poscar_lattice = poscar_lattice
        self.symbol = symbol
        self.rmax = rmax
        self.tol = tol
        self.tolerance = tolerance
        self.verbose = verbose
        
        # check error
        if not os.path.isfile(self.poscar_lattice):
            sys.exit(f"{self.poscar_lattice} is not found.")
            
        with open(self.poscar_lattice, 'r', encoding='utf-8') as f:
            contents = f.read()
            self.structure = Structure.from_str(contents, fmt='poscar')
        
        if not any(site.specie.symbol==self.symbol for site in self.structure):
            sys.exit(f"{self.symbol} is not in {self.poscar_lattice}")
        
        # attributes    
        self.path = []
        self.path_name = []
        self.site_name = None
        self.lattice_point = None
        self.lattice_parameter = self.structure.lattice.matrix
        self.find_hopping_path()
        
        # unknown path
        self.path_unknown = []
        
        # summary path
        if self.verbose:
            self.summary()
        
    def find_hopping_path(self):
        # find inequivalent sites
        sga = SpacegroupAnalyzer(self.structure)
        sym_structure = sga.get_symmetrized_structure()
        non_eq_sites = sym_structure.equivalent_sites
        non_eq_sites = [
            site_group for site_group in non_eq_sites if site_group[0].specie.symbol==self.symbol
            ]
        index = []
        for sites in non_eq_sites:
            index_sites = []
            for site in sites:
                coords = site.coords
                for i, _site in enumerate(self.structure.sites):
                    if np.linalg.norm(coords - _site.coords) < self.tolerance:
                        index_sites.append(i)
            index.append(index_sites)
        # index = np.array(index, dtype=int)
        
        # save site names
        self.site_name = [f"site{i+1}" for i in range(len(index))]
        
        # save lattice points
        self.lattice_point = []
        for i in range(min(map(min,index)), max(map(max,index))+1):
            for j, index_j in enumerate(index):
                if i in index_j:
                    site_i = j+1
                    break
            point = {}
            point['site'] = f"site{site_i}"
            point['coord'] = self.structure[i].frac_coords
            point['coord_C'] = self.structure[i].coords
            self.lattice_point.append(point)
            
        # find hopping paths
        nn_finder = VoronoiNN(tol=self.tol)
        self.path, self.path_name = [], []
        for i, idx in enumerate([index_i[0] for index_i in index]):
            paths_idx = []
            distances = np.array([], dtype=float)
            site_init = f"site{i+1}"
            neighbors = nn_finder.get_nn_info(self.structure, idx)
            neighbors = [
                neighbor for neighbor in neighbors if neighbor['site'].specie.symbol==self.symbol
                ]
            for neighbor in neighbors:
                distance = self.structure[idx].distance(neighbor['site'])
                if distance < self.rmax:
                    for j, index_j in enumerate(index):
                        if neighbor['site_index'] in index_j:
                            site_final = j+1
                            break
                    site_final = f"site{site_final}"
                    path_index = np.where(abs(distances - distance) < self.tolerance)[0]
                    if len(path_index) == 0:
                        path = {}
                        path['site_init'] = site_init
                        path['site_final'] = site_final
                        path['distance'] = float(distance)
                        path['z'] = 1
                        path['coord_init'] = self.structure[idx].frac_coords
                        path['coord_final'] = neighbor['site'].frac_coords
                        paths_idx.append(path)
                        distances = np.append(distances, distance)
                        self.path_name.append(f"{chr(i+65)}{len(paths_idx)}")
                    else:
                        paths_idx[path_index[0]]['z'] += 1
            self.path += paths_idx
        self.path = sorted(self.path, key=lambda x: (x['site_init'], x['distance']))
        self.path_name = sorted(self.path_name)
        for path, name in zip(self.path, self.path_name):
            path['name'] = name
    
    def summary(self):
        print(f"Number of inequivalent sites for {self.symbol} : {len(self.site_name)}")
        print(f"Number of inequivalent paths for {self.symbol} : {len(self.path_name)} (Rmax = {self.rmax:.2f} Å)")
        print('')
        print('Path information')
        headers = ['name', 'init', 'final', 'a(Å)', 'z', 'coord_init', 'coord_final']
        data = [
            [
                path['name'], 
                path['site_init'], 
                path['site_final'], 
                f"{path['distance']:.5f}", 
                path['z'],
                f"[{path['coord_init'][0]:.5f} {path['coord_init'][1]:.5f} {path['coord_init'][2]:.5f}]", 
                f"[{path['coord_final'][0]:.5f} {path['coord_final'][1]:.5f} {path['coord_final'][2]:.5f}]"
            ] for path in self.path
        ]
        print(tabulate(data, headers=headers, tablefmt="simple"))


class Trajectory:
    def __init__(self,
                 interval: float,
                 num_vac: int,
                 lattice,
                 pos_file: str,
                 force_file: str,
                 cond_file: str,
                 verbose: bool = True):
        
        self._validate_file(pos_file)
        self._validate_file(force_file)
        self._validate_file(cond_file)
        
        self.interval = interval
        self.num_vac = num_vac
        self.vervose = verbose
        self.cmap = self._custrom_cmap()
        
        # md conditions
        self.nsw = None  # nsw <- int(nsw / nblock)
        self.temp = None
        self.potim = None  # potim <- potim * nblock
        self.nblock = None
        self.symbol = None
        self.nsw_cut = None
        self.num_step = None
        self.num_atom = None
        self.interval_nsw = None
        self.read_cond_file(cond_file)
        
        # lattice information
        self.lattice = lattice
        self.lattice_point = None
        self.lattice_point_C = None
        self.lattice_parameter = None
        self.num_lattice_point = None
        self._init_lattice_point(lattice)
        
        self.occupation = None
        self.trace_arrows = None
        self.atomic_trajectory(pos_file, force_file)
        
        # trace arrows
        self.trace_arrows = None
        self.get_trace_arrows()
        
        # vacancy trajectory
        self.hopping_sequence = {}
        self.vacancy_trajectory_index = {}
        self.vacancy_trajectory_coord_C = {}
        self.transient_vacancy = {}
        self.get_vacancy_trajectory()
    
    def nearest_lattice_points(self, coords):  
        # coords: (num_atom, 3)
        diff = self.lattice_point[None, :, :] - coords[:, None, :]
        diff = diff - np.floor(diff)
        diff[diff > 0.5] -= 1.0
        diff[diff < -0.5] += 1.0
        diff_cart = np.tensordot(diff, self.lattice_parameter, axes=([2], [1]))
        norm = np.linalg.norm(diff_cart, axis=2)
        return np.argmin(norm, axis=1).astype(np.int16) 

    def atomic_trajectory(self, 
                          pos_file, 
                          force_file):
        pos = np.load(pos_file, mmap_mode='r')
        force = np.load(force_file, mmap_mode='r')
        
        check_init = False
        occupation = np.zeros((self.num_step, self.num_atom), dtype=np.int16)
        for i in range(self.num_step):
            # read step-wise data
            start = i * self.interval_nsw
            end = start + self.interval_nsw
            pos_chunk = np.average(pos[start:end], axis=0) # (num_atom, 3)
            force_chunk = np.average(force[start:end], axis=0) # (num_atom, 3)
            
            # proximity-based occupation
            occupation_i = self.nearest_lattice_points(pos_chunk)
            
            # initial occupation
            if not check_init:
                if len(set(occupation_i)) == self.num_atom:
                    for j in range(i+1):
                        occupation[j] = occupation_i
                    check_init = True
                continue
            
            # TS criterion
            indices_move_atom = np.where(occupation_i != occupation[i-1])[0]
            
            for index in indices_move_atom:
                site_init = occupation[i-1][index]
                site_final = occupation_i[index]
                force_atom = force_chunk[index]
                
                p_init = self.lattice_point[site_init]
                p_final = self.lattice_point[site_final]
                p_atom = pos_chunk[index]
                
                r_init = self.displacement_PBC(p_atom, p_init)
                r_final = self.displacement_PBC(p_atom, p_final)
                
                eps = 1e-12
                norm_f = np.linalg.norm(force_atom)
                norm_init = np.linalg.norm(r_init)
                norm_final = np.linalg.norm(r_final)
                
                if norm_f < eps or norm_init < eps:
                    cos_init = np.nan
                else:
                    cos_init = np.dot(force_atom, r_init) / (norm_f * norm_init)

                if norm_f < eps or norm_final < eps:
                    cos_final = np.nan
                else:
                    cos_final = np.dot(force_atom, r_final) / (norm_f * norm_final)
                    
                if np.isnan(cos_init) or np.isnan(cos_final):
                    print(f"WARNING: NaN in cos_init/final at step {i}")
                    
                if cos_init > cos_final:
                    occupation_i[index] = site_init
            occupation[i] = occupation_i
        self.occupation = occupation.T

    def get_trace_arrows(self):
        change_in_occ = np.diff(self.occupation, axis=1)
        move_atom, move_step = np.where(change_in_occ != 0)
        move_step += 1
        
        self.trace_arrows = {}
        for step, atom in zip(move_step, move_atom):
            arrow = {}
            arrow['c'] = self.cmap[(atom)%len(self.cmap)]
            arrow['lattice_point'] = [
                self.occupation[atom][step-1],
                self.occupation[atom][step]
            ]
            arrow['p'] = np.vstack((
                self.lattice_point_C[self.occupation[atom][step-1]],
                self.lattice_point_C[self.occupation[atom][step]]
            ))
            
            if step in self.trace_arrows.keys():
                self.trace_arrows[step].append(arrow)
            else:
                self.trace_arrows[step] = [arrow]
        
        for step in range(self.num_step):
            if not step in self.trace_arrows.keys():
                self.trace_arrows[step] = []              

    def get_vacancy_trajectory(self): 
        # initiallization
        self.hopping_sequence = {}
        self.vacancy_trajectory_index = {}
        self.vacancy_trajectory_coord_C = {}
        self.transient_vacancy = {0: np.array([], dtype=np.int16)}
        step_transient = {0: False}
        lattice_site = np.arange(self.num_lattice_point)
        
        # helper method to find vacancy paths
        def trace_paths(site_init, site_final, paths):
            path_map = defaultdict(list)
            for to_site, from_site in paths:
                path_map[from_site].append(to_site)

            site_final_set = set(site_final)
            candidate_routes = {s: [] for s in site_init}

            # Step 1: Collect all candidate paths for each site_init
            for s_init in site_init:
                stack = [(s_init, [s_init])]
                while stack:
                    current, route = stack.pop()
                    if current in site_final_set and not (len(route) == 1 and current == s_init):
                        candidate_routes[s_init].append(route)
                    for next_site in path_map.get(current, []):
                        if next_site not in route:
                            stack.append((next_site, route + [next_site]))

            # Step 2: Try all permutations of site_init to avoid final-site duplication
            for ordering in permutations(site_init):
                used_finals = set()
                results = []
                used_paths = set()
                for s in ordering:
                    found = False
                    for route in candidate_routes[s]:
                        if route[-1] not in used_finals:
                            results.append(route)
                            used_finals.add(route[-1])
                            # Track used path segments
                            for i in range(len(route) - 1):
                                # match input path format [to, from]
                                used_paths.add((route[i+1], route[i]))
                            found = True
                            break
                    if not found:
                        results.append(None)

                if len(results) == len(site_init) and None not in results:
                    reordered = [None] * len(site_init)
                    for i, s in enumerate(ordering):
                        idx = site_init.index(s)
                        reordered[idx] = results[i]
                    # Identify unused paths
                    path_set = set(map(tuple, paths))
                    unused_paths = list(path_set - used_paths)
                    return reordered, unused_paths

            # Fallback
            fallback = []
            used_paths = set()
            used_finals = set()
            for s in site_init:
                found = False
                for route in candidate_routes[s]:
                    if route[-1] not in used_finals:
                        fallback.append(route)
                        used_finals.add(route[-1])
                        for i in range(len(route) - 1):
                            used_paths.add((route[i+1], route[i]))
                        found = True
                        break
                if not found:
                    fallback.append(None)

            path_set = set(map(tuple, paths))
            unused_paths = list(path_set - used_paths)
            return fallback, unused_paths

        # vacancy trajectory at very first steps
        step_init = 0
        while True:
            site_vac = np.setdiff1d(lattice_site, self.occupation[:, step_init])
            if len(site_vac) == self.num_vac:
                break
            step_init += 1
        
        for step in range(step_init + 1):
            self.vacancy_trajectory_index[step] = copy.deepcopy(site_vac)
            self.vacancy_trajectory_coord_C[step] = self.lattice_point_C[site_vac]
        
        # find vacancy trajectory
        for step in range(step_init+1, self.num_step):
            site_vac_new = np.setdiff1d(lattice_site, self.occupation[:, step])
            
            # check presence of transient vacancy
            step_transient[step] = True if len(site_vac_new) > self.num_vac else False
            
            # check hopping events
            site_init = np.setdiff1d(site_vac, site_vac_new)
            site_final = np.setdiff1d(site_vac_new, site_vac)
            
            # no hop
            if len(site_init) == 0:
                self.vacancy_trajectory_index[step] = copy.deepcopy(site_vac)
                self.vacancy_trajectory_coord_C[step] = self.lattice_point_C[site_vac]
            
            # hops occur
            else:
                # check transient vacancies
                loop = 1
                site_transient = []
                paths= [arrow['lattice_point'] for arrow in self.trace_arrows[step]]
                while step_transient[step - loop]:
                    paths += [arrow['lattice_point'] for arrow in self.trace_arrows[step-loop]]
                    site_transient += list(self.transient_vacancy[step-loop])
                    loop += 1
                
                # get hoppigg paths
                site_final = np.array(list(set(list(site_final) + site_transient)))
                path_connect, unused_path = trace_paths(list(site_init), site_final, paths)
                
                for i, site in enumerate(site_init):
                    # find path index
                    path_index = None
                    for j, p in enumerate(path_connect):
                        if p is not None and p[0] == site:
                            path_index = j
                            break
                    
                    if path_index is None:
                        print(f"ERROR: Fail to find vacancy trajectory (step = {step}).")
                        print(f"site_init : {site_init}")
                        print(f"site_final : {site_final}")
                        print(f"paths : {paths}")
                        print(f"path_connect : {path_connect}")
                        sys.exit(0)
                    
                    # update vacancy site
                    site_vac[list(site_vac).index(site)] = path_connect[path_index][-1]
                    site_final = site_final[site_final != path_connect[path_index][-1]]
                        
                # compare 'updated site_vac' and 'site_vac_new' (correction for transient vac)
                site_remain = np.setdiff1d(site_vac, site_vac_new)
                if len(site_remain) > 0:
                    site_unexpect = np.setdiff1d(site_vac_new, site_vac)
                    path_unexpect, _ = trace_paths(list(site_remain), site_unexpect, unused_path)
                    
                    for i, site in enumerate(site_remain):
                        site_vac[list(site_vac).index(site)] = path_unexpect[i][-1]
                        site_final = site_final[site_final != path_unexpect[i][-1]]
                        
                        for path in path_connect:
                            if path is not None and path[-1] == site:
                                path.append(path_unexpect[i][-1])
                                break
                        
                # update vacancy trajectory
                self.hopping_sequence[step] = copy.deepcopy(path_connect)
                self.vacancy_trajectory_index[step] = copy.deepcopy(site_vac)
                self.vacancy_trajectory_coord_C[step] = self.lattice_point_C[site_vac]
                
            # update trasient vacancy
            self.transient_vacancy[step] = site_final
                   
    def _custrom_cmap(self):
        cmap = [
            'blue',
            'red',
            'teal',
            'indigo',
            'lime',
            'darkgoldenrod',
            'cyan',
            'hotpink',
            'dodgerblue',
            'dimgray',
            'forestgreen',
            'slateblue'
        ]
        return cmap  

    def _validate_file(self, file):
        if not os.path.isfile(file):
            print(f"ERROR: '{file}' is not found.")
            sys.exit(0)

    def _init_lattice_point(self, lattice):
        if lattice.symbol != self.symbol:
            print(f"Error: unmatched atomic symbol")
            sys.exit(0)
        
        self.lattice_point = np.array(
            [p['coord'] for p in lattice.lattice_point], dtype=np.float64
        )
        self.lattice_point_C = np.array(
            [p['coord_C'] for p in lattice.lattice_point], dtype=np.float64
        )
        self.num_lattice_point = len(self.lattice_point)
        self.lattice_parameter = lattice.lattice_parameter
        
    def read_cond_file(self, cond_file):
        with open(cond_file, "r") as f:
            condition = json.load(f)
        
        self.nsw = condition["nsw"]
        self.potim = condition["potim"]
        self.nblock = condition["nblock"]
        self.symbol = condition["symbol"]
        self.temp = condition["temperature"]
        self.num_atom = condition["atom_counts"][self.symbol]
        
        # effective nsw and potim
        self.nsw = int(self.nsw / self.nblock)
        self.potim *= self.nblock
        
        eps = 1e-9
        val = self.interval * 1000 / self.potim
        if math.isclose(val, round(val), abs_tol=eps):
             self.interval_nsw = int(round(val))
        else:
            print("ERROR: interval must be a multiple of 'potim * nblock'.")
            sys.exit(0)
            
        self.num_step = int(self.nsw / self.interval_nsw)
        self.nsw_cut = self.num_step * self.interval_nsw
    
    def distance_PBC(self, 
                     coord1: list, 
                     coord2: list):
        """
        Argument
        --------
        coord1 : list
            Initial coordinate in fraction (one point or multiple points)
        coord2 : list
            Final coordinate in fraction    
        """
        distance = coord1 - coord2
        distance[distance>0.5] -= 1.0
        distance[distance<-0.5] += 1.0

        if coord1.ndim == 1:
            return np.sqrt(np.sum(np.dot(distance, self.lattice_parameter)**2))
        else:
            return np.sqrt(np.sum(np.dot(distance, self.lattice_parameter)**2,axis=1))
    
    def displacement_PBC(self, 
                         r1, 
                         r2):
        disp = r2 - r1
        disp[disp > 0.5] -= 1.0
        disp[disp < -0.5] += 1.0
        return np.dot(disp, self.lattice_parameter)
    
    def plot_lattice(self, ax, label=False):
        coord_origin = np.zeros((1, 3))
        
        def plot_edge(start, end):
            edge = np.concatenate((start.reshape(1, 3), end.reshape(1, 3)), axis=0).T
            ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')

        a, b, c = self.lattice_parameter
        edges = [
            (coord_origin, a), (coord_origin, b), (coord_origin, c),
            (a + b, a), (a + b, b),
            (b + c, b), (b + c, c),
            (c + a, c), (c + a, a),
            (a + b + c, a + b), 
            (a + b + c, b + c),
            (a + b + c, c + a)
        ]

        for start, end in edges:
            plot_edge(start, end)

        ax.scatter(*self.lattice_point_C.T, facecolor='none', edgecolors='k', alpha=0.8)
        
        if label:
            for i, coord in enumerate(self.lattice_point_C):
                ax.text(*coord.T, s=f"{i}", fontsize='xx-small')

        ax.set_xlabel('x (Å)')
        ax.set_ylabel('y (Å)')
        ax.set_zlabel('z (Å)')

    def animation(self,
                  index: list = 'all',
                  step: list = 'all',
                  vac: bool = True,
                  gif: bool = True,
                  filename: str = 'traj.gif',
                  foldername: str = 'snapshot',
                  update_alpha: float = 0.75,
                  fps: int = 5,
                  loop = 0,
                  dpi: int = 100,
                  legend: bool = False,
                  label: bool = False):
        """
        Arguments
        ---------
        index : list or 'all', optional
            Index of atoms interested in (not index of lattice sites)
        step : list or 'all', optional
            Steps interested in.
        vac : bool, optional
            If true, vacancy is displayed.
        gif : bool, optional
            If true, gif file is generated.
        filename : list, optional
            Name of gif output file.
        foldername : list, optional
            Path of directory in which the snapshots save.
        update_alpha : float, optional
            Rate of increasing transparency of trace arrows.
        """
        if not os.path.isdir(foldername):
            os.mkdir(foldername)

        if index == 'all':
            index = np.arange(self.num_atom)
        
        if str(step) == 'all':
            step = np.arange(self.num_step)
        
        files = []
        for step in tqdm(step,
                         bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'%(Fore.GREEN, Fore.RESET),
                         ascii=False,
                         desc=f'{RED}{BOLD}Progress{RESET}'):
            
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # plot lattice and lattice points
            self.plot_lattice(ax, label=label)

            # plot atoms
            for i, idx in enumerate(index):
                ax.scatter(*self.lattice_point_C[self.occupation][idx][step].T,
                           facecolor=self.cmap[idx%len(self.cmap)],
                           edgecolor='none',
                           alpha=0.8,
                           label=f"{idx}")
            
            # plot trace arrows
            alpha = 1
            for i in reversed(range(step+1)):
                for arrow in self.trace_arrows[i]:
                    arrow_prop_dict = dict(mutation_scale=10,
                                           arrowstyle='->',
                                           color=arrow['c'],
                                           alpha=alpha,
                                           shrinkA=0, 
                                           shrinkB=0)
                    disp_arrow = Arrow3D(*arrow['p'].T, **arrow_prop_dict)
                    ax.add_artist(disp_arrow)
                alpha *= update_alpha

            # plot vacancy
            if vac:
                # true vacancy
                ax.plot(*self.vacancy_trajectory_coord_C[step].T,
                        color='yellow', 
                        marker='o', 
                        linestyle='none', 
                        markersize=8, 
                        alpha=0.8, 
                        zorder=1)
                
                # trasient vacancy
                ax.plot(*self.lattice_point_C[self.transient_vacancy[step]].T,
                        color='orange', 
                        marker='o', 
                        linestyle='none', 
                        markersize=8, 
                        alpha=0.8, 
                        zorder=1)

            # make snapshot
            time = step * self.interval_nsw * self.potim / 1000 # ps
            time_tot = self.nsw_cut * self.potim / 1000 # ps
            plt.title("(%.2f/%.2f) ps, (%d/%d) step"%(time, time_tot, step, self.num_step))

            if legend:
                plt.legend()

            # save snapshot
            snapshot = os.path.join(foldername, f"snapshot_{step}.png")
            files.append(snapshot)
            plt.savefig(snapshot, dpi=dpi)
            plt.close()
        
        # make gif file
        if gif:
            print(f"Merging snapshots...")
            imgs = [Image.open(file) for file in files]
            imgs[0].save(fp=filename, format='GIF', append_images=imgs[1:], 
                         save_all=True, duration=int(1000/fps), loop=loop)
            print(f"{filename} was created.")


class TrajectoryAnalyzer:
    def __init__(self,
                 lattice,
                 trajectory,
                 tolerance=1e-3,
                 verbose=True):
        """
        Arguments
        ---------
        lattice : Lattice
            Lattice object
        trajectory : Trajectory
            Trajectory object
        tolerance : float, optional
            Tolerance for numerical accuracy
        verbose : bool, optoinal
            Verbosity flag
        """
        
        # argument
        self.trajectory = trajectory
        self.num_vac = trajectory.num_vac
        self.lattice = lattice
        self.tolerance = tolerance
        self.verbose = verbose
        
        # lattice information 
        self.path = lattice.path
        self.path_name = lattice.path_name
        self.site_name = lattice.site_name
        self.path_distance = np.array([p['distance'] for p in lattice.path])
        self.lattice_point = lattice.lattice_point
        
        # unknown path
        self.path_unknown = lattice.path_unknown
        
        # hopping history
        self.hopping_history = [[] for _ in range(self.num_vac)] # self.path_vac
        self.counts = np.zeros((self.num_vac, len(self.path_name)))
        self.hopping_statistics()
        
        # counts for unknown paths
        self.unknown_name = [f"unknown{i+1}" for i in range(len(self.path_unknown))]
        self.counts_unknown = np.zeros((self.num_vac, len(self.path_unknown)))
        self.counts_unknown_path()
        
        # random walk msd
        self.msd_rand = None # unknown path들도 포함됨
        self.random_walk_msd()
        
        # vacancy residence time (ps)
        self.residence_time = np.zeros((self.num_vac, len(self.site_name)))
        self.get_residence_time()
        
        # summary
        if verbose:
            self.summary()
        
    def hopping_statistics(self):
        for step, sequence in self.trajectory.hopping_sequence.items():
            for path in sequence:
                # check vacancy index
                index_vac = list(self.trajectory.vacancy_trajectory_index[step]).index(path[-1])
                
                for i in range(len(path) - 1):
                    index_init = path[i]
                    index_final = path[i+1]
                    
                    distance = self.trajectory.distance_PBC(
                        self.trajectory.lattice_point[index_init],
                        self.trajectory.lattice_point[index_final]
                    )
                    
                    # categorize migration
                    check_normal = False
                    for i, p in enumerate(self.path):
                        check1 = True if np.abs(p['distance'] - distance) < self.tolerance else False
                        check2 = True if p['site_init'] == self.lattice_point[index_init]['site'] else False
                        check3 = True if p['site_final'] == self.lattice_point[index_final]['site'] else False
                        
                        if check1 and check2 and check3:
                            check_normal = True
                            index_path = i
                            break
                    
                    # normal path
                    if check_normal:
                        path_info = copy.deepcopy(self.path[index_path])
                        path_info['step'] = step
                        path_info['index_init'] = index_init
                        path_info['index_final'] = index_final
                        self.hopping_history[index_vac].append(path_info)
                        self.counts[index_vac][index_path] += 1
                        
                    # unknown path
                    else:
                        check_unknown = False
                        for i, p in enumerate(self.path_unknown):
                            check1 = True if np.abs(p['distance'] - distance) < self.tolerance else False
                            check2 = True if p['site_init'] == self.lattice_point[index_init]['site'] else False
                            check3 = True if p['site_final'] == self.lattice_point[index_final]['site'] else False
                            
                            if check1 and check2 and check3:
                                check_unknown = True
                                index_path = i
                        
                        if check_unknown:
                            path_info = copy.deepcopy(self.path_unknown[index_path])
                            path_info['step'] = step
                            path_info['index_init'] = index_init
                            path_info['index_final'] = index_final
                            self.hopping_history[index_vac].append(path_info)
                        
                        else:
                            unknown_new={
                                'site_init': self.lattice_point[index_init]['site'],
                                'site_final': self.lattice_point[index_final]['site'],
                                'distance': distance,
                                'coord_init': self.lattice_point[index_init]['coord'],
                                'coord_final': self.lattice_point[index_final]['coord'],
                                'name': f"unknown{len(self.path_unknown)+1}"
                            }
                            self.path_unknown.append(copy.deepcopy(unknown_new))
                            unknown_new['step'] = step
                            unknown_new['index_init'] = index_init
                            unknown_new['index_final'] = index_final
                            self.hopping_history[index_vac].append(unknown_new)
    
    def get_residence_time(self):
        for indices in self.trajectory.vacancy_trajectory_index.values():
            for i, index in enumerate(indices):
                index_site = self.site_name.index(self.lattice.lattice_point[index]['site'])
                self.residence_time[i][index_site] += 1
        self.residence_time *= self.trajectory.interval
        
    def counts_unknown_path(self):
        for index_vac in range(self.num_vac):
            for path in self.hopping_history[index_vac]:
                if 'unknown' in path['name']:
                    index_path = int(''.join(filter(str.isdigit, path['name']))) - 1
                    self.counts_unknown[index_vac][index_path] += 1
    
    def random_walk_msd(self):
        distance_all = np.array(
            list(self.path_distance) + [p['distance'] for p in self.path_unknown]
        )
        counts_all = np.hstack((self.counts, self.counts_unknown))
        self.msd_rand = np.average(
            np.sum(distance_all ** 2 * counts_all, axis=1)
        )
                
    def summary(self):
        # path counts
        name_all = self.path_name + self.unknown_name
        counts_all = np.hstack((self.counts, self.counts_unknown))
        counts_all = np.array(counts_all, dtype=np.int32)
        vacancy_name = [f"Vacancy{i+1}" for i in range(self.num_vac)]
        print("# Path counts :")
        header = ['path'] + vacancy_name
        data = np.vstack((name_all, counts_all)).T
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        # vacancy residence time
        print("# Vacancy residence time (ps) :")
        header = ['site'] + vacancy_name
        data = np.vstack((self.site_name, self.residence_time)).T
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        # hopping sequence
        print("# Hopping sequence :")
        for i in range(self.num_vac):
            print(f"# Vacancy{i+1}")
            header = ['num', 'time (ps)', 'path', 'a (Å)', 'initial site', 'final site']
            data = [
                [
                    f"{j+1}",
                    f"{path['step'] * self.trajectory.interval:.2f}",
                    f"{path['name']}",
                    f"{path['distance']:.5f}",
                    f"{path['site_init']} [{', '.join(f'{x:.5f}' for x in self.trajectory.lattice_point[path['index_init']])}]",
                    f"{path['site_final']} [{', '.join(f'{x:.5f}' for x in self.trajectory.lattice_point[path['index_final']])}]"
                ] for j, path in enumerate(self.hopping_history[i])
            ]
            print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
            print('')

        
class Encounter:
    def __init__(self,
                 analyzer,
                 use_incomplete_encounter: bool = False,
                 verbose: bool = True):
        """
        Arguments
        ---------
        analyzer : TrajectoryAnalyzer
            TrajectoryAnalyzer object
        use_incomplete_encounter : bool, optional
            If true, incomplete encounters are used together in computations.
        verbose : bool, optional
            Verbosity flag
        """
        # arguments
        self.analyzer = analyzer
        self.trajectory = analyzer.trajectory
        self.tolerance = analyzer.tolerance
        self.use_incomplete_encounter = use_incomplete_encounter
        self.verbose = verbose
        
        # path information
        self.path = analyzer.path + analyzer.path_unknown
        self.path_name = [p['name'] for p in self.path]
        
        # unwrapped vacancy trajectory
        self.vacancy_coord_unwrap = None
        self.unwrapped_vacancy_trajectory()
        
        # encounter
        self.encounter_complete = []
        self.encounter_in_process = []
        self.encounter_PBC()
        self.num_encounter_complete = len(self.encounter_complete)
        self.num_encounter_incomplete = len(self.encounter_in_process)
        
        # encounter information
        self.msd = None
        self.path_count = None
        self.f_cor = None
        self.path_distance = None
        
        if use_incomplete_encounter:
            self.encounter_complete += self.encounter_in_process
        self.num_encounter = len(self.encounter_complete)
        
        if self.num_encounter == 0:
            if verbose:
                print("All encounters are incompleted")
        else:
            self.analyze_encounter()
            self.correlation_factor()
            if verbose:
                self.summary()
               
    def find_path_name(self, 
                    index_init: int, 
                    index_final: int):
        """
        Argument
        --------
        index_init : int
            Index of initial site
        index_final : int
            Index of final site
        """
        lattice_point = self.trajectory.lattice.lattice_point
        
        site_init = lattice_point[index_init]['site']
        site_final = lattice_point[index_final]['site']
        
        coord_init = lattice_point[index_init]['coord']
        coord_final = lattice_point[index_final]['coord']
        distance = self.trajectory.distance_PBC(coord_init, coord_final)
        
        for path in self.path:
            check1 = True if path['site_init'] == site_init else False
            check2 = True if path['site_final'] == site_final else False
            check3 = True if abs(path['distance'] - distance) < self.tolerance else False
            
            if check1 and check2 and check3:
                path_name = path['name']
                break

        return path_name
        
    def unwrapped_vacancy_trajectory(self):
        coords = self.trajectory.lattice_point[
            np.array(
                list(self.trajectory.vacancy_trajectory_index.values())
            )
        ]
        displacement = np.zeros_like(coords)
        displacement[1:] = np.diff(coords, axis=0)
        displacement[displacement > 0.5] -= 1.0
        displacement[displacement < -0.5] += 1.0
        displacement = np.cumsum(displacement, axis=0)
        self.vacancy_coord_unwrap = coords[0] + displacement
        
    def encounter_PBC(self):
        for step, sequence in self.trajectory.hopping_sequence.items():
            for path_connect in sequence:
                # vacancy index
                index_vac = list(
                    self.trajectory.vacancy_trajectory_index[step]
                ).index(path_connect[-1])
                
                # decompose path
                trace_arrow = [
                    path_connect[i-1:i+1][::-1] for i in range(len(path_connect)-1, 0, -1)
                ]
                
                # encounter analysis
                coord_init = self.vacancy_coord_unwrap[step][index_vac]
                for path in trace_arrow:
                    # atom index
                    try:
                        index_atom = list(self.trajectory.occupation[:, step]).index(path[-1])
                        
                    # transient vacancy may cause VauleError
                    except ValueError:
                        loop = 1
                        while True:
                            if len(self.trajectory.transient_vacancy[step-loop]) == 0:
                                print(f"Atomic index is not defined (step : {step}).")
                                sys.exit(0)
                                
                            match = next(
                                (arrow for arrow in self.trajectory.trace_arrows[step-loop]
                                if path == arrow['lattice_point']),
                                None
                            )
                            
                            if match:
                                index_atom = list(
                                    self.trajectory.occupation[:, step-loop]
                                ).index(path[-1])
                                break
                            
                            loop += 1

                    # path name
                    path_name = self.find_path_name(*path)
                    
                    # unwrapped atomic coord
                    displacement = np.diff(self.trajectory.lattice_point[path], axis=0)
                    displacement[displacement > 0.5] -= 1.0
                    displacement[displacement < -0.5] += 1.0
                    coord_final = coord_init + displacement.flatten()
                    
                    # comparison with existing encounters
                    index_encounter = None
                    for i, enc in enumerate(self.encounter_in_process):
                        if enc['index_atom'] == index_atom:
                            index_encounter = i
                            break
                    
                    # case 1. no matching encounter
                    if index_encounter is None:
                        encounter = {
                            'index_atom': index_atom,
                            'index_vac': index_vac,
                            'coord_init': coord_init,
                            'coord_final': coord_final,
                            'hopping_history': [path_name]
                        }
                        self.encounter_in_process.append(encounter)
                        coord_init = coord_final
                        continue
                    
                    # matching encounter
                    encounter_match = self.encounter_in_process[index_encounter]
                    coord_encounter = encounter_match['coord_final']
                    
                    distance = np.linalg.norm(
                        np.dot(
                            coord_encounter - coord_init,
                            self.trajectory.lattice_parameter
                        )
                    )
                    
                    # case 2. exactly matching encounter
                    if distance < self.tolerance:
                        # exchange with the associated vacancy : update encoutner
                        if encounter_match['index_vac'] == index_vac:
                            encounter_match['coord_final'] = coord_final
                            encounter_match['hopping_history'].append(path_name)
                        
                        # exchange with a new vacancy : terminate encounter
                        else:
                            # terminate the existing encounter
                            self.encounter_complete.append(encounter_match.copy())
                            del self.encounter_in_process[index_encounter]
                            
                            # initiate a new encounter
                            encounter = {
                                'index_atom': index_atom,
                                'index_vac': index_vac,
                                'coord_init': coord_init,
                                'coord_final': coord_final,
                                'hopping_history': [path_name]
                            }
                            self.encounter_in_process.append(encounter)
                    
                    # case 3. PBC matching encounter:
                    else:
                        # terminate the existing encounter
                        self.encounter_complete.append(encounter_match.copy())
                        del self.encounter_in_process[index_encounter]
                        
                        # initiate a new encounter
                        encounter = {
                            'index_atom': index_atom,
                            'index_vac': index_vac,
                            'coord_init': coord_init,
                            'coord_final': coord_final,
                            'hopping_history': [path_name]
                        }
                        self.encounter_in_process.append(encounter)
                    
                    coord_init = coord_final
                               
    def analyze_encounter(self):
        displacement = []
        self.path_count = np.zeros_like(self.path_name, dtype=np.float64)
        
        for encounter in self.encounter_complete:
            # path count
            for name in encounter['hopping_history']:
                index_path = self.path_name.index(name)
                self.path_count[index_path] += 1
            
            # displacement
            disp = encounter['coord_final'] - encounter['coord_init']
            disp = np.dot(disp, self.trajectory.lattice_parameter)
            displacement.append(disp)
            
        displacement = np.array(displacement)
        self.msd = np.average(np.sum(displacement**2, axis=1))
            
    def correlation_factor(self):
        self.path_distance = np.array([path['distance'] for path in self.path])
        denominator = np.sum(self.path_distance ** 2 * (self.path_count / self.num_encounter))
        self.f_cor = self.msd / denominator

    def summary(self):
        print("# Encounter analysis")
        print(f"  Use_incomplete_encounter : {self.use_incomplete_encounter}")
        print(f"  Correlation factor : {self.f_cor:.5f}")
        print(f"  Mean squared displacement (MSD) : {self.msd:.5f} Å2")
        print(f"  Number of complete encounters : {self.num_encounter_complete}")
        print(f"  Number of incomplete encounters : {self.num_encounter_incomplete}")
        print(f"  Number of encounters used in computations: {self.num_encounter}")
        print(f"  Total counts of hopping events : {int(np.sum(self.path_count))}")
        print(f"  Mean counts of hopping evnets : {np.sum(self.path_count) / self.num_encounter:.5f}")
        print('') 
        print(f"# Pathwise counts in encounters : ")
        header = ['path', 'a (Å)', 'count', 'count/enc']
        data = [
            [
                name,
                f"{a:.5f}",
                f"{int(count)}",
                f"{count / self.num_encounter:.5f}"
            ]
            for name, a, count in zip(self.path_name, self.path_distance, self.path_count)
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')      
                    
