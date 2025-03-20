import os
import sys
import time
import copy   
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from PIL import Image
from tqdm import tqdm
from colorama import Fore
from tabulate import tabulate
from scipy.optimize import minimize_scalar

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
        
    def do_3d_projection(self, 
                       renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        return np.min(zs)


class Lattice:
    def __init__(self, 
                 poscar_lattice, 
                 symbol,
                 rmax=3.0,
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
        
        # contributions    
        self.path = []
        self.path_names = []
        self.site_names = None
        self.lat_points = None
        self.lattice = self.structure.lattice.matrix
        self.find_hopping_path()
        
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
        index = np.array(index, dtype=int)
        
        # save site names
        self.site_names = [f"site{i+1}" for i in range(len(index))]
        
        # save lattice points
        self.lat_points = []
        for i in range(np.min(index), np.max(index)+1):
            point = {}
            point['site'] = f"site{np.where(index==i)[0][0]+1}"
            point['coord'] = self.structure[i].frac_coords
            point['coord_C'] = self.structure[i].coords
            self.lat_points.append(point)
            
        # find hopping paths
        nn_finder = VoronoiNN(tol=self.tol)
        self.path, self.path_names = [], []
        for i, idx in enumerate(index[:,0]):
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
                    site_final = f"site{np.where(index==neighbor['site_index'])[0][0]+1}"
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
                        self.path_names.append(f"{chr(i+65)}{len(paths_idx)}")
                    else:
                        paths_idx[path_index[0]]['z'] += 1
            self.path += paths_idx
        self.path = sorted(self.path, key=lambda x: (x['site_init'], x['distance']))
        self.path_names = sorted(self.path_names)
        for path, name in zip(self.path, self.path_names):
            path['name'] = name
    
    def summary(self):
        print(f"Number of inequivalent sites for {self.symbol} : {len(self.site_names)}")
        print(f"Number of inequivalent paths for {self.symbol} : {len(self.path_names)} (Rmax = {self.rmax:.2f} Å)")
        print('')
        print('Path information')
        headers = ['name', 'init', 'final', 'a(Å)', 'z', 'coord_init', 'coord_final']
        data = [
            [path['name'], path['site_init'], path['site_final'], f"{path['distance']:.5f}", path['z'],
             path['coord_init'], path['coord_final']] 
            for path in self.path
            ]
        print(tabulate(data, headers=headers, tablefmt="simple"))


class Trajectory:
    def __init__(self,
                 xdatcar,
                 lattice,
                 force=None,
                 interval=1,
                 verbose=True):
        """
        xdatcar: (str) path for XDATCAR.
        lattice: trajectory.Lattice class
        interval: (int) step interval to be used in averaging.
        """
        if os.path.isfile(xdatcar):
            self.xdatcar = xdatcar
        else:
            print(f"'{xdatcar} is not found.")
            sys.exit(0)

        self.interval = interval
        self.verbose = verbose
        # color map for arrows
        self.cmap = ['b', 'c', 'black', 'deeppink', 'darkorange', 
                     'saddlebrown', 'red', 'lawngreen', 'grey', 'darkkhaki', 
                     'slateblue', 'purple', 'g']
        
        # lattice information
        self.target = lattice.symbol
        self.lat_points = np.array([d['coord'] for d in lattice.lat_points], dtype=float)
        self.lat_points_C = np.array([d['coord_C'] for d in lattice.lat_points], dtype=float)
        self.num_lat_points = len(self.lat_points)

        # read xdatcar
        self.lattice = None
        self.atom_species = None
        self.num_species = None
        self.num_atoms = None
        self.nsw = None
        self.num_step = None # (=nsw/interval)
        self.position = []
        self.idx_target = None
        self.read_xdatcar()

        # read force data
        self.force_file = force
        self.forces = None
        if self.force_file is not None:
            self.read_force()

        # trajectory of atom
        self.traj_on_lat_C = None
        self.occ_lat_point = None
        self.trajectory_on_lattice()

        # number of atoms preceding target
        self.count_before = 0 
        for i in range(self.idx_target):
            self.count_before += self.position[i]['num']

        # trajectory of vacancy
        self.idx_vac = {}
        self.traj_vac_C = {} 
        self.find_vacancy()

        # trace arrows
        self.trace_arrows = None
        self.get_trace_arrows()

        # check multi-vacancy issue
        self.multi_vac = None

    def read_xdatcar(self):
        # read xdatcar
        with open(self.xdatcar, 'r') as f:
            lines = np.array([s.strip() for s in f])

        self.lattice = np.array([s.split() for s in lines[2:5]], dtype=float)
        self.lattice *= float(lines[1])

        self.atom_species = np.array(lines[5].split())
        self.num_species = len(self.atom_species)

        self.num_atoms = np.array(lines[6].split(), dtype=int)
        num_atoms_tot = np.sum(self.num_atoms)

        self.nsw = int((lines.shape[0]-7) / (1+num_atoms_tot))

        if self.nsw % self.interval == 0:
            self.num_step = int(self.nsw / self.interval)  
        else:
            print("nsw is not divided by interval.")
            sys.exit(0)

        # save coordnation
        for i, spec in enumerate(self.atom_species):
            if self.target == spec:
                self.idx_target = i
            
            atom = {}
            atom['species'] = spec
            atom['num'] = self.num_atoms[i]
            
            # traj : mean coords (atom['num'], num_step, 3)
            # coords_C : original coords (atom['num'], nsw, 3)
            traj = np.zeros((atom['num'], self.num_step, 3)) 
            coords_C = np.zeros((atom['num'], self.nsw, 3)) 

            for j in range(atom['num']):
                start = np.sum(self.num_atoms[:i]) + j + 8
                end = lines.shape[0] + 1
                step = num_atoms_tot + 1
                coords = [s.split() for s in lines[start:end:step]]
                coords = np.array(coords, dtype=float)
                
                displacement = np.zeros_like(coords)
                displacement[0,:] = 0
                displacement[1:,:] = np.diff(coords, axis=0)

                # correction for periodic boundary condition
                displacement[displacement>0.5] -= 1.0
                displacement[displacement<-0.5] += 1.0
                displacement = np.cumsum(displacement, axis=0)
                coords = coords[0] + displacement

                # covert to cartesian coordination
                coords_C[j] = np.dot(coords, self.lattice)

                # averaged coordination
                coords = coords.reshape(self.num_step, self.interval, 3)
                coords = np.average(coords, axis=1)

                # wrap back into cell
                coords = coords - np.floor(coords)
                traj[j] = coords

            atom['coords_C'] = coords_C
            atom['traj'] = traj
            atom['traj_C'] = np.dot(traj, self.lattice)
            self.position += [atom]

    def read_force(self):
        # read force data
        with open(self.force_file, 'r') as f:
            lines = [s.strip() for s in f]
        
        # number of atoms
        num_tot = np.sum(self.num_atoms)
        num_pre = np.sum(self.num_atoms[:self.idx_target])
        num_tar = self.num_atoms[self.idx_target]
        
        # save forces
        self.forces = np.zeros((self.nsw, num_tar, 3))
        for i in range(self.nsw):
            start = (num_tot+1)*i + num_pre + 1
            end = start + num_tar
            self.forces[i] = np.array([s.split() for s in lines[start:end]], dtype=float)

        self.forces = self.forces.reshape(self.num_step, self.interval, num_tar, 3)
        self.forces = np.average(self.forces, axis=1)

    def distance_PBC(self, 
                     coord1, 
                     coord2):
        """
        coord1 and coord2 are direct coordinations.
        coord1 is one point or multiple points.
        coord2 is one point.
        return: cartesian distance
        """
        distance = coord1 - coord2
        distance[distance>0.5] -= 1.0
        distance[distance<-0.5] += 1.0

        if coord1.ndim == 1:
            return np.sqrt(np.sum(np.dot(distance, self.lattice)**2))
        else:
            return np.sqrt(np.sum(np.dot(distance, self.lattice)**2,axis=1))

    def displacement_PBC(self, 
                         r1, 
                         r2):
        disp = r2 - r1
        disp[disp > 0.5] -= 1.0
        disp[disp < -0.5] += 1.0
        return np.dot(disp, self.lattice)
            
    def trajectory_on_lattice(self):
        traj = self.position[self.idx_target]['traj']

        # distance from lattice points
        disp = self.lat_points[np.newaxis,np.newaxis,:,:] - traj[:,:,np.newaxis,:]
        disp[disp > 0.5] -= 1.0
        disp[disp < -0.5] += 1.0
        disp = np.linalg.norm(np.dot(disp, self.lattice), axis=3)

        # save trajectory on lattice
        self.occ_lat_point = np.argmin(disp, axis=2)
        self.traj_on_lat_C = self.lat_points_C[self.occ_lat_point]

    def save_trajectory(self,
                        interval_traj=1,
                        foldername='traj',
                        label=False):
        """
        interval_traj: trajectory is plotted with step interval of "interval_traj"
        folder: path to directory where traj files are saved
        label: if true, the lattice points are labelled
        """
        
        if not os.path.isdir(foldername):
            os.mkdir(foldername)
        
        for i in tqdm(range(self.position[self.idx_target]['num']),
                      bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'% (Fore.GREEN, Fore.RESET),
                      ascii=False,
                      desc=f'{RED}{BOLD}save traj{RESET}'):
            
            coords = self.position[self.idx_target]['coords_C'][i][0:-1:interval_traj]
            
            # plot lattice and lattice points
            fig = plt.figure()
            ax = fig.add_subplot(111, projection = '3d')
            self.plot_lattice(ax, label=label)

            # plot trajectory
            ax.plot(*coords.T, 'b-', marker=None)
            ax.scatter(*coords[0], color='red')
            ax.scatter(*coords[-1], color='red', marker='x')

            # save plot
            filename = f'traj_{self.target}{i}.png'
            plt.title(f"Atom index = {i}")
            outfile = os.path.join(foldername, filename)
            plt.savefig(outfile, format='png')
            plt.close()

    def plot_lattice(self, 
                     ax, 
                     label=False):
        coord_origin = np.zeros([1,3])

        # plot edges
        edge = np.concatenate(
            (coord_origin, self.lattice[0].reshape(1,3)), axis=0).T
        ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')
        edge = np.concatenate(
            (coord_origin, self.lattice[1].reshape(1,3)), axis=0).T
        ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')
        edge = np.concatenate(
            (coord_origin, self.lattice[2].reshape(1,3)), axis=0).T
        ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')
        edge = np.concatenate(
            ((self.lattice[0]+self.lattice[1]).reshape(1,3), 
             self.lattice[0].reshape(1,3)), axis=0).T
        ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')
        edge = np.concatenate(
            ((self.lattice[0]+self.lattice[1]).reshape(1,3), 
             self.lattice[1].reshape(1,3)), axis=0).T
        ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')
        edge = np.concatenate(
            ((self.lattice[1]+self.lattice[2]).reshape(1,3), 
             self.lattice[1].reshape(1,3)), axis=0).T
        ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')
        edge = np.concatenate(
            ((self.lattice[1]+self.lattice[2]).reshape(1,3), 
             self.lattice[2].reshape(1,3)), axis=0).T
        ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')
        edge = np.concatenate(
            ((self.lattice[2]+self.lattice[0]).reshape(1,3), 
             self.lattice[2].reshape(1,3)), axis=0).T
        ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')
        edge = np.concatenate(
            ((self.lattice[2]+self.lattice[0]).reshape(1,3), 
             self.lattice[0].reshape(1,3)), axis=0).T
        ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')
        edge = np.concatenate(
            ((self.lattice[0]+self.lattice[1]+self.lattice[2]).reshape(1,3), 
             (self.lattice[0]+self.lattice[1]).reshape(1,3)), axis=0).T
        ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')
        edge = np.concatenate(
            ((self.lattice[0]+self.lattice[1]+self.lattice[2]).reshape(1,3), 
             (self.lattice[1]+self.lattice[2]).reshape(1,3)), axis=0).T
        ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')
        edge = np.concatenate(
            ((self.lattice[0]+self.lattice[1]+self.lattice[2]).reshape(1,3), 
             (self.lattice[2]+self.lattice[0]).reshape(1,3)), axis=0).T
        ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none')

        # plot lattice points
        ax.scatter(*self.lat_points_C.T, facecolor='none', edgecolors='k', alpha=0.8)
        if label:
            for i, coord in enumerate(self.lat_points_C):
                ax.text(*coord.T, s=f"{i}", fontsize='xx-small')
        
        # axis label
        ax.set_xlabel('x (Å)')
        ax.set_ylabel('y (Å)')
        ax.set_zlabel('z (Å)')

    def find_vacancy(self):
        idx_lat = np.arange(self.num_lat_points)
        for i in range(self.num_step):
            idx_i = np.setdiff1d(idx_lat, self.occ_lat_point[:,i])
            self.idx_vac[i] = idx_i
            self.traj_vac_C[i] = self.lat_points_C[idx_i]

    def correct_transition_state(self,
                                 step_ignore=[]):
        traj = np.transpose(self.position[self.idx_target]['traj'], (1, 0, 2))
        for i in range(1, self.num_step):
            if i in step_ignore:
                continue
            # check whether vacancy moves
            try:
                idx_pre = self.idx_vac[i-1][0]
                idx_now = self.idx_vac[i][0]
            except:
                print(f"error occured in step {i}. (correction TS)")
                print(f"step{i-1} : {self.idx_vac[i-1]}")
                print(f"step{i} : {self.idx_vac[i]}")
                continue

            if idx_pre == idx_now:
                continue

            try:
                atom = np.where((self.occ_lat_point[:,i]==idx_pre)==True)[0][0]
            except:
                if self.verbose:
                    print(f"idx_pre = {idx_pre}") 
                    print(f'error occured at step {i}.')
                    print('please correct the vacancy site yourself.')
                    print('')
                sys.exit(0)

            coord = traj[i, atom]
            force = self.forces[i, atom]

            r_pre = self.displacement_PBC(coord, self.lat_points[idx_now])
            r_now = self.displacement_PBC(coord, self.lat_points[idx_pre])

            d_site = np.linalg.norm(r_now - r_pre)
            d_atom = np.linalg.norm(r_pre)

            cond1 ,cond2 = False, False
            # condition 1
            if d_atom > d_site:
                cond1 = True
            
            # condition 2
            norm_force = np.linalg.norm(force)
            norm_r_pre = np.linalg.norm(r_pre)
            norm_r_now = np.linalg.norm(r_now)

            cos_pre = np.dot(r_pre, force) / (norm_force*norm_r_pre)
            cos_now = np.dot(r_now, force) / (norm_force*norm_r_now)

            if cos_now > cos_pre:
                cond2 = True

            # atom passes TS
            if cond1 or cond2:
                continue

            # atom doesn't pass TS
            # update trajectories
            self.idx_vac[i] = self.idx_vac[i-1]
            self.traj_vac_C[i] = self.traj_vac_C[i-1]
            self.occ_lat_point[atom][i] = self.occ_lat_point[atom][i-1]
            self.traj_on_lat_C[atom][i] = self.traj_on_lat_C[atom][i-1]

        # update trace arrows
        self.get_trace_arrows()

    def get_trace_arrows(self):
        """
        displaying trajectory of moving atom at each step.
        """
        idx_diff = np.diff(self.occ_lat_point, axis=1)
        move_atom, move_step = np.where(idx_diff != 0)

        self.trace_arrows = {}
        for step, atom in zip(move_step, move_atom):
            arrow = {}
            arrow['p'] = np.vstack((self.traj_on_lat_C[atom][step], 
                                    self.traj_on_lat_C[atom][step+1]))
            arrow['c'] = self.cmap[(atom+1)%len(self.cmap)]
            arrow['lat_points'] = [self.occ_lat_point[atom][step], 
                                   self.occ_lat_point[atom][step+1]]
            
            if step in self.trace_arrows.keys():
                self.trace_arrows[step].append(arrow)
            else:
                self.trace_arrows[step] = [arrow]
        
        # for step in range(1, self.num_step):
        for step in range(self.num_step):
            if not step in self.trace_arrows.keys():
                self.trace_arrows[step] = []

    def save_gif(self, 
                 filename, 
                 files, 
                 fps=5, 
                 loop=0):
        """
        helper method to generate gif files
        """
        imgs = [Image.open(file) for file in files]
        imgs[0].save(fp=filename, format='GIF', append_images=imgs[1:], 
                     save_all=True, duration=int(1000/fps), loop=loop)
    
    def animation(self,
                  index='all',
                  step='all',
                  vac=True,
                  gif=True,
                  filename='traj.gif',
                  foldername='snapshot',
                  update_alpha=0.75,
                  potim=2,
                  fps=5,
                  loop=0,
                  dpi=100,
                  legend=False,
                  label=False):
        """
        make gif file of atom movement
        index: (list or 'all') index of atoms interested in. (Note: not index of lat_point)
        step: (list or 'all') steps interested in.
        vac: if True, vacancy is displayed.
        gif: if True, gif file is generated.
        filename: name of gif file.
        foldername: path of directory where the snapshots save.
        update_alpha: update tranparency.
        """
        if not os.path.isdir(foldername):
            os.mkdir(foldername)

        if index == 'all':
            num_target = self.num_atoms[self.idx_target]
            index = np.arange(num_target)
        
        if str(step) == 'all':
            step = np.arange(self.num_step)
        
        files = []
        for step in tqdm(step,
                         bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'%(Fore.GREEN, Fore.RESET),
                         ascii=False,
                         desc=f'{RED}{BOLD}snapshot{RESET}'):
            
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # plot lattice and lattice points
            self.plot_lattice(ax, label=label)

            # plot points
            for i, idx in enumerate(index):
                ax.scatter(*self.traj_on_lat_C[idx-1][step].T,
                           facecolor=self.cmap[i%len(self.cmap)],
                           edgecolor='none',
                           alpha=0.8,
                           label=f"{idx}")
            
            # plot trace arrows
            alpha = 1
            for i in reversed(range(step)):
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
                ax.plot(*self.traj_vac_C[step].T,
                        color='yellow', 
                        marker='o', 
                        linestyle='none', 
                        markersize=8, 
                        alpha=0.8, 
                        zorder=1)

            # make snapshot
            time = step * self.interval * potim / 1000 # ps
            time_tot = self.nsw * potim / 1000 # ps
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
            self.save_gif(filename=filename,
                              files=files,
                              fps=fps,
                              loop=loop)
            print(f"{filename} was created.")
        
    def save_poscar(self,
                    step,
                    outdir='./',
                    vac=False,
                    expression_vac='XX'):
        """
        if vac=True, vacancy is labelled by 'XX'
        """
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        
        filename = os.path.join(outdir, f"POSCAR_{step}")
        with open(filename, 'w') as f:
            f.write(f"step_{step}. generated by vachoppy.\n")
            f.write("1.0\n")

            # write dwon lattice
            for lat in self.lattice:
                f.write("%.6f %.6f %.6f\n"%(lat[0], lat[1], lat[2]))

            # write down atom species
            for atom in self.position:
                f.write(f"{atom['species']} ")

            if vac:
                f.write(expression_vac)
            f.write("\n")

            # write down number of atoms
            for atom in self.position:
                f.write(f"{atom['num']} ")

            if vac:
                f.write(f"{len(self.idx_vac[step])}")
            f.write("\n")

            # write down coordination
            f.write("Direct\n")
            for atom in self.position:
                for traj in atom['traj'][:,step,:]:
                    f.write("%.6f %.6f %.6f\n"%(traj[0], traj[1], traj[2]))
            
            if vac:
                for idx in self.idx_vac[step]:
                    coord = self.lat_points[idx]
                    f.write("%.6f %.6f %.6f\n"%(coord[0], coord[1], coord[2])) 

    def save_traj_on_lattice(self,
                             lat_point=[],
                             step=[],
                             foldername='traj_on_lat',
                             vac=True,
                             label=False,
                             potim=2,
                             dpi=300):
        """
        lat_point: label of lattice points at the first element of step array.
        step: steps interested in
        foldername: path of directory where files save.
        vac: if True, vacancy is displayed.
        label: if True, label of lattice point is displayed.
        """
        if not os.path.isdir(foldername):
            os.makedirs(foldername, exist_ok=True)

        # obtain atom numbers
        atom_idx = []
        for idx in lat_point:
            check = np.sum(self.occ_lat_point[:,step[0]]==idx)
            if check > 1:
                print(f"there are multiple atom at site {idx} in step {step[0]}.")
                sys.exit(0)
            else:
                atom_idx += [np.argmax(self.occ_lat_point[:,step[0]]==idx)]
        
        check_first = True
        points_init = []
        for s in tqdm(step,
                      bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'%(Fore.GREEN, Fore.RESET),
                      ascii=False,
                      desc=f'{RED}{BOLD}save traj_on_lat{RESET}'):
            
            # plot lattice and lattice points
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            self.plot_lattice(ax, label=label)
            
            color_atom = {}
            for i, idx in enumerate(atom_idx):
                # plot points
                ax.scatter(*self.traj_on_lat_C[idx][s].T,
                           facecolor=self.cmap[i%len(self.cmap)],
                           edgecolor='none', 
                           alpha=0.8, 
                           label=f"{lat_point[i]}")
                lat_p = self.occ_lat_point[idx][s]
                color_atom[lat_p] = self.cmap[i%len(self.cmap)]
                
                # save initial postions
                if check_first:
                    point_init = {}
                    point_init['p'] = self.traj_on_lat_C[idx][s]
                    point_init['c'] = self.cmap[i%len(self.cmap)]
                    points_init += [point_init]
            check_first = False
        
            # plot trajectory arrow
            lat_p_atoms = [self.occ_lat_point[i][s] for i in atom_idx]
            arrows = []
            color_arrows = []
            for arrow in self.trace_arrows[s-1]:
                arrow_head = arrow['lat_points'][1]
                if  arrow_head in lat_p_atoms:
                    arrows.append(arrow)
                    color_arrows.append(color_atom[arrow['lat_points'][1]])

            alpha = 1
            for i, line in enumerate(arrows):
                arrow_prop_dict = dict(mutation_scale=10, 
                                       arrowstyle='->', 
                                       color=color_arrows[(i)%len(color_arrows)],
                                       alpha=alpha, 
                                       shrinkA=0, 
                                       shrinkB=0)
                arrow = Arrow3D(*line['p'].T, **arrow_prop_dict)
                ax.add_artist(arrow)

            # show the initial positions
            for point in points_init:
                ax.plot(*point['p'].T, 
                        c=point['c'], 
                        marker='o', 
                        linestyle='none', 
                        markersize=10, 
                        alpha=0.4, 
                        zorder=0)

            if vac:
                ax.plot(*self.traj_vac_C[s].T, 
                        color='yellow', 
                        marker='o', 
                        linestyle='none', 
                        markersize=8, 
                        alpha=0.8, 
                        zorder=1)
            
            time = s * self.interval * potim / 1000
            time_tot = self.nsw * potim / 1000
            plt.title("(%.2f/%.2f) ps, (%d/%d) step"%(time, time_tot, s, self.num_step))
            outfile = os.path.join(foldername, f"traj_{s}.png")
            plt.savefig(outfile, dpi=dpi)
            plt.close()

    def check_multivacancy(self):
        check = np.array([len(i) for i in self.idx_vac.values()])
        check = np.where((check==1) == False)[0]
        
        if len(check) > 0:
            self.multi_vac = True
            if self.verbose:
                print('multi-vacancy issue occurs:')
                print('  step :', end=' ')
                for i in check:
                    print(i, end = ' ')
                print('')

        else:
            self.multi_vac = False
            if self.verbose:
                print('vacancy is unique.')

    def update_vacancy(self,
                       step,
                       lat_point):
        """
        step: step which the user want to update the vacancy site
        lat_point: label of lattice point where vacancy exist at the step
        """
        self.idx_vac[step] = [lat_point]
        self.traj_vac_C[step] = np.array([self.lat_points_C[lat_point]])

    def correct_multivacancy(self, 
                             start=1):
        """
        correction for multi-vacancy issue
        correction starts from 'start' step
        """
        trace_lines = self.trace_arrows
        vac_site = self.idx_vac[0][0]

        for step in range(start, self.num_step):
            # when only one vacancy exist
            if len(self.idx_vac[step]) == 1:
                vac_site = self.idx_vac[step][0]
                self.update_vacancy(step, vac_site)
                continue

            # when multiple vacancies exsit
            #    when vacancy is stationary
            if vac_site in self.idx_vac[step]:
                # correct fake vacancy
                idx = np.where(self.idx_vac[step]==vac_site)[0][0]
                fake_vac = np.delete(self.idx_vac[step], idx)
                for vac in fake_vac:
                    for i in range(step-1, 0, -1):
                        if vac in self.occ_lat_point[:,i]:
                            atom = np.where(self.occ_lat_point[:,i]==vac)[0][0]
                            self.occ_lat_point[atom][step] = vac
                            self.traj_on_lat_C[atom][step] = self.lat_points_C[vac]
                            break

                # update vacancy site
                self.update_vacancy(step, vac_site)
                continue

            # when vacancy moves
            #   find connected points with vacancy
            points = [vac_site]
            while True:
                check1 = len(points)
                for dic in trace_lines[step-1]:
                    if len(list(set(points) & set(dic['lat_points']))) == 1:
                        points += dic['lat_points']
                        points = list(set(points))
                
                check2 = len(points)

                # no more connected points
                if check1 == check2:
                    break

            site = list(set(points) & set(self.idx_vac[step]))
            
            if len(site) == 1:
                vac_site = site[0]

                # correct fake vacancy
                idx = np.where(self.idx_vac[step]==vac_site)[0][0]
                fake_vac = np.delete(self.idx_vac[step], idx)
                for vac in fake_vac:
                    for i in range(step-1, 0, -1):
                        if vac in self.occ_lat_point[:,i]:
                            atom = np.where(self.occ_lat_point[:,i]==vac)[0][0]
                            self.occ_lat_point[atom][step] = vac
                            self.traj_on_lat_C[atom][step] = self.lat_points_C[vac]
                            break

                # updata vacancy site
                self.update_vacancy(step, vac_site)
                continue

            elif len(site) == 0:
                if self.verbose:
                    print("there is no connected site.")       
                    print(f"find the vacancy site for your self. (step: {step})")
                break
            
            else:
                if self.verbose:
                    print("there are multiple candidates.")       
                    print(f"find the vacancy site for your self. (step: {step})")
                break

        # update trace arrows
        self.get_trace_arrows()


class TrajectoryAnalyzer:
    def __init__(self,
                 traj,
                 lattice,
                 tolerance=1e-3,
                 verbose=True):

        self.traj = traj
        self.traj_backup = traj
        self.lattice = lattice
        self.tolerance = tolerance
        self.verbose = verbose
            
        self.path = copy.deepcopy(lattice.path)
        self.path_names = copy.deepcopy(lattice.path_names)
        self.site_names = lattice.site_names
        self.lat_points = lattice.lat_points
        
        self.step_unknown = []
        self.path_unknown = {}
        self.prefix_unknown = 'unknown'
        self.path_unknown['name'] = self.prefix_unknown
        self.path_unknown['z'] = 'Nan'
        
        self.path_vac = None
        self.idx_vac = copy.deepcopy(self.traj.idx_vac)
        
        self.num_unknown = 0
        for name in self.path_names:
            if 'unknown' in name:
                self.num_unknown = int(name.replace('unknown', ''))
                
        # determine path of vacancy
        self.get_path_vacancy()
        if len(self.step_unknown) > 0:
            self.correct_multipath()
        
        # get counts
        self.hopping_sequence = [path['name'] for path in self.path_vac]
        self.counts_tot = len(self.hopping_sequence)
        self.counts = np.array(
            [self.hopping_sequence.count(name) for name in self.path_names], dtype=float
            )
        
        # random walk MSD
        self.a = np.array([path['distance'] for path in self.path], dtype=float)
        self.msd_rand = np.sum(self.a**2 * self.counts)
        
        # total steps vacancy remained at eash site
        self.total_reside_steps = None
        self.get_total_reside_step()
        
        # print results
        if self.verbose:
            self.summary()
            
    def get_path_vacancy(self):
        step_unknown = []
        self.path_vac = []
        idx = 0
        for step in range(self.traj.num_step-1):
            coord_init = self.lat_points[self.traj.idx_vac[step][0]]['coord']
            coord_final = self.lat_points[self.traj.idx_vac[step+1][0]]['coord']
            
            # check whether vacancy moves
            distance = self.traj.distance_PBC(coord_init, coord_final)
            if distance > self.tolerance:
                site_init = self.lat_points[self.traj.idx_vac[step][0]]['site']
                path = self.determine_path(site_init, distance)
                path['step'] = step+1
                if self.prefix_unknown in path['name']:
                    step_unknown += [step+1]
                self.path_vac += [copy.deepcopy(path)]
                self.path_vac[-1]['index_init'] = int(self.idx_vac[step][0])
                self.path_vac[-1]['index_final'] = int(self.idx_vac[step+1][0])
                idx += 1
        self.step_unknown = step_unknown
                
        if len(step_unknown) > 0 and self.verbose:
            print("unknown steps are detected : ", end='')
            for step in step_unknown:
                print(step, end=' ')
            print('')
        
    def determine_path(self, site_init, distance):
        candidate = []
        for p in self.path:
            err = abs(distance - p['distance'])
            if err < self.tolerance and p['site_init']==site_init:
                candidate += [p]
                
        if len(candidate) == 0:
            # add a new unknown path
            p = self.path_unknown
            p['site_init'] = site_init
            p['distance'] = distance
            return p
        
        elif len(candidate) > 1:
            print("Two path cannot be distinguished based on distance and initial site:")
            print(f"  initial site = {site_init}, distance = {distance:.6f}")
            print('please use smaller tolerance.')
            print(f"tolerance used in this calculation = {self.tolerance:.3e}")
            sys.exit(0)
            
        else:
            return candidate[0]   

    def path_tracer(self, paths, p_init, p_goal):
        """
        find sequential paths connection p_init and p_goal
        """
        answer = [p_init]
        while True:
            if answer[-1] == p_goal:
                return answer
            intersect = []
            for i, path in enumerate(paths):
                if path[0] == p_init:
                    intersect += [i]
            if len(intersect) == 1:
                p_init = paths[intersect[0]][1]
                answer += [p_init]
            elif len(intersect) == 0:
                return []
            else:
                for i in intersect:
                    answer += self.path_tracer(paths, paths[i][1], p_goal)
                if answer[-1] != p_goal:
                    return []

    def path_decomposer(self, index):
        step = self.path_vac[index]['step']
        arrows = np.zeros((len(self.traj.trace_arrows[step-1]), 2))
        for i, dic_arrow in enumerate(self.traj.trace_arrows[step-1]):
            arrows[i] = dic_arrow['lat_points']
    
        vac_now = self.traj.idx_vac[step][0]
        vac_pre = self.traj.idx_vac[step-1][0]
        path = self.path_tracer(arrows, vac_now, vac_pre)
        path = np.array(path, dtype=int)

        # update index of lattice points occupied by vacancy
        self.idx_vac[step] = path[-2::-1]
        return path

    def correct_multipath(self):
        if self.verbose:
            print('  correction for multi-path is in progress.')
        path_unwrap = []
        for idx, path in enumerate(self.path_vac):
            step = path['step']
            if self.prefix_unknown not in path['name']:
                path_unwrap += [path]
            else:
                try:
                    p = self.path_decomposer(idx)
                    p = np.flip(p)
                except:
                    print(f"error in unwrapping path_vac[{idx}].")
                    print(f"path_vac[{idx}] : ")
                    print(self.path_vac[idx])
                    return
                if len(path) == 0:
                    continue
                if len(p) == 0:
                    p_new = {}
                    p_new['step'] = step
                    p_new['index_init'] = path['index_init']
                    p_new['index_final'] = path['index_final']
                    p_new['site_init'] = self.lat_points[p_new['index_init']]['site']
                    p_new['site_final'] = self.lat_points[p_new['index_final']]['site']
                    p_new['distance'] = path['distance']
                    check_new = True
                    for _path in self.path:
                        if _path['site_init']==p_new['site_init'] and \
                            abs(_path['distance']-p_new['distance']) < self.tolerance:
                            p_new['name'] = _path['name']
                            p_new['z'] = _path['z']
                            check_new = False
                            break
                    if check_new:
                        self.num_unknown += 1
                        p_new['name'] = self.prefix_unknown + str(self.num_unknown)
                        p_new['z'] = self.path_unknown['z']
                        coord_init = self.lat_points[p_new['index_init']]['coord']
                        coord_final = self.lat_points[p_new['index_final']]['coord']
                        p_new['coord_init'] = coord_init
                        displacement = coord_final - coord_init
                        displacement[displacement>0.5] -= 1.0
                        displacement[displacement<=-0.5] += 1.0
                        p_new['coord_final'] = coord_init + displacement
                        self.path.append(copy.deepcopy(p_new))
                        self.path_names.append(p_new['name'])     
                    path_unwrap.append(copy.deepcopy(p_new))
                    continue
                
                for i in range(len(p)-1):
                    coord_init = self.lat_points[p[i]]['coord']
                    coord_final = self.lat_points[p[i+1]]['coord']
                    site_init = self.lat_points[p[i]]['site']
                    distance = self.traj.distance_PBC(coord_init, coord_final)
                    p_new = self.determine_path(site_init, distance)
                    p_new['step'] = step
                    p_new['index_init'] = p[i]
                    p_new['index_final'] = p[i+1]
                    p_new['site_init'] = self.lat_points[p[i]]['site']
                    p_new['site_final'] = self.lat_points[p[i+1]]['site']
                    check_new = True
                    for _path in self.path:
                        if _path['site_init']==p_new['site_init'] and \
                            abs(_path['distance']-p_new['distance']) < self.tolerance:
                            p_new['name'] = _path['name']
                            p_new['z'] = _path['z']
                            check_new = False
                            break
                    if check_new:
                        self.num_unknown += 1
                        p_new['name'] = self.prefix_unknown + str(self.num_unknown)
                        p_new['z'] = self.path_unknown['z']
                        p_new['coord_init'] = coord_init
                        displacement = coord_final - coord_init
                        displacement[displacement>0.5] -= 1.0
                        displacement[displacement<=-0.5] += 1.0
                        p_new['coord_final'] = coord_init + displacement
                        self.path.append(copy.deepcopy(p_new))
                        self.path_names.append(p_new['name'])            
                    path_unwrap.append(copy.deepcopy(p_new))
        self.path_vac = path_unwrap
        
        # check unknown path
        check_unknown = []
        for p_vac in self.path_vac:
            if self.prefix_unknown in p_vac['name']:
                check_unknown += [p_vac]
        if len(check_unknown) == 0:
            if self.verbose:
                print('  correction for multi-path is done.')
                print("no unknown step remains.\n")
        else:
            if self.verbose:
                print('  correction for multi-path is done.')
                print("unknown steps remain : ", end='')
                for p in check_unknown:
                    print(p['step'], end=' ')
                print('\n')
                
    def get_total_reside_step(self):
        self.total_reside_steps = np.zeros(len(self.lattice.site_names))
        step_before = 0
        for path in self.path_vac:
            index_init = self.lattice.site_names.index(path['site_init'])
            self.total_reside_steps[index_init] += path['step'] - step_before
            step_before = path['step']
        index_final = self.lattice.site_names.index(self.path_vac[-1]['site_final'])
        self.total_reside_steps[index_final] += self.traj.num_step - self.path_vac[-1]['step']   
                
    def summary(self):
        # print counts
        print('# Hopping sequence analysis')
        header = ['path', 'count', 'init', 'final', 'a(Å)', 'z']
        data = [
            [path['name'], count, path['site_init'], path['site_final'],
             path['distance'], f"{path['z']}"] for path, count in zip(self.path, self.counts)
        ]
        data.append(['Total', np.sum(self.counts)])
        print('Path information :')
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        # print hopping sequence
        header = ['num', 'path', 'step', 'init', 'final', 'a(Å)']
        data = [
            [f"{i+1}", path['name'], f"{path['step']}", f"{path['index_init']} ({path['site_init']})", 
             f"{path['index_final']} ({path['site_final']})", f"{path['distance']:.5f}"] 
            for i, path in enumerate(self.path_vac)
        ]
        print('Hopping sequence :')
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        # total steps vacancy remained at eash site
        header = ['site', 'total steps']
        data = [
            [name, step] for name, step in zip(self.lattice.site_names, self.total_reside_steps)
        ]
        data.append(['Total', self.traj.num_step])
        print('Total steps the vacancy remained at each site :')
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        # random walk msd
        print(f"MSD for random walk process = {self.msd_rand:.5f} Å2")


class Encounter:
    def __init__(self,
                 analyzer,
                 verbose=True):
        """
        Obtain information on the encounters.
        Args:
            analyzer : instance of TrajectoryAnalyzer class
            verbose  : (default: True)
        """
        self.analyzer = analyzer
        self.traj = analyzer.traj
        self.path = copy.deepcopy(self.analyzer.path)
        self.path_names = self.analyzer.path_names
        
        # check multi-vacancy
        if self.traj.multi_vac:
            print('This method is not applicable to system with multiple vacancies.')
            sys.exit(0)
        
        # trajectory of vacancy with consideration of PBC
        self.traj_vac = {}
        self.get_traj_vac()
        
        # encounters
        self.coord_i_enc = []
        self.coord_f_enc = []
        self.path_enc = []
        self.path_enc_all = []
        self.tolerance = 0.01
        self.get_encounters()
        self.num_enc = len(self.path_enc)
        
        # correlation factor
        self.msd = 0
        self.path_counts = np.zeros(len(self.path_names))
        self.path_dist = np.zeros_like(self.path_counts)
        self.get_msd()
        self.get_counts()
        
        # correlation factor
        self.f_cor = self.msd / np.sum(self.path_dist**2 * self.path_counts)
        
        # print results
        if verbose:
            self.print_summary()

    def get_traj_vac(self):
        idx_vac = np.array([site[0] for site in self.traj.idx_vac.values()])
        step_move = np.diff(idx_vac)

        # steps where vacancy moved
        step_move = np.where(step_move != 0)[0]
        step_move += 1
        
        # path of vacancy
        path_net = idx_vac[step_move]
        
        # coords considering PBC
        coords = self.traj.lat_points[path_net]
        displacement = np.zeros_like(coords)
        displacement[1:] = np.diff(coords, axis=0)
        displacement[displacement > 0.5] -= 1.0
        displacement[displacement < -0.5] += 1.0
        displacement = np.cumsum(displacement, axis=0)
        coords = coords[0] + displacement
        
        # save net path
        for step, coord in zip(step_move, coords):
            dic = {}
            dic['index'] = idx_vac[step]
            dic['coord'] = coord
            self.traj_vac[step] = dic
        
    def update_encounters(self, step):
        """
        Update encounter coordinates and paths based on the given step.
        Args:
            step (int) : Current simulation step.
        """
        # Trace arrows at the given step
        arrows = np.array([dic['lat_points'] for dic in self.traj.trace_arrows[step-1]], dtype=int)
        
        # Path of the vacancy
        path = self.analyzer.path_tracer(arrows, self.traj.idx_vac[step][0], self.traj.idx_vac[step-1][0])

        # Get the current vacancy coordinates
        coord_vac = self.traj_vac[step]['coord']

        # Check if there are any initial encounters
        if len(self.coord_i_enc) == 0:
            updated_coord_i_enc = []
            updated_coord_f_enc = []
            updated_path_enc = []
            
            # Loop through the path and update coordinates and paths
            for i in range(len(path) - 1):
                idx_i, idx_f = path[i], path[i + 1]

                coord_i = self.traj.lat_points[idx_i]
                coord_f = self.traj.lat_points[idx_f]

                displacement = coord_f - coord_i
                displacement[displacement > 0.5] -= 1.0
                displacement[displacement < -0.5] += 1.0

                coord_new = coord_vac + displacement
                site = self.analyzer.lat_points[idx_f]['site']
                distance = self.traj.distance_PBC(coord_i, coord_f)
                path_name = self.analyzer.determine_path(site, distance)['name']

                updated_coord_i_enc.append(coord_vac)
                updated_coord_f_enc.append(coord_new)
                updated_path_enc.append([path_name])

                coord_vac = coord_new

            updated_coord_i_enc = np.array(updated_coord_i_enc)
            updated_coord_f_enc = np.array(updated_coord_f_enc)

        else:
            updated_coord_i_enc = copy.deepcopy(self.coord_i_enc)
            updated_coord_f_enc = copy.deepcopy(self.coord_f_enc)
            updated_path_enc = copy.deepcopy(self.path_enc)

            # Loop through the path and update encounters
            for i in range(len(path) - 1):
                idx_i, idx_f = path[i], path[i + 1]

                coord_i = self.traj.lat_points[idx_i]
                coord_f = self.traj.lat_points[idx_f]

                displacement = coord_f - coord_i
                displacement[displacement > 0.5] -= 1.0
                displacement[displacement < -0.5] += 1.0

                coord_new = coord_vac + displacement
                site = self.analyzer.lat_points[idx_f]['site']
                distance = self.traj.distance_PBC(coord_i, coord_f)
                path_name = self.analyzer.determine_path(site, distance)['name']

                # Check if current vacancy coordinate is in the final encounter coordinates
                check = np.linalg.norm(self.coord_f_enc - coord_vac, axis=1) < self.tolerance

                if np.any(check):
                    idx = np.where(check)[0][0]
                    updated_coord_f_enc[idx] = coord_new
                    updated_path_enc[idx].append(path_name)
                else:
                    updated_coord_i_enc = np.vstack([updated_coord_i_enc, coord_vac])
                    updated_coord_f_enc = np.vstack([updated_coord_f_enc, coord_new])
                    updated_path_enc.append([path_name])
                coord_vac = coord_new
        self.coord_i_enc = updated_coord_i_enc
        self.coord_f_enc = updated_coord_f_enc
        self.path_enc = updated_path_enc
    
    def get_encounters(self):
        for step in self.traj_vac.keys():
            self.update_encounters(step)
        for path in self.path_enc:
            self.path_enc_all += path
                
    def get_msd(self):
        displacement = self.coord_f_enc - self.coord_i_enc
        displacement = np.dot(displacement, self.traj.lattice)
        self.msd = np.average(np.sum(displacement**2, axis=1))
        
    def get_counts(self):
        for i, name in enumerate(self.path_names):
            self.path_counts[i] = self.path_enc_all.count(name)
            self.path_dist[i] = self.analyzer.path[i]['distance']
        self.path_counts /= self.num_enc
            
    def print_summary(self):
        print('\n# Encounter analysis')
        print(f"Number of encounters      : {self.num_enc}")
        print(f"Mean squared displacement : {self.msd:.5f} Å2")
        count_tot = int(np.sum(self.path_counts*self.num_enc))
        print(f"Total hopping counts      : {count_tot}")
        count_mean = np.sum(self.path_counts)
        print(f"Mean hopping counts       : {count_mean:.5f}")
        print('')
        print('Counts in encounter analysis :')
        # print('Note : It can be differ from counts from vacancy path analysis, since it based on atom not vacancy.')
        header = ['path', 'a(Å)', 'count', 'count/enc']
        data = [
            [name, f"{a:.5f}", f"{round(count*self.num_enc)}", f"{count:.5f}"] 
            for name, a, count in zip(self.path_names, self.path_dist, self.path_counts)
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        print(f"Correlation factor = {self.f_cor:.5f}")
        print("")
        

# class Parameter:
#     def __init__(self,
#                  data,
#                  interval,
#                  poscar_lattice,
#                  symbol,
#                  rmax=3.0,
#                  tol=1e-3,
#                  tolerance=1e-3,
#                  verbose=False):
        
#         self.data = data
#         self.interval = interval
#         self.poscar_lattice = poscar_lattice
#         self.symbol = symbol
#         self.rmax = rmax
#         self.tol = tol
#         self.tolerance = tolerance
#         self.verbose = verbose
        
#         self.kb = 8.61733326e-5
#         self.cmap = plt.get_cmap("Set1")
#         self.temp = np.array(self.data.temp, dtype=int)
#         self.lattice = Lattice(self.poscar_lattice,
#                                self.symbol,
#                                self.rmax,
#                                self.tol,
#                                self.tolerance)
#         self.num_path_except_unknowns = len(self.lattice.path_names)
        
#         # number of paths at each site
#         self.num_path_site = np.zeros(len(self.lattice.site_names))
#         for path in self.lattice.path:
#             self.num_path_site[self.lattice.site_names.index(path['site_init'])] += 1
            
#         # total steps vacancy remained at eash site
#         self.total_reside_steps = []
        
#         # path counts
#         self.counts = []
#         self.encounter_num = []
#         self.encounter_msd = []
#         self.encounter_counts = []
#         self.labels_failed = []
#         self.get_counts()
        
#         # correlation factor
#         self.f_ind = []
#         self.f_avg = []
#         self.f_cum = []
#         self.f0 = None
#         self.Ea_f = None
#         self.get_correlation_factor()
        
#         # efffective paramters
#         self.D_rand = None
#         self.D0_rand = None
#         self.Ea = None
#         self.tau = None
#         self.tau0 = None
#         self.a_eff = None
#         self.nu_eff =None
#         self.z_mean = None
#         self.get_effective_parameters()
        
#         if self.verbose:
#             self.summary()
        
#     def get_counts(self):
#         for i, temp in enumerate(
#             tqdm(self.temp,
#                  bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'% (Fore.GREEN, Fore.GREEN),
#                  ascii=False,
#                  desc=f"{GREEN}{BOLD}Parameter")):
            
#             # check validity of interval
#             if (self.interval * 1000) % self.data.potim[i] != 0:
#                 print(f"unvalid interval : interval should be multiple of potim (T={temp}K)")
#                 print(f"potim = {self.data.potim[i]} fs")
#                 print(f"interval = {self.interval} ps")
#                 sys.exit(0)
#             else:
#                 step_interval = int(self.interval*1000 / self.data.potim[i])
            
#             counts_i = {}
#             encounter_num_i = []
#             encounter_msd_i = []
#             encounter_counts_i = []
#             total_reside_steps_i = np.zeros(len(self.lattice.site_names))
            
#             fail_i = []
#             desc = str(int(temp))+'K'
#             for j in tqdm(range(len(self.data.label[i])),
#                           bar_format='{l_bar}{bar:20}{r_bar}{bar:-10b}',
#                           ascii=True,
#                           desc=f"{RED}{BOLD}{desc:>9s}{RESET}"):
#                 if self.verbose:
#                     print(f"# -----------------( T = {temp} K, Label = {self.data.label[i][j]} )-----------------")
#                 xdatcar = self.data.xdatcar[i][j]
#                 force = self.data.force[i][j] if self.data.force is not None else None
#                 try:
#                     traj = Trajectory(xdatcar, 
#                                       self.lattice, 
#                                       force=force, 
#                                       interval=step_interval,
#                                       verbose=False)
#                     traj.correct_multivacancy(start=1)
#                     traj.check_multivacancy()
#                 except:
#                     print(f"An error occured during trajectory analysis : {temp} K, {self.data.label[i][j]}\n")
#                     fail_i.append(self.data.label[i][j])
#                     continue
                
#                 if traj.multi_vac is True:
#                     print(f"Multi-vacancy issue occured : {temp} K, {self.data.label[i][j]}\n")
#                     fail_i.append(self.data.label[i][j])
#                     continue
                
#                 if force is not None:
#                     traj.correct_transition_state()
                    
#                 anal = TrajectoryAnalyzer(traj, self.lattice, verbose=self.verbose)
#                 for name, count in zip(anal.path_names, anal.counts):
#                     if name not in counts_i.keys():
#                         counts_i[name] = count
#                     else:
#                         counts_i[name] += count
#                 total_reside_steps_i += anal.total_reside_steps
                        
#                 # encounter 
#                 enc = Encounter(anal, verbose=self.verbose)
#                 encounter_num_i.append(enc.num_enc)
#                 encounter_msd_i.append(enc.msd)
#                 _encounter_counts = {}
#                 for name, count in zip(enc.path_names, enc.path_counts):
#                     _encounter_counts[name] = count
#                 encounter_counts_i.append(_encounter_counts)
                
#                 # update path informaiton
#                 self.lattice.path = anal.path
#                 self.lattice.path_names = anal.path_names
            
#             self.counts.append(counts_i)
#             self.total_reside_steps.append(total_reside_steps_i)
#             self.encounter_num.append(encounter_num_i)
#             self.encounter_msd.append(encounter_msd_i)
#             self.encounter_counts.append(encounter_counts_i)
#             self.labels_failed.append(fail_i)
#         self.total_reside_steps = np.array(self.total_reside_steps)
        
#         # convert dictionary to numpy
#         all_keys = sorted(set().union(*self.counts))
#         counts = np.zeros((len(self.counts), len(all_keys)))
#         for i, entry in enumerate(self.counts):
#             for j, key in enumerate(all_keys):
#                 counts[i, j] = entry.get(key, 0)
#         self.counts = counts
        
#         encounter_counts = []
#         for entry in self.encounter_counts:
#             _counts = np.zeros((len(entry), len(all_keys)))
#             for i, dic in enumerate(entry):
#                 for j, key in enumerate(all_keys):
#                     _counts[i, j] = dic.get(key, 0)
#             encounter_counts.append(_counts)
#         self.encounter_counts = encounter_counts
        
#     def get_correlation_factor(self):
#         a = np.array([path['distance'] for path in self.lattice.path], dtype=float)
#         for i in range(len(self.temp)):
#             # individual correlation factor
#             msd = np.array(self.encounter_msd[i])
#             counts = np.array(self.encounter_counts[i])
#             self.f_ind.append(msd / np.sum(a**2 * counts, axis=1))
            
#             # averaged correlation factor
#             self.f_avg.append(np.average(self.f_ind[-1]))
            
#             # cumulative correlation factor
#             num_enc = np.array(self.encounter_num[i])
#             msd_cum = np.sum(msd * num_enc) / np.sum(num_enc)
#             counts_cum = np.sum(counts * num_enc.reshape(-1,1), axis=0) / np.sum(num_enc)
#             self.f_cum.append(msd_cum / np.sum(a**2 * counts_cum))
#         self.f_cum = np.array(self.f_cum)
        
#         # Arrhenius fit
#         if len(self.temp) >= 2:
#             slop, intercept = np.polyfit(1/self.temp, np.log(self.f_cum), deg=1)
#             self.f0 = np.exp(intercept)
#             self.Ea_f = -self.kb * slop
                   
#     def get_effective_parameters(self):
#         # effective z
#         z = np.array([path['z'] for path in self.lattice.path], dtype=float)[:self.num_path_except_unknowns]
#         self.z_mean = np.sum(self.counts[:,:self.num_path_except_unknowns], axis=1)
#         self.z_mean /= np.sum(self.counts[:,:self.num_path_except_unknowns] / z, axis=1)
#         self.z_mean = np.average(self.z_mean)
        
#         # D_rand
#         a = np.array([path['distance'] for path in self.lattice.path], dtype=float)
#         a = a[:self.num_path_except_unknowns] # lattice hopping only
#         num_label = np.array([len(label) for label in self.data.label], dtype=float)
#         t = self.data.nsw * self.data.potim * num_label
#         counts = self.counts[:,:self.num_path_except_unknowns] # lattice hopping only
#         self.D_rand = np.sum(a**2 * counts, axis=1) / (6*t) * 1e-5
        
#         # Ea and D0_rand
#         slop, intercept = np.polyfit(1/self.temp, np.log(self.D_rand), deg=1)
#         self.Ea = -slop * self.kb
#         self.D0_rand = np.exp(intercept)
        
#         # effective tau0
#         self.tau = t * (1e-3) / np.sum(counts, axis=1) # ps
#         error_tau = lambda tau0: np.linalg.norm(
#             self.tau - tau0 * np.exp(self.Ea / (self.kb * self.temp))
#             )
#         result = minimize_scalar(error_tau)
#         self.tau0 = result.x # ps
        
#         # effective a
#         self.a_eff = np.sqrt(6*self.D0_rand*self.tau0) * 1e4

#     def summary(self):
#         print("# -----------------( Summary )-----------------")
#         # lattice information
#         num_paths = np.zeros(len(self.lattice.site_names))
#         for path in self.lattice.path[:self.num_path_except_unknowns]:
#             num_paths[self.lattice.site_names.index(path['site_init'])] += 1
#         print('Lattice information :')
#         print(f"  Number of sites : {len(self.lattice.site_names)}")
#         print(f"  Number of paths : ", end='')
#         for num in num_paths:
#             print(int(num), end=' ')
#         print('')
#         print(f"  Number of unknown paths : {len(self.lattice.path) - self.num_path_except_unknowns}")
#         print('')
        
#         print('Vacancy hopping paths :')
#         header = ['path', 
#                   'init',
#                   'final',
#                   'a(Å)', 
#                   'z',
#                   'coord_init', 
#                   'coord_final'
#                   ]
#         data = [
#             [path['name'], 
#              path['site_init'], 
#              path['site_final'], 
#              path['distance'], 
#              path['z'],
#              f"[{path['coord_init'][0]:.5f} {path['coord_init'][1]:.5f} {path['coord_init'][2]:.5f}]", 
#              f"[{path['coord_final'][0]:.5f} {path['coord_final'][1]:.5f} {path['coord_final'][2]:.5f}]"
#              ] 
#             for path in self.lattice.path[:self.num_path_except_unknowns]
#         ]
#         print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
#         print('')
        
#         print('Unknown paths :')
#         header = ['path', 
#                   'init', 
#                   'final', 
#                   'a(Å)', 
#                   'z',
#                   'coord_init', 
#                   'coord_final'
#                   ]
#         data = [
#             [path['name'], 
#              path['site_init'], 
#              path['site_final'], 
#              path['distance'], 
#              path['z'],
#              f"[{path['coord_init'][0]:.5f} {path['coord_init'][1]:.5f} {path['coord_init'][2]:.5f}]", 
#              f"[{path['coord_final'][0]:.5f} {path['coord_final'][1]:.5f} {path['coord_final'][2]:.5f}]"
#              ]
#             for path in self.lattice.path[self.num_path_except_unknowns:]
#         ]
#         print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
#         print('')

#         # simulation temperatures
#         print('Simulation temperatures (K) :\n', end='  ')
#         for temp in self.temp:
#             print(temp, end=' ')
#         print('')

#         # residence time
#         print('\nTime vacancy remained at each site (ps) :')
#         header = ['T (K)'] + self.lattice.site_names
#         data = [
#             [temp] + step.tolist() 
#             for temp, step in zip(self.temp, self.total_reside_steps * self.interval)
#         ]
#         data.append(['Total'] + np.sum(self.total_reside_steps * self.interval, axis=0).tolist())
#         print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
#         print('')

#         # counts
#         print('Counts of occurrences for each hopping path :')
#         header = ['T (K)'] + self.lattice.path_names
#         data = [
#             [temp] + count.tolist()
#             for temp, count in zip(self.temp, self.counts)
#         ]
#         data.append(['Total'] + np.sum(self.counts, axis=0).tolist())
#         print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
#         print('')
                
#         # effective parameters
#         print('\nEffective hopping parameters : ')
#         header = ['parameter', 'value', 'description']
#         parameter = ['Drand_0 (m2/s)',
#                      'tau0 (ps)',
#                      'Ea (eV)', 
#                      'a (Å)', 
#                      'f', 
#                      '<z>'
#                      ]
#         value = [f"{self.D0_rand:.5e}", 
#                  f"{self.tau0:.5f}", 
#                  f"{self.Ea:.5f}",
#                  f"{self.a_eff:.5f}", 
#                  f"{np.average(self.f_cum):.5f}", 
#                  f"{self.z_mean:.5f}"
#                  ]
#         description = ['pre-exponential for random walk diffusivity',
#                        'pre-exponential for residence time',
#                        'hopping barrier',
#                        'hopping distance',
#                        'correlation factor',
#                        'mean number of equivalent paths per path type']
#         data = [[p, v, d] for p, v, d in zip(parameter, value, description)]
#         print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
#         print('')
        
#         # correlation factor
#         # merge f_ind
#         label_all = list(set(sum(self.data.label, [])))
#         label_all.sort()
#         f_merged = []
#         for i, f_ind_i in enumerate(self.f_ind):
#             f_merged_i = ['-'] * len(label_all)
#             for j, label in enumerate(self.data.label[i]):
#                 f_merged_i[label_all.index(label)] = f"{f_ind_i[j]:.5f}"
#             f_merged.append(f_merged_i)
#         f_merged = [list(x) for x in zip(*f_merged)] # transpose
        
#         print('Individual correlation factors : ')
#         header = ['label'] + [f"f({int(T)}K)" for T in self.temp]
#         data = [
#             [label] + f_merged_i for label, f_merged_i in zip(label_all, f_merged)
#         ]
#         data.append(['Average']+[f"{sum(f)/len(f):.5f}" for f in self.f_ind])
#         print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
#         print('')
            
#         print('Cumulative correlation factors : ')
#         print('(Note: cumulative f is more reliable than avaerage of individual f)')
#         header = ['T (K)', 'f']
#         data = [
#             [f"{temp}", f"{f:.5f}"] for temp, f in zip(self.temp, self.f_cum)
#         ]
#         data.append(['Average', f"{np.average(self.f_cum):.5f}"])
#         print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
#         print('')
        
#         # random walk diffuison coefficient
#         print('Random walk diffusion coefficient : ')
#         header = ['T (K)', 'D_rand (m2/s)']
#         data = [
#             [f"{temp}", f"{D:.5e}"] for temp, D in zip(self.temp, self.D_rand)
#         ]
#         print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        
#         # residence time
#         print('\nResidence time : ')
#         header = ['T (K)', 'tau (ps)']
#         data = [
#             [f"{temp}", f"{tau:.5f}"] for temp, tau in zip(self.temp, self.tau)
#         ]
#         print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        
#         # check errors
#         if any(labels for labels in self.labels_failed):
#             print("\nLabels where errors occurred : ")
#             header = ['T (K)', 'Label']
#             data = [
#                 [f"{temp}", f"{labels}"] for temp, labels in zip(self.temp, self.labels_failed)
#             ]
#             print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
#         print('')
#         print("# -----------------( Finish  )-----------------")

#     def save_figures(self):
#         # D_rand
#         plt.style.use('default')
#         plt.rcParams['figure.figsize'] = (3.8, 3.8)
#         plt.rcParams['font.size'] = 11
#         fig, ax = plt.subplots()
#         for axis in ['top','bottom','left','right']:
#             ax.spines[axis].set_linewidth(1.2)
        
#         for i in range(len(self.temp)):
#             plt.scatter(1/self.temp[i], np.log(self.D_rand[i]), 
#                         color=self.cmap(i), marker='s', s=50, label=str(int(self.temp[i])))
#         slop, intercept = np.polyfit(1/self.temp, np.log(self.D_rand), deg=1)
#         x = np.linspace(np.min(1/self.temp), np.max(1/self.temp), 100)
#         plt.plot(x, slop*x + intercept, 'k:', linewidth=1)
#         plt.xlabel('1/T (1/K)', fontsize=14)
#         plt.ylabel(r'ln $D_{rand}$ ($m^{2}$/s)', fontsize=14)
#         num_data = len(self.D_rand)
#         ncol = int(np.ceil(num_data / 5))
#         plt.legend(loc='best', fancybox=True, framealpha=1, edgecolor='inherit',
#                    ncol=ncol, labelspacing = 0.3, columnspacing=0.5, borderpad=0.2, handlelength=0.6,
#                    fontsize=11, title='T (K)', title_fontsize=11)
#         if num_data >= 3:
#             x = np.array([self.temp[0], self.temp[int(num_data/2)], self.temp[-1]])
#         else:
#             x = self.temp
#         x_str = [f"1/{int(T)}" for T in x]
#         x = 1/x
#         plt.xticks(x, x_str)
#         plt.savefig('D_rand.png', transparent=False, dpi=300, bbox_inches="tight")
#         plt.close()
        
#         # tau
#         plt.style.use('default')
#         plt.rcParams['figure.figsize'] = (3.8, 3.8)
#         plt.rcParams['font.size'] = 11
#         fig, ax = plt.subplots()
#         for axis in ['top','bottom','left','right']:
#             ax.spines[axis].set_linewidth(1.2)
            
#         for i, temp in enumerate(self.temp):
#             ax.bar(temp, self.tau[i], width=50, edgecolor='k', color=self.cmap(i))
#             ax.scatter(temp, self.tau[i], marker='o', edgecolors='k', color='k')
#         x = np.linspace(0.99*self.temp[0], 1.01*self.temp[-1], 1000)
#         ax.plot(x, self.tau0 * np.exp(self.Ea/(self.kb*x)), 'k:')
#         plt.xlabel('T (K)', fontsize=14)
#         plt.ylabel(r'$\tau$ (ps)', fontsize=14)
#         plt.savefig('tau.png', transparent=False, dpi=300, bbox_inches="tight")
#         plt.close()
        
#         # correlation factor
#         plt.style.use('default')
#         plt.rcParams['figure.figsize'] = (3.8, 3.8)
#         plt.rcParams['font.size'] = 11
#         fig, ax = plt.subplots()
#         for axis in ['top','bottom','left','right']:
#             ax.spines[axis].set_linewidth(1.2)
            
#         for i in range(len(self.temp)):
#             ax.scatter(self.temp[i], self.f_cum[i], color=self.cmap(i), marker='s', s=50)
#         plt.ylim([0, 1])
#         plt.xlabel('T (K)', fontsize=14)
#         plt.ylabel(r'$f$', fontsize=14)

#         # inset graph
#         axins = ax.inset_axes([1.125, 0.615, 0.35, 0.35])
#         x_ins = np.linspace(1/self.temp[-1], 1/self.temp[0], 100)
#         axins.plot(x_ins, -(self.Ea_f/self.kb) * x_ins + np.log(self.f0), 'k:')
#         for i in range(len(self.temp)):
#             axins.scatter(1/self.temp[i], np.log(self.f_cum[i]), color=self.cmap(i), marker='s')
#         axins.set_xlabel('1/T', fontsize=12)
#         axins.set_ylabel(r'ln $f$', fontsize=12)
#         axins.set_xticks([])
#         axins.set_yticks([])
#         plt.savefig('f_cor.png', transparent=False, dpi=300, bbox_inches="tight")
#         plt.close()
#         print('')
       

# class CorrelationFactor(Parameter):
#     def __init__(self,
#                  data,
#                  interval,
#                  poscar_lattice,
#                  symbol,
#                  temp,
#                  label='all',
#                  rmax=3.0,
#                  tol=1e-3,
#                  tolerance=1e-3,
#                  verbose=False):
        
#         # check validity of temp and label
#         if int(temp) in data.temp:
#             self.temp = temp
#         else:
#             print(f"{temp}K is not valid.")
#             sys.exit(0)
            
#         self.label = data.label[np.where(data.temp==self.temp)[0][0]] if label=='all' else label
#         for label in self.label:
#             path_xdatcar = os.path.join(data.prefix1,
#                                         f"{data.prefix2}.{self.temp}K",
#                                         f"XDATCAR_{label}")
#             if not os.path.isfile(path_xdatcar):
#                 print(f"{path_xdatcar} is not found.")
#                 sys.exit(0)
                
#         # refine data instance
#         self.data = self.refine_data(copy.deepcopy(data))
#         self.temp = np.array([self.temp], dtype=int)
#         self.interval = interval
#         self.poscar_lattice = poscar_lattice
#         self.symbol = symbol
#         self.rmax = rmax
#         self.tol = tol
#         self.tolerance = tolerance
#         self.verbose = verbose
        
#         self.kb = 8.61733326e-5
#         self.cmap = plt.get_cmap("Set1")
#         self.lattice = Lattice(self.poscar_lattice,
#                                self.symbol,
#                                self.rmax,
#                                self.tol,
#                                self.tolerance)
#         self.num_path_except_unknowns = len(self.lattice.path_names)
        
#         # path counts
#         self.counts = []
#         self.encounter_num = []
#         self.encounter_msd = []
#         self.encounter_counts = []
#         self.labels_failed = []
#         self.get_counts()

#         # correlation factor
#         self.f_ind = []
#         self.f_avg = []
#         self.f_cum = []
#         self.f0 = None
#         self.Ea_f = None
#         self.get_correlation_factor()

#         # summary
#         if self.verbose:
#             self.summary()
        
#     def refine_data(self, data):
#         index_temp = np.where(data.temp==self.temp)[0][0]
#         index_label = [data.label[index_temp].index(label) for label in self.label]
#         data.outcar = [data.outcar[index_temp]]
#         data.temp = [self.temp]
#         data.potim = [data.potim[index_temp]]
#         data.nblock = [data.nblock[index_temp]]
#         data.nsw = [data.nsw[index_temp]]
#         data.label = [self.label]
#         data.xdatcar = [[data.xdatcar[index_temp][i] for i in index_label]]
#         if data.force is not None:
#             data.force = [[data.force[index_temp][i] for i in index_label]]
#         return data
    
#     def summary(self):
#         print("# -----------------( Summary )-----------------")
#         print('\nIndividual correlation factors : ')
#         header = ['label', 'f']
#         for i, temp in enumerate(self.temp):
#             print(f"T = {temp} K")
#             data = [
#                 [label, f"{f:.5f}"] for label, f in zip(self.data.label[i], self.f_ind[i])
#             ]
#             data.append(['Average', f"{self.f_avg[i]:.5f}"])
#             print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
#             print('')
#         print(f"Interval : {self.interval} ps")
#         print(f"Cumulative correlation factors : {self.f_cum[0]:.5f}")
#         # check errors
#         if any(labels for labels in self.labels_failed):
#             print("\nLabels where errors occurred : ")
#             header = ['T (K)', 'Label']
#             data = [
#                 [f"{temp}", f"{labels}"] for temp, labels in zip(self.temp, self.labels_failed)
#             ]
#             print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
#         print('')
#         print("# -----------------( Finish  )-----------------")
