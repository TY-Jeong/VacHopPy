import os
import sys
import time
import copy   
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from tqdm import tqdm
#from vachoppy import einstein
# For Arrow3D
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
# For structure visualization
from ase import Atoms
from ase.visualize import view
from ase.io.vasp import read_vasp
# color map for tqdm
from colorama import Fore
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



def extract_force(file_in, 
                  file_out='force.dat'):
    """
    extract force profile from vasprun.xml
    """
    # read vasprun.xml
    with open(file_in, 'r') as f:
        lines = [s.strip() for s in f]
    
    # system info
    nsw, num_atoms = None, None
    for line in lines:
        if 'NSW' in line:
            nsw = int(line.split()[-1].split('<')[0])
        if '<atoms>' in line:
            num_atoms = int(line.split()[1])
        if nsw and num_atoms:
            break
    
    # save forces
    step = 0
    forces = np.zeros((nsw, num_atoms, 3))
    for i, line in enumerate(lines):
        if 'forces' in line:
            force = [list(map(float, s.split()[1:4])) for s in lines[i+1:i+1+num_atoms]]
            force = np.array(force)
            forces[step] = force
            step += 1

    # write out_file
    with open(file_out, 'w') as f:
        for i, force in enumerate(forces, start=1):
            f.write(f"Iteration {i}\n")
            for fx, fy, fz in force:
                fx = str(fx).rjust(12)
                fy = str(fy).rjust(12)
                fz = str(fz).rjust(12)
                f.write(f"{fx} {fy} {fz}\n")



# Helper methods for lattice
def path_TiO2(lattice, acc='high'):
    names = ['OP', 'IP1', 'IP2']
    if acc == 'high':
        d = [2.80311, 2.56299, 2.96677]
        Ea = [1.046, 0.9571, 2.1531]
    else:
        d = [2.85322, 2.24167, 2.96678]
        Ea = [0.870400, 1.044100, 1.965600] # EDIFFG = -0.02
        Ea = [0.983500, 1.045100, 2.203200] # EDIFFG = -0.05
    
z = [8, 1, 2]
    for i in range(len(names)):
        lattice.add_path(names[i], 'vac', 'vac', d[i], Ea[i], 0, z[i])
    for lat_p in lattice.lat_points:
        lat_p['site'] = 'vac'



def path_HfO2(lattice, acc='high'):
    final_A = ['cn4', 'cn3', 'cn3', 'cn3', 'cn4', 'cn4', 'cn4']
    final_B = ['cn3', 'cn4', 'cn4', 'cn4', 'cn3', 'cn3', 'cn3']
    d_A = [2.542, 2.574, 2.785, 2.837, 2.937, 2.965, 2.989]
    d_B = [2.542, 2.576, 2.662, 2.724, 2.937, 2.965, 2.989]
    if acc == 'high':
        Ea_A = [0.74, 0.84, 0.85, 1.35, 1.91, 2.07, 2.01]
        Ea_B = [0.08, 0.32, 0.86, 0.98, 1.25, 1.42, 1.36]
    else:
        Ea_A = [0.7084, 0.8449, 0.8689, 1.366, 1.8806, 2.0373, 2.0027]
        Ea_B = [0.0809, 0.2983, 0.8571, 0.9904, 1.2531, 1.4098, 1.3752]
    z_A = [1, 1, 2, 2, 1, 1, 1]
    z_B = [1, 1, 2, 1, 1, 1, 1]
    
    deltaE = 0.65 if acc=='high' else 0.63

    for i in range(7):
        dE_A = 0.0 if final_A[i]=='cn3' else deltaE
        dE_B = 0.0 if final_A[i]=='cn4' else -deltaE
        lattice.add_path(f"A{i+1}", 'cn3', final_A[i], d_A[i], Ea_A[i], dE_A, z_A[i])
        lattice.add_path(f"B{i+1}", 'cn4', final_B[i], d_B[i], Ea_B[i], dE_B, z_B[i])

    for lat_p in lattice.lat_points:
        x_coord = lat_p['coord'][0]
        if 0.13796 < x_coord < 0.36204 or 0.63796 < x_coord < 0.86204:
            lat_p['site'] = 'cn4'
        else:
            lat_p['site'] = 'cn3'



def path_Al(lattice, acc='high'):
    Ea = 0.622300 if acc=='high' else 0.456500
    lattice.add_path('hop', 'vac', 'vac', 2.85595, Ea, 0, 12)
    for lat_p in lattice.lat_points:
        lat_p['site'] = 'vac'



class Lattice:
    def __init__(self,
                 poscar_perf,
                 symbol='O'):
        """
        class containing information on hopping paths and lattice points
        """
        self.path_poscar = poscar_perf
        self.symbol = symbol

        # information on lattice points
        self.lattice = None
        self.lat_points = []
        self.read_poscar()

        # information on hopping path
        self.path = []
        self.path_names = []
        self.site_names = []


    def read_poscar(self):
        with open(self.path_poscar, 'r') as f:
            lines = np.array([line.strip() for line in f])

        # lattice
        scale = float(lines[1])
        self.lattice = np.array([line.split() for line in lines[2:5]], dtype=float)
        self.lattice *= scale

        # atom species
        atom_species = np.array(lines[5].split())
        num_atoms = np.array(lines[6].split(), dtype=int)
        idx = np.where(atom_species == self.symbol)[0][0]
        coords= lines[num_atoms[:idx].sum()+8:num_atoms[:idx+1].sum()+8]

        # lattice points
        for coord in coords:
            dic_lat = {}
            dic_lat['site'] = 'Nan'
            dic_lat['coord'] = np.array(coord.split()[:3], dtype=float)
            dic_lat['coord_C'] = np.dot(dic_lat['coord'], self.lattice)
            self.lat_points.append(dic_lat)


    def add_path(self,
                 name,
                 site_init,
                 site_final,
                 d,
                 Ea,
                 dE='Nan',
                 z='Nan'):
        
        if name in self.path_names:
            print(f"{name} already exsits.")
            sys.exit(0)
        
        dic_path = {}
        dic_path['name'] = name
        dic_path['site_init'] = site_init
        dic_path['site_final'] = site_final
        dic_path['distance'] = d
        dic_path['Ea'] = Ea
        dic_path['dE'] = dE
        dic_path['z'] = z
        
        self.path.append(dic_path)
        self.path_names.append(name)

        if not site_init in self.site_names:
            self.site_names += [site_init]
        
        if not site_final in self.site_names:
            self.site_names += [site_final]
    

    def print_path(self):
        """
        paths sorted by names
        """
        path_sorted = sorted(self.path,
                             key=lambda x:list(x.values()))

        print("name\tinit\tfinal\td (Å)\tEa (eV)\tdE (eV)")
        for path in path_sorted:
            print(f"{path['name']}", end='\t')
            print(f"{path['site_init']}", end='\t')
            print(f"{path['site_final']}", end='\t')
            print("%.3f"%path['distance'], end='\t')
            print("%.2f"%path['Ea'], end='\t')
            if path['dE'] == 'Nan':
                print(f"{path['dE']}", end='\n')
            else:
                print("%.2f"%path['dE'], end='\n')


    def print_lattice_points(self):
        print("site\tcoord(direct)\tcoord(cartesian)")
        for lat_p in self.lat_points:
            print(lat_p['site'], end='\t')
            print(f"{lat_p['coord']}", end='\t')
            print(f"{lat_p['coord_C']}", end='\n')



class RandomWalk:
    def __init__(self, T, lattice):
        """
        Calculate probability according to site and path
        Args:
            T       : temperature range (list)
            lattice : lattice (Lattice)
        """
        self.T = np.array(T, dtype=float)
        self.lattice = lattice

        # physical constant
        self.kb = 8.61733326e-5

        # energy of site
        self.energy_site = None

        # probability
        self.prob_site = None # (site, T)
        self.prob_path = None # (site, T, path)
        self.U = None         # (site, T)

        # random walk diffusion coeff
        self.D = None

        # Arrhenius fit
        self.Ea = None
        self.D0 = None


    def get_energy(self):
        energy_site = [None] * len(self.lattice.site_names)
        energy_site[0] = 0

        i, repeat = 0, 0
        while True:
            idx_i = self.lattice.site_names.index(self.lattice.path[i]['site_init'])
            idx_f = self.lattice.site_names.index(self.lattice.path[i]['site_final'])

            if (energy_site[idx_f] is None) and (energy_site[idx_i] is not None):
                energy_site[idx_f] = energy_site[idx_i] + self.lattice.path[i]['dE']
            if not(None in energy_site):
                break
            if repeat > 10:
                print('cannot define energy.')
                break

            repeat += 1
            i = (i+1) % len(self.lattice.path)
    
        self.energy_site = np.array(energy_site)


    def get_probability(self):
        if self.energy_site is None:
            self.get_energy()
        
        # probability for site
        name_site = [dic['site'] for dic in self.lattice.lat_points]
        frac_site = np.array([name_site.count(s) for s in self.lattice.site_names], dtype=float)
        frac_site /= np.sum(frac_site)

        prob = np.exp(-self.energy_site / (self.kb*self.T[:, np.newaxis]))
        prob *= frac_site
        prob /= np.sum(prob, axis=1).reshape(-1,1)
        self.prob_site = prob.T

        # probability for path
        z_path = []
        Ea_path = []
        for site in self.lattice.site_names:
            _z_path = [p['z'] for p in self.lattice.path if p['site_init'] == site]
            _Ea_path = [p['Ea'] for p in self.lattice.path if p['site_init'] == site]
            z_path.append(_z_path)
            Ea_path.append(_Ea_path)
            
        self.prob_path = []
        self.U = []
        for z, Ea in zip(z_path, Ea_path):
            z = np.array(z, dtype=float)
            Ea = np.array(Ea, dtype=float)

            prob = np.exp(-Ea / (self.kb*self.T[:, np.newaxis]))
            prob *= z
            U_i = np.sum(prob, axis=1)
            prob /= U_i.reshape(-1,1)
            self.prob_path.append(prob)
            self.U.append(U_i)
        self.U = np.array(self.U)


    def D_rand(self, nu=1e13):
        """
        Diffusion coefficient for Random walk diffusion
        Args:
            nu : jump attemp frequency in Hz (float)
        Return:
            Drand : random walk diffusion coefficient
        """
        if self.prob_site is None:
            self.get_probability()
        
        # hopping distance
        a_path = []
        for site in self.lattice.site_names:
            _a_path = [p['distance'] for p in self.lattice.path if p['site_init']==site]
            a_path.append(np.array(_a_path, dtype=float))

        num_site = len(self.lattice.site_names)
        inner_sum = np.zeros((num_site, len(self.T)))
        for i in range(num_site):
            inner_sum[i] = np.sum(self.prob_path[i] * a_path[i]**2 * 1e-20, axis=1)

        D = self.prob_site * self.U * inner_sum
        D = np.sum(D, axis=0)
        D *= nu/6
        self.D = D
        
        return D
    

    def linear_fitting(self,
                       verbose=False,
                       plot=False):
        """
        Arrhnius fitting for the random walk diffusion coefficient
        Args:
            verbose (bool) : if True, Ea and D0 are displayed.
            plot (bool)    : if True, fitting result will be displayed
        """
        if self.D is None:
            print('D is not defined.')
            sys.exit(0)
        slop, intersect = np.polyfit(1/self.T, np.log(self.D), deg=1)
        self.Ea = -self.kb * slop
        self.D0 = np.exp(intersect)

        if verbose:
            print(f"Ea (rand) = {self.Ea :.3f} eV")
            print(f"D0 (rand) = {self.D0 :.3e} m2/s")

        if plot:
            plt.plot(1/self.T, np.log(self.D), 'k-', label=r'$D_{rand}$')
            plt.plot(1/self.T, slop*(1/self.T)+intersect, 'r:', label=r'Linear fit.')
            plt.ylabel('ln D', fontsize=12)
            plt.xlabel('1/T (1/K)', fontsize=12)
            plt.legend(fontsize=11)
            plt.show()



class LatticeHopping:
    def __init__(self,
                 xdatcar,
                 lattice,
                 force=None,
                 interval=1):
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

        # color map for arrows
        self.cmap = ['b', 'c', 'g', 'deeppink', 'darkorange', 
                     'sienna', 'darkkhaki', 'lawngreen', 'grey', 'wheat', 
                     'navy', 'slateblue', 'purple', 'pink']

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
                      desc=f'{RED}save traj{RESET}'):
            
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
                atom = np.where((self.occ_lat_point[:,i]==idx_pre)==True)[0][0] # 움직인 atom 찾기
            except:
                print(f"idx_pre = {idx_pre}") ## 240912
                print(f'error occured at step {i}.')
                print('please correct the vacancy site yourself.')
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
                  foldername='gif',
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
                         desc=f'{RED}save gif{RESET}'):
            
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
            print(f"Generating {filename}...")
            self.save_gif(filename=filename,
                              files=files,
                              fps=fps,
                              loop=loop)
            print(f"{filename} was generated.")
        

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


    def show_poscar(self,
                    step=None,
                    filename=None,
                    vac=False):
        """
        recieve step or filename
        """
        if step is not None:
            self.save_poscar(step=step, vac=vac)
            filename = f"POSCAR_{step}"
        
        poscar = read_vasp(filename)
        view(poscar)


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
                      desc=f'{RED}save traj_on_lat{RESET}'):
            
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


    def check_multivacancy(self, verbose=True):
        check = np.array([len(i) for i in self.idx_vac.values()])
        check = np.where((check==1) == False)[0]
        
        if len(check) > 0:
            self.multi_vac = True
            if verbose:
                print('multi-vacancy issue occurs:')
                print('  step :', end=' ')
                for i in check:
                    print(i, end = ' ')
                print('')

        else:
            self.multi_vac = False
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
                print("there is no connected site.")       
                print(f"find the vacancy site for your self. (step: {step})")
                break
            
            else:
                print("there are multiple candidates.")       
                print(f"find the vacancy site for your self. (step: {step})")
                break

        # update trace arrows
        self.get_trace_arrows()
            


class Analyzer:
    def __init__(self,
                 traj,
                 lattice):
        """
        module to analyze trajectory.
        this module search diffusion paths in MD trajectory.
        """
        self.traj = traj
        self.traj_backup = traj
        self.lattice = lattice
        self.tolerance = 0.001
        
        # check whether hopping paths are defined
        if len(lattice.path) == 0:
            print("hopping paths are not defined.")
            sys.exit(0)
            
        self.path = lattice.path
        self.path_names = lattice.path_names
        self.site_names = lattice.site_names

        # check whether lattice points are defined
        for lat_p in lattice.lat_points:
            if lat_p['site'] == 'Nan':
                print("some lattice points are not defined")
                print(f"coord: {lat_p['coord']}")
                sys.exit(0)
        self.lat_points = lattice.lat_points
        
        self.path_unknown = {}
        self.path_unknown['name'] = 'unknown'
        self.path_vac = None # path of vacancy

        # self.idx_vac = self.traj.idx_vac # (dic) index of vacancy
        self.idx_vac = copy.deepcopy(self.traj.idx_vac)
            

    def get_path_vacancy(self, 
                         verbose=True):
        step_unknown = []
        self.path_vac = []
        idx = 0 # index in self.path_vac
        
        for step in range(self.traj.num_step-1):
            coord_init = self.lat_points[self.traj.idx_vac[step][0]]['coord']
            coord_final = self.lat_points[self.traj.idx_vac[step+1][0]]['coord']
            
            # check whether vacancy moves
            distance = self.traj.distance_PBC(coord_init, coord_final)
            
            if distance > self.tolerance:
                site_init = self.lat_points[self.traj.idx_vac[step][0]]['site']
                path = self.determine_path(site_init, distance)
                
                path['step'] = step+1
                
                if path['name'] == self.path_unknown['name']:
                    step_unknown += [step+1]
                    
                self.path_vac += [copy.deepcopy(path)]
                idx += 1
                
        if len(step_unknown) > 0 and verbose:
            print(f"unknown steps are detected.: {step_unknown}")
        

    def determine_path(self,
                       site_init,
                       distance):
        """
        path is determined based on initial site and distance.
        """
        if not site_init in self.site_names:
            print(f"{site_init} is unknown site.")
            sys.exit(0)
        
        candidate = []
        for p in self.path:
            err = abs(distance - p['distance'])
            if err < self.tolerance and p['site_init']==site_init:
                candidate += [p]
        
        if len(candidate) == 0:
            p = self.path_unknown
            p['site_init'] = site_init
            p['distance'] = distance
            return p
            
        elif len(candidate) > 1:
            print('there are many candidates.')
            print(f"initial site = {site_init}, distance = {distance}")
            print(f'please use smaller tolerance.'\
                  f'now: tolerance={self.tolerance}')
            sys.exit(0)
        
        else:
            return candidate[0]
            

    def print_path_vacancy(self):
        if self.path_vac is None:
            print("path_vac is not defines.")
            print("please run 'get_path_vacancy'.")

        else:
            print("path of vacancy :")
            for p_vac in self.path_vac:
                print(p_vac['name'], end=' ')
            print('')
    

    def plot_path_counts(self, 
                         figure='counts.png',
                         text='counts.txt',
                         disp=True,
                         save_figure=True,
                         save_text=True,
                         bar_width=0.6,
                         bar_color='c',
                         dpi=300,
                         sort=True):
        
        path_vac_names = [p_vac['name'] for p_vac in self.path_vac]
        path_type = copy.deepcopy(self.path_names)
        
        check_unknown = False
        if self.path_unknown['name'] in path_vac_names:
            if not 'U' in path_type:
                path_type.append('U')
            check_unknown = True
            num_unknown = path_vac_names.count(self.path_unknown['name'])
        
        path_count = []
        for p_type in path_type:
            path_count.append(path_vac_names.count(p_type))

        if check_unknown:
            path_count[-1] = num_unknown
        path_count = np.array(path_count)

        # sorting paths using Ea
        if sort:
            path_Ea = []
            for p_type in path_type:
                for p in self.path:
                    if p['name'] == p_type:
                        path_Ea.append(p['Ea'])
            if check_unknown:
                path_Ea.append(100)

            path_Ea = np.array(path_Ea)
            args = np.argsort(path_Ea)

            path_type_sorted = []
            for arg in args:
                path_type_sorted.append(path_type[arg])

            path_Ea = path_Ea[args]
            path_count = path_count[args]
            path_type = path_type_sorted

        # plot bar graph
        x = np.arange(len(path_count))
        plt.bar(x, path_count, color=bar_color, width=bar_width)
        plt.xticks(x, path_type)
        
        plt.xlabel('Path', fontsize=13)
        plt.ylabel('Counts', fontsize=13)
        
        if save_figure:
            plt.savefig(figure, dpi=dpi)

        if disp:
            plt.show()
        plt.close()
        
        # write counts.txt
        if save_text:
            with open(text, 'w') as f:
                f.write(f"total counts = {np.sum(path_count)}\n\n")
                f.write("path\tcounts\n")
                for name, count in zip(path_type, path_count):
                    f.write(f"{name}\t{count}\n")
                if check_unknown:
                    f.write(f"unknown\t{num_unknown}")
        

    def path_tracer(self,
                    paths,
                    p_init,
                    p_goal):
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
    

    def path_decomposer(self,
                        index):
        """
        index : index in self.path_vac
        """
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
        path_unwrap = []
        
        for idx, path in enumerate(self.path_vac):
            step = path['step']
            # known path
            if path['name'] != self.path_unknown['name']:
                path_unwrap += [path]
            
            # unknown path
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
                    # add unknown path
                    p_new = {}
                    p_new['name'] = self.path_unknown['name']
                    p_new['step'] = step
                    path_unwrap += [copy.deepcopy(p_new)]
                    continue
                
                for i in range(len(p)-1):
                    coord_init = self.lat_points[p[i]]['coord']
                    coord_final = self.lat_points[p[i+1]]['coord']
                    
                    site_init = self.lat_points[p[i]]['site']
                    distance = self.traj.distance_PBC(coord_init, coord_final)
                    
                    p_new = self.determine_path(site_init, distance)
                    p_new['step'] = step
                    
                    path_unwrap += [copy.deepcopy(p_new)]
        
        self.path_vac = path_unwrap
        
        # check unknown path
        check_unknown = []
        for p_vac in self.path_vac:
            if p_vac['name'] == self.path_unknown['name']:
                check_unknown += [p_vac]
        
        if len(check_unknown) == 0:
            print("no unknown path exist.")
        
        else:
            print("unknown path exist.", end=' ')
            print("( step:", end=' ')
            for p in check_unknown:
                print(p['step'], end=' ')
            print(")")
        

    def print_summary(self,
                      figure='counts.png',
                      text='counts.txt',
                      disp=True,
                      save_figure=False,
                      save_text=False,
                      bar_width=0.6,
                      bar_color='c',
                      dpi=300,
                      sort=True):
        
        counts_tot = len(self.path_vac)
        print(f"total counts : {counts_tot}")
        print(f"hopping sequence :" )
        
        Ea_max = 0
        for p_vac in self.path_vac:
            print(p_vac['name'], end=' ')
            if p_vac['name'] != self.path_unknown['name'] and p_vac['Ea'] > Ea_max:
                Ea_max = p_vac['Ea']
        print('')
        print(f"maximum Ea : {Ea_max} eV")
        
        if disp or save_figure or save_text:
            self.plot_path_counts(figure=figure,
                                text=text,
                                disp=disp,
                                save_figure=save_figure,
                                save_text=save_text,
                                bar_width=bar_width,
                                bar_color=bar_color,
                                dpi=dpi,
                                sort=sort)



class Encounter:
    def __init__(self,
                 analyzer,
                 verbose=True):
        """
        Obtain information on the encounters.
        
        Args:
            analyzer : instance of Analyzer class
            verbose  : (default: True)
            
        """
        
        self.analyzer = analyzer
        self.traj = analyzer.traj
        self.path = copy.deepcopy(self.analyzer.path)
        self.path_names = self.analyzer.path_names
        
        # check multi-vacancy
        if self.traj.multi_vac:
            print('multi-vacancy issue occured.')
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
        
        # unknown path
        self.path_unknown = self.analyzer.path_unknown['name']
        self.unknown = False
        self.check_unknown()
        
        # correlation factor
        self.msd = 0
        self.path_counts = np.zeros(len(self.path_names))
        self.path_dist = np.zeros_like(self.path_counts)
        self.get_msd()
        self.get_counts()
        
        # print results
        if verbose:
            self.print_summary()
            
            
    def check_unknown(self):
        self.unknown = True if self.path_unknown in self.path_enc_all else False
        
        if self.unknown:
            dic = {}
            dic['name'] = self.path_unknown
            dic['distance'] = 0
            self.path.append(dic)
        
        
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
        print(f"Number of encounters      : {self.num_enc}")
        print(f"Mean squared displacement : {self.msd:.3f}") # Å2
        count_tot = int(np.sum(self.path_counts*self.num_enc))
        print(f"Total hopping counts      : {count_tot}")
        count_mean = np.sum(self.path_counts)
        print(f"Mean hopping counts       : {count_mean:.3f}")
        print('')
        print("           ", end="")
        for name in self.path_names:
            print("{:<10}".format(name), end=" ")
        print('')
        print("Distance   ", end="")
        for dist in self.path_dist:
            print("{:<10.3f}".format(dist), end=" ")
        print('')
        print("<Counts>   ", end="")
        for count in self.path_counts:
            print("{:<10.3f}".format(count), end=" ")
        print('\n') 



class CorrelationFactor(Encounter):
    def __init__(self,
                 analyzer,
                 temp,
                 verbose=True):
        """
        Calculate correlation factor based on encounters.
        Args:
            analyzer: 
            temp    :
            verbose :
        
        """
        super().__init__(analyzer, verbose=False)
        
        self.temp = temp
        self.lattice = analyzer.lattice
        
        # hopping distance
        self.a_path = []
        self.get_a_path()
        
        # random walk probability
        rand = RandomWalk([self.temp], self.lattice)
        rand.get_probability()
        self.prob_site = rand.prob_site
        self.prob_path = rand.prob_path
        
        # correlatoin factor
        self.f_cor = None
        self.get_f_cor()
        
        # print results
        if verbose:
            self.print_summary()
        
        
    def get_a_path(self):
        for site in self.lattice.site_names:
            _a_path = [p['distance'] for p in self.lattice.path if p['site_init']==site]
            self.a_path.append(np.array(_a_path, dtype=float))
            
            
    def get_f_cor(self):
        inner_sum = np.array([p * a**2 for p, a in zip(self.prob_path, self.a_path)])
        inner_sum = np.sum(inner_sum, axis=2)
        outer_sum = np.sum(self.prob_site * inner_sum, axis=0)
        outer_sum *= np.sum(self.path_counts)
        self.f_cor = self.msd / outer_sum[0]
        
        
    def print_summary(self):
        print(f"Correlation factor        : {self.f_cor:.3f}")
        print(f"Number of encounters      : {self.num_enc}")
        count_tot = int(np.sum(self.path_counts*self.num_enc))
        print(f"Total hopping counts      : {count_tot}")
        print(f"Mean squared displacement : {self.msd:.3f}")
        count_mean = np.sum(self.path_counts)
        print(f"Mean hopping counts       : {count_mean:.3f}")
        print('')
        print("           ", end="")
        for name in self.path_names:
            print("{:<10}".format(name), end=" ")
        print('')
        print("Distance   ", end="")
        for dist in self.path_dist:
            print("{:<10.3f}".format(dist), end=" ")
        print('')
        print("<Counts>   ", end="")
        for count in self.path_counts:
            print("{:<10.3f}".format(count), end=" ")
        print('\n')
        
        
        
class CumulativeCorrelationFactor:
    def __init__(self,
                 xdatcar,
                 lattice,
                 temp,
                 label='auto',
                 force=None,
                 interval=1,
                 verbose=False,
                 multi_vac=True,
                 multi_path=True,
                 correction_TS=True,
                 update_vacancy={},
                 step_ignore={}):
        """
        Cumulative correlation function
        
        Args:
            xdatcar : directory where XDATCAR_{label} exists
            label   : label in XDATCAR_{label} (list or 'auto)
            temp    : temperature in K to calculate Boltzmann distribution
            force   : directory where force_{label}.dat exists
            
            update_vacancy : (dic) used to fix error in multi-vacancy correction
            : key = label ; item = list ex. [[2978, 11]]
            step_ignore : (dic) steps to ignore for the TS correction
            : key = label; item = list of steps
        """
        
        self.xdatcar = xdatcar
        self.force_dir = force
        self.lattice = lattice
        self.temp = temp
        self.label = self.get_label() if label=='auto' else label
        self.interval = interval
        self.verbose = verbose
        self.symbol = lattice.symbol
        self.path = self.lattice.path

        # corrections
        self.multi_vac = multi_vac
        self.multi_path = multi_path
        self.correction_TS = correction_TS
        self.update_vacancy = update_vacancy
        self.step_ignore = step_ignore
        
        # f_cor
        self.times = []
        self.path_name = []
        self.path_dist = []
        self.f_ensemble = []
        self.path_seq_cum = []
        self.msd_cum = 0
        self.num_enc_cum = 0
        self.label_err = []
        self.prob_site = None
        self.prob_path = None
        self.a_path = None
        self.label_success = []
        self.get_correlation_factors()

        # unknown path
        self.path_name.append('unknown')
        self.path_dist.append(0)
        self.path_dist = np.array(self.path_dist)

        # cumulative f_cor
        self.f_avg = np.average(self.f_ensemble) # averaged over ensmebles
        self.f_cum = None
        self.get_cumulative_correlation_factor()

        # print results
        self.print_summary()


    def get_label(self):
        label = []
        for filename in os.listdir(self.xdatcar):
            if len(filename.split('_')) == 2:
                first, second = filename.split('_')
                if first == 'XDATCAR':
                    label.append(second)
        label.sort()
        return label
    

    def get_correlation_factors(self):
        check_first = True
        for label in tqdm(self.label, 
                          bar_format='{l_bar}{bar:20}{r_bar}{bar:-10b}', 
                          ascii=True, 
                          desc=f'{RED}f_cor{RESET}'):
            start = time.time()

            print(f"# Label : {label}")
            analyzer = self.get_analyzer(label)
            if analyzer is None:
                end = time.time()
                continue
            correlation = CorrelationFactor(analyzer, self.temp, self.verbose)
            
            if check_first:
                self.prob_site = correlation.prob_site
                self.prob_path = correlation.prob_path
                self.a_path = correlation.a_path
                check_first = False
            
            self.f_ensemble.append(correlation.f_cor)
            self.msd_cum += correlation.msd * correlation.num_enc
            self.num_enc_cum += correlation.num_enc
            self.path_seq_cum += correlation.path_enc_all

            end = time.time()
            self.times.append(end-start)
            self.label_success.append(label)

        self.f_ensemble = np.array(self.f_ensemble, dtype=float)
        self.msd_cum /= self.num_enc_cum
        
        for path in self.path:
            self.path_name.append(path['name'])
            self.path_dist.append(path['distance'])
    

    def get_analyzer(self, 
                     label):
        path_xdatcar = os.path.join(self.xdatcar, f"XDATCAR_{label}")
        if self.force_dir is not None:
            path_force = os.path.join(self.force_dir, f"force_{label}.dat")
        else:
            path_force = None
            
        traj = LatticeHopping(lattice=self.lattice,
                              xdatcar=path_xdatcar,
                              force=path_force,
                              interval=self.interval)
        
        if self.multi_vac:
            traj.correct_multivacancy(start=1)
            if label in self.update_vacancy.keys():
                for [step, site] in self.update_vacancy[label]:
                    traj.update_vacancy(step, site)
                    traj.correct_multivacancy(start=step)
            traj.check_multivacancy(verbose=False)
        
        if traj.multi_vac is True:
            print(f'multi-vacancy issue in label {label}.')
            print('the calculation is skipped.')
            self.label_err.append(label)
            return None

        if self.correction_TS and (self.force_dir is not None):
            step_ignore = self.step_ignore[label] if label in self.step_ignore.keys() else []
            traj.correct_transition_state(step_ignore=step_ignore)

        analyzer = Analyzer(traj, self.lattice)
        analyzer.get_path_vacancy(verbose=self.verbose)

        if self.multi_path:
            analyzer.correct_multipath()

        if self.verbose:
            analyzer.print_summary(disp=False, 
                                   save_figure=False, 
                                   save_text=False)
        print('')

        return analyzer
    

    def get_cumulative_correlation_factor(self):
        # cumulative mean count of each path
        self.counts_cum = np.zeros(len(self.path_name))
        for i, name in enumerate(self.path_name):
            self.counts_cum[i] = self.path_seq_cum.count(name)
        self.counts_cum /= self.num_enc_cum
        
        # cumulative correlation function
        inner_sum = np.array([p * a**2 for p, a in zip(self.prob_path, self.a_path)])
        inner_sum = np.sum(inner_sum, axis=2)
        outer_sum = np.sum(self.prob_site * inner_sum, axis=0)
        outer_sum *= np.sum(self.counts_cum)
        self.f_cum = self.msd_cum / outer_sum[0]


    def print_summary(self):
        print("## Summary")
        print("#  Total counts")
        print("           ", end="")
        for name in self.path_name:
            print("{:<10}".format(name), end=" ")
        print('')
        print("Distance   ", end="")
        for dist in self.path_dist:
            print("{:<10.3f}".format(dist), end=" ")
        print('')
        print("Counts     ", end="")
        for count in self.counts_cum * self.num_enc_cum:
            count = int(count)
            print("{:<10}".format(count), end=" ")
        print('\n')
        print(f"      Mean correlation factor : {self.f_avg:.3f}")
        print(f"Cumulative correlation factor : {self.f_cum:.3f}")
        print('')

        if self.verbose:
            print("{:<10} {:<10}".format('Label', 'f_cor'))
            for label, f in zip(self.label_success, self.f_ensemble):
                print("{:<10} {:<10.3f}".format(label, f))
            print('')
            print("{:<10} {:<10}".format('Label', 'Time(s)'))
            for label, time in zip(self.label_success, self.times):
                print("{:<10} {:<10.3f}".format(label, time))
            time_tot = np.sum(np.array(self.times))
            print('')
            print(f'Total time : {time_tot:.3f} s')

        if len(self.label_err) > 0:
            print('Error occured : ', end='')
            for label in self.label_err:
                print(label, end=' ')
            print('')




# class CorrelationFactor_legacy:
#     def __init__(self,
#                  analyzer,
#                  verbose=True):
#         """
#         Calculate correlation factor using encouters.
#         Encounters are defined with consideration of PBC.
        
#         Args:
#             analyzer : instance of Analyzer class
#             verbose  : (default: True)
            
#         Return:
#             f_cor    : correlation factor
#         """
        
#         self.analyzer = analyzer
#         self.traj = analyzer.traj
#         self.path = copy.deepcopy(self.analyzer.path)
#         self.path_names = self.analyzer.path_names
        
#         # check multi-vacancy
#         if self.traj.multi_vac:
#             print('multi-vacancy issue occured.')
#             sys.exit(0)
        
#         # trajectory of vacancy with consideration of PBC
#         self.traj_vac = {}
#         self.get_traj_vac()
        
#         # encounters
#         self.coord_i_enc = []
#         self.coord_f_enc = []
#         self.path_enc = []
#         self.path_enc_all = []
#         self.tolerance = 0.01
#         self.get_encounters()
#         self.num_enc = len(self.path_enc)
        
#         # unknown path
#         self.path_unknown = self.analyzer.path_unknown['name']
#         self.unknown = False
#         self.check_unknown()
        
#         # correlation factor
#         self.msd = 0
#         self.path_counts = np.zeros(len(self.path_names))
#         self.path_dist = np.zeros_like(self.path_counts)
#         self.get_msd()
#         self.get_counts()
        
#         self.f_cor = None
#         self.get_correlation_factor()
        
#         # print results
#         if verbose:
#             self.print_summary()
            
            
#     def check_unknown(self):
#         self.unknown = True if self.path_unknown in self.path_enc_all else False
        
#         if self.unknown:
#             dic = {}
#             dic['name'] = self.path_unknown
#             dic['distance'] = 0
#             self.path.append(dic)
        
        
#     def get_traj_vac(self):
#         idx_vac = np.array([site[0] for site in self.traj.idx_vac.values()])
#         step_move = np.diff(idx_vac)
        
#         # steps where vacancy moved
#         step_move = np.where(step_move != 0)[0]
#         step_move += 1
        
#         # path of vacancy
#         path_net = idx_vac[step_move]
        
#         # coords considering PBC
#         coords = self.traj.lat_points[path_net]
#         displacement = np.zeros_like(coords)
#         displacement[1:] = np.diff(coords, axis=0)
#         displacement[displacement > 0.5] -= 1.0
#         displacement[displacement < -0.5] += 1.0
#         displacement = np.cumsum(displacement, axis=0)
#         coords = coords[0] + displacement
        
#         # save net path
#         for step, coord in zip(step_move, coords):
#             dic = {}
#             dic['index'] = idx_vac[step]
#             dic['coord'] = coord
#             self.traj_vac[step] = dic
        
            
#     def update_encounters(self, step):
#         """
#         Update encounter coordinates and paths based on the given step.

#         Args:
#             step (int) : Current simulation step.
#         """
#         # Trace arrows at the given step
#         arrows = np.array([dic['lat_points'] for dic in self.traj.trace_arrows[step-1]], dtype=int)
        
#         # Path of the vacancy
#         path = self.analyzer.path_tracer(arrows, self.traj.idx_vac[step][0], self.traj.idx_vac[step-1][0])

#         # Get the current vacancy coordinates
#         coord_vac = self.traj_vac[step]['coord']

#         # Check if there are any initial encounters
#         if len(self.coord_i_enc) == 0:
#             updated_coord_i_enc = []
#             updated_coord_f_enc = []
#             updated_path_enc = []
            
#             # Loop through the path and update coordinates and paths
#             for i in range(len(path) - 1):
#                 idx_i, idx_f = path[i], path[i + 1]

#                 coord_i = self.traj.lat_points[idx_i]
#                 coord_f = self.traj.lat_points[idx_f]

#                 displacement = coord_f - coord_i
#                 displacement[displacement > 0.5] -= 1.0
#                 displacement[displacement < -0.5] += 1.0

#                 coord_new = coord_vac + displacement
#                 site = self.analyzer.lat_points[idx_f]['site']
#                 distance = self.traj.distance_PBC(coord_i, coord_f)
#                 path_name = self.analyzer.determine_path(site, distance)['name']

#                 updated_coord_i_enc.append(coord_vac)
#                 updated_coord_f_enc.append(coord_new)
#                 updated_path_enc.append([path_name])

#                 coord_vac = coord_new

#             updated_coord_i_enc = np.array(updated_coord_i_enc)
#             updated_coord_f_enc = np.array(updated_coord_f_enc)

#         else:
#             updated_coord_i_enc = copy.deepcopy(self.coord_i_enc)
#             updated_coord_f_enc = copy.deepcopy(self.coord_f_enc)
#             updated_path_enc = copy.deepcopy(self.path_enc)

#             # Loop through the path and update encounters
#             for i in range(len(path) - 1):
#                 idx_i, idx_f = path[i], path[i + 1]

#                 coord_i = self.traj.lat_points[idx_i]
#                 coord_f = self.traj.lat_points[idx_f]

#                 displacement = coord_f - coord_i
#                 displacement[displacement > 0.5] -= 1.0
#                 displacement[displacement < -0.5] += 1.0

#                 coord_new = coord_vac + displacement
#                 site = self.analyzer.lat_points[idx_f]['site']
#                 distance = self.traj.distance_PBC(coord_i, coord_f)
#                 path_name = self.analyzer.determine_path(site, distance)['name']

#                 # Check if current vacancy coordinate is in the final encounter coordinates
                
#                 check = np.linalg.norm(self.coord_f_enc - coord_vac, axis=1) < self.tolerance

#                 if np.any(check):
#                     idx = np.where(check)[0][0]
#                     updated_coord_f_enc[idx] = coord_new
#                     updated_path_enc[idx].append(path_name)
#                 else:
#                     updated_coord_i_enc = np.vstack([updated_coord_i_enc, coord_vac])
#                     updated_coord_f_enc = np.vstack([updated_coord_f_enc, coord_new])
#                     updated_path_enc.append([path_name])

#                 coord_vac = coord_new
                
#         self.coord_i_enc = updated_coord_i_enc
#         self.coord_f_enc = updated_coord_f_enc
#         self.path_enc = updated_path_enc
    
    
#     def get_encounters(self):
#         for step in self.traj_vac.keys():
#             self.update_encounters(step)
            
#         for path in self.path_enc:
#             self.path_enc_all += path
                
                
#     def get_msd(self):
#         displacement = self.coord_f_enc - self.coord_i_enc
#         displacement = np.dot(displacement, self.traj.lattice)
#         self.msd = np.average(np.sum(displacement**2, axis=1))
        
    
#     def get_counts(self):
#         for i, name in enumerate(self.path_names):
#             self.path_counts[i] = self.path_enc_all.count(name)
#             self.path_dist[i] = self.analyzer.path[i]['distance']
#         self.path_counts /= self.num_enc
        
        
#     def get_correlation_factor(self):
#         if np.sum(self.path_counts) > 0:
#             self.f_cor = self.msd / np.sum(self.path_counts * self.path_dist**2)
#         else:
#             print('no hopping detected : f_cor is set to 0')
#             self.f_cor = 0
            
            
#     def print_summary(self):
#         print(f"Correlation factor        : {self.f_cor:.3f}")
#         print(f"Number of encounters      : {self.num_enc}")
#         count_tot = int(np.sum(self.path_counts*self.num_enc))
#         print(f"Total hopping counts      : {count_tot}")
#         print(f"Mean squared displacement : {self.msd:.3f}")
#         count_mean = np.sum(self.path_counts)
#         print(f"Mean hopping counts       : {count_mean:.3f}")
#         print('')
#         print("           ", end="")
#         for name in self.path_names:
#             print("{:<10}".format(name), end=" ")
#         print('')
#         print("Distance   ", end="")
#         for dist in self.path_dist:
#             print("{:<10.3f}".format(dist), end=" ")
#         print('')
#         print("<Counts>   ", end="")
#         for count in self.path_counts:
#             print("{:<10.3f}".format(count), end=" ")
#         print('\n')   

   

# class CumulativeCorrelationFactor_legacy:
#     def __init__(self,
#                  xdatcar,
#                  lattice,
#                  label='auto',
#                  force=None,
#                  interval=1,
#                  verbose=False,
#                  multi_vac=True,
#                  multi_path=True,
#                  correction_TS=True):
#         """
#         xdatcar : directory where XDATCAR_{label} exists
#         label   : label in XDATCAR_{label} (list or 'auto)
#         force   : directory where force_{label}.dat exists
#         """
        
#         self.xdatcar = xdatcar
#         self.force_dir = force
#         self.lattice = lattice
#         self.label = self.get_label() if label=='auto' else label
#         self.interval = interval
#         self.verbose = verbose
#         self.symbol = lattice.symbol
#         self.path = self.lattice.path

#         # corrections
#         self.multi_vac = multi_vac
#         self.multi_path = multi_path
#         self.correction_TS = correction_TS
        
#         # f_cor
#         self.times = []
#         self.path_name = []
#         self.path_dist = []
#         self.f_ensemble = []
#         self.path_seq_cum = []
#         self.msd_cum = 0
#         self.num_enc_cum = 0
#         self.label_err = []
#         self.get_correlation_factors()

#         # unknown path
#         self.path_name.append('unknown')
#         self.path_dist.append(0)
#         self.path_dist = np.array(self.path_dist)

#         # cumulative f_cor
#         self.f_avg = np.average(self.f_ensemble) # averaged over ensmebles
#         self.f_cum = None
#         self.get_cumulative_correlation_factor()

#         # print results
#         self.print_summary()


#     def get_label(self):
#         label = []
#         for filename in os.listdir(self.xdatcar):
#             if len(filename.split('_')) == 2:
#                 first, second = filename.split('_')
#                 if first == 'XDATCAR':
#                     label.append(second)
#         label.sort()
#         return label
    

#     def get_correlation_factors(self):
#         for label in tqdm(self.label, 
#                           bar_format='{l_bar}{bar:20}{r_bar}{bar:-10b}', 
#                           ascii=True, 
#                           desc=f'{RED}f_cor{RESET}'):
#             start = time.time()

#             print(f"# Label : {label}")
#             analyzer = self.get_analyzer(label)
#             if analyzer is None:
#                 end = time.time()
#                 continue
#             correlation = CorrelationFactor(analyzer, self.verbose)

#             self.f_ensemble.append(correlation.f_cor)
#             self.msd_cum += correlation.msd * correlation.num_enc
#             self.num_enc_cum += correlation.num_enc
#             self.path_seq_cum += correlation.path_enc_all

#             end = time.time()
#             self.times.append(end-start)

#         self.f_ensemble = np.array(self.f_ensemble, dtype=float)
#         self.msd_cum /= self.num_enc_cum
        
#         for path in self.path:
#             self.path_name.append(path['name'])
#             self.path_dist.append(path['distance'])
    

#     def get_analyzer(self, 
#                      label):
#         path_xdatcar = os.path.join(self.xdatcar, f"XDATCAR_{label}")
#         if self.force_dir is not None:
#             path_force = os.path.join(self.force_dir, f"force_{label}.dat")
#         else:
#             path_force = None
            
#         traj = LatticeHopping(lattice=self.lattice,
#                               xdatcar=path_xdatcar,
#                               force=path_force,
#                               interval=self.interval)
        
#         if self.multi_vac:
#             traj.correct_multivacancy()
#             traj.check_multivacancy(verbose=False)
        
#         if traj.multi_vac is True:
#             print(f'multi-vacancy issue in label {label}.')
#             print('the calculation is skipped.')
#             self.label_err.append(label)
#             return None

#         if self.correction_TS and (self.force_dir is not None):
#             traj.correct_transition_state()

#         analyzer = Analyzer(traj, self.lattice)
#         analyzer.get_path_vacancy(verbose=self.verbose)

#         if self.multi_path:
#             analyzer.correct_multipath()

#         if self.verbose:
#             analyzer.print_summary(disp=False, 
#                                    save_figure=False, 
#                                    save_text=False)
#         print('')

#         return analyzer
    

#     def get_cumulative_correlation_factor(self):
#         self.counts_cum = np.zeros(len(self.path_name))

#         for i, name in enumerate(self.path_name):
#             self.counts_cum[i] = self.path_seq_cum.count(name)

#         self.counts_cum /= self.num_enc_cum
#         self.f_cum = self.msd_cum / np.sum(self.counts_cum * self.path_dist**2)


#     def print_summary(self):
#         print("## Summary")
#         print("#  Total counts")
#         print("           ", end="")
#         for name in self.path_name:
#             print("{:<10}".format(name), end=" ")
#         print('')
#         print("Distance   ", end="")
#         for dist in self.path_dist:
#             print("{:<10.3f}".format(dist), end=" ")
#         print('')
#         print("Counts     ", end="")
#         for count in self.counts_cum * self.num_enc_cum:
#             count = int(count)
#             print("{:<10}".format(count), end=" ")
#         print('\n')
#         print(f"      Mean correlation factor : {self.f_avg:.3f}")
#         print(f"Cumulative correlation factor : {self.f_cum:.3f}")
#         print('')

#         if self.verbose:
#             print("{:<10} {:<10}".format('Label', 'f_cor'))
#             for label, f in zip(self.label, self.f_ensemble):
#                 print("{:<10} {:<10.3f}".format(label, f))
#             print('')
#             print("{:<10} {:<10}".format('Label', 'Time(s)'))
#             for label, time in zip(self.label, self.times):
#                 print("{:<10} {:<10.3f}".format(label, time))
#             time_tot = np.sum(np.array(self.times))
#             print('')
#             print(f'Total time : {time_tot:.3f} s')

#         if len(self.label_err) > 0:
#             print('Error occured : ', end='')
#             for label in self.label_err:
#                 print(label, end=' ')
#             print('')



# class CorrelationFactor_legacy:
#     def __init__(self,
#                  analyzer,
#                  verbose=True):
        
#         self.analyzer = analyzer
#         self.traj = analyzer.traj
#         self.path = self.analyzer.path
#         self.path_names = self.analyzer.path_names

#         self.encounters = []
#         self.get_encounters()
#         self.num_enc = len(self.encounters)

#         self.path_unknown = self.analyzer.path_unknown['name']
#         self.unknown = False
#         self.check_unknown()

#         # mean squared displacement
#         self.msd = 0
#         self.path_sequence = []
#         self.get_msd()

#         # mean counts
#         self.path_counts = np.zeros(len(self.path_names))
#         self.path_dist = np.zeros_like(self.path_counts)
#         self.get_counts()
        
#         # f_cor
#         self.f_cor = None
#         self.get_correlation_factor()

#         # print results
#         if verbose:
#             self.print_summary()


#     def get_encounters(self):
#         path_diff = np.zeros_like(self.traj.occ_lat_point)
#         path_diff[:,1:] = np.diff(self.traj.occ_lat_point, axis=1)

#         atom_move, step_move = np.where(path_diff != 0)
#         self.encounters = []
        
#         for atom in set(atom_move):
#             dic = {}
#             idx = np.where(atom_move==atom)[0]
#             step = np.concatenate((np.array([0]), step_move[idx]))
#             path = self.traj.occ_lat_point[atom, step]

#             dic['idx'] = atom
#             dic['path'] = path
#             dic['path_names'] = self.get_path_name(path)
#             displacement = self.displacement_PBC(path)
#             dic['squared_disp'] = np.sum(displacement**2)

#             self.encounters.append(dic)
    

#     def displacement_PBC(self, 
#                          path):
#         coords = self.traj.lat_points[path]
        
#         displacement = np.zeros_like(coords)
#         displacement[1:,:] = np.diff(coords, axis=0)

#         displacement[displacement > 0.5] -= 1.0
#         displacement[displacement < -0.5] += 1.0
#         displacement = np.sum(displacement, axis=0)
        
#         displacement_C = np.dot(displacement, self.traj.lattice)

#         return displacement_C  
    

#     def get_path_name(self, 
#                       path):
#         path_names = []
#         for i in range(1, len(path)):
#             idx_i, idx_f = path[i-1], path[i]
#             site_init = self.analyzer.lat_points[idx_i]['site']

#             coord_i = self.traj.lat_points[idx_i]
#             coord_f = self.traj.lat_points[idx_f]
#             distance = self.traj.distance_PBC(coord_i, coord_f)
            
#             name = self.analyzer.determine_path(site_init, distance)['name']
#             path_names.append(name)

#         return path_names
    

#     def check_unknown(self):
#         for dic in self.encounters:
#             if self.path_unknown in dic['path_names']:
#                 self.unknown = True
#                 self.path_names.append(self.path_unknown)
#                 break

#         if self.unknown:
#             dic = {}
#             dic['name'] = self.path_unknown
#             dic['distance'] = 0
#             self.path.append(dic)


#     def get_msd(self):
#         for dic in self.encounters:
#             self.msd += dic['squared_disp']
#             self.path_sequence += dic['path_names']
#         self.msd /= self.num_enc


#     def get_counts(self):
#         for i, name in enumerate(self.path_names):
#             self.path_counts[i] = self.path_sequence.count(name)
#             self.path_dist[i] = self.analyzer.path[i]['distance']
#         self.path_counts /= self.num_enc


#     def get_correlation_factor(self):
#         if np.sum(self.path_counts) > 0:
#             self.f_cor = self.msd / np.sum(self.path_counts * self.path_dist**2)
#         else:
#             print('no hopping detected : f_cor is set to 0')
#             self.f_cor = 0


#     def print_summary(self):
#         print(f"Correlation factor        : {self.f_cor:.3f}")
#         print(f"Number of encounters      : {self.num_enc}")
#         count_tot = int(np.sum(self.path_counts*self.num_enc))
#         print(f"Total hopping counts      : {count_tot}")
#         print(f"Mean squared displacement : {self.msd:.3f}")
#         count_mean = np.sum(self.path_counts)
#         print(f"Mean hopping counts       : {count_mean:.3f}")
#         print('')
#         print("           ", end="")
#         for name in self.path_names:
#             print("{:<10}".format(name), end=" ")
#         print('')
#         print("Distance   ", end="")
#         for dist in self.path_dist:
#             print("{:<10.3f}".format(dist), end=" ")
#         print('')
#         print("<Counts>   ", end="")
#         for count in self.path_counts:
#             print("{:<10.3f}".format(count), end=" ")
#         print('\n')
             
    
    
# class CumulativeCorrelationFactor_legacy:
#     def __init__(self,
#                  xdatcar,
#                  lattice,
#                  label='auto',
#                  force=None,
#                  interval=1,
#                  verbose=False,
#                  multi_vac=True,
#                  multi_path=True,
#                  correction_TS=True):
#         """
#         xdatcar : directory where XDATCAR_{label} exists
#         label   : label in XDATCAR_{label} (list or 'auto)
#         force   : directory where force_{label}.dat exists
#         """
        
#         self.xdatcar = xdatcar
#         self.force_dir = force
#         self.lattice = lattice
#         self.label = self.get_label() if label=='auto' else label
#         self.interval = interval
#         self.verbose = verbose
#         self.symbol = lattice.symbol
#         self.path = self.lattice.path

#         # corrections
#         self.multi_vac = multi_vac
#         self.multi_path = multi_path
#         self.correction_TS = correction_TS
        
#         # f_cor
#         self.times = []
#         self.path_name = []
#         self.path_dist = []
#         self.f_ensemble = []
#         self.path_seq_cum = []
#         self.msd_cum = 0
#         self.num_enc_cum = 0
#         self.label_err = []
#         self.get_correlation_factors()

#         # unknown path
#         self.path_name.append('unknown')
#         self.path_dist.append(0)
#         self.path_dist = np.array(self.path_dist)

#         # cumulative f_cor
#         self.f_avg = np.average(self.f_ensemble) # averaged over ensmebles
#         self.f_cum = None
#         self.get_cumulative_correlation_factor()

#         # print results
#         self.print_summary()


#     def get_label(self):
#         label = []
#         for filename in os.listdir(self.xdatcar):
#             if len(filename.split('_')) == 2:
#                 first, second = filename.split('_')
#                 if first == 'XDATCAR':
#                     label.append(second)
#         label.sort()
#         return label
    

#     def get_correlation_factors(self):
#         for label in tqdm(self.label, 
#                           bar_format='{l_bar}{bar:20}{r_bar}{bar:-10b}', 
#                           ascii=True, 
#                           desc=f'{RED}f_cor{RESET}'):
#             start = time.time()

#             print(f"# Label : {label}")
#             analyzer = self.get_analyzer(label)
#             if analyzer is None:
#                 end = time.time()
#                 continue
#             correlation = CorrelationFactor(analyzer, self.verbose)

#             self.f_ensemble.append(correlation.f_cor)
#             self.msd_cum += correlation.msd * correlation.num_enc
#             self.num_enc_cum += correlation.num_enc
#             self.path_seq_cum += correlation.path_sequence

#             end = time.time()
#             self.times.append(end-start)

#         self.f_ensemble = np.array(self.f_ensemble, dtype=float)
#         self.msd_cum /= self.num_enc_cum
        
#         for path in self.path:
#             self.path_name.append(path['name'])
#             self.path_dist.append(path['distance'])
    

#     def get_analyzer(self, 
#                      label):
#         path_xdatcar = os.path.join(self.xdatcar, f"XDATCAR_{label}")
#         if self.force_dir is not None:
#             path_force = os.path.join(self.force_dir, f"force_{label}.dat")
#         else:
#             path_force = None
            
#         traj = LatticeHopping(lattice=self.lattice,
#                               xdatcar=path_xdatcar,
#                               force=path_force,
#                               interval=self.interval)
        
#         if self.multi_vac:
#             traj.correct_multivacancy()
#             traj.check_multivacancy(verbose=False)
        
#         if traj.multi_vac is True:
#             print(f'multi-vacancy issue in label {label}.')
#             print('the calculation is skipped.')
#             self.label_err.append(label)
#             return None

#         if self.correction_TS and (self.force_dir is not None):
#             traj.correct_transition_state()

#         analyzer = Analyzer(traj, self.lattice)
#         analyzer.get_path_vacancy(verbose=self.verbose)

#         if self.multi_path:
#             analyzer.correct_multipath()

#         if self.verbose:
#             analyzer.print_summary(disp=False, 
#                                    save_figure=False, 
#                                    save_text=False)
#         print('')

#         return analyzer
    

#     def get_cumulative_correlation_factor(self):
#         self.counts_cum = np.zeros(len(self.path_name))

#         for i, name in enumerate(self.path_name):
#             self.counts_cum[i] = self.path_seq_cum.count(name)

#         self.counts_cum /= self.num_enc_cum
#         self.f_cum = self.msd_cum / np.sum(self.counts_cum * self.path_dist**2)


#     def print_summary(self):
#         print("## Summary")
#         print("#  Total counts")
#         print("           ", end="")
#         for name in self.path_name:
#             print("{:<10}".format(name), end=" ")
#         print('')
#         print("Distance   ", end="")
#         for dist in self.path_dist:
#             print("{:<10.3f}".format(dist), end=" ")
#         print('')
#         print("Counts     ", end="")
#         for count in self.counts_cum * self.num_enc_cum:
#             count = int(count)
#             print("{:<10}".format(count), end=" ")
#         print('\n')
#         print(f"      Mean correlation factor : {self.f_avg:.3f}")
#         print(f"Cumulative correlation factor : {self.f_cum:.3f}")
#         print('')

#         if self.verbose:
#             print("{:<10} {:<10}".format('Label', 'f_cor'))
#             for label, f in zip(self.label, self.f_ensemble):
#                 print("{:<10} {:<10.3f}".format(label, f))
#             print('')
#             print("{:<10} {:<10}".format('Label', 'Time(s)'))
#             for label, time in zip(self.label, self.times):
#                 print("{:<10} {:<10.3f}".format(label, time))
#             time_tot = np.sum(np.array(self.times))
#             print('')
#             print(f'Total time : {time_tot:.3f} s')

#         if len(self.label_err) > 0:
#             print('Error occured : ', end='')
#             for label in self.label_err:
#                 print(label, end=' ')
#             print('')
