import os
import sys
import copy   
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from tqdm import tqdm
from vachoppy import einstein
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
    def __init__(self, xs, ys, zs, *args, **kwargs):
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



def extractForce(file_in, 
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


    def addPath(self,
                name,
                site_init,
                site_final,
                d,
                Ea,
                dE='Nan'):
        
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
        
        self.path.append(dic_path)
        self.path_names.append(name)

        if not site_init in self.site_names:
            self.site_names += [site_init]
        
        if not site_final in self.site_names:
            self.site_names += [site_final]
    

    def printPath(self):
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


    def printLatticePoints(self):
        print("site\tcoord(direct)\tcoord(cartesian)")
        for lat_p in self.lat_points:
            print(lat_p['site'], end='\t')
            print(f"{lat_p['coord']}", end='\t')
            print(f"{lat_p['coord_C']}", end='\n')



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
        self.traj_on_lattice()

        # number of atoms preceding target
        self.count_before = 0 
        for i in range(self.idx_target):
            self.count_before += self.position[i]['num']

        # trajectory of vacancy
        self.idx_vac = {}
        self.traj_vac_C = {} 
        self.find_vacancy()

        # trace arrows
        self.trace_arrows = {}
        self.get_trace_lines()

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


    def distance_pbc(self, coord1, coord2):
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


    def getDisplacement(self, r1, r2):
        disp = r2 - r1
        disp[disp > 0.5] -= 1.0
        disp[disp < -0.5] += 1.0
        return np.dot(disp, self.lattice)
            

    def traj_on_lattice(self):
        traj = self.position[self.idx_target]['traj']

        # distance from lattice points
        disp = self.lat_points[np.newaxis,np.newaxis,:,:] - traj[:,:,np.newaxis,:]
        disp[disp > 0.5] -= 1.0
        disp[disp < -0.5] += 1.0
        disp = np.linalg.norm(np.dot(disp, self.lattice), axis=3)

        # save trajectory on lattice
        self.occ_lat_point = np.argmin(disp, axis=2)
        self.traj_on_lat_C = self.lat_points_C[self.occ_lat_point]


    def save_traj(self,
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


    def plot_lattice(self, ax, label=False):
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


    def correctionTS(self):
        traj = np.transpose(self.position[self.idx_target]['traj'], (1, 0, 2))
        for i in range(1, self.num_step):
            # check whether vacancy moves
            try:
                idx_pre = self.idx_vac[i-1][0]
                idx_now = self.idx_vac[i][0]
            except:
                print(f"error raised in step {i}. (correction TS)")
                print(f"step{i-1} : {self.idx_vac[i-1]}")
                print(f"step{i} : {self.idx_vac[i]}")
                continue

            if idx_pre == idx_now:
                continue

            try:
                atom = np.where((self.occ_lat_point[:,i]==idx_pre)==True)[0][0]
            except:
                print(f'error raised at step {i}.')
                print('please correct the vacancy-site your-self.')
                sys.exit(0)

            coord = traj[i, atom]
            force = self.forces[i, atom]

            r_pre = self.getDisplacement(coord, self.lat_points[idx_now])
            r_now = self.getDisplacement(coord, self.lat_points[idx_pre])

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
        self.trace_arrows = {}
        self.get_trace_lines()


    def get_trace_lines(self):
        """
        displaying trajectory of moving atom at each step.
        """
        num_target = self.num_atoms[self.idx_target]
        idx_atoms = np.arange(num_target)

        for step in range(1, self.num_step):
            arrows = []
            for idx in idx_atoms:
                traj_now = self.lat_points[self.occ_lat_point[idx][step]]
                traj_pre = self.lat_points[self.occ_lat_point[idx][step-1]]
                
                # check whether atom moves
                distance = self.distance_pbc(traj_now, traj_pre)
                if distance > 0.001:
                    arrow = {}
                    arrow['p'] = np.vstack((self.traj_on_lat_C[idx][step-1],
                                            self.traj_on_lat_C[idx][step]))
                    arrow['c'] = self.cmap[(idx+1)%len(self.cmap)]
                    arrow['lat_points'] = [self.occ_lat_point[idx][step-1],
                                           self.occ_lat_point[idx][step]]
                    arrows += [arrow]
                
            self.trace_arrows[step-1] = arrows


    def save_gif_PIL(self, 
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
                  dpi=300,
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
            self.save_gif_PIL(filename=filename,
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


    def save_traj_on_lat(self,
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


    def check_unique_vac(self):
        check = np.array([len(i) for i in self.idx_vac.values()])
        check = np.where((check==1) == False)[0]
        
        if len(check) > 0:
            self.multi_vac = True
            print('multi-vacancy issue raised:')
            print('  step :', end=' ')

            for i in check:
                print(i, end = ' ')
            print('')

        else:
            self.multi_vac = False
            print('vacancy is unique.')


    def update_vac(self,
                   step,
                   lat_point):
        """
        step: step which the user want to update the vacancy site
        lat_point: label of lattice point where vacancy exist at the step
        """
        self.idx_vac[step] = [lat_point]
        self.traj_vac_C[step] = np.array([self.lat_points_C[lat_point]])


    def check_connectivity(self, start=1):
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
                self.update_vac(step, vac_site)
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
                self.update_vac(step, vac_site)
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
                self.update_vac(step, vac_site)
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
        self.trace_arrows = {}
        self.get_trace_lines()
    

    # def check_connectivity_backup(self, start=1):
    #     """
    #     correction for multi-vacancy issue
    #     correction starts from 'start' step
    #     """
    #     trace_lines = self.trace_arrows
    #     vac_site = self.idx_vac[0][0]

    #     for step in range(start, self.num_step):
    #         # when only one vacancy exist
    #         if len(self.idx_vac[step]) == 1:
    #             vac_site = self.idx_vac[step][0]
    #             self.update_vac(step, vac_site)
    #             continue

    #         # when multiple vacancies exsit
    #         #    when vacancy is stationary
    #         if vac_site in self.idx_vac[step]:
    #             self.update_vac(step, vac_site)
    #             continue

    #         # when vacancy moves
    #         #   find connected points with vacancy
    #         points = [vac_site]
    #         while True:
    #             check1 = len(points)
    #             for dic in trace_lines[step-1]:
    #                 if len(list(set(points) & set(dic['lat_points']))) == 1:
    #                     points += dic['lat_points']
    #                     points = list(set(points))
                
    #             check2 = len(points)

    #             # no more connected points
    #             if check1 == check2:
    #                 break

    #         site = list(set(points) & set(self.idx_vac[step]))
            
    #         if len(site) == 1:
    #             vac_site = site[0]
    #             self.update_vac(step, vac_site)
            
    #         elif len(site) == 0:
    #             print("there is no connected site.")       
    #             print(f"find the vacancy site for your self. (step: {step})")
    #             break
            
    #         else:
    #             print("there are multiple candidates.")       
    #             print(f"find the vacancy site for your self. (step: {step})")
    #             break
            


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

        self.idx_vac = self.traj.idx_vac # (dic) index of vacancy
            

    def search_path_vac(self, verbose=True):
        step_unknown = []
        self.path_vac = []
        idx = 0 # index in self.path_vac
        
        for step in range(self.traj.num_step-1):
            coord_init = self.lat_points[self.traj.idx_vac[step][0]]['coord']
            coord_final = self.lat_points[self.traj.idx_vac[step+1][0]]['coord']
            
            # check whether vacancy moves
            distance = self.traj.distance_pbc(coord_init, coord_final)
            
            if distance > self.tolerance:
                site_init = self.lat_points[self.traj.idx_vac[step][0]]['site']
                path = self.get_path(site_init, distance)
                
                path['step'] = step+1
                
                if path['name'] == self.path_unknown['name']:
                    step_unknown += [step+1]
                    
                self.path_vac += [copy.deepcopy(path)]
                idx += 1
                
        if len(step_unknown) > 0 and verbose:
            print(f"unknown steps are detected.: {step_unknown}")
        

    def get_path(self,
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
            

    def print_path_vac(self):
        if self.path_vac is None:
            print("path_vac is not defines.")
            print("please run 'search_path_vac'.")
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
                         sort = True):
        
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
        
        arrows = np.zeros((len(self.traj.trace_arrows[step-1]),2 ))
        
        for i, dic_arrow in enumerate(self.traj.trace_arrows[step-1]):
            arrows[i] = dic_arrow['lat_points']
        
        vac_now = self.traj.idx_vac[step][0]
        vac_pre = self.traj.idx_vac[step-1][0]
        
        path = self.path_tracer(arrows, vac_now, vac_pre)
        path = np.array(path, dtype=int)

        # update index of lattice points occupied by vacancy
        self.idx_vac[step] = path[-2::-1]
        
        return path
    

    def unwrap_path(self):
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
                    distance = self.traj.distance_pbc(coord_init, coord_final)
                    
                    p_new = self.get_path(site_init, distance)
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
                      save_figure=True,
                      save_text=True,
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
        
        self.plot_path_counts(figure=figure,
                              text=text,
                              disp=disp,
                              save_figure=save_figure,
                              save_text=save_text,
                              bar_width=bar_width,
                              bar_color=bar_color,
                              dpi=dpi,
                              sort=sort)


class CorrelationFactor:
    def __init__(self, 
                 lattice,
                 xdatcar,
                 label='auto',
                 interval=1,
                 verbose=True):
        """
        encounter MSD is calculated using projected vacancy path
        """
        # save arguments
        self.lattice = lattice
        self.xdatcar = xdatcar
        self.interval = interval
        self.label = self.getLabel() if label=='auto' else label
        self.verbose = verbose
        
        # lattice information
        self.symbol = lattice.symbol
        self.path = lattice.path
        for p in self.path:
            p['counts'] = 0
        
        # unknown path
        self.path_unknown = {}
        self.path_unknown['name'] = 'unknown'
        self.path_unknown['counts'] = 0

        # (test)
        self.cum_displacement = 0
        
        # correlation factor of each ensemble 
        self.f_ensemble = []
        self.getCorrelationFactor()
        self.f_ensemble = np.array(self.f_ensemble)

        # correlation factor
        self.f_cor = self.f_ensemble[~np.isnan(self.f_ensemble)]
        self.f_cor = np.average(self.f_cor)
        
        if self.verbose:
            self.print_fcor()

    def getLabel(self):
        label = []
        for filename in os.listdir(self.xdatcar):
                if len(filename.split('_')) == 2:
                    first, second = filename.split('_')
                    if first == 'XDATCAR':
                        label.append(second)
        label.sort()
        return label

    def encounterMSD(self, analyzer):
        path_vac_idx = []  
        for i in range(analyzer.traj.num_step):
            path_vac_idx += list(analyzer.idx_vac[i])

        coord_vac = [analyzer.traj.lat_points[idx] for idx in path_vac_idx]
        coord_vac = np.array(coord_vac)
        
        displacement = np.zeros_like(coord_vac)
        displacement[1:,:] = np.diff(coord_vac, axis=0)

        displacement[displacement > 0.5] -= 1.0
        displacement[displacement < -0.5] += 1.0
        
        displacement = np.sum(displacement, axis=0)
        self.cum_displacement += displacement # (test)

        lattice = analyzer.traj.lattice
        return np.sum(np.dot(displacement, lattice)**2)
    
    def randomwalkMSD(self, analyzer):
        path_vac = []
        for p_vac in analyzer.path_vac:
            print(p_vac['name'], end=' ')
            path_vac.append(p_vac['name'])
        print('')
        
        msd_rand = 0
        for p in self.path:
            num_p = path_vac.count(p['name'])
            p['counts'] += num_p
            msd_rand += num_p * p['distance']**2
        self.path_unknown['counts'] += path_vac.count(self.path_unknown['name'])
        return msd_rand
    
    def makeAnalyzer(self, label):
        path_xdatcar = os.path.join(self.xdatcar, f"XDATCAR_{label}")
        traj = LatticeHopping(lattice=self.lattice,
                              xdatcar=path_xdatcar,
                              interval=self.interval)
        traj.check_connectivity()
        traj.check_unique_vac()

        analyzer = Analyzer(traj=traj,
                            lattice=self.lattice)
        analyzer.search_path_vac(verbose=False)
        analyzer.unwrap_path()
        return analyzer
    
    def getCorrelationFactor(self):
        for label in tqdm(self.label,
                          bar_format='{l_bar}{bar:20}{r_bar}{bar:-10b}',
                          ascii=True,
                          desc=f'{RED}f_cor{RESET}'):
            print(label)
            analyzer = self.makeAnalyzer(label)
            
            msd_enc = self.encounterMSD(analyzer)
            msd_rand = self.randomwalkMSD(analyzer)

            if msd_rand > 1e-5:
                f_label = msd_enc/msd_rand
                self.f_ensemble.append(f_label)
                print(f"f_cor = {f_label}")
            else:
                self.f_ensemble.append(np.NaN)
                print("no hopping detected.")
            print('')
    
    def plotCounts(self,
                   title='',
                   disp=True,
                   save=True,
                   outfile='counts.png',
                   dpi=300):
        
        name, Ea, counts = [], [], []
        for p in self.path:
            name.append(p['name'])
            Ea.append(p['Ea'])
            counts.append(p['counts'])
        
        # unknown path
        name.append('U')
        Ea.append(np.inf)
        counts.append(self.path_unknown['counts'])

        num_path = len(name)
        Ea, counts = np.array(Ea), np.array(counts)
        
        # sort by activation energy
        idx_sort = Ea.argsort()
        counts_sort = counts[idx_sort]
        name_sort = []
        for idx in idx_sort:
            name_sort.append(name[idx])

        # plot counts
        x = np.arange(num_path)
        plt.bar(x, counts_sort)
        plt.xticks(x, name_sort)
        plt.ylabel('Counts')
        plt.title(title + f" (total counts: {np.sum(counts)})")
        if save:
            plt.savefig(outfile, dpi=dpi)
        if disp:
            plt.show()

    def printCounts(self):
        name, Ea, counts = [], [], []
        for p in self.path:
            name.append(p['name'])
            Ea.append(p['Ea'])
            counts.append(p['counts'])
        
        # unknown path
        name.append('U')
        Ea.append(np.inf)
        counts.append(self.path_unknown['counts'])
        Ea, counts = np.array(Ea), np.array(counts)
        
        # sort by activation energy
        idx_sort = Ea.argsort()
        counts_sort = counts[idx_sort]
        name_sort = []
        for idx in idx_sort:
            name_sort.append(name[idx])
        
        print(f"total counts : {np.sum(counts)} ")
        for n, c in zip(name_sort, counts_sort):
            print(f"{n}({c})", end=' ')

    def print_fcor(self):
        print("label\tf_cor")
        for label, f in zip(self.label, self.f_ensemble):
            print(f"{label}\t{f}")
        print('')
        print(f"averaged f_cor = {self.f_cor}")


class CorrelationFactor_legacy:
    def __init__(self,
                 lattice,
                 fraction,
                 xdatcar,
                 label='auto',
                 interval=1):
        """
        encounter MSD is calculated by Einstein relation
        """
        
        self.lattice = lattice
        self.xdatcar = xdatcar
        self.fraction = fraction
        self.interval = interval
        self.label = self.getLabel() if label=='auto' else label
        
        # lattice information
        self.symbol = lattice.symbol
        self.path = lattice.path

        # unknown path
        self.path_unknown = {}
        self.path_unknown['name'] = 'unknown'
        self.path_unknown['counts'] = 0

        # encounter MSD
        # print("Calculating encounter MSD...")
        msd_enc = self.encounterMSD()

        # random walk MSD
        # print("\nCalculating random walk MSD...")
        msd_random = self.randomMSD()

        # correlation factor
        self.f = msd_enc / msd_random
        print(f"correlation factor = {self.f}")

    def getLabel(self):
        label = []
        for filename in os.listdir(self.xdatcar):
                if len(filename.split('_')) == 2:
                    first, second = filename.split('_')
                    if first == 'XDATCAR':
                        label.append(second)
        label.sort()
        return label

    def encounterMSD(self):
        ensembleEin = einstein.EnsembleEinstein(symbol=self.symbol,
                                                prefix=self.xdatcar,
                                                labels=self.label,
                                                segments=1,
                                                skip=0,
                                                start=None)
        return ensembleEin.msd[-1] / self.fraction
    
    def makeAnalyzer(self, label):
        path_xdatcar = os.path.join(self.xdatcar, f"XDATCAR_{label}")

        traj = LatticeHopping(lattice=self.lattice,
                              xdatcar=path_xdatcar,
                              interval=self.interval)
        traj.check_connectivity()
        traj.check_unique_vac()

        analyzer = Analyzer(traj=traj,
                            lattice=self.lattice)
        analyzer.search_path_vac(verbose=False)
        analyzer.unwrap_path()
        return analyzer
    
    def randomMSD(self):
        # path of vacancy
        path_vac = []
        for l in tqdm(self.label,
                      bar_format='{l_bar}{bar:20}{r_bar}{bar:-10b}',
                      ascii=True,
                      desc=f'{RED}MSD_rand_walk{RESET}'):
            print(l)
            analyzer = self.makeAnalyzer(label=l)
            for p_vac in analyzer.path_vac:
                print(p_vac['name'], end=' ')
                path_vac.append(p_vac['name'])
            print('\n')

        msd_random = 0
        self.path_unknown['counts'] = path_vac.count(self.path_unknown['name'])
        for p in self.path:
            p['counts'] = path_vac.count(p['name'])
            msd_random += p['counts'] * p['distance']**2
        msd_random /= len(self.label)
        return msd_random
    
    def plotCounts(self,
                   title='',
                   disp=True,
                   save=True,
                   outfile='counts.png',
                   dpi=300):
        
        name, Ea, counts = [], [], []
        for p in self.path:
            name.append(p['name'])
            Ea.append(p['Ea'])
            counts.append(p['counts'])
        
        # unknown path
        name.append('U')
        Ea.append(np.inf)
        counts.append(self.path_unknown['counts'])

        num_path = len(name)
        Ea, counts = np.array(Ea), np.array(counts)
        
        # sort by activation energy
        idx_sort = Ea.argsort()
        counts_sort = counts[idx_sort]
        name_sort = []
        for idx in idx_sort:
            name_sort.append(name[idx])

        # plot counts
        x = np.arange(num_path)
        plt.bar(x, counts_sort)
        plt.xticks(x, name_sort)
        plt.ylabel('Counts')
        plt.title(title + f" (total counts: {np.sum(counts)})")
        if save:
            plt.savefig(outfile, dpi=dpi)
        if disp:
            plt.show()

    def printCounts(self):
        name, Ea, counts = [], [], []
        for p in self.path:
            name.append(p['name'])
            Ea.append(p['Ea'])
            counts.append(p['counts'])
        
        # unknown path
        name.append('U')
        Ea.append(np.inf)
        counts.append(self.path_unknown['counts'])
        Ea, counts = np.array(Ea), np.array(counts)
        
        # sort by activation energy
        idx_sort = Ea.argsort()
        counts_sort = counts[idx_sort]
        name_sort = []
        for idx in idx_sort:
            name_sort.append(name[idx])
        
        print(f"total counts : {np.sum(counts)}")
        for n, c in zip(name_sort, counts_sort):
            print(f"{n}({c})", end=' ')




        



                
    


