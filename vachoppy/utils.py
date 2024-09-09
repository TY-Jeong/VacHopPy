import os
import sys
import shutil
import numpy as np
from .amorphous import genInput


class getMDset:
    def __init__(self,
                 path_poscar,
                 temp,
                 label=None,
                 potcar='pbe',
                 nsw=10000,
                 potim=2,
                 charge=0,
                 ncore=4):
        """
        Arg 1: (str) path_poscar; path of directory containing poscar files
        Arg 2: (list) temp; temperature in K 
        Arg 3: (list; opt) label; labels of poscar. poscar format should be POSCAR_{label}
        """
        self.path_poscar = path_poscar
        self.label = label
        self.temp = temp
        self.potcar = potcar
        self.nsw = nsw
        self.potim = potim
        self.charge = charge
        self.ncore = ncore

        if self.label is None:
            self.label = []
            path_now = os.getcwd()
            for name in os.listdir(self.path_poscar):
                if len(name.split('_')) == 2:
                    poscar, label = name.split('_')
                    if poscar=='POSCAR':
                        self.label += [label]
        
        self.foldername=[]
        self.make_input_set()


    def make_input_set(self):
        path_now = os.getcwd()
        for t in self.temp:
            outer_folder = f"{t}K"
            for l in self.label:
                # make folder
                path_dir = os.path.join(outer_folder, f"{l}")
                self.foldername += [path_dir]
                os.makedirs(path_dir, exist_ok=True)
                
                # copy poscar
                from_pos_path = os.path.join(self.path_poscar, f"POSCAR_{l}")
                to_pos_path = os.path.join(path_dir, 'POSCAR')
                shutil.copyfile(from_pos_path, to_pos_path)
                
                # make inpur files
                os.chdir(path_dir)
                _ = genInput(potcar=self.potcar,
                             nsw=self.nsw,
                             potim=self.potim,
                             temp=t,
                             charge=self.charge,
                             ncore=self.ncore)
                os.chdir(path_now)                



class getMDresult:
    def __init__(self,
                 temp=None,
                 label=None,
                 outdir='xdatcar'):
        """
        Arg 1: (list; opt) temp;
        Arg 2: (list; opt) list;
        """
        
        self.temp = temp
        self.label = label
        self.outdir = outdir

        # folders where MD was conducted
        self.foldername=[]
        if self.temp is None and self.label is None:
            self.auto_search()
        else:
            for t in self.temp:
                for l in self.label:
                    foldername = f"{t}K/{l}"
                    self.foldername += [foldername]
        
        # copy XDATCAR files
        if os.path.isdir(self.outdir):
            os.makedirs(self.outdir, exist_ok=True)
            
        check_temp = []   
        for path in self.foldername:
            temp, label = path.split('/')
            path_save = os.path.join(self.outdir, f"xdatcar.{temp}")
            os.makedirs(path_save, exist_ok=True)
            
            # copy xdatcar
            from_xdat_path = os.path.join(path, 'XDATCAR')
            to_xdat_path = os.path.join(path_save, f'XDATCAR_{label}')
            if not os.path.isfile(from_xdat_path):
                print(f"no XDATCAR in {path}.")
            else:
                shutil.copyfile(from_xdat_path, to_xdat_path)    

            # copy outcar
            if not temp in check_temp:
                from_out_path = os.path.join(path, 'OUTCAR')
                to_out_path = os.path.join(path_save, 'OUTCAR')
                if not os.path.isfile(from_out_path):
                    print(f"no OUTCAR in {path}.")
                else:
                    check_temp += [temp]
                    shutil.copyfile(from_out_path, to_out_path)


    def auto_search(self):
        path_now = os.getcwd()
        for name_out in os.listdir(path_now):
            if os.path.isdir(name_out) and name_out[-1]=='K':
                outer_folder = name_out
                os.chdir(name_out)
                for name_in in os.listdir(os.getcwd()):
                    if os.path.isdir(name_in):
                        inner_folder = name_in
                        foldername = os.path.join(outer_folder, inner_folder)
                        self.foldername += [foldername]
                os.chdir(path_now)



def concat_xdatcar(xdatcar1, 
                  xdatcar2, 
                  outcar1, 
                  outcar2, 
                  xdatcar_out, 
                  outcar_out):

    if not os.path.isfile(xdatcar1):
        print(f"{xdatcar1} does not exist.")
        sys.exit(0)

    if not os.path.isfile(xdatcar2):
        print(f"{xdatcar2} does not exist.")
        sys.exit(0)

    if not os.path.isfile(outcar1):
        print(f"{outcar1} does not exist.")
        sys.exit(0)

    if not os.path.isfile(outcar2):
        print(f"{outcar2} does not exist.")
        sys.exit(0)

    # read outcars
    with open(outcar1, 'r') as f:
        for line in f:
            if 'NSW' in line:
                nsw1 = int(line.split()[2])
                break

    with open(outcar2, 'r') as f:
        for line in f:
            if 'NSW' in line:
                nsw2 = int(line.split()[2])
                break

    # write outcar_out
    if not os.path.isfile(outcar_out):
        with open(outcar1, 'r') as f1:
            with open(outcar_out, 'w') as f_out:
                for line in f1:
                    if 'NSW' in line:
                        line = f"   NSW    =  {nsw1+nsw2}    number of steps for IOM\n"
                    f_out.write(line)
                
    # copy xdatcar1
    shutil.copyfile(xdatcar1, xdatcar_out)

    # concatnate xdatcars
    with open(xdatcar2, 'r') as f:
        lines = [line.strip() for line in f]

    num_atoms = np.array(lines[6].split(), dtype=int)
    num_atoms = np.sum(num_atoms)

    coords = np.zeros((nsw2, num_atoms, 3))
    for i in range(nsw2):
        idx_i = 8 + (num_atoms+1)*i
        idx_f = idx_i + num_atoms

        coord = lines[idx_i:idx_f]
        coord = [list(map(float, s.split())) for s in coord]
        coord = np.array(coord)
        coords[i] = coord

    with open(xdatcar_out, 'a') as f:
        for idx, coord in enumerate(coords, start=1):
            step = str(nsw1 + idx).rjust(6)
            f.write(f'Direct configuration={step}\n')
            for x, y, z in coord:
                f.write("   %.8f  %.8f  %.8f\n"%(x,y,z))


def concat_force(force1, 
                force2, 
                force_out):
    if not os.path.isfile(force1):
        print(f'there is no {force1}.')
        sys.exit(0)

    if not os.path.isfile(force2):
        print(f'there is no {force2}.')
        sys.exit(0)
        
    with open(force1, 'r') as f:
        lines = [s.split()[0].isalpha() for s in f]
        lines = np.array(lines)
    nsw1 = len(np.where(lines==True)[0])

    with open(force2, 'r') as f:
        lines = [s for s in f]
        
    shutil.copyfile(force1, force_new)
    with open(force_new, 'a') as f:
        for s in lines:
            if s.split()[0].isalpha():
                step = int(s.strip().split()[1])
                step += nsw1
                f.write(f'Iteration {step}\n')
            else:
                f.write(s)
                

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
            nsw = int(line.split('>')[-2].split('<')[0])
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