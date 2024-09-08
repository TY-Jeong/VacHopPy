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
                
                # make input files
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
