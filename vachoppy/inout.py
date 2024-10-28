import os
import sys
import numpy as np

class DataInfo:
    def __init__(self, 
                 prefix1='traj',
                 prefix2='traj',
                 verbose=False):
        self.prefix1 = prefix1
        self.prefix2 = prefix2
        if prefix1 not in os.listdir(os.getcwd()):
            print(f'{self.prefix1} directory is not found.')
            sys.exit(0)
            
        self.outcar = []
        self.xdatcar = []
        self.force = []
        
        # temperature
        self.temp = None
        self.get_temperature()
        
        # read outcar
        self.potim = np.zeros_like(self.temp)
        self.nblock = np.zeros_like(self.temp)
        self.nsw = np.zeros_like(self.temp)
        self.read_outcar()
        
        # read label
        self.label = []
        self.read_labels()
        
        # read force
        self.force = None
        self.read_force()
        
        if verbose:
            self.summary()
        
    def get_temperature(self):
        dirList = os.listdir(os.path.join(os.getcwd(), self.prefix1))
        dirList.sort()
        self.temp = []
        for dir in dirList:
            if dir.split('.')[0] == self.prefix2:
                self.temp.append(int(dir.split('.')[-1].split('K')[0]))
        self.temp = np.array(self.temp)
        
    def read_outcar(self):
        for i, T in enumerate(self.temp):
            outcar = [self.prefix1, f'{self.prefix2}.{T}K', 'OUTCAR']
            outcar = os.path.join(*outcar)
            self.outcar.append(outcar)
            check_potim, check_nblock, check_nsw = False, False, False
            with open(outcar, 'r') as f:
                for line in f:
                    if 'POTIM ' in line:
                        self.potim[i] = float(line.split()[2])
                        check_potim = True
                    if 'NBLOCK' in line:
                        self.nblock[i] = int(line.split()[2].replace(';',''))
                        check_nblock = True
                    if 'NSW' in line:
                        self.nsw[i] = int(line.split()[2])
                        check_nsw = True
                    if check_potim and check_nblock and check_nsw:
                        break
    
    def read_labels(self):
        for T in self.temp:
            list_xdat = []
            path_dir = os.path.join(self.prefix1, f'{self.prefix2}.{T}K')
            fileList = os.listdir(path_dir)
            fileList.sort()
            label = []
            for file in fileList:
                _file = file.split('_')
                if _file[0] == 'XDATCAR' and len(_file) > 1:
                    label.append(file.split('_')[1])
                    list_xdat.append(os.path.join(path_dir, file))
            self.xdatcar.append(list_xdat)
            self.label.append(label)
            
    def read_force(self):
        self.force = []
        check = True
        for i, T in enumerate(self.temp):
            list_force = []
            path_dir = os.path.join(self.prefix1, f'{self.prefix2}.{T}K')
            for label in self.label[i]:
                if f'FORCE_{label}' not in os.listdir(path_dir):
                    check = False
                    break
                else:
                    list_force.append(os.path.join(path_dir, f'FORCE_{label}'))
            self.force.append(list_force)
        if check is False:
            self.force = None
            
    def summary(self):
        print('Temperature (K) :', end=' ')
        for T in self.temp:
            print(T, end=' ')
        print('')
        use_force = False if self.force is None else True
        print(f'Use FORCE data : {use_force}')
        print('')
        for i in range(len(self.temp)):
            print(f'T = {self.temp[i]} K :')
            print(f'  number of labels = {len(self.label[i])}')
            print(f'  nsw = {self.nsw[i]}')
            print(f'  potim = {self.potim[i]}')
            print(f'  nblock = {self.nblock[i]}')
            print(f'  total md time = {self.nsw[i]*self.potim[i]/1000} ps')
            print('')   