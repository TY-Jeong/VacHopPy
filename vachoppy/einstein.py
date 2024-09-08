import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from colorama import Fore   # color map for tqdm

GREEN = '\033[92m' # Green color
RED = '\033[91m'   # Red color
RESET = '\033[0m'  # Reset to default color

class EinsteinRelation:
    def __init__(self, 
                 xdatcar="",
                 outcar="",
                 segments=1, 
                 skip=0,
                 skip2=0, 
                 verbose=True,
                 getD=True):
        
        if os.path.isfile(xdatcar):
            self.xdatcar = xdatcar
        else:
            print(f"{xdatcar} is not found.")
            sys.exit(0)

        if os.path.isfile(outcar):
            self.outcar = outcar
        else:
            print(f"{outcar} is not found.")
            sys.exit()
            
        self.verbose = verbose
        
        # read outcar
        self.potim = None
        self.nblock = None
        self.readOUTCAR()

        # read XDATCAR
        self.TypeName = None    # List of atom species
        self.Ntype = None       # Number of atom species
        self.Nions = None       # Total number of atom
        self.Nelem = None       # Number of each atoms
        self.Niter = None       # Number of steps
        self.position = None    # Direct coordination of atoms
        self.cell = None        # Lattice of cell
        self.readXDATCAR()

        self.ChemSymb = None
        self.symbolList()

        # Parameters for skip
        self.segments = segments    # Dividing whole MD steps into segments
        self.skip = skip            # Skipping initial steps for whole MD steps
        self.skip2 = skip2          # Skipping initial steps for each segments

        Time = self.Niter * self.potim * self.nblock
        if self.verbose:
            print("\t--- MD information ---")
            print(f"\ttotal time = {Time/1000} ps")

        # Mean square distance
        self.msd = None
        self.getMSD()
        
        # Diffusivity.
        # For ensemble calculation, please set getD = False to save time
        self.diffcoeff = None
        self.intercept = None
        self.ddiffcoeff = None  # Deviation of D

        if getD is True:
            self.getDiff()

        if self.verbose:
            self.plotMSD()
    
    
    def readOUTCAR(self):
        if self.verbose:
            print(f"reading {self.outcar}...")
        with open(self.outcar, 'r') as file:
            lines_outcar = [line.strip() for line in file]

        self.potim, self.nblock, lm = 0, 0, 0
        for num, line in enumerate(lines_outcar):
            if 'POTIM ' in line:
                self.potim = float(line.split()[2])
            if 'NBLOCK' in line:
                self.nblock = int(line.split()[2].replace(';',''))
            if 'Mass of Ions in am' in line:
                lm = num + 1
            if self.potim and self.nblock and lm:
                break

        if self.verbose:
            print(f"\tpotim = {self.potim}")
            print(f"\tnblock = {self.nblock}\n")
    
    
    def readXDATCAR(self):
        if self.verbose:
            print(f"reading {self.xdatcar}...")
        with open(self.xdatcar,'r') as file:
            inp = [line.strip() for line in file]
        scale = float(inp[1])
        self.cell = np.array([line.split() for line in inp[2:5]], dtype=float)
        self.cell *= scale

        ta = inp[5].split()
        tb = inp[6].split()
        if ta[0].isalpha():
            self.TypeName = ta
            self.Ntype = len(ta)
            self.Nelem = np.array(tb, dtype=int)
            self.Nions = self.Nelem.sum()
        if self.verbose:
            print(f"\tatom type = {self.TypeName}")
            print(f"\tnumber of type = {self.Ntype}")
            print(f"\tnumber of atom = {self.Nelem}")
            print(f"\ttotal number of atom = {self.Nions}")

        pos = np.array(
            [line.split() for line in inp[7:] if not line.split()[0].isalpha()],
            dtype=float)
        self.position = pos.flatten().reshape((-1,self.Nions,3))
        self.Niter = self.position.shape[0]
        if self.verbose:
            print(f"\tshape of position = {self.position.shape} #(nsw, number of atom, xyz)\n")


    def symbolList(self):
        self.ChemSymb = []
        for i in range(self.Ntype):
            self.ChemSymb += [np.tile(self.TypeName[i], self.Nelem[i])]
        self.ChemSymb = np.concatenate(self.ChemSymb)
    
    
    def getMSD(self):
        # length of segment
        seglength = int(np.floor((self.Niter-self.skip)/self.segments))
        if self.verbose:
            print(f"\tnumber of segments = {self.segments}")
            print(f"\tsegment length = {seglength}")

        # msd for (each segment), (each type of element), (each step), and, (x, y, z)
        self.msd=np.zeros(shape=(self.segments, self.Ntype, seglength, 3))

        for i in range(self.segments):
            # displacement [step,atom_idx,xyz]
            displacement = np.zeros_like(
                self.position[self.skip+seglength*i:self.skip+seglength*(i+1),:])
            displacement[0,:,:] = 0
            displacement[1:,:,:] = np.diff(
                self.position[self.skip+seglength*i:self.skip+seglength*(i+1),:], axis=0)

            # wrap back into box
            displacement[displacement > 0.5] -= 1.0
            displacement[displacement < -0.5] += 1.0

            # vector of accumulated displacement
            displacement = np.cumsum(displacement, axis=0)
        
            displacementC = np.zeros_like(displacement)
            displacementC = np.dot(displacement[:,:,:], self.cell)  # Direct to Cartesian

            # squared displacement in x, y, z
            displacementC[:,:,:] =  displacementC[:,:,:]**2

            for j in range(self.Ntype):
                if j==0:
                    labelst=0
                else:
                    labelst = np.sum(self.Nelem[:j])
                labeled = np.sum(self.Nelem[:j+1])
                self.msd[i,j,:,:] = np.average(displacementC[:,labelst:labeled,:], axis=1) 

        if self.verbose:
            print(f"\tshape of msd = {self.msd.shape} #(segment, number of type, nsw - skip, xyz)\n")
                
                
    def getDiff(self):
        self.diffcoeff = np.zeros(shape=(self.segments, self.Ntype, 3))
        self.intercept = np.zeros(shape=(self.segments, self.Ntype, 3))
        self.ddiffcoeff = np.zeros(shape=(self.Ntype, 3))

        for i in range(self.segments):           # loop over segments
            for j in range(self.Ntype):          # loop over type of elements
                for k in range(3):               # loop over x.y,z
                    # time unit in fs, distance unit in A, unit in m^2/s
                    if np.floor((self.Niter-self.skip)/self.segments) < 2*self.skip2:   
                        # too large skip2... skip2 is automatically set to 0
                        self.diffcoeff[i,j,k], self.intercept[i,j,k]=np.polyfit(
                                self.potim*self.nblock*
                                np.arange(np.floor((self.Niter-self.skip)/self.segments)),
                                1E-5*self.msd[i,j,:,k], 
                                deg=1)
                    else:
                        self.diffcoeff[i,j,k], self.intercept[i,j,k] = np.polyfit(
                                self.potim*self.nblock*
                                np.arange(np.floor((self.Niter-self.skip)/self.segments)-self.skip2), 
                                1E-5*self.msd[i,j,self.skip2:,k],
                                deg=1)
                        
        # standard deviation of D over segments for x,y,z directions
        self.diffcoeff /= 2      # x, y, z 각각에 대해 계산되었음 -> 1차원임 -> 2로 나누어짐
        self.ddiffcoeff = np.std(self.diffcoeff, axis=0) / 2
    
    
    def getMSD_Atom(self, 
                    symbol):
        if symbol in self.TypeName:
            idx = np.where(np.array(self.TypeName) == symbol)[0][0]
            return self.msd[:,idx,:,:]  # shape: [1, steps, 3:xyz]
        else:
            print("No matched atom!")
            return None
        
        
    def plotMSD(self):
        time = np.arange(self.skip, self.Niter) * self.potim / 1000
        for type in self.TypeName:
            msd = np.sum(self.getMSD_Atom(symbol=type).squeeze(), axis=1)
            plt.plot(time, msd, label=type)
        plt.xlabel('t (ps)', fontsize=13)
        plt.ylabel(r'MSD ($Å^2$)', fontsize=13)
        # plt.ylabel(r'$<\Delta r^{2}> (Å^2)$', fontsize=13)
        plt.legend(fontsize=13)
        plt.show()



class EnsembleEinstein:
    def __init__(self, 
                 symbol, 
                 prefix, 
                 labels, 
                 segments, 
                 skip, 
                 start=None, 
                 end=None):
        ## XDATCAR format : XDATCAR_{label}
        ## OUTCAR format : OUTCAR

        ## Arguments
        self.symbol = symbol        # ex: 'O'
        self.prefix = prefix        # ex: 'xdatcar/xdatcar.2000K/'
        if self.prefix[-1] != '/':
            self.prefix += '/'
        self.labels = labels        # ex: np.arange(1000, 20000+1, 1000)
        self.segments = segments    # See MD2D for description about segments
        self.skip = skip
        ## args : 'start' and 'end' : Fitting range

        # Properties of NVT ensembles
        self.msd_xyz = None         # shape: (steps, 3;xyz)
        self.msd = None             # shape: (steps, )
        self.msd_x = None
        self.msd_y = None
        self.msd_z = None

        self.potim = None
        self.nblock = None

        self.diffcoeff = None       # m2/s
        self.intercept = None
        
        self.diffcoeff_x = None     
        self.intercept_x = None
        
        self.diffcoeff_y = None     
        self.intercept_y = None
        
        self.diffcoeff_z = None     
        self.intercept_z = None

        # For linear fitting
        self.start = start          # Step number in MD simulation
        self.end = end
        self.timestep = None        # [ps]

        self.getEnsembleMSD()

        if start is not None and end is not None:
            _, __ = self.getEnsembleD()
            _, __ = self.getEnsembleD_x()
            _, __ = self.getEnsembleD_y()
            _, __ = self.getEnsembleD_z()


    def getEnsembleMSD(self):
        desc = self.prefix.split('/')[-1] if len(self.prefix.split('/')[-1]) > 0 \
            else  self.prefix.split('/')[-2]
        for label in tqdm(self.labels,
                          bar_format='{l_bar}{bar:20}{r_bar}{bar:-10b}',
                          ascii=True,
                          desc=f'{RED}{desc}{RESET}'):
            xdatcar = self.prefix + "XDATCAR_" + str(label)
            outcar = self.prefix + "OUTCAR"
            ensemble = EinsteinRelation(xdatcar=xdatcar,
                                        outcar=outcar,
                                        segments=self.segments,
                                        skip=self.skip,
                                        verbose=False,
                                        getD=False)
            if self.msd_xyz is None:
                self.msd_xyz = ensemble.getMSD_Atom(symbol=self.symbol)
                self.potim = ensemble.potim
                self.nblock = ensemble.nblock
            else:
                self.msd_xyz = np.append(self.msd_xyz, 
                                         ensemble.getMSD_Atom(symbol=self.symbol), axis=0)
        
        self.msd_xyz = np.mean(self.msd_xyz, axis=0)
        self.msd = np.sum(self.msd_xyz, axis=1)
        self.timestep = np.arange(self.msd[:].shape[0]) * self.potim * self.nblock / 1000


    def getEnsembleD(self, 
                     *args):
        if len(args) == 2:
            self.start, self.end = args[0], args[1]
        if self.start is None or self.end is None:
            sys.exit("start and end are needed to determine fitting range.")

        self.diffcoeff, self.intercept = np.polyfit(6*self.potim*self.nblock*np.arange(self.end-self.start),
                                                    1E-5*self.msd[self.start:self.end], 
                                                    deg=1)
        return self.diffcoeff, self.intercept
    
    def plotEnsembleMSD(self):
        temp = self.prefix.split('/')[-2].split('.')[-1]
        plt.plot(self.timestep, self.msd[:], label=f"{temp}")

        if self.diffcoeff is not None:
            # visualize fitting line
            plt.plot(np.arange(self.start, self.end)*self.potim*self.nblock/1000, 
                     self.diffcoeff*(6*self.potim*self.nblock*1.0E5)*
                     np.arange(self.end-self.start)+self.intercept*1.0E5, 
                     'k:')
            # visualize fitting range
            plt.axvline(self.start*self.potim*self.nblock/1000, 0, 1, color='black', linewidth=0.5)
            plt.axvline(self.end*self.potim*self.nblock/1000, 0, 1, color='black', linewidth=0.5)

        plt.legend(loc='upper left')
        plt.xlabel('t (ps)', fontsize=13)
        plt.ylabel(r'MSD ($Å^2$)', fontsize=13)


    def getEnsembleD_x(self, 
                       *args):
        if len(args) == 2:
            self.start, self.end = args[0], args[1]
        if self.start is None or self.end is None:
            sys.exit("start and end are needed to determine fitting range.")

        self.msd_x = self.msd_xyz[:,0]
        self.diffcoeff_x, self.intercept_x = np.polyfit(2*self.potim*self.nblock*
                                                        np.arange(self.end-self.start),
                                                        1E-5*self.msd_x[self.start:self.end],
                                                        deg=1)
        return self.diffcoeff_x, self.intercept_x
    
    
    def plotEnsembleMSD_x(self):
        temp = self.prefix.split('/')[-2].split('.')[-1]
        plt.plot(self.timestep, self.msd_x[:], label=f"{temp}")

        if self.diffcoeff_x is not None:
            # visualize fitting line
            plt.plot(np.arange(self.start, self.end) * self.potim * self.nblock / 1000, 
                     self.diffcoeff_x*(2*self.potim*self.nblock*1.0E5)*
                     np.arange(self.end-self.start)+self.intercept_x*1.0E5, 
                     'k:')
            # visualize fitting range
            plt.axvline(self.start* self.potim * self.nblock / 1000, 0, 1, color='black', linewidth=0.5)
            plt.axvline(self.end* self.potim * self.nblock / 1000, 0, 1, color='black', linewidth=0.5)

        plt.title("x component", fontsize=14)
        plt.legend(loc='upper left')
        plt.xlabel('t (ps)', fontsize=13)
        plt.ylabel(r'MSD ($Å^2$)', fontsize=13)


    def getEnsembleD_y(self, 
                       *args):
        if len(args) == 2:
            self.start, self.end = args[0], args[1]
        if self.start is None or self.end is None:
            sys.exit("start and end are needed to determine fitting range.")

        self.msd_y = self.msd_xyz[:,1]
        self.diffcoeff_y, self.intercept_y = np.polyfit(2*self.potim*self.nblock*
                                                        np.arange(self.end - self.start), 
                                                        1E-5*self.msd_y[self.start:self.end],
                                                        deg=1)
        return self.diffcoeff_y, self.intercept_y
    
    
    def plotEnsembleMSD_y(self):
        temp = self.prefix.split('/')[-2].split('.')[-1]
        plt.plot(self.timestep, self.msd_y[:], label=f"{temp}")

        if self.diffcoeff_y is not None:
            # visualize fitting line
            plt.plot(np.arange(self.start, self.end) * self.potim * self.nblock / 1000, 
                     self.diffcoeff_y*(2*self.potim*self.nblock*1.0E5)*
                     np.arange(self.end-self.start)+self.intercept_y*1.0E5, 
                     'k:')
            # visualize fitting range
            plt.axvline(self.start* self.potim * self.nblock / 1000, 0, 1, color='black', linewidth=0.5)
            plt.axvline(self.end* self.potim * self.nblock / 1000, 0, 1, color='black', linewidth=0.5)

        plt.title("y component", fontsize=14)
        plt.legend(loc='upper left')
        plt.xlabel('t (ps)', fontsize=13)
        plt.ylabel(r'MSD ($Å^2$)', fontsize=13)


    def getEnsembleD_z(self, 
                       *args):
        if len(args) == 2:
            self.start, self.end = args[0], args[1]
        if self.start is None or self.end is None:
            sys.exit("start and end are needed to determine fitting range.")

        self.msd_z = self.msd_xyz[:,2]
        self.diffcoeff_z, self.intercept_z = np.polyfit(2 * self.potim * self.nblock * 
                                                        np.arange(self.end - self.start), 
                                                        1E-5 * self.msd_z[self.start:self.end],
                                                        deg=1)
        return self.diffcoeff_z, self.intercept_z
    
    
    def plotEnsembleMSD_z(self):
        temp = self.prefix.split('/')[-2].split('.')[-1]
        plt.plot(self.timestep, self.msd_z[:], label=f"{temp}")

        if self.diffcoeff_z is not None:
            # visualize fitting line
            plt.plot(np.arange(self.start, self.end) * self.potim * self.nblock / 1000, 
                     self.diffcoeff_z*(2*self.potim*self.nblock*1.0E5)*
                     np.arange(self.end-self.start)+self.intercept_z*1.0E5, 
                     'k:')
            # visualize fitting range
            plt.axvline(self.start* self.potim * self.nblock / 1000, 0, 1, color='black', linewidth=0.5)
            plt.axvline(self.end* self.potim * self.nblock / 1000, 0, 1, color='black', linewidth=0.5)

        plt.title("z component", fontsize=14)
        plt.legend(loc='upper left')
        plt.xlabel('t (ps)', fontsize=13)
        plt.ylabel(r'MSD ($Å^2$)', fontsize=13)
    
    
    def saveMSD(self):
        if not os.path.isdir('msd'):
            os.mkdir('msd')
        prefix = './msd/'
        timestep = self.timestep.reshape(self.timestep.shape[0], 1)
        msd = self.msd.reshape(self.msd.shape[0], 1)
        data_save = np.concatenate((timestep, msd), axis=1)

        temp = self.prefix.split('/')[-2].split('.')[-1]
        filename = prefix+'msd_'+temp+'.txt'
        np.savetxt(filename, data_save)



class getDiffusivity:
    def __init__(self, 
                 symbol,
                 path_xdatcar='./xdatcar',
                 skip=500,
                 segment=1,
                 start=500,
                 end='auto',
                 xyz=False,
                 label='auto',
                 temp='auto'):
        """
        Arg 1: (str) symbol; target atom (ex. 'O')
        Arg 2: (str) path_xdatcar; location of xdatcar folder
        Arg 3: (int; opt) skip; initial steps to be skipped (default=500)
        Arg 4: (int; opt) segment; number of segment (default=1)
        Arg 5: (int; opt) start; start step for fitting (default=500)
        Arg 6: (int) end; end; end step for fitting (default='auto')
        Arg 7: (bool; opt) xyz; if True, x, y, z component of D is calculated. (default=False)
        Arg 8: (list; opt) temp; range of temperature (default='auto')
        Arg 9: (list; opt) label; label of XDATCAR files (default='auto')
        """
        self.symbol = symbol
        self.prefix = path_xdatcar
        self.skip = skip
        self.segment = segment
        self.start = start
        
        # get temp
        if temp == 'auto':            
            self.temp = []
            for filename in os.listdir(self.prefix):
                if len(filename.split('.')) == 2:
                    first, second = filename.split('.')
                    if first == 'xdatcar':
                        t = second.split('K')[0]
                        self.temp += [t]
            self.temp = np.sort(np.array(self.temp, dtype=int))
        else:
            self.temp = np.array(temp, dtype=int)

        # get label
        self.label = {}
        if label == 'auto':
            for t in self.temp:
                label_t = []
                path_t = os.path.join(self.prefix, f"xdatcar.{t}K")
                for filename in os.listdir(path_t):
                    if len(filename.split('_')) == 2:
                        first, second = filename.split('_')
                        if first == 'XDATCAR':
                            label_t += [second]
                self.label[t] = label_t
        else:
            for t in self.temp:
                self.label[t] = label
        
        # automately set end
        if end == 'auto':
            path_example = os.path.join(self.prefix,
                                        f"xdatcar.{self.temp[0]}K",
                                        f"XDATCAR_{self.label[self.temp[0]][0]}",)
            with open(path_example, 'r') as file:
                lines = [line.strip() for line in file]
            for line in reversed(lines):
                contents = line.split()
                if contents[0] == 'Direct':
                    nsw = int(contents[-1])
                    break
            self.end = int((nsw - self.skip)/self.segment)
        else:
            self.end = end

        self.ensembles = []
        self.getEnsembles()

        self.diffcoeffs = []
        self.getD()
        self.diffcoeffs = np.array(self.diffcoeffs)

        # get Ea and D0
        self.kb = 8.617332478E-5 # eV/K
        slop, intercept = np.polyfit(1/self.temp, np.log(self.diffcoeffs), deg=1)

        self.Ea = -slop * self.kb
        self.D0 = np.exp(intercept)

        self.plotMSD()
        self.plotArrhenius(self.diffcoeffs)
        self.saveD()

        print("")
        print("Ea = %.3f eV"%(self.Ea))
        print("D0 = %.3e m2/s"%(self.D0))

        if xyz:
            self.diffcoeffs_x = []
            self.diffcoeffs_y = []
            self.diffcoeffs_z = []
            self.getDx()
            self.getDy()
            self.getDz()
            self.plot_xyz()
            self.saveDxyz()


    def getEnsembles(self):
        for t in tqdm(self.temp, 
                      bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'% (Fore.GREEN, Fore.RESET),
                      ascii=False,
                      desc=f'{GREEN}TOTAL{RESET}'):
            ensemble = EnsembleEinstein(symbol=self.symbol,
                                        prefix=os.path.join(self.prefix,f"xdatcar.{t}K"),
                                        labels=self.label[t],
                                        segments=self.segment,
                                        skip=self.skip,
                                        start=self.start,
                                        end=self.end)
            self.ensembles += [ensemble]
            ensemble.saveMSD()


    def getD(self):
        self.diffcoeffs = [ensemble.diffcoeff for ensemble in self.ensembles]


    def getDx(self):
        self.diffcoeffs_x = [ensemble.diffcoeff_x for ensemble in self.ensembles]


    def getDy(self):
        self.diffcoeffs_y = [ensemble.diffcoeff_y for ensemble in self.ensembles]


    def getDz(self):
        self.diffcoeffs_z = [ensemble.diffcoeff_z for ensemble in self.ensembles]


    def plotMSD(self):
        for ensemble in self.ensembles:
            ensemble.plotEnsembleMSD()
            plt.legend(loc='upper left')
            plt.xlabel('t (ps)', fontsize=13)
            plt.ylabel(r'MSD ($Å^2$)', fontsize=13)
        plt.savefig('msd.png', dpi=300, transparent=False)
        plt.show()
        plt.close()


    def plotArrhenius(self, 
                      diffcoeffs, 
                      disp=True):
        kb = 8.617332478E-5 # eV/K
        slop, intercept = np.polyfit(1/self.temp, np.log(diffcoeffs), deg=1)
        tick = [r'$\frac{{1}}{{{}}}$'.format(t) for t in self.temp]
        Ea = -slop * kb
        D0 = np.exp(intercept)

        plt.plot(1/self.temp, np.log(diffcoeffs), 'bo')
        plt.plot(1/self.temp, slop*(1/self.temp)+intercept, 'k:')
        
        plt.xlabel('1/T (1/K)', fontsize=13)
        plt.xticks(1/self.temp, tick)
        plt.ylabel('ln D', fontsize=13)

        mid = int(len(self.temp)/2)
        plt.text(1/self.temp[mid]*1.05, 
                 slop*(1/self.temp[mid])+intercept, 
                 r"$E_a$=%.2f eV"%(Ea), 
                 fontsize=11)
        if disp:
            plt.savefig('Arrhenius.png', dpi=300, transparent=False)
            plt.show()
            plt.close()


    def saveD(self):
        with open('D.txt', 'w') as f:
            f.write(f"Ea = {self.Ea} eV\n")
            f.write(f"D0 = {self.D0} m2/s\n\n")
            f.write("T (K) \tD (m2/s)\n")
            for t, D in zip(self.temp, self.diffcoeffs):
                f.write(f"{t}\t{D}\n")


    def plot_xyz(self):
        figsize = (10,15)
        wspace = 0.4
        hspace = 0.5
        plt.figure(figsize=figsize)
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        
        # x
        plt.subplot(321)
        for ensemble in self.ensembles:
            ensemble.plotEnsembleMSD_x()
        
        plt.subplot(322)
        self.plotArrhenius(self.diffcoeffs_x, disp=False)

        # y
        plt.subplot(323)
        for ensemble in self.ensembles:
            ensemble.plotEnsembleMSD_y()
        
        plt.subplot(324)
        self.plotArrhenius(self.diffcoeffs_y, disp=False)

        # z
        plt.subplot(325)
        for ensemble in self.ensembles:
            ensemble.plotEnsembleMSD_z()
        
        plt.subplot(326)
        self.plotArrhenius(self.diffcoeffs_z, disp=False)

        plt.savefig('Dxyz.png', dpi=300, transparent=False)
        plt.show()
        plt.close()


    def saveDxyz(self):
        with open('Dxyz.txt', 'w') as f:
            f.write("x-component\n")
            f.write("T (K) \tD (m2/s)\n")
            for t, D in zip(self.temp, self.diffcoeffs_x):
                f.write(f"{t}\t{D}\n")
            f.write("\n")

            f.write("y-component\n")
            f.write("T (K) \tD (m2/s)\n")
            for t, D in zip(self.temp, self.diffcoeffs_y):
                f.write(f"{t}\t{D}\n")
            f.write("\n")

            f.write("z-component\n")
            f.write("T (K) \tD (m2/s)\n")
            for t, D in zip(self.temp, self.diffcoeffs_z):
                f.write(f"{t}\t{D}\n")
            f.write("\n")



class CompareNEB:
    def __init__(self):
        self.kb = 8.617332478e-5 # eV/K

        # data from Einstein relation
        self.D_ER = []

        # site information
        self.site = []
        self.site_name = []
        self.site_E = []


    def addEinstein(self, temp, D):
        if len(temp) != len(D):
            print("length of temp and D should be the same.")
            sys.exit(0)
        
        for t, d in zip(temp, D):
            dic_D = {}
            dic_D['T'] = t
            dic_D['D'] = d

            self.D_ER.append(dic_D)


    def printEinstein(self):
        D_ER_sorted = sorted(self.D_ER,
                             key=lambda x:list(x.values()))
        
        print("T(K)\tD(m2/s)")
        for dic_D in D_ER_sorted:
            print("%d"%dic_D['T'], end='\t')
            print("%.2e"%dic_D['D'])


    def addSite(self, 
                site, 
                E):
        if len(site) != len(E):
            print("length of site and E should be the same.")
            sys.exit(0)
        
        for name, energy in zip(site, E):
            if name in self.site_name:
                print(f"{name} already exsits.")
                sys.exit(0)
            dic_site = {}
            dic_site['name'] = name
            dic_site['E'] = energy
            dic_site['d'] = []
            dic_site['Ea'] = []
            dic_site['z'] = []

            self.site_name.append(name)
            self.site_E.append(energy)
            self.site.append(dic_site)


    def addPath(self, 
                init_site, 
                distance, 
                Ea, 
                z):
        num_path = len(init_site)

        for i in range(num_path):
            if not init_site[i] in self.site_name:
                print(f"{init_site[i]} is not defined.")
                sys.exit(0)
            else:
                idx = self.site_name.index(init_site[i])
            
            self.site[idx]['d'].append(distance[i])
            self.site[idx]['Ea'].append(Ea[i])
            self.site[idx]['z'].append(z[i])
    
    
    def printPath(self):
        print("site\td(Å)\tEa(eV)\tE(eV)\tz")
        for dic_site in self.site:
            num_path = len(dic_site['d'])

            for i in range(num_path):
                print(f"{dic_site['name']}", end='\t')
                print("%.3f"%dic_site['d'][i], end='\t')
                print("%.2f"%dic_site['Ea'][i], end='\t')
                print("%.2f"%dic_site['E'], end='\t')
                print("%d"%dic_site['z'][i], end='\n')


    def getD_NEB(self, 
                 T, 
                 f=1e13):
        T = np.array(T, dtype=float)
        E = np.array(self.site_E, dtype=float).reshape(-1, 1)
        
        prob = np.exp(-E/(self.kb*T))
        partition = np.sum(prob, axis=0)
        prob /= partition

        D_site = np.zeros((len(self.site), len(T)))
        for i, dic_site in enumerate(self.site):
            num_path = len(dic_site['d'])
            D_i = 0
            for j in range(num_path):
                a = dic_site['d'][j] * 1e-10
                z = dic_site['z'][j]
                D_i += (1/6)*z*a**2*f * np.exp(-dic_site['Ea'][j]/(self.kb*T))
            D_site[i] = D_i

        return np.sum(prob * D_site, axis=0)
    
    
    def plotNEB(self, 
                T, 
                f=1e13, 
                label=None):
        T = np.array(T)
        D_NEB = self.getD_NEB(T, f)
        plt.plot(1/T, np.log(D_NEB), 'k:', label=label)


    def plotER(self):
        # plot D from ER
        T_ER = np.array([dic['T'] for dic in self.D_ER])
        D_ER = np.array([dic['D'] for dic in self.D_ER])
        plt.scatter(1/T_ER, np.log(D_ER), 
                    facecolor='None', edgecolor='r', s=50, linewidth=1.5, 
                    label='Einstein')
        
        # x ticks for ER
        if len(T_ER) >= 3:
            ticks_idx = np.array([0, len(T_ER)/2, len(T_ER)-1], dtype=int)
            T_ticks = T_ER[ticks_idx]
            ticks = [f"1/{temp}" for temp in T_ticks]
            plt.xticks(1/T_ticks, ticks)
        else:
            ticks = [f"1/{temp}" for temp in T_ER]
            plt.xticks(1/T_ER, ticks)
        
        
    def plotArrhenius(self,
                      T,
                      f=[],
                      save=True,
                      outfig='Arrhenius.png',
                      dpi=300,
                      fmin=1e12,
                      fmax=1e14,
                      step=10000):
        # plot NEB
        label_NEB = None
        if len(f) == 0:
            f = [self.getFrequency(fmin, fmax, step)]
            label_NEB = "NEB " + "(" + "%.2e"%f[0] + " Hz)"

        if label_NEB is None:
            label_NEB = "NEB"

        for freq in f:
            self.plotNEB(T, freq, label=label_NEB)
            label_NEB = None
        
        # plot ER
        self.plotER()

        plt.xlabel('1/T (1/K)', fontsize=12)
        plt.ylabel('ln D', fontsize=12)
        plt.legend(fontsize=10.5)

        if save:
            plt.savefig(outfig, dpi=dpi)
        plt.show()


    def getFrequency(self,
                     fmin=1e12,
                     fmax=1e14,
                     step=10000):
        
        freq = np.linspace(fmin, fmax, step)
        T_ER = np.array([dic['T'] for dic in self.D_ER])
        D_ER = np.array([dic['D'] for dic in self.D_ER])

        err_min = np.inf
        f_opt = None
        for f in freq:
            D_NEB = self.getD_NEB(T_ER, f)
            err = np.sum((D_NEB - D_ER)**2)
            if err < err_min:
                err_min = err
                f_opt = f
                
        return f_opt
