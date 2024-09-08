import os
import sys
import numpy as np
import matplotlib.pyplot as plt


def CosineDistance(fp1, fp2):
    """
    fp1 : fingerprint of structure 1
    fp2 : fingerprint of structure 2
    """
    dot = np.dot(fp1, fp2)
    norm1 = np.linalg.norm(fp1, ord=2)
    norm2 = np.linalg.norm(fp2, ord=2)

    return 0.5 * (1 - dot/(norm1*norm2))



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
        self.readPOSCAR()
        
        # Search indices of A and B
        self.idx_A = self.searchAtomInex(self.A)
        self.idx_B = self.searchAtomInex(self.B)
        
        # volume of unit cell
        self.V_unit = np.abs(
            np.dot(np.cross(self.lattice[0], self.lattice[1]), 
                   self.lattice[2]))
        
        # extended coordinaitons of B
        self.B_ext = None
        self.getExtendedCoords_B()

        # calculate fingerprint
        self.getFingerprint()
    
    
    def readPOSCAR(self):
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
            
            
    def searchAtomInex(self, atom_name):
        for i, atom in enumerate(self.position):
            if atom['name'] == atom_name:
                return i
        raise ValueError(f"Atom {atom_name} not found in POSCAR file.")
    
    
    def getExtendedCoords_B(self):
        # supercells within Rmax
        l_lat = np.linalg.norm(self.lattice, axis=1)
        m = np.floor(self.Rmax / l_lat) + 1
        
        # make 3D grid of supercells
        mx, my, mz = [np.arange(-mi, mi+1) for mi in m]
        shifts = np.array(np.meshgrid(mx, my, mz)).T.reshape(-1, 3)

        # save candidated of B atoms
        coords_B = self.position[self.idx_B]['coords']
        self.B_ext = np.vstack([shifts + coord for coord in coords_B])
             
             
    def getFingerprint(self):
        for coord_A in self.position[self.idx_A]['coords']:
            self.fingerprint += self.getFingerprint_i(coord_A)
        
        self.fingerprint *= self.V_unit / self.position[self.idx_A]['num']
        self.fingerprint -= 1
        
        
    def getFingerprint_i(self, 
                         coord_A_i):
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
    
    
    def plotFingerprint(self, 
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