import os
import numpy as np
import matplotlib.pyplot as plt

from typing import Optional

from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor

class FingerPrint:
    """
    Calculates the atomic environment fingerprint between two atom types.

    This class reads a crystallographic structure file, identifies two specified
    atom types (A and B), and computes the partial pair correlation function g(r)
    between them. The calculation is vectorized using NumPy for performance.

    The primary workflow is to initialize the class and then call the .calculate()
    method to perform the computation. Results can be visualized using .plot_fingerprint().

    Args:
        A (str): 
            Chemical symbol of the central atom type.
        B (str): 
            Chemical symbol of the neighboring atom type.
        structure_file (str): 
            Path to the crystallographic structure file (e.g., POSCAR, cif).
        Rmax (float, optional): 
            Cutoff radius in Angstroms for the calculation. Defaults to 10.
        delta (float, optional): 
            Discretization step for the distance axis (r). Defaults to 0.08.
        sigma (float, optional): 
            Gaussian broadening width for interatomic distances. Defaults to 0.03.
        dirac (str, optional): 
            Type of Dirac delta function approximation. 'g' for Gaussian or 's' for a square function. 
            Defaults to 'g'
        verbose (bool, optional): 
            Verbosity flag. Defaults to True.

    Attributes:
        fingerprint (np.ndarray): 
            The calculated fingerprint g(r) - 1 as a 1D NumPy array.
        R (np.ndarray): 
            The distance values (r) for which the fingerprint is calculated.
    """
    def __init__(self,
                 A: str,
                 B: str,
                 structure_file: str,
                 Rmax: float = 10.0,
                 delta: float = 0.08,
                 sigma: float = 0.03,
                 dirac: str = 'g',
                 verbose: bool = True):
        
        self.A = A
        self.B = B
        self.structure_file = structure_file
        self.Rmax = Rmax
        self.delta = delta
        self.sigma = sigma
        self.dirac = dirac
        self.verbose = verbose
        
        # --- Internal attributes ---
        self.R = np.arange(0, self.Rmax, self.delta)
        self.structure = self._read_structure(structure_file)
        self.lattice = self.structure.lattice.matrix
        self.volume = self.structure.volume
        
        # Get indices of atom types A and B
        self.indices_A = [i for i, site in enumerate(self.structure) if site.specie.symbol == self.A]
        self.indices_B = [i for i, site in enumerate(self.structure) if site.specie.symbol == self.B]
        if not self.indices_A: raise ValueError(f"Atom type '{A}' not found in structure.")
        if not self.indices_B: raise ValueError(f"Atom type '{B}' not found in structure.")
            
        self.num_A = len(self.indices_A)
        self.num_B = len(self.indices_B)

        # --- Result attribute ---
        self.fingerprint: Optional[np.ndarray] = None

    def _read_structure(self, file_path: str):
        """Reads a structure file using ASE and converts it to a Pymatgen object."""
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Error: Input structure file '{file_path}' not found.")
        try:
            atoms = read(file_path)
            return AseAtomsAdaptor.get_structure(atoms)
        except Exception as e:
            raise IOError(f"Failed to read or convert '{file_path}'. ASE/Pymatgen error: {e}")

    def _gaussian_func(self, x: np.ndarray) -> np.ndarray:
        return (1 / np.sqrt(2 * np.pi * self.sigma**2)) * np.exp(-x**2 / (2 * self.sigma**2))

    def _square_func(self, x: np.ndarray) -> np.ndarray:
        return (np.abs(x) <= self.sigma) / (2 * self.sigma)

    def _dirac_func(self, x: np.ndarray) -> np.ndarray:
        if self.dirac.lower().startswith('g'):
            return self._gaussian_func(x)
        elif self.dirac.lower().startswith('s'):
            return self._square_func(x)
        else:
            raise ValueError(f"Dirac function type '{self.dirac}' is not defined. Use 'g' or 's'.")

    def _get_extended_coords(self, indices: list) -> np.ndarray:
        """Creates a supercell of atom coordinates to handle periodic boundaries."""
        l_lat = np.linalg.norm(self.lattice, axis=1)
        m = np.floor(self.Rmax / l_lat) + 1
        
        mx, my, mz = (np.arange(-mi, mi + 1) for mi in m)
        shifts = np.array(np.meshgrid(mx, my, mz)).T.reshape(-1, 3)
        
        coords_frac = np.array([self.structure[i].frac_coords for i in indices])
        
        # Use broadcasting to create all periodic images
        extended_coords = coords_frac[:, np.newaxis, :] + shifts[np.newaxis, :, :]
        return extended_coords.reshape(-1, 3)

    def _calculate_fingerprint_for_atom_A(self, index_A: int, extended_coords_B: np.ndarray) -> np.ndarray:
        """Calculates the fingerprint for a single central atom A."""
        coord_A = self.structure[index_A].frac_coords
        
        disp = extended_coords_B - coord_A
        disp_cart = np.dot(disp, self.lattice)
        R_ij = np.linalg.norm(disp_cart, axis=1)
        
        # If A and B are the same element, exclude the self-distance (R_ij=0)
        if self.A == self.B:
            R_ij[R_ij < 1e-4] = np.inf
        
        # Number density of B atoms
        rho_B = self.num_B / self.volume
        
        diff_matrix = self.R[:, np.newaxis] - R_ij[np.newaxis, :]
        dirac_matrix = self._dirac_func(diff_matrix)
        fingerprint_i = np.sum(dirac_matrix / (4 * np.pi * rho_B * R_ij**2 + 1e-12), axis=1)
        
        return fingerprint_i
        
    def calculate(self) -> None:
        """
        Runs the main fingerprint calculation.

        This method computes the extended coordinates for neighbor atoms and then
        iterates through each central atom to calculate its partial fingerprint,
        averaging the results. The final g(r) - 1 is stored in `self.fingerprint`.
        """
        extended_coords_B = self._get_extended_coords(self.indices_B)
        
        total_fingerprint = np.zeros_like(self.R)
        for idx_A in self.indices_A:
            total_fingerprint += self._calculate_fingerprint_for_atom_A(idx_A, extended_coords_B)
        
        self.fingerprint = (total_fingerprint / self.num_A) - 1
        
        if self.verbose:
            self.summary()

    def summary(self):
        """Prints a summary of the fingerprint analysis settings."""
        print("\n" + "="*50)
        print(f"           Atomic Fingerprint Summary")
        print("="*50)
        print(f"  - Structure File : {self.structure_file}")
        print(f"  - Central Atom A : {self.A} (Found {self.num_A})")
        print(f"  - Neighbor Atom B: {self.B} (Found {self.num_B})")
        print(f"  - Rmax           : {self.Rmax} Ang")
        print(f"  - Delta          : {self.delta} Ang")
        print(f"  - Sigma          : {self.sigma} Ang")
        print(f"  - Dirac Type     : {'Gaussian' if self.dirac == 'g' else 'Square'}")
        print("="*50 + "\n")

    def plot_fingerprint(self, 
                         title: Optional[str] = None,
                         save: bool = True,
                         filename: str = None,
                         dpi: int = 300) -> None:
        """ Plots the calculated fingerprint g(r) - 1.
        
        This method generates a 2D plot of the atomic fingerprint, showing
        g(r) - 1 as a function of distance (r). The plot is always displayed,
        and can optionally be saved to a file.

        Args:
            title (str, optional): 
                A custom title for the plot. If None, no title is set.
                Defaults to None.
            save (bool, optional): 
                If True, saves the plot to a file. Defaults to True.
            filename (str, optional): 
                The filename for the saved plot. If None, a default filename
                is automatically generated based on the atom types
                (e.g., 'FP_A-B.png'). Defaults to None.
            dpi (int, optional): 
                The resolution (dots per inch) for the saved figure.
                Defaults to 300.
        """

        if self.fingerprint is None:
            raise RuntimeError("Fingerprint has not been calculated. Please run the .calculate() method first.")

        fig, ax = plt.subplots(figsize=(8, 5))
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.2)

        ax.plot(self.R, self.fingerprint, label=f"{self.A}-{self.B}")
        ax.axhline(0, 0, 1, color='k', linestyle='--', linewidth=1)
        
        ax.set_xlabel("Distance (Å)", fontsize=13)
        ax.set_ylabel('g(r) - 1', fontsize=13)
        ax.set_title(title)
        ax.legend(fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.6)
        
        fig.tight_layout()

        if save:
            if filename is None:
                filename = f"FP_{self.A}-{self.B}.png"
            plt.savefig(filename, dpi=dpi)
            
        plt.show()
        plt.close(fig)