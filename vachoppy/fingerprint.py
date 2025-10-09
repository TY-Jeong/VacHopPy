import os
import sys
import tempfile
import numpy as np
import matplotlib.pyplot as plt

from ase import Atoms
from ase.io import read
from tqdm.auto import tqdm
from joblib import Parallel, delayed
from typing import List, Tuple, Optional, Union
from itertools import combinations_with_replacement

from vachoppy.utils import Snapshots


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
        path_structure (str): 
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
                 path_structure: str,
                 Rmax: float = 10.0,
                 delta: float = 0.08,
                 sigma: float = 0.03,
                 dirac: str = 'g',
                 verbose: bool = True):
        
        self.A = A
        self.B = B
        self.path_structure = path_structure
        self.Rmax = Rmax
        self.delta = delta
        self.sigma = sigma
        self.dirac = dirac
        self.verbose = verbose
        
        self.R = np.arange(0, self.Rmax, self.delta)
        self.structure = self._read_structure(path_structure)
        self.lattice = self.structure.get_cell()
        self.volume = self.structure.get_volume()
        
        symbols = np.array(self.structure.get_chemical_symbols())
        self.indices_A = np.where(symbols == self.A)[0]
        self.indices_B = np.where(symbols == self.B)[0]
        
        if self.indices_A.size == 0: 
            raise ValueError(f"Atom type '{A}' not found in structure.")
        if self.indices_B.size == 0: 
            raise ValueError(f"Atom type '{B}' not found in structure.")

        self.num_A = len(self.indices_A)
        self.num_B = len(self.indices_B)

        self.fingerprint: Optional[np.ndarray] = None

    def _read_structure(self, file_path: str) -> 'Atoms':
        """"Reads a structure file using ASE and returns an Atoms object."""
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Error: Input structure file '{file_path}' not found.")
        try:
            return read(file_path)
        except Exception as e:
            raise IOError(f"Failed to read '{file_path}' with ASE. Error: {e}")

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
        
        all_scaled_positions = self.structure.get_scaled_positions()
        coords_frac = all_scaled_positions[indices]
        
        extended_coords = coords_frac[:, np.newaxis, :] + shifts[np.newaxis, :, :]
        return extended_coords.reshape(-1, 3)

    def _calculate_fingerprint_for_atom_A(self, index_A: int, extended_coords_B: np.ndarray) -> np.ndarray:
        """Calculates the fingerprint for a single central atom A."""
        coord_A = self.structure.get_scaled_positions()[index_A]

        disp = extended_coords_B - coord_A
        disp_cart = np.dot(disp, self.lattice)
        R_ij = np.linalg.norm(disp_cart, axis=1)
        
        if self.A == self.B:
            R_ij[R_ij < 1e-4] = np.inf

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
        print(f"  - Structure File : {self.path_structure}")
        print(f"  - Central Atom A : {self.A} (Found {self.num_A})")
        print(f"  - Neighbor Atom B: {self.B} (Found {self.num_B})")
        print(f"  - Rmax           : {self.Rmax} Ang")
        print(f"  - Delta          : {self.delta} Ang")
        print(f"  - Sigma          : {self.sigma} Ang")
        print(f"  - Dirac Type     : {'Gaussian' if self.dirac == 'g' else 'Square'}")
        print("="*50 + "\n")

    def plot_fingerprint(self, 
                         title: Optional[str] = None,
                         disp: bool = True,
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
            if filename is None: filename = f"FP_{self.A}-{self.B}.png"
            plt.savefig(filename, dpi=dpi)
        if disp: plt.show()
        plt.close(fig)
        
        
def cosine_distance(fp1: np.ndarray, fp2: np.ndarray) -> float:
    """Calculates a scaled cosine distance between two fingerprint vectors.

    This function computes the cosine similarity between two vectors and
    transforms it into a distance metric ranging from 0 to 1. A distance of 0
    means the vectors are identical, while 1 means they are opposite.

    Note:
        The formula used is `0.5 * (1 - cos_similarity)`, where `cos_similarity`
        is the dot product of the unit vectors. This scales the standard
        cosine distance `(1 - cos_similarity)` to the [0, 1] range.

    Args:
        fp1 (np.ndarray): 
            The first fingerprint vector (1D NumPy array).
        fp2 (np.ndarray): 
            The second fingerprint vector (1D NumPy array).

    Returns:
        float: The scaled cosine distance between the two vectors.
    """
    # Numerically stable dot product of unit vectors (cosine similarity)
    fp1_norm = np.linalg.norm(fp1)
    fp2_norm = np.linalg.norm(fp2)
    
    # Add a small epsilon to prevent division by zero if a vector has zero length
    epsilon = 1e-9
    
    similarity = np.dot(fp1, fp2) / (fp1_norm * fp2_norm + epsilon)
    
    return 0.5 * (1.0 - similarity)


def get_fingerprint(path_structure: str, 
                    filename: str, 
                    atom_pairs: List[Tuple[str, str]], 
                    Rmax: float = 10.0,
                    delta: float = 0.08,
                    sigma: float = 0.03,
                    dirac: str = 'g',
                    disp: bool = True,
                    verbose: bool = True) -> np.ndarray:
    """
    Calculates and concatenates fingerprints for multiple atom pairs, saves them to
    a file, and optionally plots the results.

    This function iterates through a list of specified atom pairs (A, B),
    calculates the atomic fingerprint for each using the FingerPrint class,
    and concatenates them into a single 1D array.

    The final result is saved to a two-column text file. The first column is a
    composite distance axis, where each pair's distance range [0, Rmax) is
    shifted by multiples of Rmax. The second column is the concatenated
    fingerprint data. The optional plot displays each pair's fingerprint
    concatenated along the x-axis.

    Args:
        path_structure (str):
            Path to the crystallographic structure file (e.g., POSCAR, cif).
        filename (str):
            The name of the output file to save the fingerprint data.
        atom_pairs (List[Tuple[str, str]]):
            A list of tuples, where each tuple contains the chemical symbols
            for an atom pair, e.g., [('Ti', 'O'), ('O', 'O')].
        Rmax (float):
            Cutoff radius in Angstroms for the fingerprint calculation.
        delta (float):
            Discretization step for the distance axis (r).
        sigma (float):
            Gaussian broadening width for interatomic distances.
        disp (bool, optional): 
            If True, displays a plot of the fingerprints for all atom pairs. 
            Defaults to True.
        verbose (bool, optional):
            Verbosity flag. Defaults to True
        
    Returns:
        np.ndarray: A single 1D NumPy array containing the concatenated
                    fingerprints of all specified atom pairs.
    """
    all_fingerprints = []
    all_R_values = []

    if verbose: print(f"Calculating fingerprints for {len(atom_pairs)} pairs...")
    for i, (A, B) in enumerate(atom_pairs):
        fp_instance = FingerPrint(A, B, path_structure, Rmax, delta, sigma, dirac, verbose=False)
        fp_instance.calculate()
        
        all_fingerprints.append(fp_instance.fingerprint)
        shifted_R = fp_instance.R + i * Rmax
        all_R_values.append(shifted_R)

    final_fingerprint = np.concatenate(all_fingerprints)
    final_R = np.concatenate(all_R_values)
    
    with open(filename, 'w') as f:
        f.write(f'# Rmax, delta, sigma, dirac = {Rmax}, {delta}, {sigma}, {dirac}\n')
        f.write('# pair : ')
        f.write(', '.join([f'{A}-{B}' for A, B in atom_pairs]))
        f.write('\n')
        
        data_to_save = np.vstack((final_R, final_fingerprint)).T
        np.savetxt(f, data_to_save, fmt='%.6f', delimiter='\t')

    if verbose: print(f"Fingerprint data successfully saved to '{filename}'")
    
    if disp:
        fig, ax = plt.subplots(figsize=(10, 5))
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.2)
        
        for i, (A, B) in enumerate(atom_pairs):
            label = f"{A}-{B}"
            ax.plot(all_R_values[i], all_fingerprints[i], label=label)

        ax.axhline(0, 0, 1, color='k', linestyle='--', linewidth=1)
        
        for i in range(1, len(atom_pairs)):
            ax.axvline(x=i * Rmax, color='gray', linestyle=':', linewidth=1.2)
        
        ax.set_xlabel("Distance (Å)", fontsize=13)
        ax.set_ylabel('g(r) - 1', fontsize=13)
        ax.legend(fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.6)
        
        fig.tight_layout()
        plt.show()
        plt.close(fig)

    return final_fingerprint


def _worker_calculate_distance(args: tuple) -> List[float]:
    """[Parallel Worker] Calculates a fingerprint and its cosine distance to a reference."""
    snapshot_index, snapshot_path, ref_fingerprint, atom_pairs, Rmax, delta, sigma, dirac, output_dir = args
    
    snapshot_fingerprint = get_fingerprint(
        path_structure=snapshot_path,
        filename=os.path.join(output_dir, f"fingerprint_{snapshot_index:04d}.txt"),
        atom_pairs=atom_pairs,
        Rmax=Rmax, delta=delta, sigma=sigma, dirac=dirac,
        disp=False,
        verbose=False
    )

    dist = cosine_distance(ref_fingerprint, snapshot_fingerprint)
    
    return [snapshot_index, dist]

def plot_cosine_distance(path_traj: Union[str, List[str]],
                         t_interval: float,
                         reference_structure: str,
                         atom_pairs: Optional[List[Tuple[str, str]]] = None,
                         Rmax: float = 10.0,
                         delta: float = 0.08,
                         sigma: float = 0.03,
                         dirac: str = 'g',
                         prefix: str = 'cosine_distance_trace',
                         dpi: int = 300,
                         path_dir: str = 'fingerprint_trace',
                         n_jobs: int = -1,
                         find_fluctuations: bool = True,
                         window_size: int = 50,
                         threshold_std: float = None,
                         verbose: bool = True) -> None:
    """
    Traces the change in atomic fingerprint over time against a reference structure.

    This function internally generates snapshots from the given trajectory files,
    then calculates the cosine distance between each snapshot's fingerprint and
    a reference structure's fingerprint. The process is parallelized using joblib.
    Results (cosine distance vs. time) are saved to a text file and plotted.

    Args:
        path_traj (Union[str, List[str]]):
            Path to a single HDF5 trajectory file or a list of paths.
        t_interval (float):
            The time interval in picoseconds (ps) for averaging snapshots.
        reference_structure (str):
            Path to the reference structure file (e.g., POSCAR of the initial phase).
        atom_pairs (Optional[List[Tuple[str, str]]], optional):
            A list of atom pairs to analyze. If None, all unique pair
            combinations from the reference_structure will be generated
            automatically. Defaults to None.
        Rmax (float, optional): 
            Cutoff radius for the fingerprint. Defaults to 10.0.
        delta (float, optional): 
            Discretization step for the fingerprint. Defaults to 0.08.
        sigma (float, optional): 
            Gaussian broadening for the fingerprint. Defaults to 0.03.
        dirac (str, optional): 
            Dirac function type ('g' or 's'). Defaults to 'g'.
        prefix (str, optional):
            A prefix for the output plot and data files (e.g., 'cosine_distance'). 
            Defaults to 'cosine_distance_trace'.
        dpi (int, optional):
            The resolution (dots per inch) for the saved plot figure.
            Defaults to 300.
        path_dir (str, optional): 
            Directory to save final output files. Defaults to 'fingerprint_trace'.
        n_jobs (int, optional): 
            Number of CPU cores for parallel processing. -1 uses all. Defaults to -1.
        find_fluctuations (bool, optional): 
            If True, performs an analysis to find intervals where the cosine distance 
            deviates significantly from the mean. Defaults to True.
        window_size (int, optional): 
            The number of data points to include in the moving average window. 
            A larger window provides more smoothing. Defaults to 50.
        threshold_std (float, optional): 
            The number of standard deviations (sigma) from the global mean 
            to define the fluctuation threshold. If `None`, the threshold line
            is not drawn and fluctuation intervals are not detected.
            Common values are:
                - 1.0 (flags ~31.8% of data as potential outliers)
                - 2.0 (flags ~4.6% of data as outliers)
                - 3.0 (flags ~0.3% of data as extreme outliers)
            Defaults to None.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.
    """
    if not os.path.isdir(path_dir):
        os.makedirs(path_dir)
        if verbose: print(f"Created output directory: '{path_dir}'")

    if atom_pairs is None:
        if verbose: print("Argument 'atom_pairs' not provided: Auto-generating all unique pairs...")
        try:
            atoms = read(reference_structure)
            atom_species = sorted(list(set(atoms.get_chemical_symbols())))
            atom_pairs = list(combinations_with_replacement(atom_species, 2))
            if verbose: print(f"-> Generated pairs: {atom_pairs}\n")
        except Exception as e:
            raise IOError(f"Failed to read reference structure '{reference_structure}' to auto-generate pairs. Error: {e}")
    with tempfile.TemporaryDirectory() as temp_dir:
        snapshots = Snapshots(
            path_traj=path_traj,
            t_interval=t_interval,
            verbose=False
        )
        snapshots.save_snapshots(
            path_dir=temp_dir,
            format='vasp',
            prefix='POSCAR'
        )
        ref_fingerprint = get_fingerprint(
            path_structure=reference_structure,
            filename=os.path.join(path_dir, "fingerprint_ref.txt"),
            atom_pairs=atom_pairs,
            Rmax=Rmax, 
            delta=delta, 
            sigma=sigma,
            dirac=dirac,
            disp=False,
            verbose=False
        )
        
        tasks = []
        for i in range(snapshots.num_steps):
            snapshot_path = os.path.join(temp_dir, f"POSCAR_{i:0{snapshots.digit}d}")
            tasks.append((i, snapshot_path, ref_fingerprint, atom_pairs, Rmax, delta, sigma, dirac, path_dir))
        results = Parallel(n_jobs=n_jobs, verbose=0)(
            delayed(_worker_calculate_distance)(task)
            for task in tqdm(tasks,
                             desc=f'Compute Fingerprint',
                             bar_format='{l_bar}{bar:30}{r_bar}',
                             ascii=True)
        )
        
    results.sort(key=lambda x: x[0])
    results = np.array(results, dtype=np.float64)
    results[:, 0] = (results[:, 0] + 1) * snapshots.t_interval
    time_data, distance_data = results[:, 0], results[:, 1]
    
    fluctuation_intervals = []
    if find_fluctuations:
        global_mean = np.mean(distance_data)
        global_std = np.std(distance_data)
        if threshold_std is not None:
            upper_threshold = global_mean + threshold_std * global_std
        
        moving_avg = np.convolve(distance_data, np.ones(window_size)/window_size, mode='valid')
        time_avg = time_data[(window_size-1)//2 : -(window_size-1)//2][:len(moving_avg)]
        if threshold_std is not None:
            suspicious_indices = np.where(moving_avg > upper_threshold)[0]
            if suspicious_indices.size > 0:
                for group in np.split(suspicious_indices, np.where(np.diff(suspicious_indices) != 1)[0] + 1):
                    start_time = time_avg[group[0]]
                    end_time = time_avg[group[-1]]
                    fluctuation_intervals.append((start_time, end_time))
                
    plt.figure(figsize=(12, 5))
    ax = plt.gca()
    ax.scatter(time_data, distance_data, alpha=0.5, label='Cosine Distance', s=15, zorder=2)
    
    if find_fluctuations:
        ax.plot(time_avg, moving_avg, color='crimson', lw=2, label=f'Moving Avg (win={window_size})', zorder=3)
        ax.axhline(global_mean, color='k', linestyle='--', label='Global Mean', zorder=1)
        if threshold_std is not None:
            ax.axhline(upper_threshold, color='k', linestyle=':', label=f'Threshold ({threshold_std}σ)', zorder=1)
        for start, end in fluctuation_intervals:
            ax.axvspan(start, end, color='orange', alpha=0.3, label='Detected Fluctuation')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), fontsize=9)
    
    out_png = f'{prefix}.png'
    ax.set_xlabel("Time (ps)", fontsize=13)
    ax.set_ylabel('Cosine Distance to Reference', fontsize=13)
    ax.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)
    plt.show()
    plt.close()
    print(f"\n'{out_png}' created successfully.")
    
    output_txt = f'{prefix}.txt'
    with open(output_txt, 'w') as f:
        f.write(f'# Rmax, delta, sigma, dirac = {Rmax}, {delta}, {sigma}, {dirac}\n')
        f.write('# pair : ' + ', '.join([f'{A}-{B}' for A, B in atom_pairs]) + '\n')
        np.savetxt(f, results, fmt='%.6f', delimiter='\t')
    print(f"'{output_txt}' created successfully.")
    
    if (threshold_std is not None) and find_fluctuations and fluctuation_intervals:
        print("\nDetected fluctuation intervals:")
        for start, end in fluctuation_intervals:
            print(f"  - From {start:.2f} ps to {end:.2f} ps")