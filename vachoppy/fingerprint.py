"""
vachoppy.fingerprint
====================

Provides tools for calculating and analyzing atomic environment fingerprints,
which are based on the partial pair correlation function, g(r).

This module offers a powerful way to characterize crystal structures and track
their evolution over time. A "fingerprint" serves as a quantitative measure of
the local atomic environment. The module supports two primary workflows:

1.  **Static Analysis**: Calculating the fingerprint for a single, static
    crystal structure to characterize its atomic arrangement.
2.  **Dynamic Analysis**: Tracking how a system's fingerprint changes over the
    course of a molecular dynamics trajectory by comparing it to a reference
    structure.

Main Components
---------------
- **FingerPrint**: The core class for calculating the g(r) fingerprint for a
  single pair of atom types in a given crystal structure.
- **get_fingerprint**: A convenience function that uses `FingerPrint` to
  calculate and concatenate fingerprints for multiple atom pairs.
- **cosine_distance**: A utility function to compute a scaled distance metric
  (from 0 to 1) between two fingerprint vectors.
- **plot_cosine_distance**: A high-level analysis function that takes an MD
  trajectory, generates snapshots, and plots the cosine distance of each
  snapshot's fingerprint to a reference, revealing structural changes over time.

Typical Usage
-------------
**1. Static Fingerprint Calculation:**

.. code-block:: python

    from vachoppy.fingerprint import get_fingerprint

    # Define atom pairs to analyze
    pairs = [('Hf', 'Hf'), ('Hf', 'O'), ('O', 'O')]

    # Calculate and save the combined fingerprint for a structure
    get_fingerprint(
        path_structure='path/to/hfo2.cif',
        filename='hfo2_fingerprint.txt',
        atom_pairs=pairs,
        disp=True  # Also display a plot
    )

**2. Dynamic Trajectory Analysis:**

.. code-block:: python

    from vachoppy.fingerprint import plot_cosine_distance

    # Trace the cosine distance of a trajectory relative to a reference
    plot_cosine_distance(
        path_traj='path/to/trajectory.h5',
        t_interval=0.1,  # Generate a snapshot every 0.1 ps
        reference_structure='path/to/initial_structure.cif',
        # Atom pairs will be auto-generated if not specified
    )
"""

from __future__ import annotations

__all__ =['FingerPrint', 'cosine_distance', 'get_fingerprint', 'plot_cosine_distance']

import os
import tempfile
import numpy as np
import matplotlib.pyplot as plt

from ase import Atoms
from ase.io import read
from tqdm.auto import tqdm
from joblib import Parallel, delayed
from typing import List, Tuple, Optional
from itertools import combinations_with_replacement

from vachoppy.utils import Snapshots


class FingerPrint:
    """Calculates the atomic environment fingerprint between two atom types.

    This class computes the partial pair correlation function, g(r), between two
    specified atom types (A and B) within a crystal structure. It serves as a
    "fingerprint" to characterize the local atomic environment. The calculation
    is vectorized using NumPy for performance.

    The primary workflow is to initialize the class, which automatically triggers
    the calculation. Results can then be visualized using `.plot_fingerprint()`.

    Args:
        A (str):
            Chemical symbol of the central atom type.
        B (str):
            Chemical symbol of the neighboring atom type.
        path_structure (str):
            Path to the crystallographic structure file (e.g., POSCAR, cif).
        Rmax (float, optional):
            Cutoff radius in Angstroms for the calculation. Defaults to 10.0.
        delta (float, optional):
            Discretization step (bin size) for the distance axis (r).
            Defaults to 0.08.
        sigma (float, optional):
            Gaussian broadening width applied to interatomic distances.
            Defaults to 0.03.
        dirac (str, optional):
            Type of Dirac delta function approximation: 'g' for Gaussian or
            's' for a square function. Defaults to 'g'.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Attributes:
        fingerprint (numpy.ndarray):
            The calculated fingerprint, g(r) - 1, as a 1D NumPy array.
        R (numpy.ndarray):
            The distance values (r) in Angstroms for which the fingerprint is calculated.
        num_A (int):
            The number of atoms of type A found in the structure.
        num_B (int):
            The number of atoms of type B found in the structure.

    Raises:
        FileNotFoundError:
            If the specified structure file does not exist.
        IOError:
            If the structure file cannot be read by ASE.
        ValueError:
            If atom types A or B are not found in the structure, or if an
            invalid `dirac` type is specified.

    Examples:
        >>> fp = FingerPrint(
        ...     A='Ti',
        ...     B='O',
        ...     path_structure='POSCAR',
        ...     rmax=10.0,
        ...     verbose=True
        ... )
        >>> fp.plot_fingerprint()
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
        """Runs the main fingerprint calculation.

        This method computes the extended coordinates for neighbor atoms to handle
        periodic boundaries, then iterates through each central atom to calculate
        its partial fingerprint, and finally averages the results. The final
        g(r) - 1 is stored in the `self.fingerprint` attribute.
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
                         title: str | None = None,
                         disp: bool = True,
                         save: bool = True,
                         filename: str | None = None,
                         dpi: int = 300) -> None:
        """Plots the calculated fingerprint, g(r) - 1.

        This method generates a 2D plot of the atomic fingerprint, showing
        g(r) - 1 as a function of distance (r). The plot can be displayed
        interactively and optionally saved to a file.

        Args:
            title (str | None, optional):
                A custom title for the plot. If None, no title is set.
                Defaults to None.
            disp (bool, optional):
                If True, displays the plot interactively (`plt.show()`).
                Defaults to True.
            save (bool, optional):
                If True, saves the plot to a file. Defaults to True.
            filename (str | None, optional):
                The filename for the saved plot. If None, a default filename
                is automatically generated (e.g., 'FP_A-B.png').
                Defaults to None.
            dpi (int, optional):
                The resolution (dots per inch) for the saved figure.
                Defaults to 300.

        Returns:
            None: This method does not return any value.

        Raises:
            RuntimeError: If the fingerprint has not been calculated. Call the
                `.calculate()` method first.

        Examples:
            >>> fp = FingerPrint(A='Ti', B='O', path_structure='POSCAR')
            >>> fp.calculate()
            >>> # Display the plot and save it with a default name
            >>> fp.plot_fingerprint()
            
            >>> # Save the plot with a custom name without displaying it
            >>> fp.plot_fingerprint(
            ...     title="Ti-O Fingerprint",
            ...     disp=False,
            ...     save=True,
            ...     filename="tio2_fp.png"
            ... )
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
    transforms it into a distance metric scaled to the range [0, 1]. A distance
    of 0 indicates identical vectors, while 1 indicates opposite vectors.

    Note:
        The formula used is `0.5 * (1 - cos_similarity)`, where `cos_similarity`
        is the dot product of the unit vectors.

    Args:
        fp1 (np.ndarray):
            The first fingerprint vector (1D NumPy array).
        fp2 (np.ndarray):
            The second fingerprint vector (1D NumPy array).

    Returns:
        float:
            The scaled cosine distance, a value between 0.0 and 1.0.

    Raises:
        ValueError:
            If the input arrays have mismatched shapes, which would prevent
            the dot product calculation.

    Examples:
        >>> import numpy as np
        >>> vec1 = np.array([1.0, 0.0, 0.0])
        >>> vec2 = np.array([0.0, 1.0, 0.0])
        >>> # Identical vectors
        >>> cosine_distance(vec1, vec1)
        0.0
        
        >>> # Orthogonal vectors
        >>> cosine_distance(vec1, vec2)
        0.5

        >>> # Opposite vectors
        >>> cosine_distance(vec1, -vec1)
        1.0
    """
    fp1_norm = np.linalg.norm(fp1)
    fp2_norm = np.linalg.norm(fp2)
    
    # Add a small epsilon to prevent division by zero
    epsilon = 1e-9
    
    similarity = np.dot(fp1, fp2) / (fp1_norm * fp2_norm + epsilon)
    
    return 0.5 * (1.0 - similarity)


def get_fingerprint(path_structure: str, 
                    filename: str | None, 
                    atom_pairs: List[Tuple[str, str]], 
                    Rmax: float = 10.0,
                    delta: float = 0.08,
                    sigma: float = 0.03,
                    dirac: str = 'g',
                    disp: bool = True,
                    verbose: bool = True) -> np.ndarray:
    """Calculates, concatenates, and saves fingerprints for multiple atom pairs.

    This function serves as a convenient wrapper around the `FingerPrint` class.
    It iterates through a list of specified atom pairs (A, B), calculates the
    fingerprint for each, and concatenates them into a single 1D array.

    The final result is saved to a two-column text file, where the first column
    is a composite distance axis and the second is the concatenated fingerprint
    data. An optional plot visualizes all fingerprints sequentially.

    Args:
        path_structure (str):
            Path to the crystallographic structure file (e.g., POSCAR, cif).
        filename (str | None):
            The name of the output file to save the concatenated fingerprint data.
            If None, output file is not saved.
        atom_pairs (List[Tuple[str, str]]):
            A list of tuples, each containing the chemical symbols for an atom
            pair, e.g., `[('Ti', 'O'), ('O', 'O')]`.
        Rmax (float, optional):
            Cutoff radius in Angstroms for the calculation. Defaults to 10.0.
        delta (float, optional):
            Discretization step (bin size) for the distance axis (r).
            Defaults to 0.08.
        sigma (float, optional):
            Gaussian broadening width for interatomic distances. Defaults to 0.03.
        dirac (str, optional):
            Type of Dirac delta function approximation: 'g' for Gaussian or
            's' for a square function. Defaults to 'g'.
        disp (bool, optional):
            If True, displays a plot of the concatenated fingerprints.
            Defaults to True.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Returns:
        np.ndarray:
            A single 1D NumPy array containing the concatenated fingerprints of
            all specified atom pairs.

    Raises:
        FileNotFoundError:
            If the `path_structure` file does not exist.
        ValueError:
            If an atom type specified in `atom_pairs` is not found in the
            structure file.

    Examples:
        >>> pairs = [('Ti', 'Ti'), ('Ti', 'O'), ('O', 'O')]
        >>> combined_fp = get_fingerprint(
        ...     path_structure='POSCAR',
        ...     filename='tio2_full_fp.txt',
        ...     atom_pairs=pairs,
        ...     Rmax=10.0,
        ...     disp=True
        ... )
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
    
    if filename is not None:
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


def plot_cosine_distance(path_traj: str | list[str],
                         t_interval: float,
                         reference_structure: str,
                         atom_pairs: list[tuple[str, str]] | None = None,
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
                         threshold_std: float | None = None,
                         disp: bool = True,
                         verbose: bool = True) -> None:
    """Traces structural evolution by plotting fingerprint cosine distance over time.

    This function provides a comprehensive workflow to analyze how a system's
    atomic structure deviates from a reference state over a trajectory. It
    generates structural snapshots, calculates the fingerprint for each in
    parallel, and computes the cosine distance to a reference fingerprint.
    Optionally, it can analyze the resulting time-series data to detect
    significant fluctuations. The final data is saved to a text file and plotted.

    Args:
        path_traj (str | list[str]):
            Path to a single HDF5 trajectory file or a list of such paths.
        t_interval (float):
            The time interval in picoseconds (ps) for generating snapshots.
        reference_structure (str):
            Path to the reference structure file (e.g., POSCAR of the initial phase).
        atom_pairs (list[tuple[str, str]] | None, optional):
            A list of atom pairs to include in the fingerprint calculation. If None,
            all unique pair combinations are auto-generated. Defaults to None.
        Rmax (float, optional):
            Cutoff radius (Å) for the fingerprint calculation. Defaults to 10.0.
        delta (float, optional):
            Discretization step (Å) for the fingerprint. Defaults to 0.08.
        sigma (float, optional):
            Gaussian broadening width (Å) for the fingerprint. Defaults to 0.03.
        dirac (str, optional):
            Dirac function type ('g' for Gaussian or 's' for square). Defaults to 'g'.
        prefix (str, optional):
            A prefix for the output plot and data files.
            Defaults to 'cosine_distance_trace'.
        dpi (int, optional):
            Resolution in dots per inch for the saved plot. Defaults to 300.
        path_dir (str, optional):
            Directory to save the final output files. Defaults to 'fingerprint_trace'.
        n_jobs (int, optional):
            Number of CPU cores for parallel processing. -1 uses all available cores.
            Defaults to -1.
        find_fluctuations (bool, optional):
            If True, analyzes the data for significant deviations from the mean.
            Defaults to True.
        window_size (int, optional):
            The window size (number of data points) for the moving average filter.
            Defaults to 50.
        threshold_std (float | None, optional):
            The threshold in standard deviations (σ) from the global mean to define
            a fluctuation. If None, fluctuation intervals are not detected.
            Defaults to None.
        disp (bool, optional):
            If True, displays a plot of the cosine dinstance vs time.
            Defaults to True.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Returns:
        None:
            This function does not return any value. It saves a plot (`.png`) and
            a data file (`.txt`) to the specified `path_dir`.

    Raises:
        IOError:
            If the `reference_structure` cannot be read to auto-generate atom pairs.
        FileNotFoundError:
            If any of the input trajectory files are not found.
    
    Examples:
        >>> plot_cosine_distance(
        ...     path_traj='path/to/trajectory.h5',
        ...     t_interval=0.1,
        ...     reference_structure='path/to/initial_POSCAR',
        ...     atom_pairs=[('Ti', 'O'), ('O', 'O')],
        ...     threshold_std=2.5
        ... )
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

        ref_fingerprint_file = os.path.join(path_dir, "fingerprint_ref.txt")
        ref_fingerprint = get_fingerprint(
            path_structure=reference_structure,
            filename=ref_fingerprint_file,
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
                             ascii=True,
                             disable=not verbose))

    results.sort(key=lambda x: x[0])
    results = np.array(results, dtype=np.float64)
    results[:, 0] = (results[:, 0]) * snapshots.t_interval
    time_data, distance_data = results[:, 0], results[:, 1]

    # Fluctuation analysis (optional)
    fluctuation_intervals = []
    moving_avg = None
    global_mean = np.mean(distance_data)
    if find_fluctuations:
        if len(distance_data) < window_size:
            if verbose:
                print(f"\n[Warning] Data points ({len(distance_data)}) are fewer than window size ({window_size})."
                      "\n          Skipping moving average and fluctuation analysis.")
        else:
            global_std = np.std(distance_data)
            upper_threshold = np.inf
            if threshold_std is not None:
                upper_threshold = global_mean + threshold_std * global_std

            moving_avg = np.convolve(distance_data, np.ones(window_size)/window_size, mode='valid')
            time_avg_start_index = (window_size - 1) // 2
            time_avg_end_index = len(time_data) - (window_size // 2)
            time_avg = time_data[time_avg_start_index:time_avg_end_index]


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

    if find_fluctuations and moving_avg is not None:
        ax.plot(time_avg, moving_avg, color='crimson', lw=2, label=f'Moving Avg (win={window_size})', zorder=3)
        ax.axhline(global_mean, color='k', linestyle='--', label='Global Mean', zorder=1)
        if threshold_std is not None:
            ax.axhline(upper_threshold, color='k', linestyle=':', label=f'Threshold ({threshold_std:.1f}σ)', zorder=1) # Formatted label

        fluctuation_label_added = False
        for start, end in fluctuation_intervals:
            label = 'Detected Fluctuation' if not fluctuation_label_added else ""
            ax.axvspan(start, end, color='orange', alpha=0.3, label=label)
            fluctuation_label_added = True

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), fontsize=9)

    out_png_filename = f'{prefix}.png'
    full_png_path = os.path.join(path_dir, out_png_filename)
    
    ax.set_xlabel("Time (ps)", fontsize=13)
    ax.set_ylabel('Cosine Distance to Reference', fontsize=13)
    ax.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(full_png_path, dpi=dpi)
    if disp: plt.show()
    plt.close()
    if verbose: print(f"\n'{full_png_path}' created successfully.")
    
    output_txt_filename = f'{prefix}.txt'
    full_txt_path = os.path.join(path_dir, output_txt_filename)

    with open(full_txt_path, 'w') as f:
        f.write(f'# Rmax={Rmax}, delta={delta}, sigma={sigma}, dirac={dirac}\n')
        f.write('# pair: ' + ', '.join([f'{A}-{B}' for A, B in atom_pairs]) + '\n')
        f.write('# Time (ps)\tCosine Distance\n')
        np.savetxt(f, results, fmt='%.6f', delimiter='\t')
    if verbose: print(f"'{full_txt_path}' created successfully.")

    if (threshold_std is not None) and find_fluctuations and fluctuation_intervals:
        if verbose: print("\nDetected fluctuation intervals:")
        for start, end in fluctuation_intervals:
            if verbose: print(f"  - From {start:.2f} ps to {end:.2f} ps")