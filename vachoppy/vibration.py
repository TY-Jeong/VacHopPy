"""
vachoppy.vibration
==================

Provides the `Vibration` class for calculating the characteristic atomic
vibrational frequency from a molecular dynamics trajectory.

The primary purpose of this module is to determine a suitable time scale for
coarse-graining the main trajectory analysis. The calculated mean vibrational
frequency is the inverse of the `t_interval` parameter used throughout `vachoppy`.
This allows for a clear, data-driven distinction between rapid atomic vibrations
and slower, diffusive hopping events.

Main Components
---------------
- **Vibration**: A class that analyzes atomic displacements within a short
  trajectory segment. It uses a statistical approach and a Fast Fourier Transform
  (FFT) to determine the mean vibrational frequency.

Typical Usage
-------------
This class is often used to automatically estimate the `t_interval` parameter
before running a full diffusion analysis.

.. code-block:: python

    from vachoppy.core import Site
    from vachoppy.vibration import Vibration

    # 1. First, set up the site information
    site_info = Site(path_structure="path/to/POSCAR", symbol="O")

    # 2. Initialize the Vibration class with a trajectory
    vib_analyzer = Vibration(
        path_traj="path/to/TRAJ_O.h5",
        site=site_info
    )

    # 3. Run the frequency calculation
    vib_analyzer.calculate()

    # 4. Access the result to determine a suitable t_interval
    if vib_analyzer.mean_frequency > 0:
        estimated_t_interval = 1 / vib_analyzer.mean_frequency
        print(f"Estimated t_interval: {estimated_t_interval:.3f} ps")

    # This estimated_t_interval can then be passed to a Calculator object.
"""

from __future__ import annotations

__all__ =['Vibration']

import os
import h5py
import json
import itertools
import numpy as np
import matplotlib.pyplot as plt

from tqdm.auto import tqdm
from scipy.stats import norm
from scipy.spatial.distance import cdist
from joblib import Parallel, delayed

from vachoppy.utils import monitor_performance


# ============================================
#   Helper functions for parallel processing
# ============================================

def _helper_distance_pbc(coord1, coord2, lattice):
    """(Helper) Calculates PBC-aware distance between two fractional coordinates."""
    displacement_frac = coord1 - coord2
    displacement_frac -= np.round(displacement_frac)
    return np.linalg.norm(np.dot(displacement_frac, lattice))

def _helper_get_segment_trajectory(trajectory, jump_detection_radius, lattice):
    """(Helper) Segments a single atom's trajectory into vibrational periods."""
    segments = []; start_index = 0
    if len(trajectory) < 10: return segments
    for i in range(1, len(trajectory)):
        segment_center = np.mean(trajectory[start_index:i], axis=0)
        distance = _helper_distance_pbc(trajectory[i], segment_center, lattice)
        if distance > jump_detection_radius:
            segment = trajectory[start_index:i]
            if len(segment) > 10: segments.append(segment)
            start_index = i
    segment_final = trajectory[start_index:]
    if len(segment_final) > 10: segments.append(segment_final)
    return segments

def _worker_get_displacements(args):
    """[Parallel Worker] Calculates displacements for a single atom's trajectory."""
    atom_traj_frac, jump_detection_radius, lattice = args
    segments = _helper_get_segment_trajectory(atom_traj_frac, jump_detection_radius, lattice)
    displacements_cart = []
    for seg in segments:
        center = np.mean(seg, axis=0)
        displacement = seg - center
        displacement -= np.round(displacement)
        displacements_cart.append(np.dot(displacement, lattice))
    return displacements_cart

def _worker_assign_sites(args):
    """[Parallel Worker] Assigns atoms to the nearest lattice site for a single timestep."""
    positions_at_t, site_positions_frac, site_radius, lattice = args
    positions_at_t, site_positions_frac, site_radius, lattice = args
    distance_matrix = cdist(positions_at_t, site_positions_frac,
                            lambda u, v: _helper_distance_pbc(u, v, lattice))
    closest_site_indices = np.argmin(distance_matrix, axis=1)
    min_distances = np.min(distance_matrix, axis=1)
    assignments_at_t = np.full(positions_at_t.shape[0], -1, dtype=int)
    assigned_mask = min_distances < site_radius
    assignments_at_t[assigned_mask] = closest_site_indices[assigned_mask]
    return assignments_at_t

def _worker_get_frequencies(args):
    """[Parallel Worker] Calculates vibrational frequencies for a single atom."""
    atom_traj_frac, atom_assignments, site_positions_frac, dt_s, lattice = args
    frequencies = []
    
    first_valid_step_arr = np.where(atom_assignments > -1)[0]
    if len(first_valid_step_arr) == 0: return []
    first_valid_step = first_valid_step_arr[0]
    
    filtered_traj = atom_traj_frac[first_valid_step:]
    filtered_assign = atom_assignments[first_valid_step:]
    
    jump_indices = np.where(filtered_assign[:-1] != filtered_assign[1:])[0] + 1
    seg_starts = np.insert(jump_indices, 0, 0)
    seg_ends = np.append(jump_indices, len(filtered_assign))
    
    for start, end in zip(seg_starts, seg_ends):
        assigned_site_id = filtered_assign[start]
        if assigned_site_id == -1: continue
        segment = filtered_traj[start:end]
        if len(segment) < 20 : continue
        
        site_center = site_positions_frac[assigned_site_id]
        displacement = segment - site_center
        displacement -= np.round(displacement)
        displacement_cart = np.dot(displacement, lattice)
        
        for axis in range(3):
            disp_axis = displacement_cart[:, axis]
            n = len(disp_axis); hann = np.hanning(n)
            power = np.abs(np.fft.fft(disp_axis * hann))**2
            freq_hz = np.fft.fftfreq(n, d=dt_s)
            mask = freq_hz > 0
            if not np.any(mask): continue
            freqs, power = freq_hz[mask], power[mask]
            frequencies.append(freqs[np.argmax(power)] / 1e12)
    return frequencies

# ============================================

class Vibration:
    """Calculates atomic vibrational frequencies from a molecular dynamics trajectory.

    This class analyzes a short segment of a trajectory to determine the
    characteristic vibrational frequency of atoms. The workflow involves:
    1.  Determining a data-driven "site radius" based on atomic displacements.
    2.  Assigning atoms to their nearest lattice sites.
    3.  Segmenting individual atom trajectories into periods of stable vibration.
    4.  Calculating the frequency spectrum for all segments using FFT.

    The main entry point is the `.calculate()` method. The results can be
    visualized with the `.plot_*()` methods or inspected via attributes.

    Args:
        path_traj (str):
            Path to the HDF5 trajectory file.
        site (Site):
            An initialized `Site` object containing lattice site information.
        sampling_size (int, optional):
            Number of initial trajectory frames to use for the analysis.
            Defaults to 5000.
        filter_high_freq (bool, optional):
            If True, filters out high-frequency outliers using the IQR method.
            Defaults to True.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Attributes:
        mean_frequency (float):
            The mean of all calculated vibrational frequencies in THz.
        frequencies (list[float]):
            A list containing all individual vibrational frequencies calculated
            from the trajectory segments.
        displacements (numpy.ndarray):
            A flattened array of all measured atomic displacements (Å) during
            vibrational periods.
        site_radius (float):
            The calculated site radius (2σ of displacements) in Å, used for
            atom-to-site assignment.

    Raises:
        FileNotFoundError: If the `path_traj` file does not exist.
        ValueError: If the HDF5 file is missing required data or metadata.
        IOError: If the HDF5 file cannot be read.
    """
    def __init__(self,
                 path_traj: str,
                 site: Site,
                 sampling_size: int = 5000,
                 filter_high_freq : bool = True,
                 verbose: bool = True):

        self.path_traj = path_traj
        self.site = site
        self.filter_high_freq = filter_high_freq
        self.verbose = verbose
        
        self._validate_traj(self.path_traj)
        self.site_positions = np.array([s['coord'] for s in self.site.lattice_sites])
        
        self.dt = None
        self.symbol = None
        self.lattice = None
        self.sampling_size = None
        self.positions = None   # fractional coordinates
        self._read_traj(sampling_size)
        
        self.displacements = None
        self.mu_displacements = None
        self.sigma_displacements = None
        self.site_radius = None
        
        self.frequencies = None
        self.mean_frequency = None
        
    def _validate_traj(self, path_traj: str) -> None:
        """
        Validates the structure and content of the HDF5 trajectory file.

        This method checks for the correct file extension, file existence, and
        the presence of required datasets ('positions', 'forces') and metadata
        attributes ('symbol', 'nsw', 'dt', 'temperature', 'atom_counts', 'lattice').

        Args:
            path_traj (str): The file path to validate.

        Raises:
            ValueError: If the file extension is not '.h5' or if required
                        metadata or datasets are missing.
            FileNotFoundError: If the trajectory file does not exist.
            IOError: If the file cannot be read as an HDF5 file.
        """
        if not path_traj.endswith('.h5'):
            raise ValueError(f"Error: Trajectory file must have a .h5 extension, but got '{path_traj}'.")

        if not os.path.isfile(path_traj):
            raise FileNotFoundError(f"Error: Input file '{path_traj}' not found.")
        
        try:
            with h5py.File(path_traj, "r") as f:
                required_datasets = ["positions", "forces"]
                for dataset in required_datasets:
                    if dataset not in f:
                        raise ValueError(f"Error: Required dataset '{dataset}' not found in '{path_traj}'.")

                metadata_str = f.attrs.get("metadata")
                if not metadata_str:
                    raise ValueError(f"Error: Required attribute 'metadata' not found in '{path_traj}'.")
                
                cond = json.loads(metadata_str)
                required_keys = ["symbol", "nsw", "dt", "temperature", "atom_counts", "lattice"]
                for key in required_keys:
                    if key not in cond:
                        raise ValueError(f"Error: Required key '{key}' not found in metadata of '{path_traj}'.")

        except (IOError, OSError) as e:
            raise IOError(f"Error: Failed to read '{path_traj}' as an HDF5 file. Reason: {e}")
        
    def _read_traj(self, sampling_size: int) -> None:
        """Reads metadata and a chunk of trajectory data from the HDF5 file."""
        with h5py.File(self.path_traj, 'r') as f:
            cond = json.loads(f.attrs['metadata'])
            self.dt = cond.get('dt')
            self.symbol = cond.get('symbol')
            self.lattice = np.array(cond.get('lattice'), dtype=np.float64)
            
            num_frames = cond.get('nsw')
            self.sampling_size = min(sampling_size, num_frames)
            self.positions = f['positions'][:self.sampling_size].astype(np.float64)
            
    def _filter_frequencies_iqr(self, frequencies: list) -> list:
        """Filters high-frequency outliers from a list of frequencies using the IQR method."""
        freq_array = np.array(frequencies)
        q1 = np.percentile(freq_array, 25)
        q3 = np.percentile(freq_array, 75)
        iqr = q3 - q1
        upper_bound = q3 + 1.5 * iqr
        
        filtered_frequencies = freq_array[freq_array < upper_bound]
        removed_count = len(frequencies) - len(filtered_frequencies)
        
        if self.verbose:
            print("="*52)
            print(f"       High-Frequency Filtering Results (IQR)")
            print("="*52)
            print(f"  - Cutoff Frequency              : {upper_bound:.2f} THz")
            print(f"  - Removed Outlier Frequencies   : {removed_count} (out of {len(frequencies)})")

        return filtered_frequencies.tolist()
    
    def _get_site_radius(self, 
                         n_jobs: int = -1,
                         jump_detection_radius: float = 1.0) -> None:
        """Calculates the vibrational amplitude and determines the site radius."""
        n_atoms = self.positions.shape[1]
        tasks = [(self.positions[:, i, :], jump_detection_radius, self.lattice) for i in range(n_atoms)]
        
        results = Parallel(n_jobs=n_jobs, verbose=0)(
            delayed(_worker_get_displacements)(task) 
            for task in tqdm(tasks, 
                             desc=f"Compute Displacement",
                             bar_format='{l_bar}{bar:30}{r_bar}',
                             ascii=True,
                             disable=not self.verbose)
        )
        
        all_displacements_cart = list(itertools.chain.from_iterable(results))
        
        if not all_displacements_cart:
            raise ValueError("Could not find any valid vibrational segments to analyze.")
        
        self.displacements = np.concatenate(all_displacements_cart).flatten()
        self.mu_displacements, self.sigma_displacements = norm.fit(self.displacements)
        self.site_radius = 2 * self.sigma_displacements
        
    def plot_displacements(self, 
                           bins: int = 50,
                           disp: bool = True,
                           save: bool = True,
                           title: str | None = "Displacement Distribution",
                           filename: str = "displacement.png",
                           dpi : int = 300) -> None:
        """Plots a histogram of the atomic displacements with a Gaussian fit.

        Args:
            bins (int, optional): Number of bins for the histogram. Defaults to 50.
            disp (bool, optional): If True, displays the plot. Defaults to True.
            save (bool, optional): If True, saves the plot to a file. Defaults to True.
            title (str | None, optional): A custom title for the plot.
            filename (str, optional): Filename for the saved plot.
            dpi (int, optional): Resolution for the saved figure.

        Raises:
            AttributeError: If `.calculate()` has not been run yet.
        """
        if self.displacements is None:
            raise AttributeError("Displacement data not available. Please run the .calculate() method first.")
        
        plt.figure(figsize=(10, 6))
        plt.hist(
            self.displacements,
            bins=bins,
            density=True,
            color='skyblue',
            alpha=0.7,
            edgecolor='black',
            label=f"{self.symbol} Displacements"
        )
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 200)
        p = norm.pdf(x, self.mu_displacements, self.sigma_displacements)
        plt.plot(x, p, 'r-', linewidth=2, label="Gaussian Fit")
        plt.title(title)
        plt.xlabel("Displacement (Ang)", fontsize=12)
        plt.ylabel("Probability Density", fontsize=12)
        plt.legend()
        plt.grid(True, linestyle='--')
        
        if save: plt.savefig(filename, dpi=dpi)
        if disp: plt.show()
        
    def plot_frequencies(self,
                         bins: int = 50,
                         disp: bool = True,
                         save: bool = True,
                         title: str | None = "Frequency Distribution",
                         filename: str = "frequency.png",
                         dpi: int = 300) -> None:
        """Plots a histogram of the calculated vibrational frequencies.

        Args:
            bins (int, optional): Number of bins for the histogram. Defaults to 50.
            disp (bool, optional): If True, displays the plot. Defaults to True.
            save (bool, optional): If True, saves the plot to a file. Defaults to True.
            title (str | None, optional): A custom title for the plot.
            filename (str, optional): Filename for the saved plot.
            dpi (int,optional): Resolution for the saved figure.

        Raises:
            AttributeError: If `.calculate()` has not been run yet.
        """
        if self.frequencies is None:
            raise AttributeError("Frequency data not available. Please run the .calculate() method first.")
        
        plt.figure(figsize=(10, 6))
        plt.hist(
            self.frequencies,
            bins=bins,
            density=True,
            color='mediumpurple',
            alpha=0.7,
            edgecolor='black',
            label=f'{self.symbol} Frequencies'
        )
        plt.axvline(
            self.mean_frequency, 
            color='r', 
            linestyle='--', 
            linewidth=2, 
            label=f'Mean: {self.mean_frequency:.2f} THz'
        )
        plt.title(title)
        plt.xlabel('Frequency (THz)', fontsize=12) 
        plt.ylabel('Probability Density', fontsize=12)
        plt.legend()
        plt.grid(True, linestyle='--')
        if save: plt.savefig(filename, dpi=dpi)
        if disp: plt.show()
    
    @monitor_performance
    def calculate(self, 
                  n_jobs: int = -1,
                  jump_detection_radius: float = 1.0,
                  verbose: bool | None = None) -> None:
        """Executes the full vibrational frequency analysis workflow.

        This is the main method to run the analysis. It performs the following steps:
        1. Calculates the average vibrational amplitude to determine a site radius.
        2. Assigns each atom to its nearest lattice site at each timestep.
        3. Segments the trajectory for each atom based on these site assignments.
        4. Calculates the vibrational frequencies for all segments using FFT.
        5. Optionally filters high-frequency outliers.

        The results are stored in the object's attributes (e.g., `self.mean_frequency`).

        Args:
            n_jobs (int, optional):
                The number of CPU cores for parallel processing. -1 uses all
                available cores. Defaults to -1.
            jump_detection_radius (float, optional):
                The radius (Å) used to distinguish between vibrations and jumps
                during the initial amplitude estimation. Defaults to 1.0.
            verbose (bool | None, optional):
                Overrides the class-level verbosity for this method run.
                If None, the class-level setting is used. Defaults to None.
        """
        if verbose is None:
            verbose = self.verbose
        
        # Site radius estimation
        self._get_site_radius(n_jobs=n_jobs, 
                              jump_detection_radius=jump_detection_radius)

        # Site assignment
        n_steps, n_atoms, _ = self.positions.shape
        dt_s = self.dt * 1e-15
        site_assign_tasks = [(self.positions[i], 
                              self.site_positions, 
                              self.site_radius, 
                              self.lattice) for i in range(n_steps)]
        site_assignments_list = Parallel(n_jobs=n_jobs, verbose=0)(
            delayed(_worker_assign_sites)(task) 
            for task in tqdm(site_assign_tasks, 
                             desc=f'Capture Vibrations  ',
                             bar_format='{l_bar}{bar:30}{r_bar}',
                             ascii=True,
                             disable=not self.verbose)
        )
        site_assignments = np.array(site_assignments_list)
        
        # Frequency calculation
        freq_tasks = [(self.positions[:, i, :], 
                       site_assignments[:, i], 
                       self.site_positions, 
                       dt_s, 
                       self.lattice) for i in range(n_atoms)]
        results = Parallel(n_jobs=n_jobs, verbose=0)(
            delayed(_worker_get_frequencies)(task) 
            for task in tqdm(freq_tasks, 
                             desc=f"Compute Frequenciy  ",
                             bar_format='{l_bar}{bar:30}{r_bar}',
                             ascii=True,
                             disable=not self.verbose)
        )
        frequencies = list(itertools.chain.from_iterable(results))
        if self.verbose: print("")
        
        if self.filter_high_freq:
            frequencies = self._filter_frequencies_iqr(frequencies)
        self.frequencies = frequencies
        
        if self.frequencies:
            self.mean_frequency = np.mean(self.frequencies)
        else:
            self.mean_frequency = 0
        
        if self.verbose:
            self.summary()
            
    def summary(self):
        """Prints a formatted summary of the vibrational analysis results.
        
        The summary includes the mean displacement, determined site radius,
        and the final mean vibrational frequency.
        """
        print("="*52)
        print(f"       Vibrational Analysis Results Summary")
        print("="*52)
        print(f"  - Mean Vibrational Amplitude (σ) : {self.sigma_displacements:.3f} Ang")
        print(f"  - Determined Site Radius (2 x σ) : {self.site_radius:.3f} Ang")
        print(f"  - Total Vibrational Frequencies  : {len(self.frequencies)} found")
        print(f"  - Mean Vibrational Frequency     : {self.mean_frequency:.3f} THz")
        print("="*52 + "\n")