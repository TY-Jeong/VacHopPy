import h5py
import json
import itertools
import numpy as np
from tqdm import tqdm
from pathlib import Path

from vachoppy.utils import monitor_performance

BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m' 
RED = '\033[91m' 
RESET = '\033[0m'


class TrajBundle:
    """Manages a collection of HDF5 trajectory files from simulations.

    This class automatically finds, validates, and groups trajectory files based
    on simulation parameters. It searches for files in a given path, groups them
    by temperature, and ensures that all grouped files are from consistent
    simulations (e.g., same timestep, atom counts, lattice).

    Args:
        path (str): 
            The root directory to search for trajectory files.
        symbol (str): 
            The chemical symbol of the target element to filter files.
        prefix (str, optional): 
            The prefix of the trajectory files to search for.
            Defaults to "TRAJ".
        depth (int, optional): 
            The maximum directory depth to search.
            Defaults to 2.
        eps (float, optional): 
            The tolerance for comparing floating-point values
            like temperature, dt, and lattice vectors. Defaults to 1.0e-3.
        verbose (bool, optional): 
            Verbosity tag. Defaults to True.

    Attributes:
        temperatures (list[float]): 
            A sorted list of unique temperatures found.
        traj (list[list[str]]): 
            A 2D list where each sublist contains file paths
            for a corresponding temperature in `self.temperatures`.
        atom_count (int): 
            The number of atoms for the specified symbol, confirmed
            to be consistent across all files.
        lattice (np.ndarray): 
            The 3x3 lattice vectors as a NumPy array, confirmed
            to be consistent across all files.
    """
    def __init__(self,
                 path: str,
                 symbol: str,
                 prefix: str = "TRAJ",
                 depth: int = 2,
                 eps: float = 1.0e-3,
                 verbose: bool = True):
        
        self.eps = eps
        self.path = Path(path)
        self.depth = depth
        self.symbol = symbol
        self.prefix = prefix
        self.verbose = verbose
        
        if not self.path.is_dir():
            raise NotADirectoryError(f"Error: Invalid directory: '{path}'")
        if self.depth < 1:
            raise ValueError("Error: depth must be 1 or greater.")

        # List of TRAJ_*_.h5 files 
        self.temperatures = []
        self.traj = []
        self.search_traj()
        
        self.atom_count = None
        self.lattice = None
        self._validate_consistency()
        
        if self.verbose:
            self.summary()
        
    def _validate_consistency(self) -> None:
        """
        Validates that all found trajectory files share consistent simulation parameters
        and sets the class attributes for lattice and atom_count.
        """
        all_paths = list(itertools.chain.from_iterable(self.traj))
        
        if not all_paths:
            return

        ref_path = all_paths[0]
        try:
            with h5py.File(ref_path, "r") as f:
                cond = json.loads(f.attrs["metadata"])
                ref_dt = cond.get("dt")
                ref_atom_counts = cond.get("atom_counts")
                ref_lattice = np.array(cond.get("lattice"))
                
                self.lattice = ref_lattice
                self.atom_count = ref_atom_counts.get(self.symbol)

        except Exception as e:
            raise IOError(f"Could not read reference metadata from '{ref_path}'. Reason: {e}")
        
        if len(all_paths) > 1:
            for path in all_paths[1:]:
                try:
                    with h5py.File(path, "r") as f:
                        cond = json.loads(f.attrs['metadata'])
                        
                        current_dt = cond.get("dt")
                        if abs(current_dt - ref_dt) > self.eps:
                            raise ValueError(
                                f"Inconsistent 'dt' parameter found.\n"
                                f"  - Reference '{ref_path}': {ref_dt}\n"
                                f"  - Conflicting '{path}': {current_dt}"
                            )
                        
                        current_atom_counts = cond.get('atom_counts')
                        if current_atom_counts != ref_atom_counts:
                            raise ValueError(
                                f"Inconsistent 'atom_counts' parameter found.\n"
                                f"  - Reference '{ref_path}': {ref_atom_counts}\n"
                                f"  - Conflicting '{path}': {current_atom_counts}"
                            )
                            
                        current_lattice = np.array(cond.get('lattice'))
                        if not np.all(np.abs(current_lattice - ref_lattice) <= self.eps):
                            raise ValueError(
                                f"Inconsistent 'lattice' parameter found.\n"
                                f"  - Reference '{ref_path}':\n{ref_lattice}\n"
                                f"  - Conflicting '{path}':\n{current_lattice}"
                            )
                except Exception as e:
                    if isinstance(e, ValueError):
                        raise
                    raise IOError(f"Could not read or validate metadata from '{path}'. Reason: {e}")
        
    def search_traj(self) -> None:
        """
        Searches for HDF5 trajectory files, groups them by temperature,
        and populates self.temperatures and self.traj lists.
        """
        glob_iterators = []
        for i in range(self.depth):
            dir_prefix = '*/' * i
            pattern = f"{dir_prefix}{self.prefix}*.h5"
            glob_iterators.append(self.path.glob(pattern))    
        candidate_files = itertools.chain.from_iterable(glob_iterators)
        
        found_paths = []
        for file_path in sorted(candidate_files):
            try:
                with h5py.File(file_path, "r") as f:
                    metadata_str = f.attrs.get("metadata")
                    if not metadata_str:
                        if self.verbose:
                            print(f"Warning: File '{file_path.name}' is "+
                                  "missing 'metadata' attribute. Skipping.")
                        continue
                    
                    cond = json.loads(metadata_str)
                    file_symbol = cond.get("symbol")
                    file_temp = cond.get("temperature")
                    
                    if file_symbol == self.symbol:
                        if file_temp is None:
                            if self.verbose:
                                print(f"Warning: File '{file_path.name}' is missing "+
                                      "'temperature' in metadata. Skipping.")
                            continue
                        
                        if "positions" in f and "forces" in f:
                            found_paths.append((str(file_path.resolve()), float(file_temp)))
                        else:
                            if self.verbose:
                                print(f"Warning: File '{file_path.name}' is missing " +
                                    "required datasets ('positions', 'forces'). Skipping.")
            except Exception as e:
                if self.verbose:
                    print(f"Warning: Could not read or verify '{file_path.name}'. "
                          f"Skipping file. Reason: {e}")
        
        if not found_paths:
            raise FileNotFoundError(
                f"Error: No valid trajectory files found for symbol '{self.symbol}' "
                f"in path '{self.path}' with depth {self.depth}."
            )
            
        sorted_files = sorted(found_paths, key=lambda item: item[1])
        temp_groups, traj_groups = [], []
        for path, temp in sorted_files:
            if not traj_groups or abs(temp - temp_groups[-1]) > self.eps:
                temp_groups.append(temp)
                traj_groups.append([path])
            else:
                traj_groups[-1].append(path)
        
        self.temperatures = temp_groups
        self.traj = traj_groups
    
    def summary(self) -> None:
        """
        Creates a 'bundle.txt' file, summarizing the found trajectories.
        Includes system information at the top.
        """
        output_filename = "bundle.txt"
        with open(output_filename, 'w', encoding='utf-8') as f:
            f.write("--- System Information ---\n")
            f.write(f"Symbol        : {self.symbol}\n")
            f.write(f"Atom Count    : {self.atom_count}\n")
            f.write(f"Lattice (Ang) :\n{self.lattice}\n\n")
            f.write("--- Trajectory Summary ---\n\n")
            
            for temp, paths in zip(self.temperatures, self.traj):
                f.write(f"--- Temperature: {temp:.1f} K (Found {len(paths)} files) ---\n")
                
                total_sim_time = 0.0
                
                for path in paths:
                    try:
                        with h5py.File(path, 'r') as h5f:
                            cond = json.loads(h5f.attrs['metadata'])
                            nsw = cond.get('nsw')
                            dt = cond.get('dt')
                            
                            if nsw is not None and dt is not None:
                                sim_time = (nsw * dt) / 1000.0  # in ps
                                total_sim_time += sim_time
                                f.write(f"  - {path}: {sim_time:.2f} ps\n")
                            else:
                                f.write(f"  - {path}: [nsw/dt not found in metadata]\n")
                    except Exception as e:
                        f.write(f"  - {path}: [Error reading metadata: {e}]\n")
                        
                f.write(f"\n  Total simulation time: {total_sim_time:.2f} ps\n\n")
        print(f"Summary written to '{output_filename}'")