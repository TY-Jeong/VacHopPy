import os
import h5py
import json
import time
import inspect
import functools
import tracemalloc
import numpy as np

from tqdm.auto import tqdm
from tabulate import tabulate
from colorama import Fore, Style
from typing import List, Union

from ase import Atoms
from ase.io import write

BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m' # Green color
RED = '\033[91m'   # Red color
RESET = '\033[0m'  # Reset to default color

def monitor_performance(func):
    """
    A decorator that measures and prints the execution time 
    and peak memory usage of a function.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        is_verbose = False
        try:
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()
            is_verbose = bound_args.arguments.get('verbose', False)
        except TypeError:
            pass

        if not is_verbose:
            return func(*args, **kwargs)

        tracemalloc.start()
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        _, peak_mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        print(f"Execution Time: {elapsed_time:.3f} seconds")
        print(f"Peak RAM Usage: {peak_mem / 1024**3:.3f} GB")
        
        return result
    return wrapper


@monitor_performance
def concat_traj(traj1:str,
                traj2:str,
                label:str="CONCAT",
                chunk_size:int=10000,
                eps:float=1.0e-4,
                verbose:bool=True) -> None:
    """
    Concatenates two HDF5 trajectory files after checking for consistency.

    Args:
        traj1 (str):
            Path to the first HDF5 trajectory file.
        traj2 (str):
            Path to the second HDF5 trajectory file.
        label (str, optional):
            A label for the output concatenated file, resulting in a filename
            like 'TRAJ_SYMBOL_LABEL.h5'. Defaults to "CONCAT".
        chunk_size (int, optional):
            The number of frames to process at once during the copy operation.
            Defaults to 10000.
        eps (float, optional):
            Tolerance for comparing floating-point metadata values such as
            temperature, time step, and lattice vectors. Defaults to 1.0e-4.
        verbose (bool, optional):
            Verbosity tag. Defaults to True.

    Returns:
        None: The function saves a new, concatenated HDF5 file to disk.

    Raises:
        FileNotFoundError: If `traj1` or `traj2` is not found.
        ValueError: If the metadata of the two files is inconsistent.
    """
    if not os.path.exists(traj1):
        raise FileNotFoundError(f"Input file not found: {traj1}")
    if not os.path.exists(traj2):
        raise FileNotFoundError(f"Input file not found: {traj2}")
    
    # Check consistency
    with h5py.File(traj1, "r") as f1, h5py.File(traj2, "r") as f2:
        cond1 = json.loads(f1.attrs["metadata"])
        cond2 = json.loads(f2.attrs["metadata"])
        
        # Check symbol
        symbol1 = cond1["symbol"]
        symbol2 = cond2["symbol"]
        if symbol1 == symbol2:
            symbol = symbol1
        else:
            raise ValueError(
                f"Mismatch in chemical symbol: Cannot concatenate. "
                f"'{traj1}' has '{symbol1}', but '{traj2}' has '{symbol2}'."
            )
            
        # Check composition (atom_counts)
        comp1 = cond1["atom_counts"]
        comp2 = cond2["atom_counts"]
        if comp1 == comp2:
            atom_counts = comp1
        else:
            expr1 = "".join(f"{k}{v}" for k, v in sorted(comp1.items()))
            expr2 = "".join(f"{k}{v}" for k, v in sorted(comp2.items()))
            raise ValueError(
                f"Mismatch in composition: Cannot concatenate. "
                f"'{traj1}' has '{expr1}', but '{traj2}' has '{expr2}'."
            )
        
        # Check temperature
        temp1 = cond1["temperature"]
        temp2 = cond2["temperature"]
        if abs(temp1 - temp2) < eps:
            temp = temp1
        else:
            raise ValueError(
                f"Mismatch in temperature: Cannot concatenate. "
                f"'{traj1}' is at {temp1} K, but '{traj2}' is at {temp2} K."
            )
        
        # Check time step (assuming 'potim' is the time step)
        dt1 = cond1["dt"]
        dt2 = cond2["dt"]
        if abs(dt1 - dt2) < eps:
            dt = dt1
        else:
            raise ValueError(
                f"Mismatch in time step: Cannot concatenate. "
                f"'{traj1}' has {dt1} fs/step, but '{traj2}' has {dt2} fs/step."
            )
        
        # Check lattice
        lat1 = np.array(cond1["lattice"])
        lat2 = np.array(cond2["lattice"])
        if np.allclose(lat1, lat2, atol=eps, rtol=0):
            lattice = lat1
        else:
            raise ValueError(
                f"Mismatch in lattice vectors: Cannot concatenate. "
                f"The lattice parameters of '{traj1}' and '{traj2}' differ."
            )
            
        # Concatenate two traj files
        num_frames1 = cond1["nsw"]
        num_frames2 = cond2["nsw"]
        total_frames = num_frames1 + num_frames2
        out_file = f"TRAJ_{symbol}_{label}.h5"
        with h5py.File(out_file, "w") as f_out:
            cond = {
                "symbol": symbol, 
                "nsw": total_frames,
                "temperature": temp,  
                "dt": dt, 
                "atom_counts": atom_counts, 
                "lattice": lattice.tolist()
            } 
            f_out.attrs['metadata'] = json.dumps(cond)

            pos_out = f_out.create_dataset(
                "positions", 
                shape=(total_frames, atom_counts[symbol], 3), 
                dtype=np.float64
            )
            force_out = f_out.create_dataset(
                "forces", 
                shape=(total_frames, atom_counts[symbol], 3), 
                dtype=np.float64
            )

            pbar = tqdm(
                desc=f"{RED}{BOLD}Progress{RESET}", 
                unit=" frames", 
                total=total_frames,
                ascii=False,
                bar_format='{l_bar}%s{bar:20}%s{r_bar}'%(Fore.GREEN, Fore.RESET)
            )

            # Copy from the first file
            for i in range(0, num_frames1, chunk_size):
                end = min(i + chunk_size, num_frames1)
                pos_out[i:end] = f1['positions'][i:end]
                force_out[i:end] = f1['forces'][i:end]
                pbar.update(end - i)

            # Copy from the second file
            for i in range(0, num_frames2, chunk_size):
                end = min(i + chunk_size, num_frames2)
                pos_out[num_frames1 + i : num_frames1 + end] = f2['positions'][i:end]
                force_out[num_frames1 + i : num_frames1 + end] = f2['forces'][i:end]
                pbar.update(end - i)
            pbar.close()
            
    if verbose:
        print(f"Successfully created concatenated file: '{out_file}'")


def cut_traj(traj) -> None:
    # will be updated
    if not os.path.exists(traj):
        raise FileNotFoundError(f"Input file not found: {traj}")
    pass


def show_traj(traj: str) -> None:
    """Displays metadata and dataset info from a trajectory HDF5 file.
    
    Args:
        traj (str):
            Path to the HDF5 trajectory file to inspect.

    Returns:
        None: The function prints information to the console and does not
              return any value.

    Raises:
        FileNotFoundError: If the input HDF5 file is not found.
    """
    if not os.path.exists(traj):
        raise FileNotFoundError(f"Input file not found: {traj}")
    
    with h5py.File(traj, "r") as f:
        try:
            cond = json.loads(f.attrs["metadata"])
        except KeyError:
            print(f"Error: Metadata attribute not found in '{traj}'.")
            return
        
        print("="*40)
        print(f"Trajectory File: {os.path.basename(traj)}")
        print("="*40)
        
        print("\n[Simulation Parameters]")
        print(f"  - Atomic Symbol:      {cond.get('symbol', 'N/A')}")
        print(f"  - Number of Frames:   {cond.get('nsw', 'N/A')}")
        print(f"  - Temperature:        {cond.get('temperature', 'N/A')} K")
        print(f"  - Time Step:          {cond.get('dt', 'N/A')} fs")
        
        atom_counts = cond.get('atom_counts', {})
        if atom_counts:
            print("\n[Composition]")
            composition_str = ", ".join(f"{k}: {v}" for k, v in sorted(atom_counts.items()))
            total_atoms = sum(atom_counts.values())
            print(f"  - Counts:             {composition_str}")
            print(f"  - Total Atoms:        {total_atoms}")
            
        lattice = cond.get('lattice', [])
        if lattice:
            print("\n[Lattice Vectors (Å)]")
            for vector in lattice:
                print(f"  [{vector[0]:>9.5f}, {vector[1]:>9.5f}, {vector[2]:>9.5f}]")
        
        print("\n[Stored Datasets]")
        if 'positions' in f:
            pos_shape = f['positions'].shape
            print(f"  - positions:          Shape = {pos_shape}")
        else:
            print("  - positions:          Not found")
            
        if 'forces' in f:
            force_shape = f['forces'].shape
            print(f"  - forces:             Shape = {force_shape}")
        else:
            print("  - forces:             Not found")    
            
        print("="*40)


class Snapshots:
    """
    Generates step-wise structure files from a set of HDF5 trajectory files.

    This class reads HDF5 trajectory files, reconstructs the full trajectory,
    and calculates averaged atomic positions for specified time intervals.
    The results are stored in memory and can be saved to files using the
    `save_snapshots` method.

    Args:
        traj_files (Union[str, List[str]]):
            A path to a single HDF5 trajectory file or a list of paths.
        t_interval (float):
            The time interval in picoseconds (ps) for averaging snapshots.
            Must be a multiple of the simulation's `dt`.
        eps (float, optional):
            A tolerance for floating-point comparisons (e.g., for temperature).
            Defaults to 1.0e-3.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.
    """
    def __init__(self,
                 traj_files: Union[str, List[str]],
                 t_interval: float,
                 eps: float = 1.0e-3,
                 verbose: bool = True):

        if isinstance(traj_files, str):
            self.traj_files = [traj_files]
        elif isinstance(traj_files, list) and traj_files:
            self.traj_files = traj_files
        else:
            raise ValueError("Input 'traj_files' must be a non-empty string or a list of strings.")

        self.t_interval = t_interval
        self.eps = eps
        self.verbose = verbose
        
        self.total_frames: int = None
        self.dt: float = None
        self.lattice: np.ndarray = None
        self.atom_counts: dict = None
        self.num_atoms: int = None
        self.pos: np.ndarray = None
        self.num_steps: int = None
        self.digit: str = None
        self.frame_interval: int = None
        self._process_trajectories()

    def _process_trajectories(self):
        """Main workflow to validate, load, average, and wrap trajectory data."""
        full_pos_unwrapped = self._validate_and_load_trajectories()

        val = self.t_interval * 1000 / self.dt
        if not np.isclose(val, round(val), atol=self.eps):
            raise ValueError(f"The t_interval ({self.t_interval} ps) must be a multiple "
                             f"of the simulation timestep ({self.dt / 1000.0} ps).")
        self.frame_interval = round(val)

        if self.frame_interval == 0:
            raise ValueError("The t_interval is too small, resulting in zero frames per snapshot.")
            
        self.num_steps = self.total_frames // self.frame_interval
        num_frames_to_use = self.num_steps * self.frame_interval
        self.digit = len(str(self.num_steps - 1)) if self.num_steps > 0 else 1

        pos_sliced = full_pos_unwrapped[:num_frames_to_use]
        
        pos_reshaped = pos_sliced.reshape(self.num_steps, self.frame_interval, self.num_atoms, 3)
        pos_averaged = np.average(pos_reshaped, axis=1)
        
        self.pos = pos_averaged - np.floor(pos_averaged)
        if self.verbose: print(f"Trajectory processed into {self.num_steps} snapshots.")

    def _validate_and_load_trajectories(self) -> np.ndarray:
        if self.verbose: print("Validating and loading HDF5 trajectory files...")
        ref_meta = None
        positions_by_symbol = {}
        
        for traj_file in self.traj_files:
            if not os.path.isfile(traj_file):
                raise FileNotFoundError(f"Input file '{traj_file}' not found.")
            
            with h5py.File(traj_file, 'r') as f:
                meta_str = f.attrs.get('metadata')
                if not meta_str:
                    raise ValueError(f"File '{traj_file}' is missing 'metadata' attribute.")
                
                meta = json.loads(meta_str)
                symbol = meta.get('symbol')
                
                if ref_meta is None:
                    ref_meta = meta
                    self.dt = ref_meta['dt']
                    self.total_frames = ref_meta['nsw']
                    self.atom_counts = ref_meta['atom_counts']
                    self.lattice = np.array(ref_meta['lattice'], dtype=np.float64)
                    self.temperature = ref_meta.get('temperature')
                else:
                    for key in ['nsw', 'atom_counts']:
                        if meta[key] != ref_meta[key]:
                            raise ValueError(f"Metadata mismatch in '{traj_file}': '{key}' differs from reference file.")
                    for key in ['dt', 'temperature']:
                        if key in meta and key in ref_meta and not np.isclose(meta.get(key), ref_meta.get(key), atol=self.eps):
                             raise ValueError(f"Metadata mismatch in '{traj_file}': '{key}' differs from reference file.")
                    if not np.allclose(meta['lattice'], ref_meta['lattice'], atol=self.eps):
                        raise ValueError(f"Metadata mismatch in '{traj_file}': 'lattice' differs from reference file.")

                positions_by_symbol[symbol] = f['positions'][:].astype(np.float64)
                
                if positions_by_symbol[symbol].shape[1] != self.atom_counts[symbol]:
                     raise ValueError(f"Atom count for '{symbol}' in '{traj_file}' does not match metadata.")

        if set(self.atom_counts.keys()) != set(positions_by_symbol.keys()):
            missing = set(self.atom_counts.keys()) - set(positions_by_symbol.keys())
            raise ValueError(f"Trajectory files for the following symbols are missing: {missing}")

        full_pos_list = []
        for symbol in sorted(self.atom_counts.keys()):
            full_pos_list.append(positions_by_symbol[symbol])
        
        full_pos_unwrapped = np.concatenate(full_pos_list, axis=1)
        self.num_atoms = full_pos_unwrapped.shape[1]
        
        if self.verbose: print("All files validated and loaded successfully.")
        
        return full_pos_unwrapped

    def save_snapshots(self,
                       path_dir: str = 'snapshots',
                       format: str = 'vasp',
                       prefix: str = 'POSCAR'):
        """Saves the averaged snapshots as a series of structure files using ASE.

        This method uses ASE (Atomic Simulation Environment) to write the
        structure of each snapshot to a separate file. The output format,
        directory, and filename prefix can be specified.

        Args:
            path_dir (str, optional):
                The directory where output files will be saved. It will be
                created if it does not exist. Defaults to 'snapshots'.
            format (str, optional):
                The output file format supported by `ase.io.write`.
                Defaults to 'vasp'.
            prefix (str, optional):
                The prefix for the output filenames (e.g., 'POSCAR', 'snapshot').
                Defaults to 'POSCAR'.
        """
        if not os.path.isdir(path_dir):
            os.makedirs(path_dir, exist_ok=True)
            if self.verbose: print(f"Created output directory: '{path_dir}'")
            
        sorted_symbols = sorted(self.atom_counts.keys())
        full_symbol_list = [sym for sym in sorted_symbols for _ in range(self.atom_counts[sym])]
            
        if self.verbose: print(f"Saving {self.num_steps} snapshot files to '{path_dir}' (format: {format})...")
        
        for i in tqdm(range(self.num_steps), desc="Saving Snapshots", disable=not self.verbose):
            atoms = Atoms(symbols=full_symbol_list, scaled_positions=self.pos[i], cell=self.lattice, pbc=True)
            filename = f"{prefix}_{i:0{self.digit}d}"
            snapshot_path = os.path.join(path_dir, filename)
            write(snapshot_path, atoms, format=format)

        desc_path = os.path.join(path_dir, "description.txt")
        table_data = []
        headers = ["Filename", "Time (ps)", "Original Frame Range"]
        
        for i in range(self.num_steps):
            filename = f"{prefix}_{i:0{self.digit}d}"
            time_ps = (i + 1) * self.t_interval
            frame_start = i * self.frame_interval
            frame_end = frame_start + self.frame_interval - 1
            table_data.append([filename, f"{time_ps:.2f}", f"{frame_start} - {frame_end}"])
            
        with open(desc_path, 'w') as f:
            f.write("="*60 + "\n")
            f.write("          Snapshot Analysis Description\n")
            f.write("="*60 + "\n\n")
            
            f.write("-- Simulation Parameters --\n")
            f.write(f"  - Source Files        : {', '.join(self.traj_files)}\n")
            f.write(f"  - Temperature         : {self.temperature:.1f} K\n")
            f.write(f"  - Timestep (dt)       : {self.dt:.3f} fs\n")
            f.write(f"  - Total Frames (NSW)  : {self.total_frames}\n")
            f.write(f"  - Atom Counts         : {json.dumps(self.atom_counts)}\n\n")

            f.write("-- Snapshot Parameters --\n")
            f.write(f"  - Snapshot Interval   : {self.t_interval:.3f} ps\n")
            f.write(f"  - Frames per Snapshot : {self.frame_interval}\n")
            f.write(f"  - Total Snapshots     : {self.num_steps}\n\n")

            f.write("-- File Details --\n")
            table = tabulate(table_data, headers=headers, tablefmt="simple", numalign="center")
            f.write(table)
            
        if self.verbose: print(f"Snapshot descriptions saved to '{desc_path}'")