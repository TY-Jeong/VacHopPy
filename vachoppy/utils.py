"""
vachoppy.utils
==============

Provides utility functions and classes for managing, inspecting, and preparing
trajectory and structural data for the `VacHopPy` analysis workflow.

This module contains a collection of helper tools for common data-handling
tasks. These tools operate on the standardized HDF5 trajectory files used
by the package or help in generating data for other analyses.

Main Components
---------------
- **Data Management Functions**:
    - `concat_traj`: Joins two HDF5 trajectory files into a single, continuous file.
    - `cut_traj`: Extracts a specific range of frames from a trajectory file.
    - `show_traj`: Prints a human-readable summary of a trajectory file's metadata.
- **Data Processing Classes**:
    - `Snapshots`: A class that reads one or more trajectories and generates
      time-averaged structural snapshots at regular intervals.
- **Performance Decorators**:
    - `monitor_performance`: A decorator to measure the execution time and peak
      memory usage of a function.

Typical Usage
-------------
**1. Inspecting and Combining Files:**

.. code-block:: python

    from vachoppy.utils import show_traj, concat_traj

    # Inspect two trajectory files
    show_traj('run1.h5')
    show_traj('run2.h5')

    # Concatenate them into a new file
    concat_traj('run1.h5', 'run2.h5', label='concat')

**2. Generating Structural Snapshots:**

.. code-block:: python

    from vachoppy.utils import Snapshots

    # Process a trajectory to get snapshots every 0.1 ps
    snapshot_generator = Snapshots(path_traj='full_run.h5', t_interval=0.1)

    # Save the snapshots as POSCAR files
    snapshot_generator.save_snapshots(path_dir='snapshots_for_analysis')
"""

__all__ =['concat_traj', 'cut_traj', 'show_traj', 'Snapshots', 'monitor_performance']

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
from typing import List, Union
from numpy.typing import DTypeLike

from ase import Atoms
from ase.io import write


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
def concat_traj(path_traj1: str,
                path_traj2: str,
                label:str = "CONCAT",
                chunk_size: int = 10000,
                eps: float = 1.0e-3,
                dtype: DTypeLike = np.float64,
                verbose:bool=True) -> None:
    """Concatenates two HDF5 trajectory files after checking for consistency.

    This function joins two sequential trajectory simulations. It performs a
    thorough check to ensure that critical metadata (symbol, composition,
    temperature, dt, lattice) are consistent between the two files before
    proceeding. It also handles periodic boundary conditions by calculating an
    offset to create a continuous, unwrapped trajectory.

    Args:
        path_traj1 (str):
            Path to the first HDF5 trajectory file.
        path_traj2 (str):
            Path to the second HDF5 trajectory file.
        label (str, optional):
            A label for the output concatenated file, creating a filename
            like 'TRAJ_SYMBOL_LABEL.h5'. Defaults to "CONCAT".
        chunk_size (int, optional):
            The number of frames to process at once during the copy operation.
            Defaults to 10000.
        eps (float, optional):
            Tolerance for comparing floating-point metadata values.
            Defaults to 1.0e-3.
        dtype (DTypeLike, optional):
            NumPy data type for the output arrays. Defaults to np.float64.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Returns:
        None: This function saves a new, concatenated HDF5 file to disk.

    Raises:
        FileNotFoundError: If `path_traj1` or `path_traj2` is not found.
        ValueError: If the metadata of the two files is inconsistent or if
            either file contains no frames.
    """
    if not os.path.exists(path_traj1):
        raise FileNotFoundError(f"Input file not found: {path_traj1}")
    if not os.path.exists(path_traj2):
        raise FileNotFoundError(f"Input file not found: {path_traj2}")
    
    # Check consistency
    with h5py.File(path_traj1, "r") as f1, h5py.File(path_traj2, "r") as f2:
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
                f"'{path_traj1}' has '{symbol1}', but '{path_traj2}' has '{symbol2}'."
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
                f"'{path_traj1}' has '{expr1}', but '{path_traj2}' has '{expr2}'."
            )
        
        # Check temperature
        temp1 = cond1["temperature"]
        temp2 = cond2["temperature"]
        if abs(temp1 - temp2) < eps:
            temp = temp1
        else:
            raise ValueError(
                f"Mismatch in temperature: Cannot concatenate. "
                f"'{path_traj1}' is at {temp1} K, but '{path_traj2}' is at {temp2} K."
            )
        
        # Check time step
        dt1 = cond1["dt"]
        dt2 = cond2["dt"]
        if abs(dt1 - dt2) < eps:
            dt = dt1
        else:
            raise ValueError(
                f"Mismatch in time step: Cannot concatenate. "
                f"'{path_traj1}' has {dt1} fs/step, but '{path_traj2}' has {dt2} fs/step."
            )
        
        # Check lattice
        lat1 = np.array(cond1["lattice"])
        lat2 = np.array(cond2["lattice"])
        if np.allclose(lat1, lat2, atol=eps, rtol=0):
            lattice = lat1
        else:
            raise ValueError(
                f"Mismatch in lattice vectors: Cannot concatenate. "
                f"The lattice parameters of '{path_traj1}' and '{path_traj2}' differ."
            )
            
        # Concatenate two traj files
        num_frames1 = cond1["nsw"]
        num_frames2 = cond2["nsw"]
        if num_frames1 == 0 or num_frames2 == 0:
            raise ValueError("One of the trajectory files has no frames.")
        
        # PBC-unwrapped displacements
        displacement = f2['positions'][0] - f1['positions'][-1]
        offset = -np.round(displacement)
        
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
                dtype=dtype
            )
            force_out = f_out.create_dataset(
                "forces", 
                shape=(total_frames, atom_counts[symbol], 3), 
                dtype=dtype
            )

            pbar = tqdm(
                desc=f"Concatenate Files", 
                unit=" frames", 
                total=total_frames,
                ascii=True,
                bar_format='{l_bar}{bar:30}{r_bar}'
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
                pos_out[num_frames1 + i : num_frames1 + end] = f2['positions'][i:end] + offset
                force_out[num_frames1 + i : num_frames1 + end] = f2['forces'][i:end]
                pbar.update(end - i)
            pbar.close()
            
    if verbose: print(f"Successfully created concatenated file: '{out_file}'")


def cut_traj(path_traj: str,
             start_frame: int | None = None,
             end_frame: int | None = None,
             label: str = "CUT",
             chunk_size: int = 5000) -> None:
    """Cuts a portion of a trajectory file and saves it as a new file.

    This function extracts a specific range of frames (from `start_frame` to
    `end_frame`) from a source HDF5 trajectory and saves it into a new HDF5
    file with updated metadata.

    Args:
        path_traj (str):
            Path to the source HDF5 trajectory file.
        start_frame (int):
            The starting frame number to include in the new trajectory (inclusive).
        end_frame (int):
            The ending frame number to include in the new trajectory (exclusive).
        label (str, optional):
            A label for the output cut file, creating a filename like
            'TRAJ_SYMBOL_LABEL.h5'. Defaults to "CUT".
        chunk_size (int, optional):
            The number of frames to process in each chunk to conserve memory.
            Defaults to 5000.

    Returns:
        None: This function saves a new, shorter HDF5 file to disk.

    Raises:
        FileNotFoundError: If `path_traj` is not found.
        ValueError: If the specified frame range is invalid.
    """
    if not os.path.exists(path_traj):
        raise FileNotFoundError(f"{path_traj} not found.")
    
    base, ext = os.path.splitext(path_traj)
    output_path = f"{base}_{label}{ext}"

    with h5py.File(path_traj, 'r') as f_in:
        if 'metadata' not in f_in.attrs:
            raise ValueError("HDF5 file is missing the 'metadata' attribute.")
            
        metadata = json.loads(f_in.attrs['metadata'])
        total_frames = metadata.get('nsw')
        
        if total_frames is None:
            raise ValueError("Key 'nsw' (number of frames) not found in metadata.")

        start = start_frame if start_frame is not None else 0
        end = end_frame if end_frame is not None else total_frames

        if not (0 <= start < end <= total_frames):
            raise ValueError(
                f"Invalid frame range. Must satisfy 0 <= start < end <= {total_frames}. "
                f"Got start={start}, end={end}."
            )
        
        new_num_frames = end - start
        
        with h5py.File(output_path, 'w') as f_out:
            new_metadata = metadata.copy()
            new_metadata['nsw'] = new_num_frames
            f_out.attrs['metadata'] = json.dumps(new_metadata)

            datasets_to_copy = ['positions', 'forces']
            for dset_name in datasets_to_copy:
                if dset_name not in f_in:
                    raise ValueError(f"Required dataset '{dset_name}' not found in source file.")

                source_dset = f_in[dset_name]
                
                if new_num_frames <= chunk_size:
                    data_slice = source_dset[start:end]
                    f_out.create_dataset(dset_name, data=data_slice)
                else:
                    new_shape = (new_num_frames,) + source_dset.shape[1:]
                    out_dset = f_out.create_dataset(dset_name, shape=new_shape, dtype=source_dset.dtype)
                    
                    for i in tqdm(range(0, new_num_frames, chunk_size),
                                  desc=f'Cut Trajectory',
                                  bar_format='{l_bar}{bar:30}{r_bar}',
                                  ascii=True):
                        
                        read_start = start + i
                        read_end = min(start + i + chunk_size, end)
                        
                        write_start = i
                        write_end = i + (read_end - read_start)
                        
                        out_dset[write_start:write_end] = source_dset[read_start:read_end]

    print(f"\nSuccessfully created cut trajectory file: {output_path}")

def show_traj(path_traj: str) -> None:
    """Displays metadata and dataset info from a trajectory HDF5 file.

    This utility function reads the metadata and dataset shapes from a given
    HDF5 trajectory file and prints a formatted, human-readable summary to the
    console.

    Args:
        path_traj (str):
            Path to the HDF5 trajectory file to inspect.

    Returns:
        None: This function prints information to the console.

    Raises:
        FileNotFoundError: If the input HDF5 file is not found.

    Examples:
        >>> show_traj('path/to/my_trajectory.h5')
    """
    if not os.path.exists(path_traj):
        raise FileNotFoundError(f"Input file not found: {path_traj}")
    
    with h5py.File(path_traj, "r") as f:
        try:
            cond = json.loads(f.attrs["metadata"])
        except KeyError:
            print(f"Error: Metadata attribute not found in '{path_traj}'.")
            return
        
        print("="*50)
        print(f"  Trajectory File: {os.path.basename(path_traj)}")
        print("="*50)
        
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
            print("\n[Lattice Vectors (Ang)]")
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
            
        print("="*50)


class Snapshots:
    """Generates step-wise structure files from a set of HDF5 trajectory files.

    This class reads one or more HDF5 trajectory files, validates their
    consistency, and reconstructs the full, continuous atomic trajectory. It then
    calculates averaged atomic positions for user-specified time intervals. The
    resulting structures ("snapshots") are stored in memory and can be saved to
    individual files using the `save_snapshots` method.

    Args:
        path_traj (str | list[str]):
            A path to a single HDF5 trajectory file or a list of such paths.
        t_interval (float):
            The time interval in picoseconds (ps) for averaging snapshots. This
            must be a multiple of the simulation's `dt`.
        eps (float, optional):
            A tolerance for floating-point comparisons of metadata.
            Defaults to 1.0e-3.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Attributes:
        pos (numpy.ndarray):
            A 3D array of shape (num_steps, num_atoms, 3) containing the
            time-averaged, wrapped fractional coordinates for each snapshot.
        num_steps (int):
            The total number of snapshots generated.
        dt (float):
            The simulation timestep in femtoseconds (fs), read from metadata.
        lattice (numpy.ndarray):
            The 3x3 lattice matrix of the simulation cell.
        atom_counts (dict):
            A dictionary of atom counts for each chemical symbol.
        num_atoms (int):
            The total number of atoms in the system.

    Raises:
        FileNotFoundError: If any of the input trajectory files are not found.
        ValueError: If `path_traj` is empty, if `t_interval` is invalid, or if
            metadata is inconsistent across multiple trajectory files.

    Examples:
        >>> # Create snapshots every 10 ps from a list of files
        >>> snapshot_generator = Snapshots(
        ...     path_traj=['TRAJ_Hf_run1.h5', 'TRAJ_O_run1.h5'],
        ...     t_interval=0.1
        ... )
        >>> # Save the snapshots as VASP POSCAR files
        >>> snapshot_generator.save_snapshots(
        ...     path_dir='poscar_snapshots',
        ...     file_format='vasp',
        ...     prefix='POSCAR'
        ... )
    """
    def __init__(self,
                 path_traj: str | list[str],
                 t_interval: float,
                 eps: float = 1.0e-3,
                 verbose: bool = True):

        if isinstance(path_traj, str):
            self.path_traj = [path_traj]
        elif isinstance(path_traj, list) and path_traj:
            self.path_traj = path_traj
        else:
            raise ValueError("Input 'path_traj' must be a non-empty string or a list of strings.")

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
        
        for traj_file in self.path_traj:
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

        This method uses the Atomic Simulation Environment (ASE) to write the
        atomic structure of each generated snapshot to a separate file. It also
        creates a `description.txt` file in the output directory summarizing the
        snapshot parameters and the mapping from filename to simulation time.

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

        Returns:
            None: This method does not return a value; it writes files to disk.
        """
        if not os.path.isdir(path_dir):
            os.makedirs(path_dir, exist_ok=True)
            if self.verbose: print(f"Created output directory: '{path_dir}'")
            
        sorted_symbols = sorted(self.atom_counts.keys())
        full_symbol_list = [sym for sym in sorted_symbols for _ in range(self.atom_counts[sym])]
            
        if self.verbose: print(f"Saving {self.num_steps} snapshot files to '{path_dir}' (format: {format})...")
        
        for i in tqdm(range(self.num_steps), 
                      desc=f'Generate Snapshots',
                      bar_format='{l_bar}{bar:30}{r_bar}',
                      ascii=True,
                      disable=not self.verbose):
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
            f.write(f"  - Source Files        : {', '.join(self.path_traj)}\n")
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
