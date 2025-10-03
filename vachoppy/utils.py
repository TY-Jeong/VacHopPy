import os
import h5py
import json
import time
import inspect
import functools
import tracemalloc
import numpy as np

from tqdm.auto import tqdm
from colorama import Fore, Style

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