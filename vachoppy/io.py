import os
import time
import h5py
import json
import itertools
import numpy as np
import MDAnalysis as mda
from tqdm import tqdm
from ase.io import read, iread
from pathlib import Path
from colorama import Fore
from tabulate import tabulate

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from vachoppy.utils import monitor_performance

BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m' 
RED = '\033[91m' 
RESET = '\033[0m'

@monitor_performance
def parse_md(filename: str,
             format: str,
             temp: float,
             dt: float = 1.0,
             label: str = None,
             chunk_size: int = 5000,
             verbose: bool = True) -> None:
    """Parses a generic MD trajectory using ASE and saves data to HDF5 files.

    This function processes large trajectory files (e.g., VASP outputs, extxyz)
    supported by ASE. It reads the trajectory in chunks to manage memory usage
    efficiently. The core calculation involves converting atomic positions to
    unwrapped fractional coordinates to correctly handle periodic boundaries.

    For each chemical symbol, a separate HDF5 file (`TRAJ_SYMBOL_LABEL.h5`) is
    created. Each file contains 'positions' and 'forces' datasets, along with
    simulation metadata (lattice, temperature, etc.) stored as an attribute.
    This method avoids creating intermediate files for improved I/O performance.

    Args:
        filename (str):
            Path to the MD trajectory file (e.g., 'vasprun.xml', 'OUTCAR').
        format (str):
            File format string supported by ASE (e.g., 'vasp-xml', 'extxyz').
        temp (float):
            Simulation temperature in Kelvin.
        dt (float, optional):
            Timestep in femtoseconds. Defaults to 1.0.
        label (str, optional):
            A custom suffix for the output HDF5 filenames. Defaults to None.
        chunk_size (int, optional):
            The number of frames to process in memory at once. A larger size
            can be faster but uses more RAM. Defaults to 5000.
        verbose (bool, optional):
            Verbosity tag. Defaults to True.

    Returns:
        None: The function saves one or more HDF5 files to disk.

    Raises:
        FileNotFoundError: If the input trajectory file is not found.
        ValueError: If force data is missing from the trajectory.
    """
    
    temp_files = {'pos': {}, 'force': {}}
    
    try:
        # Read the last frame
        atoms = read(filename, format=format)
        lattice = atoms.cell.tolist()
        inv_lattice = np.linalg.inv(lattice)
        symbols = atoms.get_chemical_symbols()
        unique_symbols = np.unique(symbols)
        atom_indices = {sym: np.where(np.array(symbols) == sym)[0] for sym in unique_symbols}
        atom_counts = {sym: len(indices) for sym, indices in atom_indices.items()}
        
        for sym in unique_symbols:
            temp_files['pos'][sym] = []
            temp_files['force'][sym] = []
        
        num_frames = 0
        prev_positions = None
        generator = iread(filename, index=':', format=format)
        
        pbar = tqdm(desc=f"{RED}{BOLD}progress{RESET}", unit=" frames")
        
        while True:
            chunk_atoms = list(itertools.islice(generator, chunk_size))
            if not chunk_atoms:
                break
            
            current_chunk_size = len(chunk_atoms)
            num_frames += current_chunk_size
            pbar.update(current_chunk_size)
            
            try:
                chunk_forces = np.array(
                    [atoms.get_forces() for atoms in chunk_atoms], dtype=np.float64
                )
            except AttributeError:
                raise ValueError("Force data not found in the trajectory.")
            
            chunk_positions = np.array(
                [atoms.get_positions() for atoms in chunk_atoms], dtype=np.float64
            )
            chunk_positions = chunk_positions @ inv_lattice # convert to fractional coords
            
            displacement = np.zeros_like(chunk_positions)
            if current_chunk_size > 1:
                displacement[1:, :] = np.diff(chunk_positions, axis=0)
                displacement -= np.round(displacement)
                displacement = np.cumsum(displacement, axis=0)
            
            if prev_positions is None:
                positions_init = chunk_positions[0]
            else:
                displacement_init = chunk_positions[0] - prev_positions
                displacement_init -= np.round(displacement_init)
                positions_init = prev_positions + displacement_init
            
            chunk_positions = positions_init + displacement
            prev_positions = chunk_positions[-1].copy()
            
            for sym, indices in atom_indices.items():
                temp_pos_path = f"temp_pos_{sym}_{pbar.n}.npy"
                temp_force_path = f"temp_force_{sym}_{pbar.n}.npy"
                temp_files['pos'][sym].append(temp_pos_path)
                temp_files['force'][sym].append(temp_force_path)
                
                temp_pos_memmap = np.lib.format.open_memmap(
                    temp_pos_path,
                    mode="w+",
                    dtype=np.float64,
                    shape=(len(chunk_atoms), len(indices), 3)
                )
                temp_force_memmap = np.lib.format.open_memmap(
                    temp_force_path,
                    mode="w+",
                    dtype=np.float64,
                    shape=(len(chunk_atoms), len(indices), 3)
                )
                
                temp_pos_memmap[:] = chunk_positions[:, indices, :]
                temp_force_memmap[:] = chunk_forces[:, indices, :]       
        pbar.close()
                
        for sym in unique_symbols:
            base_name = f"{sym}" if label is None else f"{sym}_{label}"
            out_file = f"TRAJ_{base_name}.h5"
            
            with h5py.File(out_file, "w") as f_h5:
                pos_dataset = f_h5.create_dataset(
                    "positions", 
                    shape=(num_frames, atom_counts[sym], 3), 
                    dtype=np.float64
                )
                force_dataset = f_h5.create_dataset(
                    "forces", 
                    shape=(num_frames, atom_counts[sym], 3), 
                    dtype=np.float64
                )
                
                current_pos = 0
                for temp_path in temp_files['pos'][sym]:
                    temp_data = np.load(temp_path)
                    chunk_len = len(temp_data)
                    pos_dataset[current_pos : current_pos + chunk_len] = temp_data
                    current_pos += chunk_len
                    os.remove(temp_path) # Clean up

                current_pos = 0
                for temp_path in temp_files['force'][sym]:
                    temp_data = np.load(temp_path)
                    chunk_len = len(temp_data)
                    force_dataset[current_pos : current_pos + chunk_len] = temp_data
                    current_pos += chunk_len
                    os.remove(temp_path)
 
                cond = {
                    "symbol": sym, 
                    "nsw": num_frames,
                    "temperature": temp,  
                    "dt": dt, 
                    "atom_counts": atom_counts, 
                    "lattice": lattice
                }   
                f_h5.attrs['metadata'] = json.dumps(cond)
                
                if verbose:
                    print(f"'{out_file}' created successfully.")
                    
        if 'num_frames' in locals() and verbose:
            print(f"Successfully processed {num_frames} frames.")
               
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        
    except ValueError as e:
        print(f"Data Error: {e}")
        
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


@monitor_performance
def parse_lammps(lammps_data:str,
                 lammps_dump:str,
                 atom_style_data:str,
                 atom_style_dump:str,
                 atom_symbols:dict,
                 temp:float,
                 dt:float=1.0,
                 label:str=None,
                 chunk_size=5000,
                 verbose:bool=True) -> None:
    """Parses a LAMMPS trajectory using MDAnalysis and saves data to HDF5 files.

    This function processes large LAMMPS trajectories by leveraging MDAnalysis
    for efficient file indexing and reading. The trajectory is read in chunks
    to manage memory usage. The core calculation involves converting atomic
    positions to unwrapped fractional coordinates to correctly handle
    periodic boundaries.

    For each chemical symbol, a separate HDF5 file (`TRAJ_SYMBOL_LABEL.h5`) is
    created. Each file contains 'positions' and 'forces' datasets, along with
    simulation metadata (lattice, temperature, etc.) stored as an attribute.
    This method avoids creating intermediate files for improved I/O performance.

    Args:
        lammps_data (str):
            Path to the LAMMPS data file (topology).
        lammps_dump (str):
            Path to the LAMMPS dump file (trajectory).
        atom_style_data (str):
            Atom style string for the LAMMPS data file.
        atom_style_dump (str):
            Atom style string for the LAMMPS dump file.
        atom_symbols (dict):
            A dictionary mapping atom type IDs (int) to chemical symbols (str),
            e.g., {1: 'Ti', 2: 'O'}.
        temp (float):
            Simulation temperature in Kelvin.
        dt (float, optional):
            Timestep in femtoseconds. Defaults to 1.0.
        label (str, optional):
            A custom suffix for the output HDF5 filenames. Defaults to None.
        chunk_size (int, optional):
            The number of frames to process in memory at once. A larger size
            can be faster but uses more RAM. Defaults to 5000.
        verbose (bool, optional):
            Verbosity tag. Defaults to True.

    Returns:
        None: The function saves one or more HDF5 files to disk.

    Raises:
        FileNotFoundError: If a LAMMPS input file is not found.
        ValueError: If an atom type is missing from `atom_symbols` or if
                    force data is not found in the trajectory.
    """
    h5_files = {}
    h5_datasets = {}

    try:
        u = mda.Universe(
            lammps_data,
            topology_format="DATA",
            atom_style=atom_style_data
        )
        u.load_new(
            lammps_dump,
            format="LAMMPSDUMP",
            atom_style=atom_style_dump,
            dt=dt/1000    # fs to ps
        )

        lattice = u.trajectory[0].triclinic_dimensions.tolist()
        inv_lattice = np.linalg.inv(lattice)
        num_frames = len(u.trajectory)
        atom_types = np.array(u.atoms.types, dtype=int)
        unique_types = np.unique(atom_types)

        if chunk_size < 0:
            chunk_size = num_frames
        
        for type_id in unique_types:
            if not type_id in atom_symbols:
                raise ValueError(f"Atomic symbol for type {type_id} is not specified.")
            
        atom_indices = {atom_symbols[k]: np.where(atom_types == k)[0] for k in unique_types}
        atom_counts = {sym: len(indices) for sym, indices in atom_indices.items()}
            
        for type_id in unique_types:
            sym = atom_symbols[type_id]
            base_name = f"{sym}" if label is None else f"{sym}_{label}"
            out_file_name = f"TRAJ_{base_name}.h5"
            
            h5_file = h5py.File(out_file_name, "w")
            h5_files[sym] = h5_file # Store file handle
            
            h5_datasets[sym] = {
                'positions': h5_file.create_dataset(
                    "positions",
                    shape=(num_frames, atom_counts[sym], 3),
                    dtype=np.float64
                ),
                'forces': h5_file.create_dataset(
                    "forces",
                    shape=(num_frames, atom_counts[sym], 3),
                    dtype=np.float64
                )
            }
            
            cond = {
                "symbol": sym, 
                "nsw": num_frames,
                "temperature": temp,  
                "dt": dt, 
                "atom_counts": atom_counts, 
                "lattice": lattice
            }   
            h5_file.attrs['metadata'] = json.dumps(cond)

        pbar = tqdm(
            desc=f"{RED}{BOLD}progress{RESET}", 
            unit=" frames", 
            total=num_frames,
            ascii=False,
            bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'%(Fore.GREEN, Fore.RESET)
        )
        prev_positions = None
        for start_frame in range(0, num_frames, chunk_size):
            end_frame = min(start_frame + chunk_size, num_frames)
            chunk_iterator = u.trajectory[start_frame:end_frame]
            
            chunk_positions = []
            chunk_forces = []
            for atoms in chunk_iterator:
                chunk_positions.append(atoms.positions)
                try:
                    chunk_forces.append(atoms.forces)
                except AttributeError:
                    raise ValueError("Force data not found in the trajectory.")
                
            chunk_positions = np.array(chunk_positions, dtype=np.float64)
            chunk_forces = np.array(chunk_forces, dtype=np.float64)
            chunk_positions = chunk_positions @ inv_lattice # convert to fractional coords

            displacement = np.zeros_like(chunk_positions)
            if (end_frame - start_frame) > 1:
                displacement[1:, :] = np.diff(chunk_positions, axis=0)
                displacement -= np.round(displacement)
                displacement = np.cumsum(displacement, axis=0)
            
            if prev_positions is None:
                positions_init = chunk_positions[0]
            else:
                displacement_init = chunk_positions[0] - prev_positions
                displacement_init -= np.round(displacement_init)
                positions_init = prev_positions + displacement_init
            
            chunk_positions = positions_init + displacement
            prev_positions = chunk_positions[-1].copy()
            
            for sym, indices in atom_indices.items():
                pos_dset = h5_datasets[sym]['positions']
                force_dset = h5_datasets[sym]['forces']
                
                pos_dset[start_frame:end_frame] = chunk_positions[:, indices, :]
                force_dset[start_frame:end_frame] = chunk_forces[:, indices, :]
            
            pbar.update(end_frame - start_frame)
        pbar.close()
        
        if verbose:
            for sym in atom_counts.keys():
                base_name = f"{sym}" if label is None else f"{sym}_{label}"
                print(f"'TRAJ_{base_name}.h5' created successfully.")
            print(f"Successfully processed {num_frames} frames.")
            
    except FileNotFoundError as e:
        print(f"Error: The file {e.filename} was not found.")
        
    except ValueError as e:
        print(f"Data Error: {e}")
        
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        
    finally:
        for sym, h5_file in h5_files.items():
            if h5_file:
                h5_file.close()
                

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
                desc=f"{RED}{BOLD}progress{RESET}", 
                unit=" frames", 
                total=total_frames,
                ascii=False,
                bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'%(Fore.GREEN, Fore.RESET)
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


class Site:
    """
    Analyzes a crystal structure to find inequivalent atomic sites and hopping paths.

    This class reads a standard crystallographic structure file (e.g., POSCAR, cif),
    identifies the symmetrically inequivalent sites for a specified chemical
    element, and calculates all unique nearest-neighbor hopping paths up to a
    given cutoff radius. It is useful for setting up kinetic Monte Carlo or
    diffusion simulations.

    Args:
        structure (str): 
            Path to the crystallographic structure file.
        symbol (str): 
            The atomic symbol of the diffusing species to analyze.
        format (str, optional): 
            The format of the structure file. If None, ASE
            will attempt to automatically determine the format. Defaults to None.
        rmax (float, optional): 
            The maximum distance (in Angstroms) to consider
            when searching for neighbor sites for hopping paths. Defaults to 3.25.
        eps (float, optional): 
            A small tolerance value used for distance
            comparisons and finding atoms at specific coordinates. Defaults to 1.0e-3.
        verbose (bool, optional): 
            If True, a summary of the analysis will be
            printed to the console upon initialization. Defaults to False.

    Attributes:
        structure (pymatgen.core.structure.Structure): 
            The crystal structure as a pymatgen Structure object.
        symbol (str):
            The atomic symbol of the diffusing species to analyze.
        path (list[dict]): 
            A list of dictionaries, where each dictionary describes a unique hopping path.
        path_name (list[str]): 
            A list of unique names generated for each path (e.g., 'A1', 'A2', 'B1').
        site_name (list[str]): 
            A list of names for the inequivalent sites (e.g., 'site1', 'site2').
        lattice_sites (list[dict]): 
            A list containing information about each atomic site of the specified symbol, 
            including its fractional and Cartesian coordinates.
        lattice_parameter (np.ndarray): 
            The 3x3 lattice matrix of the crystal structure.
    """
    def __init__(self,
                 structure: str,
                 symbol: str,
                 format: str = None,
                 rmax: float = 3.25,
                 eps: float = 1.0e-3,
                 verbose: bool = False):
        
        self.structure_file = structure
        self.symbol = symbol
        self.format = format
        self.rmax = rmax
        self.eps = eps
        self.verbose = verbose
        
        if not os.path.isfile(self.structure_file):
            raise FileNotFoundError(f"Error: Input file '{self.structure_file}' not found.")
        
        try:
            if self.format is None:
                atoms = read(self.structure_file)
            else:
                atoms = read(self.structure_file, format=self.format)
            self.structure = AseAtomsAdaptor.get_structure(atoms)
        except Exception as e:
            raise IOError(f"Failed to read or convert '{self.structure_file}'. Error: {e}")
        
        if not any(site.specie.symbol == self.symbol for site in self.structure):
            raise ValueError(f"Error: Symbol '{self.symbol}' not found in '{self.structure_file}'.")
        
        self.path = []
        self.path_name = []
        self.site_name = None
        self.lattice_sites = None
        self.lattice_parameter = self.structure.lattice.matrix
        self.path_unknown = []
        
        self.find_hopping_path()
        
        if self.verbose:
            self.summary()
            
    def find_hopping_path(self):
        """
        Identifies inequivalent sites and calculates all unique hopping paths.

        This method is the core of the analysis. It first uses space group
        symmetry to find all symmetrically distinct sites for the given chemical
        symbol. Then, for a representative of each inequivalent site, it finds
        all neighbors within the `rmax` cutoff distance and consolidates them
        into a unique set of hopping paths, counting the coordination number (z)
        for each path.
        """
        sga = SpacegroupAnalyzer(self.structure)
        sym_structure = sga.get_symmetrized_structure()
        non_eq_sites = [
            site_group for site_group in sym_structure.equivalent_sites
            if site_group[0].specie.symbol == self.symbol
        ]
        
        index = []
        for sites in non_eq_sites:
            index_sites = []
            for site in sites:
                coords = site.coords
                for i, _site in enumerate(self.structure.sites):
                    if np.linalg.norm(coords - _site.coords) < self.eps:
                        index_sites.append(i)
            index.append(index_sites)
        
        self.site_name = [f"site{i+1}" for i in range(len(index))]
        
        self.lattice_sites = []
        min_idx = min(map(min, index)) if index else 0
        max_idx = max(map(max, index)) if index else -1
        for i in range(min_idx, max_idx + 1):
            for j, index_j in enumerate(index):
                if i in index_j:
                    site_i = j + 1
                    break
            point = {
                'site': f"site{site_i}",
                'coord': self.structure[i].frac_coords,
                'coord_cart': self.structure[i].coords
            }
            self.lattice_sites.append(point)
            
        nn_finder = VoronoiNN(tol=self.eps)
        self.path, self.path_name = [], []
        for i, idx in enumerate([index_i[0] for index_i in index]):
            paths_idx = []
            distances = np.array([], dtype=float)
            site_init = f"site{i+1}"
            neighbors = nn_finder.get_nn_info(self.structure, idx)
            neighbors = [
                n for n in neighbors if n['site'].specie.symbol == self.symbol
            ]
            
            for neighbor in neighbors:
                distance = self.structure[idx].distance(neighbor['site'])
                if distance < self.rmax:
                    for j, index_j in enumerate(index):
                        if neighbor['site_index'] in index_j:
                            site_final = j + 1
                            break
                    site_final = f"site{site_final}"
                    path_index = np.where(abs(distances - distance) < self.eps)[0]
                    if len(path_index) == 0:
                        path = {
                            'site_init': site_init,
                            'site_final': site_final,
                            'distance': float(distance),
                            'z': 1,
                            'coord_init': self.structure[idx].frac_coords,
                            'coord_final': neighbor['site'].frac_coords
                        }
                        paths_idx.append(path)
                        distances = np.append(distances, distance)
                        self.path_name.append(f"{chr(i+65)}{len(paths_idx)}")
                    else:
                        paths_idx[path_index[0]]['z'] += 1
            self.path += paths_idx
            
        self.path = sorted(self.path, key=lambda x: (x['site_init'], x['distance']))
        self.path_name = sorted(self.path_name)
        for path, name in zip(self.path, self.path_name):
            path['name'] = name
    
    def summary(self):
        """
        Prints a formatted summary of the site and path analysis to the console.

        The summary includes the number of inequivalent sites and paths found,
        followed by a detailed table of each unique path, including its name,
        initial and final sites, distance, coordination number, and fractional
        coordinates.
        """
        print(f"Number of inequivalent sites for {self.symbol} : {len(self.site_name)}")
        print(f"Number of inequivalent paths for {self.symbol} : {len(self.path_name)} (Rmax = {self.rmax:.2f} Å)")
        print('')
        print('Path information')
        headers = ['name', 'init', 'final', 'a(Ang)', 'z', 'coord_init', 'coord_final']
        data = [
            [
                path['name'], 
                path['site_init'], 
                path['site_final'], 
                f"{path['distance']:.5f}", 
                path['z'],
                f"[{path['coord_init'][0]:.5f} {path['coord_init'][1]:.5f} {path['coord_init'][2]:.5f}]", 
                f"[{path['coord_final'][0]:.5f} {path['coord_final'][1]:.5f} {path['coord_final'][2]:.5f}]"
            ] for path in self.path
        ]
        print(tabulate(data, headers=headers, tablefmt="simple"))


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