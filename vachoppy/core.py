"""
vachoppy.core
=============

Provides the primary user-facing classes and functions for setting up and
running a complete vacancy diffusion analysis with `VacHopPy`.

Core Components
---------------
- **parse_md / parse_lammps**: Functions to convert raw MD trajectories from
  common simulation packages (e.g., VASP, LAMMPS) into the standardized HDF5
  format required for analysis.
- **Site**: A fundamental class for analyzing a crystal structure to identify
  symmetrically inequivalent sites and potential hopping paths. This is the
  first object to create in an analysis workflow.
- **Calculator**: A unified function that simplifies the setup of an analysis.
  It intelligently handles both single and multiple trajectory files and returns
  a configured `CalculatorEnsemble` instance, ready for computation.

Typical Workflow
----------------
A standard analysis involves two main stages: data preparation and calculation.

**1. Data Preparation (CLI)**

First, convert your raw MD trajectory into the required HDF5 format. This is
typically done once using the command-line interface:

.. code-block:: bash

    vachoppy convert path/to/vasprun.xml 2000 1.0 --label 2000K

**2. Analysis Workflow (Python)**

Once the HDF5 files are ready, the analysis is performed in Python:

.. code-block:: python

    from vachoppy.core import Site, Calculator

    # a. Define sites and paths from the crystal structure
    site_info = Site(path_structure="path/to/POSCAR", symbol="O")

    # b. Set up the analysis for the HDF5 data
    #    (This handles both single files and directories automatically)
    calc = Calculator(
        path_traj="path/to/hdf5_files/",
        site=site_info,
        t_interval=0.1  # Coarse-graining time in ps
    )

    # c. Run the full analysis pipeline in parallel
    calc.calculate()

    # d. View results and generate plots
    calc.summary()
    calc.plot_D()
"""

from __future__ import annotations

__all__ =['parse_md', 'parse_lammps', 'Site', 'Calculator']

import os
import h5py
import json
import itertools
import numpy as np

from pathlib import Path
from typing import Union
from tqdm.auto import tqdm
from tabulate import tabulate
from ase.io import read, iread
from numpy.typing import DTypeLike

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from vachoppy.utils import monitor_performance
from vachoppy.vibration import *
from vachoppy.trajectory import *


@monitor_performance
def parse_md(filename: str,
             format: str,
             temperature: float,
             dt: float = 1.0,
             label: str = None,
             chunk_size: int = 5000,
             dtype: DTypeLike = np.float64,
             verbose: bool = True) -> None:
    """Parses a molecular dynamics (MD) trajectory and saves the data to HDF5 files.

    This function provides a memory-efficient way to process large MD trajectory
    files supported by ASE (e.g., VASP outputs, extxyz). It reads and processes
    the trajectory in chunks, converting atomic positions to unwrapped fractional
    coordinates to correctly handle periodic boundary conditions.

    Key Features:
    - Processes large files in memory-efficient chunks using ASE's iread.
    - Unwraps atomic coordinates across periodic boundaries for accurate analysis.
    - Separates data by chemical symbol into distinct HDF5 files.
    - Includes essential simulation metadata (lattice, temp, etc.) in each file.

    Args:
        filename (str):
            Path to the input MD trajectory file (e.g., 'vasprun.xml').
        file_format (str):
            The file format string recognized by ASE (e.g., 'vasp-xml', 'extxyz').
        temperature (float):
            Simulation temperature in Kelvin, to be stored as metadata.
        dt (float, optional):
            Timestep in femtoseconds (fs). Defaults to 1.0.
        label (str | None, optional):
            A custom suffix for output filenames (e.g., 'TRAJ_SYMBOL_LABEL.h5').
            Defaults to None.
        chunk_size (int, optional):
            Number of frames to read into memory per chunk. Larger values may
            improve speed but increase RAM usage. Defaults to 5000.
        dtype (DTypeLike, optional):
            NumPy data type for storing positions and forces, affecting
            precision and file size. Defaults to np.float64.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Returns:
        None:
            This function does not return a value; it writes one or more HDF5
            files to disk in the current directory.

    Raises:
        FileNotFoundError:
            If the input trajectory file is not found.
        ValueError:
            If force data is missing from the trajectory frames.

    Examples:
        >>> parse_md(
        ...     filename='path/to/vasprun.xml',
        ...     file_format='vasp-xml',
        ...     temperature=2000.0,
        ...     dt=2.0,
        ...     label='2000K'
        ... )
        # This will create files like 'TRAJ_Ti_2000K.h5' and 'TRAJ_O_2000K.h5'.
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
        
        pbar = tqdm(desc=f"Extract Trejaectory", unit=" frames")
        
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
                    dtype=dtype,
                    shape=(len(chunk_atoms), len(indices), 3)
                )
                temp_force_memmap = np.lib.format.open_memmap(
                    temp_force_path,
                    mode="w+",
                    dtype=dtype,
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
                    dtype=dtype
                )
                force_dataset = f_h5.create_dataset(
                    "forces", 
                    shape=(num_frames, atom_counts[sym], 3), 
                    dtype=dtype
                )
                
                current_pos = 0
                for temp_path in temp_files['pos'][sym]:
                    temp_data = np.load(temp_path)
                    chunk_len = len(temp_data)
                    pos_dataset[current_pos : current_pos + chunk_len] = temp_data
                    current_pos += chunk_len
                    os.remove(temp_path)

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
                    "temperature": temperature,  
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
def parse_lammps(lammps_data: str,
                 lammps_dump: str,
                 atom_style_data: str,
                 atom_style_dump: str,
                 atom_symbols: dict[int, str],
                 temperature: float,
                 dt: float = 1.0,
                 label: str | None = None,
                 chunk_size: int = 5000,
                 dtype: DTypeLike = np.float64,
                 verbose: bool = True) -> None:
    """Parses a LAMMPS trajectory using MDAnalysis and saves data to HDF5 files.

    This function leverages the MDAnalysis library to efficiently process large
    LAMMPS data and dump files. It reads the trajectory in chunks to maintain
    low memory usage while converting atomic positions to unwrapped fractional
    coordinates for accurate periodic boundary handling.

    Key Features:
    - Efficiently processes large LAMMPS dump files using MDAnalysis.
    - Unwraps atomic coordinates across periodic boundaries for accurate analysis.
    - Separates atom data by chemical symbol into distinct HDF5 files.
    - Stores simulation parameters as metadata in each output file.

    Args:
        lammps_data (str):
            Path to the LAMMPS data file containing topology and box info.
        lammps_dump (str):
            Path to the LAMMPS dump file containing the trajectory.
        atom_style_data (str):
            Atom style string for the LAMMPS data file (e.g., 'full', 'atomic').
        atom_style_dump (str):
            Atom style string for the LAMMPS dump file (e.g., 'atomic').
        atom_symbols (dict[int, str]):
            A dictionary mapping atom type IDs (int) to chemical symbols (str),
            e.g., {1: 'Ti', 2: 'O'}.
        temperature (float):
            Simulation temperature in Kelvin, to be stored as metadata.
        dt (float, optional):
            Timestep in femtoseconds (fs). Defaults to 1.0.
        label (str | None, optional):
            A custom suffix for output filenames (e.g., 'TRAJ_SYMBOL_LABEL.h5').
            Defaults to None.
        chunk_size (int, optional):
            Number of frames to read into memory per chunk. Larger values may
            improve speed but increase RAM usage. Defaults to 5000.
        dtype (DTypeLike, optional):
            NumPy data type for storing positions and forces, affecting
            precision and file size. Defaults to np.float64.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Returns:
        None:
            This function does not return a value; it writes one or more HDF5
            files to disk in the current directory.

    Raises:
        ImportError:
            If the MDAnalysis library is not installed.
        FileNotFoundError:
            If a specified LAMMPS input file is not found.
        ValueError:
            If an atom type in the trajectory is missing from `atom_symbols`,
            or if force data is not found in the dump file.

    Examples:
        >>> parse_lammps(
        ...     lammps_data='tio2.data',
        ...     lammps_dump='tio2.dump',
        ...     atom_style_data='id type x y z',
        ...     atom_style_dump='id type x y z fx fy fz',
        ...     atom_symbols={1: 'Ti', 2: 'O'},
        ...     temperature=1200.0,
        ...     dt=1.0,
        ...     label='1200K'
        ... )
        # This will create files like 'TRAJ_Ti_1200K.h5' and 'TRAJ_O_1200K.h5'.
    """
    
    try:
        import MDAnalysis as mda
    except ImportError:
        raise ImportError(
            "This feature requires the MDAnalysis library. "
            "Please install it using 'pip install MDAnalysis'."
        )
    
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
                    dtype=dtype
                ),
                'forces': h5_file.create_dataset(
                    "forces",
                    shape=(num_frames, atom_counts[sym], 3),
                    dtype=dtype
                )
            }
            
            cond = {
                "symbol": sym, 
                "nsw": num_frames,
                "temperature": temperature,  
                "dt": dt, 
                "atom_counts": atom_counts, 
                "lattice": lattice
            }   
            h5_file.attrs['metadata'] = json.dumps(cond)

        pbar = tqdm(
            desc=f"Extract Trajectory", 
            unit=" frames", 
            total=num_frames,
            ascii=True,
            bar_format='{l_bar}{bar:30}{r_bar}'
        )
        prev_positions = None
        for start_frame in range(0, num_frames, chunk_size):
            end_frame = min(start_frame + chunk_size, num_frames)
            chunk_iterator = u.trajectory[start_frame:end_frame]
            
            chunk_positions = []
            chunk_forces = []
            for atoms in chunk_iterator:
                chunk_positions.append(atoms.positions.astype(dtype))
                try:
                    chunk_forces.append(atoms.forces.astype(dtype))
                except AttributeError:
                    raise ValueError("Force data not found in the trajectory.")
                
            chunk_positions = np.array(chunk_positions, dtype=dtype)
            chunk_forces = np.array(chunk_forces, dtype=dtype)
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
                

class Site:
    """Analyzes a crystal structure to find inequivalent sites and hopping paths.

    This class reads a standard crystallographic structure file (e.g., POSCAR,
    cif), identifies symmetrically inequivalent sites for a specified element,
    and calculates all unique nearest-neighbor hopping paths up to a given
    cutoff radius. It is a foundational tool for setting up kinetic Monte Carlo
    or diffusion analyses.

    Args:
        path_structure (str):
            Path to the crystallographic structure file.
        symbol (str):
            The atomic symbol of the diffusing species to analyze.
        structure_format (str | None, optional):
            Format of the structure file, as recognized by ASE. If None, ASE
            will attempt to determine the format automatically. Defaults to None.
        rmax (float, optional):
            The maximum distance (in Angstroms) to search for neighbors when
            defining hopping paths. Defaults to 3.25.
        eps (float, optional):
            A small tolerance value for distance comparisons and identifying
            atomic coordinates. Defaults to 1.0e-3.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Attributes:
        structure (pymatgen.core.structure.Structure):
            The crystal structure represented as a pymatgen Structure object.
        symbol (str):
            The atomic symbol of the diffusing species being analyzed.
        path (list[dict]):
            A list of dictionaries, each describing a unique hopping path with
            details like distance, coordination number (z), and coordinates.
        path_name (list[str]):
            A list of unique names for each hopping path (e.g., 'A1', 'B1').
        site_name (list[str]):
            A list of names for the inequivalent sites (e.g., 'site1', 'site2').
        lattice_sites (list[dict]):
            A list detailing each atomic site of the specified symbol,
            including its fractional and Cartesian coordinates.
        lattice_parameter (numpy.ndarray):
            The 3x3 lattice matrix of the crystal structure.

    Raises:
        FileNotFoundError:
            If the specified structure file does not exist.
        IOError:
            If the structure file cannot be read or converted by ASE/pymatgen.
        ValueError:
            If the specified `symbol` is not found in the structure.

    Examples:
        >>> tio2_site = Site(
        ...     path_structure='TiO2_rutile.cif',
        ...     symbol='Ti',
        ...     rmax=3.5,
        ...     verbose=True
        ... )
        # This will create the Site object and print a summary of the analysis.
        # Hopping path data can then be accessed via tio2_site.path
    """
    def __init__(self,
                 path_structure: str,
                 symbol: str,
                 format: str = None,
                 rmax: float = 3.25,
                 eps: float = 1.0e-3,
                 verbose: bool = False):
        
        self.path_structure = path_structure
        self.symbol = symbol
        self.format = format
        self.rmax = rmax
        self.eps = eps
        self.verbose = verbose
        
        if not os.path.isfile(self.path_structure):
            raise FileNotFoundError(f"Error: Input file '{self.path_structure}' not found.")
        
        try:
            if self.format is None:
                atoms = read(self.path_structure)
            else:
                atoms = read(self.path_structure, format=self.format)
            self.structure = AseAtomsAdaptor.get_structure(atoms)
        except Exception as e:
            raise IOError(f"Failed to read or convert '{self.path_structure}'. Error: {e}")
        
        if not any(site.specie.symbol == self.symbol for site in self.structure):
            raise ValueError(f"Error: Symbol '{self.symbol}' not found in '{self.path_structure}'.")
        
        self.path = []
        self.path_name = []
        self.site_name = None
        self.lattice_sites = None
        self.lattice_parameter = self.structure.lattice.matrix
        self.path_unknown = []
        
        self._find_hopping_path()
        
        if self.verbose:
            self.summary()
            
    def _find_hopping_path(self) -> None:
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
    
    def summary(self) -> None:
        """
        Prints a formatted summary of the site and path analysis to the console.

        The summary includes the number of inequivalent sites and paths found,
        followed by a detailed table of each unique path, including its name,
        initial and final sites, distance, coordination number, and fractional
        coordinates.
        """

        print("\n" + "=" * 100)
        # print(f"{' ' * 40}Site Analysis Summary")
        print(f'  Structure File: {self.path_structure}')
        print("=" * 100)
        print("[Structure Information]")
        # print(f"    - Structure File        : {self.path_structure}")
        print(f"    - Structure Composition : {str(self.structure.composition)}")
        print(f"    - Lattice Vectors (Ang) :")
        for vector in self.structure.lattice.matrix:
            print(" "*8 + f"[{vector[0]:>9.5f}, {vector[1]:>9.5f}, {vector[2]:>9.5f}]")
        
        print("\n" + "[Hopping Path Information]")
        headers = ['Name', 'Init Site', 'Final Site', 'a (Å)', 'z',
                   'Initial Coord (Frac)', 'Final Coord (Frac)']
        data = [
            [
                path['name'], 
                path['site_init'], 
                path['site_final'], 
                f"{path['distance']:.4f}", 
                path['z'],
                f"[{path['coord_init'][0]:.4f}, {path['coord_init'][1]:.4f}, {path['coord_init'][2]:.4f}]", 
                f"[{path['coord_final'][0]:.4f}, {path['coord_final'][1]:.4f}, {path['coord_final'][2]:.4f}]"
            ] for path in self.path
        ]
        
        print(f"    - Diffusing Symbol   : {self.symbol}")
        print(f"    - Inequivalent Sites : {len(self.site_name)} found")
        print(f"    - Inequivalent Paths : {len(self.path_name)} found (with Rmax = {self.rmax:.2f} Å)\n")
        
        if not data:
            print("No hopping paths were found within the specified rmax.")
        else:
            print(tabulate(data, headers=headers, tablefmt="simple", stralign='left', numalign='left'))
        print("=" * 100 + "\n")
        

def Calculator(path_traj: str,
               site: Site,
               *,
               t_interval: float | None = None,
               **kwargs) -> CalculatorEnsemble:
    """Initializes and configures a CalculatorEnsemble for analysis.

    This function serves as the primary user entry point for setting up a
    calculation. It intelligently handles both single HDF5 trajectory files and
    directories containing multiple files by using the `TrajectoryBundle` class.

    If `t_interval` is not provided, this function will automatically estimate
    an optimal value based on the mean vibration frequency of the atoms in a
    representative trajectory.

    Args:
        path_traj (str):
            Path to a single HDF5 file or a root directory to search for files.
        site (Site):
            An initialized `Site` object containing lattice and hopping path data.
        t_interval (float | None, optional):
            The time interval in picoseconds (ps) for analysis. If None, the
            interval is automatically estimated. Defaults to None.
        **kwargs:
            Additional keyword arguments passed to underlying classes.
            Accepted arguments include:
            - `prefix` (str, optional): File prefix for directory scans.
            Defaults to "TRAJ".
            - `depth` (int, optional): Directory search depth. Defaults to 2.
            - `sampling_size` (int, optional): Frames for t_interval
            estimation. Defaults to 5000.
            - `use_incomplete_encounter` (bool, optional): Flag for Encounter
            analysis. Defaults to True.
            - `eps` (float, optional): Tolerance for float comparisons.
            Defaults to 1.0e-3.
            - `verbose` (bool, optional): Verbosity flag. Defaults to True.

    Returns:
        CalculatorEnsemble:
            An initialized and ready-to-use CalculatorEnsemble instance.

    Raises:
        FileNotFoundError:
            If `path_traj` does not exist or no valid files are found.
        ValueError:
            If `t_interval` cannot be estimated due to zero mean frequency.

    Examples:
        >>> # Analyze a single trajectory file
        >>> site_info = Site("POSCAR", symbol="O")
        >>> calc = Calculator("TRAJ_O.h5", site=site_info)

        >>> # Analyze a directory of trajectories
        >>> calc = Calculator("trajectories/", site=site_info)
    """
    
    p = Path(path_traj)
    if not p.exists():
        raise FileNotFoundError(f"Error: The path '{path_traj}' was not found.")
    
    # Helper function for t_interval estimation
    def _get_t_interval(path_traj: str) -> float:
        vib_init_keys = ['sampling_size', 'filter_high_freq']
        vib_init_kwargs = {key: kwargs.get(key) for key in vib_init_keys if key in kwargs}
        vib_params = {
            'path_traj': path_traj, 
            'site': site,
            'verbose': False
        }

        vib_params.update(vib_init_kwargs)
        vib = Vibration(**vib_params)
        vib.calculate()
        if vib.mean_frequency > 0:
            estimated_interval = 1 / vib.mean_frequency
            print(" "*13 + f"-> t_interval : {estimated_interval:.3f} ps")
            return estimated_interval
        else:
            raise ValueError(f"Could not estimate t_interval from '{path_traj}' as mean frequency is zero.")
        
    # Helper function to get dt
    def _get_dt(traj: str) -> float:
        with h5py.File(traj, 'r') as f:
            return json.loads(f.attrs['metadata']).get('dt')
    
    bundle = None   
    representative_traj = ""
    bundle_init_keys = ['prefix', 'depth', 'eps', 'verbose']
    bundle_init_kwargs = {key: kwargs.get(key) for key in bundle_init_keys if key in kwargs}

    if p.is_file():
        representative_traj = str(p.resolve())
    elif p.is_dir():
        bundle = TrajectoryBundle(path_traj=path_traj, symbol=site.symbol, **bundle_init_kwargs)
        if not bundle.traj or not bundle.traj[0]:
            raise FileNotFoundError(f"No valid trajectory files found in '{path_traj}' to use for analysis.")
        representative_traj = bundle.traj[0][0]
        
    dt_fs = _get_dt(representative_traj)
    dt_ps = dt_fs / 1000.0
    
    if t_interval is None:
        print("="*68)
        print(" "*19 + "Automatic t_interval Estimation")
        print("="*68)
        if p.is_dir():
            t_interval_list = []
            for i, temp in enumerate(bundle.temperatures):
                file_path = Path(bundle.traj[i][0])
                short_path = os.path.join(*file_path.parts[-bundle.depth:])
                print(f"  [{temp} K] Estimating from {short_path}")
                t_interval_list.append(_get_t_interval(bundle.traj[i][0]))
            t_interval = np.mean(t_interval_list) 
        elif p.is_file():
            print(f"  Estimating from {p.name}")
            t_interval = _get_t_interval(representative_traj)
        print("="*68)
    
        original_t_interval = t_interval            
        num_dt_steps = round(original_t_interval / dt_ps)
        adjusted_t_interval = num_dt_steps * dt_ps

        print(" "*8 + "Adjusting t_interval to the nearest multiple of dt")
        print("="*68)
        print(f"    - dt                  : {dt_ps:.4f} ps")
        print(f"    - Original t_interval : {original_t_interval:.4f} ps")
        print(f"    - Adjusted t_interval : {adjusted_t_interval:.4f} ps ({num_dt_steps} frames)")
        t_interval = adjusted_t_interval
        print("="*68)
    
    calc_keys = ['prefix', 'depth', 'use_incomplete_encounter', 'eps', 'verbose']
    calc_kwargs = {key: kwargs.get(key) for key in calc_keys if key in kwargs}
    
    return CalculatorEnsemble(path_traj=path_traj, site=site, t_interval=t_interval, **calc_kwargs)


