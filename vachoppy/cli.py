import os
import numpy as np
from pathlib import Path
from typing import List, Union, Optional

from vachoppy.core import *
from vachoppy.einstein import *
from vachoppy.vibration import *
from vachoppy.frequency import *
from vachoppy.trajectory import *
from vachoppy.utils import *
from vachoppy.fingerprint import *

   
@ monitor_performance
def cli_trajectory(path_traj: str,
                   path_structure: str,
                   symbol: str,
                   verbose: bool = True,
                   **kwargs) -> None:
    """[CLI Method] This method is for vacancy trajectory analysis."""
    STEP_FLAG = 1
    p = Path(path_traj)
    if not p.is_file():
        raise ValueError(f"Error: '{path_traj}' is not a regular file.")
    
    site_keys = ['format', 'rmax', 'eps']
    site_kwargs = {key: kwargs[key] for key in site_keys if key in kwargs}
    site_kwargs['verbose'] = False
    site = Site(path_structure, symbol, **site_kwargs)
    
    calc_keys = ['t_interval', 'sampling_size', 'eps', 'use_incomplete_encounter']
    calc_kwargs = {key: kwargs[key] for key in calc_keys if key in kwargs}
    calc_kwargs['verbose'] = False
    
    t_interval = kwargs.get('t_interval', None)
    if t_interval is None:
        print(f"[STEP{STEP_FLAG}] Automatic t_interval Estimation:"); STEP_FLAG += 1
        
    calc = Calculator(path_traj, site, **calc_kwargs)
    
    print(f"\n\n[STEP{STEP_FLAG}] Identifying Vacancy Trajectory:"); STEP_FLAG += 1
    calc.calculate(verbose=False)
    
    if not calc.calculators[0].hopping_sequence:
        print("\n[INFO] No hopping events found. Stopping analysis.\n")
        return
    
    print(f"\n\n[STEP{STEP_FLAG}] Summary of Hopping Paths:"); STEP_FLAG += 1
    calc.calculators[0].show_hopping_paths()
    
    print(f"\n[STEP{STEP_FLAG}] Summary of Hopping Histories:"); STEP_FLAG += 1
    calc.calculators[0].show_hopping_history()
    filename = 'trajectory.json'
    calc.calculators[0].save_trajectory(filename=filename)
    
    plot_keys = ['vacancy_indices', 'filename']
    plot_kwargs = {key: kwargs[key] for key in plot_keys if key in kwargs}
    plot_kwargs['unwrap'] = kwargs.get('unwrap', True)
    plot_kwargs['save'] = kwargs.get('save_plot', True)
    plot_kwargs['disp'] = False
    calc.calculators[0].plot_vacancy_trajectory(np.arange(calc.num_vacancies), **plot_kwargs)
    
    print(f"Results are saved in '{filename}'.")
    print(f"Trajectory is saved in 'trajectory.html'.\n")
    
    try:
        import platform
        import subprocess
        
        if platform.system() == 'Windows':
            os.startfile('trajectory.html')
        elif platform.system() == 'Darwin':  # macOS
            subprocess.run(['open', 'trajectory.html'])
        else:  # Linux
            subprocess.run(['xdg-open', 'trajectory.html'])
            
    except Exception as e:
        print("Could not open the image automatically. " + 
                f"Please open '{'trajectory.png'}'")
        

@ monitor_performance
def cli_analyze(path_traj: str,
                path_structure: str,
                symbol: str,
                neb_csv: str = None,
                dir_imgs: str = 'imgs',
                xyz: bool = False,
                disp: bool = True,
                verbose: bool = True,
                **kwargs) -> None:
    """[CLI Method] This method is for hopping parameter extraction."""
    STEP_FLAG = 1
    p = Path(path_traj)
    
    site_keys = ['format', 'rmax', 'eps']
    site_kwargs = {key: kwargs[key] for key in site_keys if key in kwargs}
    site_kwargs['verbose'] = False
    site = Site(path_structure, symbol, **site_kwargs)
    
    calc_keys = ['depth', 't_interval', 'sampling_size', 'eps', 'use_incomplete_encounter']
    calc_kwargs = {key: kwargs[key] for key in calc_keys if key in kwargs}
    calc_kwargs['verbose'] = False
    
    t_interval = kwargs.get('t_interval', None)
    if t_interval is None:
        print(f"[STEP{STEP_FLAG}] Automatic t_interval Estimation:"); STEP_FLAG += 1
    calc = Calculator(path_traj, site, **calc_kwargs)
    if t_interval is None: print('\n')
    
    print(f"[STEP{STEP_FLAG}] Vacancy Trajectory Identification:"); STEP_FLAG += 1
    n_jobs = kwargs.get('n_jobs', -1)
    calc.calculate(n_jobs=n_jobs, verbose=False)
    print('\n')
    filename: str = "parameters.json"
    calc.save_parameters(filename)
    print(f"[STEP{STEP_FLAG}] Hopping Parameter Calculation:"); STEP_FLAG += 1
    calc.summary()
    
    if neb_csv is not None:
        print(f"\n\n[STEP{STEP_FLAG}] Attempt Frequency Calculation:"); STEP_FLAG += 1
        calc.calculate_attempt_frequency(neb_csv=neb_csv, filename=filename)
        calc.attempt_frequency.summary()
    
    if xyz:
        print(f"\n[STEP{STEP_FLAG}] Decomposing Diffusivity into xyz-Components:"); STEP_FLAG += 1
        calc.decompose_diffusivity(verbose=True)
        
    if  dir_imgs is not None:
        if not os.path.isdir(dir_imgs): os.makedirs(dir_imgs)
        calc.plot_counts(disp=disp, filename=os.path.join(dir_imgs, 'counts.png'))
        
        if len(calc.temperatures) > 1:
            calc.plot_D_rand(disp=disp, filename=os.path.join(dir_imgs, 'D_rand.png'))
            calc.plot_f(disp=disp, filename=os.path.join(dir_imgs, 'f.png'))
            calc.plot_D(disp=disp, filename=os.path.join(dir_imgs, 'D.png'))
            calc.plot_tau(disp=disp, filename=os.path.join(dir_imgs, 'tau.png'))
            calc.plot_a(disp=disp, filename=os.path.join(dir_imgs, 'a.png'))
            
            if neb_csv is not None:
                calc.plot_nu(disp=disp, filename=os.path.join(dir_imgs, 'nu.png'))
                calc.plot_z(disp=disp, filename=os.path.join(dir_imgs, 'z.png'))
            
        if xyz:
            calc.plot_msd_xyz(disp=disp, filename=os.path.join(dir_imgs, 'msd_xyz.png'))
            if len(calc.temperatures) > 1:
                calc.plot_D_xyz(disp=disp, filename=os.path.join(dir_imgs, 'D_xyz.png'))    
                
        else:
            print("\n[INFO] Skipping plots (e.g., Arrhenius plots), " + 
                "as they require data from more than one temperature.")
    
    print(f"\nResults are saved in '{filename}'.")      
    if dir_imgs is not None: print(f"Images are saved in '{dir_imgs}'.")
    print('')
    
    
@ monitor_performance
def cli_vibration(path_traj:str,
                   path_structure: str,
                   symbol: str,
                   dir_imgs: str = 'imgs',
                   verbose: bool = True,
                   **kwargs):
    """[CLI Method] This method is for atomic vibration analysis"""
    STEP_FLAG = 1
    p = Path(path_traj)
    if not p.is_file():
        raise ValueError(f"Error: '{path_traj}' is not a regular file.")
    
    site_keys = ['format', 'rmax', 'eps']
    site_kwargs = {key: kwargs[key] for key in site_keys if key in kwargs}
    site_kwargs['verbose'] = False
    site = Site(path_structure, symbol, **site_kwargs)
    
    vib_keys = ['sampling_size', 'filter_high_freq', 'verbose']
    vib_kwargs = {key: kwargs[key] for key in vib_keys if key in kwargs}
    sampling_size = kwargs.get('sampling_size', 5000)
    print(f"[STEP{STEP_FLAG}] Atomic Vibration Analysis (Using initial {sampling_size} frames):")
    vib = Vibration(path_traj, site, **vib_kwargs)
    
    cal_keys = ['n_jobs', 'jump_detection_radius']
    cal_kwargs = {key: kwargs[key] for key in cal_keys if key in kwargs}
    vib.calculate(**cal_kwargs)
    
    if dir_imgs is not None:
        if not os.path.isdir(dir_imgs): os.makedirs(dir_imgs)
        vib.plot_displacements(disp=True, filename=os.path.join(dir_imgs, 'displacement.png'))
        vib.plot_frequencies(disp=True, filename=os.path.join(dir_imgs, 'frequency.png'))
        print(f"Images are saved in '{dir_imgs}'.\n")
        
        
@ monitor_performance
def cli_msd(path_traj: str,
            symbol: str,
            segment_length: Optional[Union[float, List[float]]] = None,
            dir_imgs: str = 'imgs',
            verbose=True,
            **kwargs):
    """[CLI method] This method runs einstein.Einstein. """
    STEP_FLAG = 1
    p = Path(path_traj)
    
    if p.is_dir():
        ein_keys = ['skip', 'start', 'end', 'n_jobs', 'prefix', 'depth', 'eps']
        
        ein_kwargs = {key: kwargs[key] for key in ein_keys if key in kwargs}
        ein_kwargs['verbose'] = False
        ein = Einstein(path_traj, symbol, segment_length=segment_length, **ein_kwargs)
        
        print(f"[STEP{STEP_FLAG}] MSD Calculation Based on Einstein Relation:"); STEP_FLAG += 1
        n_jobs_val = kwargs.get('n_jobs', -1)
        ein.calculate(n_jobs=n_jobs_val)
        
        print(f"\n\n[STEP{STEP_FLAG}] Summary of MSD Analysis:"); STEP_FLAG += 1
        ein.summary()
        ein.save_parameters()
        
        if not os.path.isdir(dir_imgs): os.makedirs(dir_imgs)
        ein.plot_msd(disp=True, filename=os.path.join(dir_imgs, 'msd.png'))
        ein.plot_D(disp=True, filename=os.path.join(dir_imgs, 'D_atom.png'))
        ein.save_parameters()
        
        print(f"\nResults are saved in 'einstein.json'.")
        print(f"Images are saved in {dir_imgs}.\n")
        
    if p.is_file():
        ein_keys = ['skip', 'start', 'end']
        ein_kwargs = {key: kwargs[key] for key in ein_keys if key in kwargs}
        ein_kwargs['verbose'] = False
        ein = Einstein(path_traj, symbol, segment_length=segment_length, **ein_kwargs)
        ein.calculate()
        
        ein.summary()
        
        if not os.path.isdir(dir_imgs): os.makedirs(dir_imgs)
        _ = ein.plot_msd(disp=True, filename=os.path.join(dir_imgs, 'msd.png'))
        
        print(f"Images are saved in {dir_imgs}.\n")


@ monitor_performance
def cli_distance(path_traj: Union[str, List[str]],
                 t_interval: float,
                 reference_structure: str,
                 verbose: bool = True,
                 **kwargs):
    """[CLI method] This method displays changes in cosine distance along time."""
    dist_keys = ['Rmax', 'delta', 'sigma', 'dirac', 'atom_pairs',
                'n_jobs', 'window_size', 'threshold_std']
    dist_kwargs = {key: kwargs[key] for key in dist_keys if key in kwargs}
    dist_kwargs['verbose'] = True
    
    plot_cosine_distance(
        path_traj,
        t_interval,
        reference_structure,
        **dist_kwargs
    )
    print('')
    
    
@ monitor_performance
def cli_fingerprint(path_structure: str, 
                    verbose: bool = True,
                    **kwargs):
    """[CLI method] This method displays a fingerprint plot."""
    
    from ase.io import read
    from itertools import combinations_with_replacement

    atom_pairs = kwargs.get('atom_pairs', None)
    if atom_pairs is None:
        if verbose: print("Argument 'atom_pairs' not provided: Auto-generating all unique pairs...")
        try:
            atoms = read(path_structure)
            atom_species = sorted(list(set(atoms.get_chemical_symbols())))
            atom_pairs = list(combinations_with_replacement(atom_species, 2))
            if verbose: print(f"-> Generated pairs: {atom_pairs}\n")
        except Exception as e:
            raise IOError(f"Failed to read reference structure '{path_structure}' to auto-generate pairs. Error: {e}")

    fp_keys = ['Rmax', 'delta', 'sigma', 'dirac']
    fp_kwargs = {key: kwargs[key] for key in fp_keys if key in kwargs}
    fp_kwargs['verbose'] = True
    
    _ = get_fingerprint(path_structure=path_structure,
                        filename='fingerprint.txt',
                        atom_pairs=atom_pairs,
                        disp=True,
                        **fp_kwargs)
    print('')
    
    
@monitor_performance
def cli_convert(filename: str,
                format: str,
                temperature: float,
                dt: float = 1.0,
                **kwargs) -> None:
    """[CLI method] This method converts a MD result to HDF5 files"""
    
    if format != "lammps-dump-text":
        md_keys = ['label', 'chunk_size', 'dtype', 'verbose']
        md_kwargs = {key: kwargs[key] for key in md_keys if key in kwargs}
        parse_md(filename, format, temperature, dt, **md_kwargs)      
    else:
        required_keys = [
            'lammps_data',
            'atom_style_dump',
            'atom_style_data',
            'atom_symbols'
        ]
        if not all(key in kwargs for key in required_keys):
            raise ValueError(f"For lammps format, these args are required: {required_keys}")
        
        for key in required_keys:
            if key not in kwargs:
                raise ValueError(f"This command requires the '{key}' argument.")
        
        lammps_keys = required_keys + ['label', 'chunk_size', 'dtype', 'verbose']
        lammps_kwargs = {key: kwargs.get(key) for key in lammps_keys}
        parse_lammps(lammps_dump=filename, 
                     temperature=temperature, dt=dt, **lammps_kwargs)
    print('')
    
    
@monitor_performance
def cli_concat(path_traj1: str,
               path_traj2: str,
               **kwargs):
    """[CLI method] This method concatenates two successive HDFT trajectory files"""
    con_keys = ['label', 'chunk_size', 'eps', 'dtype', 'verbose']
    con_kwargs = {key: kwargs[key] for key in con_keys if key in kwargs}
    concat_traj(path_traj1, path_traj2, **con_kwargs)
    print('')
    
    
@ monitor_performance
def cli_show(path_traj:str):
    """[CLI method] This method displays metadata of a HDF5 file"""
    show_traj(path_traj)
    print('')

@ monitor_performance
def cli_cut(path_traj:str,
            start_frame: int | None = None,
            end_frame: int | None = None,
            label: str = "CUT",
            chunk_size: int = 5000,
            verbose=True):
    """[CLI method] This method cut HDF5 file into user-defined size"""
    cut_traj(path_traj=path_traj,
             start_frame=start_frame,
             end_frame=end_frame,
             label=label,
             chunk_size=chunk_size) 
    print('')   
    