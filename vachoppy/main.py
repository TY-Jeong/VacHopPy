import argparse
from typing import Tuple
from vachoppy.cli import *


def _atom_pair_type(s: str) -> Tuple[str, str]:
    """Custom type for argparse to parse atom pairs like 'A-B'."""
    try:
        A, B = s.split('-')
        return (A, B)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Atom pair must be in 'A-B' format, but got '{s}'")
    
def main():
    """
    Main function to parse command-line arguments and run the appropriate
    Vachoppy CLI command.
    """
    # ==============================================================================
    # Main Parser
    # ==============================================================================
    parser = argparse.ArgumentParser(
        description='VacHopPy: A Python package for vacancy hopping analysis from MD simulations.',
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='command', 
                                       required=True,
                                       help='Available commands')

    # ==============================================================================
    # Sub-parser for the 'trajectory' command
    # ==============================================================================
    p_traj = subparsers.add_parser(
        'trajectory', 
        help='Identify and visualize vacancy trajectories from a single trajectory file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_traj.add_argument(
        'path_traj', 
        type=str, 
        help='Path to the HDF5 trajectory file (e.g., TRAJ_O.h5).')
    p_traj.add_argument(
        'path_structure', 
        type=str, 
        help='Path to the crystal structure file (e.g., POSCAR, cif).')
    p_traj.add_argument(
        'symbol', 
        type=str, 
        help='Chemical symbol of the diffusing species (e.g., O, Li).')
    p_traj.add_argument(
        '--t_interval', 
        type=float, 
        default=None, 
        help='Time interval in picoseconds (ps) for time-averaging the trajectory.')
    p_traj.add_argument(
        '--unwrap',
        action='store_true',
        help='Flag to plot the unwrapped trajectory across periodic boundaries.')
    p_traj.add_argument(
        '--structure_format',
        dest='format',
        default=None,
        type=str, 
        help="File format supported by ASE, e.g., 'vasp', 'cif', 'xyz'.")
    p_traj.add_argument(
        '--sampling_size',
        type=int,
        default=5000,
        help='Number of initial frames to use for `t_interval` auto-estimation.')
    p_traj.add_argument(
        '--a_max',
        dest='rmax', 
        type=float, 
        default=3.25, 
        help='Maximum distance (Angstrom) for neighbor search in path finding.')
    p_traj.add_argument(
        '--exclude-incomplete',
        dest='use_incomplete_encounter',
        action='store_false',
        help="Exclude encounters that are not completed by the end of the simulation.")
    p_traj.add_argument(
        '--eps',
        type=float,
        default=1e-3,
        help='Tolerance for floating-point comparisons.')
    
    # ==============================================================================
    # Sub-parser for the 'analyze' command
    # ==============================================================================
    p_anal = subparsers.add_parser(
        'analyze', 
        help='Extract hopping parameters from an ensemble of trajectories.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_anal.add_argument(
        'path_traj', 
        type=str, 
        help='Path to a single HDF5 file or a root directory to search for files.')
    p_anal.add_argument(
        'path_structure', 
        type=str, 
        help='Path to the crystal structure file (e.g., POSCAR, cif).')
    p_anal.add_argument(
        'symbol', 
        type=str, 
        help='Chemical symbol of the diffusing species (e.g., O, Li).')
    p_anal.add_argument(
        '--neb_csv', 
        type=str, 
        default=None, 
        help='Path to the NEB results CSV file for attempt frequency calculation.')
    p_anal.add_argument(
        '--xyz', 
        action='store_true',
        help='Decompose the diffusivity into x, y, z components.')
    p_anal.add_argument(
        '--no_disp',
        dest='disp',
        action='store_false',
        help='Do not display plots in a popup window upon creation.')
    p_anal.add_argument(
        '--t_interval', 
        type=float, 
        default=None, 
        help='Time interval in picoseconds (ps) for time-averaging the trajectory.')
    p_anal.add_argument(
        '--depth',
        type=int,
        default=2,
        help='Maximum directory depth to search for trajectory files.')
    p_anal.add_argument(
        '--prefix', 
        type=str, 
        default='TRAJ', 
        help='Prefix of trajectory files.')
    p_anal.add_argument(
        '--structure_format',
        dest='format',
        default=None,
        type=str, 
        help="File format supported by ASE, e.g., 'vasp', 'cif', 'xyz'.")
    p_anal.add_argument(
        '--n_jobs',
        type=int,
        default=-1,
        help='Number of CPU cores for parallel processing (-1 means all cores).')
    p_anal.add_argument(
        '--dir_imgs', 
        type=str, 
        default='imgs', 
        help="Directory to save output images.")
    p_anal.add_argument(
        '--sampling_size',
        type=int,
        default=5000,
        help='Number of initial frames to use for `t_interval` auto-estimation.')
    p_anal.add_argument(
        '--a_max',
        dest='rmax', 
        type=float, 
        default=3.25, 
        help='Maximum distance (Angstrom) for neighbor search in path finding.')
    p_anal.add_argument(
        '--exclude-incomplete',
        dest='use_incomplete_encounter',
        action='store_false',
        help="Exclude encounters that are not completed by the end of the simulation.")
    p_anal.add_argument(
        '--eps',
        type=float,
        default=1e-3,
        help='Tolerance for floating-point comparisons.')

    # ==============================================================================
    # Sub-parser for the 'vibration' command
    # ==============================================================================
    p_vib = subparsers.add_parser(
        'vibration', 
        help='Extract atomic vibration frequency from a single trajectory',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_vib.add_argument(
        'path_traj', 
        type=str, 
        help='Path to the HDF5 trajectory file (e.g., TRAJ_O.h5).')
    p_vib.add_argument(
        'path_structure', 
        type=str, 
        help='Path to the crystal structure file (e.g., POSCAR, cif).')
    p_vib.add_argument(
        'symbol', 
        type=str, 
        help='Chemical symbol of the diffusing species (e.g., O, Li).')  
    p_vib.add_argument(
        '--structure_format',
        dest='format',
        default=None,
        type=str, 
        help="File format supported by ASE, e.g., 'vasp', 'cif', 'xyz'.")
    p_vib.add_argument(
        '--n_jobs',
        type=int,
        default=-1,
        help='Number of CPU cores for parallel processing (-1 means all cores).')
    p_vib.add_argument(
        '--dir_imgs', 
        type=str, 
        default='imgs', 
        help="Directory to save output images.")
    p_vib.add_argument(
        '--sampling_size',
        type=int,
        default=5000,
        help='Number of initial frames to use for `t_interval` auto-estimation.')
    p_vib.add_argument(
        '-jump_detection_radius',
        type=float,
        default=1.0,
        help='The radius in Angstroms used in the initial step to distinguish')
    p_vib.add_argument(
        '--no-filter', 
        dest='filter_high_freq', 
        action='store_false', 
        help='Disable high-frequency outlier filtering.')
    p_vib.add_argument(
        '--a_max',
        dest='rmax', 
        type=float, 
        default=3.25, 
        help='Maximum distance (Angstrom) for neighbor search in path finding.')
    p_vib.add_argument(
        '--eps',
        type=float,
        default=1e-3,
        help='Tolerance for floating-point comparisons.')
    
# ==============================================================================
# Sub-parser for the 'distance' command
# ==============================================================================
    p_dist = subparsers.add_parser(
        'distance', 
        help='Trace change in cosine distance from reference structure over time.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_dist.add_argument(
        'path_traj',
        nargs='+',
        help='Path(s) to one or more HDF5 trajectory files (e.g., TRAJ_Hf.h5 TRAJ_O.h5).')
    p_dist.add_argument(
        't_interval',
        type=float,
        help='Time interval in picoseconds (ps) for time-averaging the trajectory.')
    p_dist.add_argument(
        'reference_structure',
        type=str,
        help='Path to the reference structure file (e.g., POSCAR).')
    p_dist.add_argument(
        '--atom_pairs',
        type=_atom_pair_type,
        nargs='+',
        default=None,
        help="Space-separated atom pairs to analyze (e.g., --atom_pairs Ti-O O-O). "
             "If omitted, all unique pairs are generated automatically.")
    p_dist.add_argument(
        '--Rmax', 
        type=float, 
        default=10.0, 
        help='Cutoff radius (Angstrom) for the fingerprint.')
    p_dist.add_argument(
        '--delta', 
        type=float, 
        default=0.08, 
        help='Discretization step (Angstrom) for the fingerprint.')
    p_dist.add_argument(
        '--sigma', 
        type=float, 
        default=0.03, 
        help='Gaussian broadening (Angstrom) for the fingerprint.')
    p_dist.add_argument(
        '--dirac', 
        type=str, 
        default='g', 
        choices=['g', 's'], 
        help="Dirac function type: 'g' for Gaussian, 's' for square.")
    p_dist.add_argument(
        '--n_jobs', 
        type=int, 
        default=-1, 
        help='Number of CPU cores for parallel processing (-1 means all).')
    p_dist.add_argument(
        '--window_size', 
        type=int, 
        default=50, 
        help='Window size for the moving average in fluctuation analysis.')
    p_dist.add_argument(
        '--threshold_std', 
        type=float, 
        default=None, 
        help='Sigma threshold for detecting fluctuations.')
    
    # ==============================================================================
    # Sub-parser for the 'fingerprint' command
    # ==============================================================================
    p_fp = subparsers.add_parser(
        'fingerprint',
        help='Calculate and plot the fingerprint for a single static structure.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_fp.add_argument(
        'path_structure', 
        type=str, 
        help='Path to the crystallographic structure file (e.g., POSCAR, cif).')
    p_fp.add_argument(
        '--atom_pairs',
        type=_atom_pair_type,
        nargs='+',
        default=None,
        help="Space-separated atom pairs to analyze (e.g., --atom_pairs Ti-O O-O). "
             "If omitted, all unique pairs are generated automatically.")
    p_fp.add_argument(
        '--Rmax', 
        type=float, 
        default=10.0, 
        help='Cutoff radius (Angstrom).')
    p_fp.add_argument(
        '--delta', 
        type=float, 
        default=0.08, 
        help='Discretization step (Angstrom).')
    p_fp.add_argument(
        '--sigma', 
        type=float, 
        default=0.03, 
        help='Gaussian broadening (Angstrom).')
    p_fp.add_argument(
        '--dirac', 
        type=str, 
        default='g', 
        choices=['g', 's'], 
        help="Dirac function type ('g' or 's').")

    # ==============================================================================
    # Sub-parser for the 'msd' command
    # ==============================================================================
    p_msd = subparsers.add_parser(
        'msd', 
        help='Calculate diffusivity from Mean Squared Displacement (Einstein relation).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_msd.add_argument(
        'path_traj', 
        type=str, 
        help='Path to a single HDF5 file or a root directory to search for files.')
    p_msd.add_argument(
        'symbol', 
        type=str, 
        help='Chemical symbol of the diffusing species to analyze.')
    p_msd.add_argument(
        '--segment_length',
        type=float,
        nargs='*', # 단일 값 사용가능한지 디버깅 해봐야함
        default=None,
        help='Time length in ps for each segment. Can be a single float or a list.')
    p_msd.add_argument(
        '--skip', 
        type=float, 
        default=0.0, 
        help='Initial time in ps to skip for equilibration.')
    p_msd.add_argument(
        '--start', 
        type=float, 
        default=1.0, 
        help='Start time in ps for the linear fitting range of MSD.')
    p_msd.add_argument(
        '--end', 
        type=float, 
        default=None, 
        help='End time in ps for the linear fitting range.')
    p_msd.add_argument(
        '--n_jobs', 
        type=int, 
        default=-1, 
        help='Number of CPU cores for parallel processing (-1 means all).')
    p_msd.add_argument(
        '--dir_imgs',
        type=str,
        default='imgs',
        help='Directory to save output images.')
    p_msd.add_argument(
        '--prefix', 
        type=str, 
        default='TRAJ', 
        help='Prefix of trajectory files.')
    p_msd.add_argument(
        '--depth',
        type=int,
        default=2,
        help='Maximum directory depth to search for trajectory files.')
    p_msd.add_argument(
        '--eps',
        type=float,
        default=1e-3,
        help='Tolerance for floating-point comparisons.')
    
    # ==============================================================================
    # Sub-parser for the 'convert' command
    # ==============================================================================
    p_conv = subparsers.add_parser(
        'convert', 
        help='Convert various MD trajectory formats to the standard HDF5 format.',
        formatter_class=argparse.RawTextHelpFormatter)
    p_conv.add_argument(
        'filename', 
        type=str, 
        help='Path to the source MD trajectory file (e.g., vasprun.xml, lammps.dump).')
    p_conv.add_argument(
        'temperature', 
        type=float, 
        help='Simulation temperature in Kelvin.')
    p_conv.add_argument(
        'dt', 
        type=float, 
        help='Timestep in femtoseconds.')
    p_conv.add_argument(
        '--format', 
        type=str,
        default=None,
        help="Format of the source file. Use 'lammps-dump-text' for LAMMPS, \n"
             "or any other ASE-supported format (e.g., 'vasp-xml', 'extxyz').")
    p_conv.add_argument(
        '--label', 
        type=str, 
        default=None, 
        help='Custom label for the output TRAJ_*.h5 files.')
    p_conv.add_argument(
        '--chunk_size', 
        type=int, 
        default=5000, 
        help='Number of frames to process in memory at once.')
    
    lammps_group = p_conv.add_argument_group(
        'LAMMPS Specific Arguments', 
        "These are REQUIRED when format is 'lammps-dump-text'")
    lammps_group.add_argument(
        '--lammps_data', 
        type=str, 
        help="Path to the LAMMPS data file (topology).")
    lammps_group.add_argument(
        '--atom_style_data', 
        type=str, 
        help="Atom style for the LAMMPS data file.")
    lammps_group.add_argument(
        '--atom_style_dump', 
        type=str, 
        help="Atom style for the LAMMPS dump file.")
    lammps_group.add_argument(
        '--atom_symbols',
        nargs='+',
        metavar='ID=SYMBOL',
        help="Space-separated key=value pairs mapping type IDs to symbols. "
             "Example: --atom_symbols 1=Hf 2=O")

    # ==============================================================================
    # Sub-parser for the 'concat' command
    # ==============================================================================
    p_concat = subparsers.add_parser(
        'concat', 
        help='Concatenate two successive HDF5 trajectory files into a new one.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p_concat.add_argument(
        'path_traj1', 
        type=str, 
        help='Path to the first HDF5 trajectory file (chronologically).')
    p_concat.add_argument(
        'path_traj2', 
        type=str, 
        help='Path to the second HDF5 trajectory file to be appended.')
    p_concat.add_argument(
        '--label', 
        type=str, 
        default='CONCAT', 
        help='A custom label for the output concatenated HDF5 file.')            
    p_concat.add_argument(
        '--chunk_size', 
        type=int, 
        default=5000, 
        help='Number of frames to process in memory at once during concatenation.')
    p_concat.add_argument(
        '--eps', 
        type=float, 
        default=1e-3, 
        help='Tolerance for floating-point comparisons during metadata validation.')
    p_concat.add_argument(
        '--no-verbose', 
        dest='verbose', 
        action='store_false', 
        help='Disable verbose output.')
    
    # ==============================================================================
    # Sub-parser for the 'cut' command
    # ==============================================================================
    p_cut = subparsers.add_parser(
        'cut', 
        help='Cut a portion of a HDF5 file and saves it as a new file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p_cut.add_argument(
        'path_traj', 
        type=str, 
        help='Path to the source HDF5 trajectory file.')
    p_cut.add_argument(
        '--start', 
        dest='start_frame',
        type=int,
        default=None,
        help='The starting frame number (inclusive). Defaults to the beginning.')
    p_cut.add_argument(
        '--end',
        dest='end_frame', 
        type=int,
        default=None,
        help='The ending frame number (exclusive). Defaults to the end.')
    p_cut.add_argument(
        '--label', 
        type=str,
        default='CUT',
        help='A label for the output cut file.')
    p_cut.add_argument(
        '--chunk_size', 
        type=int,
        default=5000,
        help='The number of frames to process in each chunk to conserve memory.')
    
    # ==============================================================================
    # Sub-parser for the 'show' command
    # ==============================================================================
    p_show = subparsers.add_parser(
        'show', 
        help='Display metadata summary of a HDF5 trajectory file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p_show.add_argument(
        'path_traj', 
        type=str, 
        help='Path to the HDF5 trajectory file to inspect.')

    args = parser.parse_args()
    kwargs = vars(args)
    command = kwargs.pop('command')
    
    if command == 'trajectory': cli_trajectory(**kwargs)
    elif command == 'analyze': cli_analyze(**kwargs)
    elif command == 'vibration': cli_vibration(**kwargs)
    elif command == 'distance': cli_distance(**kwargs)
    elif command == 'fingerprint': cli_fingerprint(**kwargs)
    elif command == 'msd': 
        if kwargs['segment_length'] and len(kwargs['segment_length']) == 1:
            kwargs['segment_length'] = kwargs['segment_length'][0]
        cli_msd(**kwargs)
    elif command == 'convert':
        if args.format == 'lammps-dump-text':
            required_for_lammps = ['lammps_data', 'atom_style_data', 'atom_style_dump', 'atom_symbols']
            if not all(getattr(args, key) for key in required_for_lammps):
                parser.error(
                    "For 'lammps-dump-text' format, the following arguments are required:\n"
                    "  --lammps_data, --atom_style_data, --atom_style_dump, --atom_symbols"
                )
            if 'atom_symbols' in kwargs and kwargs['atom_symbols'] is not None:
                try:
                    symbols_dict = {
                        int(pair.split('=')[0]): pair.split('=')[1]
                        for pair in kwargs['atom_symbols']
                    }
                    kwargs['atom_symbols'] = symbols_dict
                except (ValueError, IndexError):
                    parser.error("Argument --atom_symbols must be in 'ID=SYMBOL' format (e.g., '1=Hf').")
        cli_convert(**kwargs)
    elif command == 'concat': cli_concat(**kwargs)
    elif command == 'cut': cli_cut(**kwargs)
    elif command == 'show': cli_show(**kwargs)
    else:
        print(f"Unknown command: {command}")
        parser.print_help()

if __name__ == '__main__':
    main()