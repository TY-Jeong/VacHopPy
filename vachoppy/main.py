import os
import sys
import argparse
from fractions import Fraction
from vachoppy.inout import DataInfo
from vachoppy.core import *
from vachoppy.utils import *
from colorama import Fore

BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m'
RED = '\033[91m'
RESET = '\033[0m'

parser = argparse.ArgumentParser(
    description='VacHopPy: A Python Code to Analyze Vacancy Hopping',
    formatter_class=argparse.RawTextHelpFormatter
    )

group = parser.add_mutually_exclusive_group(required=True)

# key functionalities
group.add_argument(
    '-m', '--mode', 
    choices=['e', 'p', 'ep', 't', 'f'], 
    help=(
        """Choose mode:
        'e'  - For only diffusion coefficient (Einstein relation)
        'p'  - For effective diffusion parameter set
        'ep' - Both 'e' and 'p'
        't'  - For animation of trajectory
        'f'  - For fingerprint analysis
        """
        )
    )

# utilities
group.add_argument(
    '-u',
    choices=['extract_force', 'concat_xdatcar', 'concat_force', 'update_outcar'],
    dest='util',
    help=(
        """Choose mode:
        'extract_force'  - Extract FORCE from vasprun.xml
        'concat_xdatcar' - Concatenate two XDATCAR files
        'concat_force'   - Concatenate two FORCE files
        'update_outcar'  - Update OUTCAR : nsw = nsw1 + nsw2
        """
        )
)

check_mode = True if '-m' in sys.argv or '--mode' in sys.argv else False
check_util = True if '-u' in sys.argv else False

if check_mode:
    if '-m' in sys.argv:
        mode_index = sys.argv.index('-m') + 1
    elif '--mode' in sys.argv:
        mode_index = sys.argv.index('--mode') + 1

if check_util:
    mode_index = sys.argv.index('-u') + 1
    if mode_index < len(sys.argv):
        mode_value = sys.argv[mode_index]
        
        if mode_value == 'extract_force':
            parser.add_argument('-in', '--file_in',
                                default='vasprun.xml',
                                help='input vasprun.xml file (default: vasprun.xml)')
            parser.add_argument('-out', '--file_out',
                                default='FORCE',
                                help='output FORCE file (default: FORCE)')
            
        if mode_value == 'concat_xdatcar':
            parser.add_argument('-in1', '--xdatcar_in1',
                                required=True,
                                help='first XDATCAR file')
            parser.add_argument('-in2', '--xdatcar_in2',
                                required=True,
                                help='second XDATCAR file')
            parser.add_argument('-out', '--xdatcar_out',
                                default='XDATCAR_NEW',
                                help='new XDATCAR file (default: XDATCAR_NEW)')
        
        if mode_value == 'concat_force':
            parser.add_argument('-in1', '--force_in1',
                                required=True,
                                help='first FORCE file')
            parser.add_argument('-in2', '--force_in2',
                                required=True,
                                help='second FORCE file')
            parser.add_argument('-out', '--force_out',
                                default='FORCE_NEW',
                                help='new Force file (default: FORCE_NEW)')
        
        if mode_value == 'update_outcar':
            parser.add_argument('-in1', '--outcar_in1',
                                required=True,
                                help='first OUTCAR file')
            parser.add_argument('-in2', '--outcar_in2',
                                required=True,
                                help='second OUTCAR file')
            parser.add_argument('-out', '--outcar_out',
                                default='OUTCAR_NEW',
                                help='new OUTCAR file (default: OUTCAR_NEW)')
        

if check_mode:
    if mode_index < len(sys.argv):
        mode_value = sys.argv[mode_index]

        parser.add_argument('symbol',
                            type=str,
                            help='symbol of moving atom')
        
        # Arguments for DataInfo
        if 'e' in mode_value or 'p' in mode_value:
            parser.add_argument('-p1', '--prefix1', 
                                default='traj', 
                                help='name of outer directory')
            parser.add_argument('-p2', '--prefix2', 
                                default='traj', 
                                help='prefix of inner directories, ex.{prefix2}.{temp}K')
            
        if 'e' in mode_value:
            parser.add_argument('x_vac',
                                type=str,
                                help='fraction of vacancy (use 1 for atom)')
            parser.add_argument('t_width',
                                type=float,
                                help='x-range in msd plot in ps')
            parser.add_argument('--skip',
                                type=float,
                                default=0,
                                help='steps to be skipped in ps (default: 0)')
            parser.add_argument('--start',
                                type=float,
                                default=1,
                                help='initial time used for linear fit (default: 1)')
            
        if 'p' in mode_value:
            parser.add_argument('interval',
                                type=float,
                                help='time interval for averaging in ps')
            parser.add_argument('-l', '--lattice',
                                type=str,
                                default='POSCAR_LATTICE',
                                help='lattice file in POSCAR format (default: POSCAR_LATICE)')
            parser.add_argument('-n', '--neb',
                                type=str,
                                default='NEB.csv',
                                help='neb data in csv format (default: NEB.csv)')
            parser.add_argument('-e', '--einstein',
                                type=str,
                                default='Einstein.txt',
                                help='Einstein relation results (default: Einstein.txt)')
            parser.add_argument('-v', '--verbose',
                                action='store_true',
                                help='verbosity for parameter calculation')
            # experiments
            # parser.add_argument('--release_Ea',
            #                     action='store_true',
            #                     help='release Ea for residence time')
            # parser.add_argument('--tolerance',
            #                     type=float,
            #                     default=0,
            #                     help='tolerance for Ea of correlation factor in eV (default: 0)')
            
        if 't' in mode_value:
            parser.add_argument('interval',
                                type=float,
                                help='time interval for averaging in ps')
            parser.add_argument('-x', '--xdatcar',
                                type=str,
                                default='XDATCAR',
                                help='path to XDATCAR file (default: XDATCAR)')
            parser.add_argument('-f', '--force',
                                type=str,
                                default='FORCE',
                                help='path to FORCE file (default: FORCE)')
            parser.add_argument('-l', '--lattice',
                                type=str,
                                default='POSCAR_LATTICE',
                                help='lattice file in POSCAR format (default: POSCAR_LATTICE)')
            parser.add_argument('-o', '--outcar',
                                type=str,
                                default='OUTCAR',
                                help='outcar file (default: OUTCAR)')
            parser.add_argument('--no_correction',
                                action='store_true',
                                help='if use, corrections will be skipped')
            parser.add_argument('--label',
                                action='store_true',
                                help='if use, label of each atom will be shown')
            parser.add_argument('-v', '--verbose',
                                action='store_true',
                                help='verbosity for parameter calculation')


args = parser.parse_args()

def main():
    if check_mode:
        print(f'{CYAN}{BOLD}VacHopPy is in progress{RESET}')
        
        # save arguments
        with open('arg.txt', 'w') as f:
            print(f'{GREEN}{BOLD}Arguments and Values :{RESET}')
            f.write('Arguments and Values :\n')
            for arg, value in vars(args).items():
                print(f'    {arg} = {value}')
                f.write(f'    {arg} = {value}\n')
        print('')
        
        # functionalities
        if 'e' in mode_value or 'p' in mode_value:
            data = DataInfo(prefix1=args.prefix1,
                            prefix2=args.prefix2,
                            verbose=False)

        if 'e' in mode_value:
            msd = MSD(data=data, 
                      tmax=args.t_width, 
                      skip=args.skip,
                      start=args.start,
                      symbol=args.symbol, 
                      x_vac=float(Fraction(args.x_vac)))
            
        if 'p' in mode_value:
            params = Parameter(data=data,
                               symbol=args.symbol,
                               interval=args.interval,
                               poscar=args.lattice,
                               neb=args.neb,
                               einstein=args.einstein,
                               verbose=args.verbose,
                               fix_Ea_t_res=False,
                               tolerance_Ea_f=0)
            
        if mode_value == 't':
            traj = Trajectory(xdatcar=args.xdatcar,
                              force=args.force,
                              poscar=args.lattice,
                              symbol=args.symbol,
                              outcar=args.outcar,
                              interval=args.interval,
                              correction=not(args.no_correction),
                              label=args.label,
                              verbose=args.verbose)
            
    if check_util:
        if mode_value == 'extract_force':
            extract_force(args.file_in, args.file_out)
            print(f'{args.file_out} is created')
            
        if mode_value == 'concat_xdatcar':
            concat_xdatcar(args.xdatcar_in1, args.xdatcar_in2, args.xdatcar_out)
            print(f'{args.xdatcar_out} is created')
        
        if mode_value == 'concat_force':
            concat_force(args.force_in1, args.force_in2, args.force_out)
            print(f'{args.force_out} is created')
        
        if mode_value == 'update_outcar':
            update_outcar(args.outcar_in1, args.outcar_in2, args.outcar_out)
            print(f'{args.outcar_out} is created')

