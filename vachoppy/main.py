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
    choices=['p', 'pp', 'e', 't', 'a', 'd'], 
    help=(
        """Choose mode:
        'p'  - For effective diffusion parameter set
        'pp'  - Post-process for effective parameter calculation (neb.csv file is required)
        'e'  - Diffusion coefficient calculation using Einstein relation
        't'  - For animation of trajectory
        'a'  - For path analysis
        'd'  - For cosine distance analysis
        """
        )
    )

# utilities
group.add_argument(
    '-u',
    choices=['extract_force', 'concat_xdatcar', 'concat_force', 
             'update_outcar', 'fingerprint', 'cosine_distance'],
    dest='util',
    help=(
        """Choose mode:
        'extract_force'   - Extract FORCE from vasprun.xml
        'concat_xdatcar'  - Concatenate two XDATCAR files
        'concat_force'    - Concatenate two FORCE files
        'update_outcar'   - Update OUTCAR : nsw = nsw1 + nsw2
        'fingerprint'     - Extract fingerprint
        'cosine_distance' - Calculate cosine distance
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
            
        if mode_value == 'cosine_distance':
            parser.add_argument('-in1', '--fingerprint_in1',
                                required=True,
                                help='first fingerprint file')
            parser.add_argument('-in2', '--fingerprint_in2',
                                required=True,
                                help='second fingerprint file')
            
        if mode_value == 'fingerprint':
            parser.add_argument('-p','--poscar',
                                required=True,
                                type=str,
                                help='POSCAR to be used for fingerprint extraction')
            parser.add_argument('--prefix',
                                type=str,
                                default='fingerprint',
                                help='prefix for output files (default: fingerprint)')
            parser.add_argument('-d', '--disp',
                                action='store_true',
                                help='if use, fingerprint plot will be displayed')
        

if check_mode:
    if mode_index < len(sys.argv):
        mode_value = sys.argv[mode_index]

        # Arguments for DataInfo
        if mode_value in ['p', 'e', 't', 'a']:
            parser.add_argument('symbol',
                                type=str,
                                help='symbol of moving atom')
        
        if mode_value in ['e', 'p', 't', 'a']:
            parser.add_argument('-p1', '--prefix1', 
                                default='traj', 
                                help='name of outer directory (default: traj)')
            parser.add_argument('-p2', '--prefix2', 
                                default='traj', 
                                help='prefix of inner directories, ex.{prefix2}.{temp}K (default: traj)')
            
        if mode_value == 'e':
            parser.add_argument('t_width',
                                type=float,
                                help='x-range in msd plot in ps')
            parser.add_argument('--skip',
                                type=float,
                                default=0,
                                help='steps to be skipped in ps (default: 0)')
            parser.add_argument('--x_vac',
                                type=str,
                                default=1,
                                help='fraction of vacancy (use 1 for atom) (default: 1)')
            parser.add_argument('--start',
                                type=float,
                                default=1,
                                help='initial time used for linear fit (default: 1)')
            
        if mode_value == 'p':
            parser.add_argument('interval',
                                type=float,
                                help='time interval for averaging in ps')
            parser.add_argument('-l', '--lattice',
                                type=str,
                                default='POSCAR_LATTICE',
                                help='lattice file in POSCAR format (default: POSCAR_LATICE)')
            parser.add_argument('--rmax',
                                type=float,
                                default=3.0,
                                help='maximum distance for hopping path identification (default: 3.0)')
            parser.add_argument('--tol',
                                type=float,
                                default=1e-3,
                                help='tolerance for VoronoiNN (default: 1e-3)')
            parser.add_argument('--tolerance',
                                type=float,
                                default=1e-3,
                                help='tolerance for distance comparison (default: 1e-3)')
        if mode_value == 'pp':
            parser.add_argument('-p', '--parameter',
                                type=str,
                                default='parameter.txt',
                                help='parameter.txt file (default: parameter.txt)')
            parser.add_argument('-n', '--neb',
                                type=str,
                                default='neb.csv',
                                help='neb.csv file containing hopping barriers for each path (default: neb.csv)')
            
        if mode_value == 't':
            parser.add_argument('interval',
                                type=float,
                                help='time interval for averaging in ps')
            parser.add_argument('temp',
                                type=int,
                                help='temperatue in K')
            parser.add_argument('label',
                                type=str,
                                help='labels')
            parser.add_argument('-l', '--lattice',
                                type=str,
                                default='POSCAR_LATTICE',
                                help='lattice file in POSCAR format (default: POSCAR_LATTICE)')
            parser.add_argument('--update_alpha',
                                type=float,
                                default=0.75,
                                help='adjust the remaining time of the arrow (default: 0.75)')
            parser.add_argument('--no_correction',
                                action='store_true',
                                help='if use, corrections will be skipped')
            parser.add_argument('--show_index',
                                action='store_true',
                                help='if use, index of each atom will be shown')
            parser.add_argument('--dpi',
                                type=int,
                                default=300,
                                help='adjust dpi of snapshots (default: 300)')
            parser.add_argument('-v', '--verbose',
                                action='store_true',
                                help='verbosity for parameter calculation')
            
        if mode_value == 'a':
            parser.add_argument('interval',
                                type=float,
                                help='time interval for averaging in ps')
            parser.add_argument('temp',
                                type=int,
                                help='temperatue in K')
            parser.add_argument('--label',
                                nargs="+",
                                type=str,
                                help='labels')
            parser.add_argument('-l', '--lattice',
                                type=str,
                                default='POSCAR_LATTICE',
                                help='lattice file in POSCAR format (default: POSCAR_LATICE)')
            parser.add_argument('--rmax',
                                type=float,
                                default=3.0,
                                help='maximum distance for hopping path identification (default: 3.0)')
            parser.add_argument('--tol',
                                type=float,
                                default=1e-3,
                                help='tolerance for VoronoiNN (default: 1e-3)')
            parser.add_argument('--tolerance',
                                type=float,
                                default=1e-3,
                                help='tolerance for distance comparison (default: 1e-3)')
            
        if mode_value == 'd':
            parser.add_argument('interval',
                                type=float,
                                help='time interval for averaging in ps')
            parser.add_argument('-x','--xdatcar',
                                type=str,
                                default='XDATCAR',
                                help='path to XDATCAR file (default: XDATCAR)')
            parser.add_argument('-o','--outcar',
                                type=str,
                                default='OUTCAR',
                                help='path to OUTCAR file (default: OUTCAR)')
            parser.add_argument('-p','--poscar_mother',
                                default='POSCAR_MOTHER',
                                type=str,
                                help='path to POSCAR of mother phase (default: POSCAR_MOTHER)')
            parser.add_argument('--prefix1',
                                type=str,
                                default='snapshots',
                                help='directory to save xdatcar snapshots (default: snapshots)')
            parser.add_argument('--prefix2',
                                type=str,
                                default='fingerprints',
                                help='directory to save fingerprints (default: fingerprints)')


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
        if mode_value in ['e', 'p', 'a', 't']:
            data = DataInfo(prefix1=args.prefix1,
                            prefix2=args.prefix2,
                            verbose=True)

        if mode_value == 'e':
            msd = MSD(data=data, 
                      tmax=args.t_width, 
                      skip=args.skip,
                      start=args.start,
                      symbol=args.symbol, 
                      x_vac=float(Fraction(args.x_vac)))
            
        if mode_value == 'p':
            effective_params = EffectiveDiffusionParameter(data=data,
                                                           interval=args.interval,
                                                           poscar_lattice=args.lattice,
                                                           symbol=args.symbol,
                                                           file_out='parameter.txt',
                                                           rmax=args.rmax,
                                                           tol=args.tol,
                                                           tolerance=args.tolerance,
                                                           verbose=True)
        
        if mode_value == 'pp':
            post = Post_EffectiveDiffusionParameter(file_params=args.parameter,
                                                    file_neb=args.neb)
            
        if mode_value == 't':
            traj = Trajectory(data=data,
                              temp=args.temp,
                              label=args.label,
                              interval=args.interval,
                              poscar_lattice=args.lattice,
                              symbol=args.symbol,
                              correlation=not(args.no_correction),
                              update_alpha=args.update_alpha,
                              show_index=args.show_index,
                              dpi=args.dpi,
                              verbose=args.verbose)
        
        if mode_value == 'a':
            anal = PathAnalyzer(data=data,
                                interval=args.interval,
                                poscar_lattice=args.lattice,
                                symbol=args.symbol,
                                temp=args.temp,
                                label=args.label,
                                rmax=args.rmax,
                                tol=args.tol,
                                tolerance=args.tolerance,
                                verbose=True)
            
        if mode_value == 'd':
            phase = PhaseTransition(xdatcar=args.xdatcar,
                                    outcar=args.outcar,
                                    interval=args.interval,
                                    poscar_mother=args.poscar_mother,
                                    prefix1=args.prefix1,
                                    prefix2=args.prefix2)
            
                
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
        
        if mode_value == 'cosine_distance':
            GetCosineDistance(args.fingerprint_in1, args.fingerprint_in2)
            
        if mode_value == 'fingerprint':
            finger = GetFingerPrint(poscar=args.poscar,
                                    prefix=args.prefix,
                                    disp=args.disp)
            

