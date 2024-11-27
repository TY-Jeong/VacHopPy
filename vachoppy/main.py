import os
import sys
import argparse
from fractions import Fraction
from vachoppy.inout import DataInfo
from vachoppy.core import MSD, Parameter, Trajectory
from colorama import Fore

BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m' # Green color
RED = '\033[91m'   # Red color
RESET = '\033[0m'  # Reset to default color

parser = argparse.ArgumentParser(
    description='VacHopPy: A Python Code to Analyze Vacancy Hopping',
    formatter_class=argparse.RawTextHelpFormatter
    )

# Arguments for DataInfo
parser.add_argument('-p1', '--prefix1', default='traj', 
                    help='name of outer directory')
parser.add_argument('-p2', '--prefix2', default='traj', 
                    help='prefix of inner directory : ex.(prefix2).(temp)K')

parser.add_argument(
    '-m', '--mode', 
    choices=['e', 'p', 'ep', 't', 'f'], 
    required=True, 
    help=(
        """Choose mode:
        'e' - For only diffusion coefficient (Einstein relation)
        'p' - For effective diffusion parameter set
        'ep' - Both 'e' and 'p'
        't' - For animation of trajectory
        'f' - For fingerprint analysis
        """
        )
    )

check_mode = True
if '-m' in sys.argv:
    mode_index = sys.argv.index('-m') + 1
elif '--mode' in sys.argv:
    mode_index = sys.argv.index('--mode') + 1
else:
    check_mode = False
    
if check_mode:
    if mode_index < len(sys.argv):
        mode_value = sys.argv[mode_index]

        parser.add_argument('symbol',
                            type=str,
                            help='symbol of moving atom')
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
            parser.add_argument('-f', '--fix_Ea',
                                action='store_true',
                                help='fix Ea for residence time to <Ea_hop>_vhp (experiment)')
            parser.add_argument('-t', '--tolerance',
                                type=float,
                                default=0.1,
                                help='tolerance for Ea of correlation factor in eV (default: 0.1)')
            
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
print(f'{CYAN}{BOLD}VacHopPy is in progress{RESET}')

with open('arg.txt', 'w') as f:
    print(f'{GREEN}{BOLD}Arguments and Values :{RESET}')
    f.write('Arguments and Values :\n')
    for arg, value in vars(args).items():
        print(f'    {arg} = {value}')
        f.write(f'    {arg} = {value}\n')
print('')

def main():
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
                           fix_Ea_t_res=args.fix_Ea,
                           tolerance_Ea_f=args.tolerance)
        
    if mode_value == 't':
        traj = Trajectory(xdatcar=args.xdatcar,
                          force=args.force,
                          poscar=args.lattice,
                          symbol=args.symbol,
                          outcar=args.outcar,
                          interval=args.interval,
                          correction=not(args.no_correction),
                          label=args.label,
                          verbose=args.verbose
                          )
            

