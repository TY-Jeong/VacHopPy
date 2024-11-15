import os
import sys
import argparse
from fractions import Fraction
from vachoppy.inout import DataInfo
from vachoppy.core import MSD, Parameter
from colorama import Fore
BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m' # Green color
RED = '\033[91m'   # Red color
RESET = '\033[0m'  # Reset to default color

parser = argparse.ArgumentParser(
    description='VacHopPy: A Python Code to Analyze Vacancy Hopping')

# Arguments for DataInfo
parser.add_argument('-p1', '--prefix1', default='traj', 
                    help='name of outer directory')
parser.add_argument('-p2', '--prefix2', default='traj', 
                    help='prefix of inner directory : ex.(prefix2).(temp)K')

parser.add_argument(
    '-m', '--mode', 
    choices=['e', 'p', 'ep'], 
    required=True, 
    help="'e' for Einstein relation\n'p' for Parameter calculation\n'ep' for both"
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
                                help='interval for trajectory module')
            parser.add_argument('-l', '--lattice',
                                type=str,
                                default='POSCAR',
                                help='lattice file in POSCAR format (default: POSCAR)')
            parser.add_argument('-n', '--neb',
                                type=str,
                                default='neb.csv',
                                help='neb data in csv format (default: neb.csv)')
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
            

