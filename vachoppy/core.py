import os
import sys
import csv
import copy
import numpy as np
from tqdm import tqdm
from colorama import Fore
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from itertools import combinations_with_replacement

from vachoppy.inout import *
from vachoppy.einstein import *
from vachoppy.trajectory import *
from vachoppy.fingerprint import *
from vachoppy.utils import *

BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m' # Green color
RED = '\033[91m'   # Red color
RESET = '\033[0m'  # Reset to default color
    
      
class MSD:
    def __init__(self,
                 data,
                 tmax,
                 skip=0,
                 start=1,
                 symbol='O',
                 x_vac=1):
        '''
        Calculate MSD at temperatures in data.
        Arguments:
            data : instance of dataInfo
            tmax : x-axis of msd plot (ps)
            skip : steps to be skipped (ps)
            start :  (ps)
        '''
        print(f'{CYAN}{BOLD}Diffusion coefficient from Einstein relation.{RESET}')
        self.data = data
        self.tmax = tmax
        self.skip = skip
        self.start = start
        self.symbol = symbol
        self.x_vac = x_vac
        
        self.msd = []
        self.D = []
        self.runEnsembleEinstein()
        self.saveMSD()
        self.D /= x_vac
        
        self.plot_results()
        
    def runEnsembleEinstein(self):
        desc = 'Einstein'
        for i, T in enumerate(
            tqdm(self.data.temp, 
                 bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'% (Fore.GREEN, Fore.GREEN),
                 ascii=False,
                 desc=f'{GREEN}{BOLD}{desc:>9s}')):
            path_dir = os.path.join(self.data.prefix1, f'{self.data.prefix2}.{T}K')
            step_skip = int(self.skip * 1000 / self.data.potim[i])
            step_tmax = int(self.tmax * 1000 / self.data.potim[i])
            step_start = int(self.start * 1000 / self.data.potim[i])
            step_end = int(self.tmax * 1000 / self.data.potim[i])
            
            if (self.data.nsw[i] - step_skip) % step_tmax != 0:
                print(f'The MD time is not divided by the tmax. (T={T} K)')
                time_tot = self.data.nsw[i] * self.data.potim[i] / 1000
                print(f'  total MD time = {time_tot} ps ({self.data.nsw[i]} step)')
                print(f'  skip = {self.skip} ps ({step_skip} step)')
                print(f'  tmax = {self.tmax} ps ({step_tmax} step)')
                print(f'  ({self.data.nsw[i]} - {step_skip}) % {step_tmax} != 0')
                sys.exit(0)
            else:
                segment = int((self.data.nsw[i] - step_skip) / step_tmax)
                
            msd = EnsembleEinstein(symbol=self.symbol,
                                   prefix=path_dir,
                                   labels=self.data.label[i],
                                   segments=segment,
                                   skip=step_skip,
                                   start=step_start,
                                   end=step_end)
            self.msd.append(msd)
    
    def saveMSD(self):
        for msd in self.msd:
            msd.save_msd()
            self.D.append(msd.diffcoeff)
        self.D = np.array(self.D)
        print('msd directory is created.')
        print('msd raw data is written in msd directory.')

    def plot_results(self):
        kb = 8.61733326e-5
        cmap = plt.get_cmap("Set1")
        
        # msd plot
        plt.rcParams['figure.figsize'] = (3.8, 3.8)
        plt.rcParams['font.size'] = 11
        fig, ax = plt.subplots()
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.2)
            
        for i, msd in enumerate(self.msd):
            x, y = msd.timestep, msd.msd
            x_fit = np.linspace(self.start, self.tmax, 1000)
            slop, intercept = np.polyfit(x, y, deg=1)
            ax.plot(x, y, 
                    c=cmap(i), linestyle='-', linewidth=2.5, alpha=0.35)
            ax.plot(x_fit, slop*x_fit+intercept, 
                    c=cmap(i), linestyle='-', linewidth=2.5, label=str(self.data.temp[i]))
        
        ax.axvline(self.start, 0, 1, c='k', linewidth=1, linestyle=':')
        plt.xlim([-1, self.tmax])
        plt.xlabel('t (ps)', fontsize=14)
        plt.ylabel(r'MSD ($Å^2$)', fontsize=14)
        plt.legend(loc='upper left', fontsize=12,
                   title='T (K)', title_fontsize=13, 
                   fancybox=True, framealpha=1, edgecolor='inherit')
        plt.savefig('msd.svg', transparent=True, dpi=600, bbox_inches="tight")
        # plt.show()
        print('msd.svg is created.')
            
        # Arrhenius plot
        plt.rcParams['figure.figsize'] = (3.8, 3.8)
        plt.rcParams['font.size'] = 11
        fig, ax = plt.subplots()
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.2)
            
        slop, intercept = np.polyfit(1/self.data.temp, np.log(self.D), deg=1)
        x_fit = np.linspace(1/self.data.temp[-1], 1/self.data.temp[0], 100)
        plt.plot(x_fit, slop*x_fit+intercept, 'k:')
        
        for i, D in enumerate(self.D):
            plt.scatter(1/self.data.temp[i], np.log(D), 
                        color=cmap(i), marker='s', s=50, label=str(self.data.temp[i]))
        plt.xlabel('1/T (1/K)', fontsize=14)
        plt.ylabel(r'ln $D$', fontsize=14)
        num_data = len(self.D)
        ncol = int(np.ceil(num_data / 5))
        plt.legend(loc='best', fancybox=True, framealpha=1, edgecolor='inherit',
                   ncol=ncol, labelspacing = 0.3, columnspacing=0.5, borderpad=0.2, handlelength=0.6,
                   fontsize=11, title='T (K)', title_fontsize=11)
        
        if num_data >= 3:
            x = np.array([self.data.temp[0], self.data.temp[int(num_data/2)] ,self.data.temp[-1]])
        else:
            x = self.data.temp
            
        x_str = [f'1/{T}' for T in x]
        x = 1/x
        plt.xticks(x, x_str)
        plt.savefig('Arrhenius.svg', transparent=True, dpi=600, bbox_inches="tight")
        # plt.show()
        print('Arrhenius.svg is created.')
        
        # write results
        with open('Einstein.txt', 'w') as f:
            f.write(f'symbol of moving atom = {self.symbol}\n')
            f.write(f'x_vac = {self.x_vac :.5f}\n')
            if self.x_vac == 1:
                f.write('target : atom\n')
            else:
                f.write('target : vacancy\n')
            f.write('\nParameters for diffusion coefficient : \n')
            f.write(f'  D0   = {np.exp(intercept):.5e} m2/s\n')
            f.write(f'  Ea_D = {-kb * slop:.5f} eV\n')
            f.write('\n')
            f.write('Raw data\n')
            f.write('T(K) \tD(m2/s)\n')
            for temp, D in zip(self.data.temp, self.D):
                f.write(f'{temp} \t{D :.5e}\n')
                   
        print('Einstein.txt is created.')
        print('parameters for diffusion coefficient : ')
        print(f'  D0   = {np.exp(intercept) :.5e} m2/s')
        print(f'  Ea_D = {-kb * slop :.5f} eV')
        print('')
        
        
        
class Parameter:
    def __init__(self,
                 data,
                 interval,  # ps
                 poscar='POSCAR',
                 neb='neb.csv',
                 einstein='Einstein.txt',
                 symbol='O',
                 verbose=False,
                 fix_Ea_t_res=True,
                 tolerance_Ea_f=0.0):
        
        self.data = data
        self.interval = interval
        self.poscar = poscar
        self.neb = neb
        self.einstein = einstein
        self.symbol = symbol
        self.temp = self.data.temp
        self.verbose = verbose
        self.cmap = plt.get_cmap("Set1")
        self.kb = 8.61733326e-5
        self.fix_Ea_t_res = fix_Ea_t_res
        self.tolerance_Ea_f = tolerance_Ea_f
        
        # check file
        if not os.path.isfile(self.poscar):
            print(f'{self.poscar} is not found.')
            sys.exit(0)
        if not os.path.isfile(self.neb):
            print(f'{self.neb} is not found.')
            sys.exit(0)
        if not os.path.isfile(self.einstein):
            print(f'{self.einstein} is not found.')
            print('you should calculate Einstein relation in advance.')
            sys.exit(0)
        
        print(f'{CYAN}{BOLD}Effective diffusion parameters.{RESET}')
        
        # lattice
        self.lattice = Lattice(self.poscar, symbol=self.symbol)
        self.read_neb()
        
        # correlation factor
        self.cor = []
        self.correlation_factor()
        self.f_cum = [cor.f_cum for cor in self.cor]
        self.count = np.array([cor.counts_cum * cor.num_enc_cum for cor in self.cor])
        self.label_err = [cor.label_err for cor in self.cor]
        self.save_correlation_factor()
        
        # Arrhenius fitting (f_cor)
        slop_f, intercept_f = np.polyfit(1/self.temp, np.log(self.f_cum), deg=1)
        self.f0 = np.exp(intercept_f)
        self.Ea_f = -self.kb * slop_f
        if abs(self.Ea_f) < self.tolerance_Ea_f:
            self.Ea_f = 0
            self.f0 = np.average(self.f_cum)
            slop_f, intercept_f = 0, self.f0
        self.plot_correlation_factor(slop_f, intercept_f)
        
        # read einstein
        self.D=[]
        self.D0 = None
        self.Ea_D = None
        self.read_einstein()
        
        # random walk
        temp_rand = np.linspace(self.temp[0], self.temp[-1], 1000)
        self.rand = RandomWalk(temp_rand, self.lattice)
        
        # representative nu
        result = minimize_scalar(self.error_nu)
        self.nu_rep = result.x
        self.plot_nu()
        
        # mean hopping barrier from Boltzmann
        self.rand.D_rand(nu=self.nu_rep)
        self.rand.linear_fitting()
        self.Ea_hop_rand = self.rand.Ea
        
        # <Ea_hop>_vhp
        self.Ea_hop_vhp = None
        self.mean_Ea_hop_vhp()
        
        # representative z
        self.z_rep = None
        if self.fix_Ea_t_res:
            self.plot_z_fixed_Ea()
        else:
            self.plot_z() # Ea_hop_vhp will be updated
            
        # activations
        self.Ea_f1 = self.Ea_hop_vhp - self.Ea_hop_rand
        self.Ea_f2 = self.Ea_f - self.Ea_f1
        
        # representative a
        self.a_rep = np.sqrt(6*self.D0/(self.z_rep*self.nu_rep*self.f0))*1e10
        
        # save representative parameters
        self.save_parameter()
        
        # save probability
        self.plot_prob()
        
    def read_neb(self):
        with open(self.neb, 'r', encoding='UTF-8') as f:
            reader = csv.reader(f)
            next(reader) # pass header
            path_neb = [row for row in reader]
            
        for path in path_neb:
            if len(path[0]) > 0:
                path[3] = float(path[3])
                path[4] = float(path[4])
                path[5] = float(path[5])
                path[6] = int(path[6])
                self.lattice.add_path(*path)
                
    def correlation_factor(self):
        for i, temp in enumerate(
                tqdm(self.data.temp, 
                     bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'% (Fore.GREEN, Fore.GREEN),
                     ascii=False,
                     desc=f'{GREEN}{BOLD}Parameter')):
            
            if (self.interval * 1000) % self.data.potim[i] != 0:
                print(f'unvalid interval : interval should be multiple of potim (T={temp}K)')
                print(f'potim = {self.data.potim[i]} fs')
                print(f'interval = {self.interval} ps')
                sys.exit(0)
            else:
                step_interval = int(self.interval*1000 / self.data.potim[i])
                
            path_dir = os.path.join(self.data.prefix1, f'{self.data.prefix2}.{temp}K')
            path_force = None if self.data.force is None else path_dir
            cor = CumulativeCorrelationFactor(xdatcar=path_dir,
                                              force=path_force,
                                              lattice=self.lattice,
                                              temp=temp,
                                              label=self.data.label[i],
                                              interval=step_interval,
                                              verbose=self.verbose)
            self.cor.append(cor)
            
    def save_correlation_factor(self):
        with open('f_cor.txt', 'w', encoding='UTF-8') as f:
            f.write(f'Data        : {self.data.prefix1}\n')
            temp_str = ' '.join(list(map(str, self.temp)))
            f.write(f'Temperature : {temp_str} (in K)\n')
            use_force = False if self.data.force is None else True
            f.write(f'Use FORCE   : {use_force}\n\n')
            
            for i, temp in enumerate(self.temp):
                f.write(f'# T = {temp} K\n')
                f.write("Total counts\n")
                f.write("           ")
                for name in self.cor[i].path_name:
                    f.write("{:<10} ".format(name))
                f.write('\n')
                f.write("Distance   ")
                for dist in self.cor[i].path_dist:
                    f.write("{:<10.3f} ".format(dist))
                f.write('\n')
                f.write("Counts     ")
                for count in self.cor[i].counts_cum * self.cor[i].num_enc_cum:
                    count = int(count)
                    f.write("{:<10} ".format(count))
                f.write('\n\n')
                f.write(f"      Mean correlation factor : {self.cor[i].f_avg:.5f}\n")
                f.write(f"Cumulative correlation factor : {self.cor[i].f_cum:.5f}\n")
                f.write('\n')
                f.write("{:<10} {:<10}\n".format('Label', 'f_cor'))
                for label, f_cor in zip(self.cor[i].label_success, self.cor[i].f_ensemble):
                    f.write("{:<10} {:<10.5f}\n".format(label, f_cor))
                f.write('\n')
                time_tot = np.sum(np.array(self.cor[i].times))
                f.write(f'Total time : {time_tot:.3f} s\n')
                
                if len(self.cor[i].label_err) > 0:
                    f.write('Error occured : ')
                    for label in self.cor[i].label_err:
                        f.write(f'{label} ')
                f.write('\n\n')
                
        print('f_cor.txt is created.')
        
    def plot_correlation_factor(self, slop, intercept):
        plt.style.use('default')
        plt.rcParams['figure.figsize'] = (3.8, 3.8)
        plt.rcParams['font.size'] = 11
        fig, ax = plt.subplots()
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.2)
            
        for i in range(len(self.temp)):
            ax.scatter(self.temp[i], self.f_cum[i], color=self.cmap(i), marker='s', s=50)
        plt.ylim([0, 1])
        plt.xlabel('T (K)', fontsize=14)
        plt.ylabel(r'$f$', fontsize=14)

        # inset graph
        axins = ax.inset_axes([1.125, 0.615, 0.35, 0.35])
        x_ins = np.linspace(1/self.temp[-1], 1/self.temp[0], 100)
        axins.plot(x_ins, slop * x_ins + intercept, 'k:')
        for i in range(len(self.temp)):
            axins.scatter(1/self.temp[i], np.log(self.f_cum[i]), color=self.cmap(i), marker='s')
        axins.set_xlabel('1/T', fontsize=12)
        axins.set_ylabel('ln f', fontsize=12)
        axins.set_xticks([])
        axins.set_yticks([])

        plt.savefig('f_cor.svg', transparent=True, dpi=600, bbox_inches="tight")
        print('f_cor.svg is created.')
    
    def read_einstein(self):
        with open(self.einstein, 'r', encoding='UTF-8') as f:
            lines = [line.strip() for line in f]
            
        for i, line in enumerate(lines):
            if 'D0' in line:
                self.D0 = float(line.split()[2])
            if 'Ea_D' in line:
                self.Ea_D = float(line.split()[2])
            if 'Raw' in line:
                num_line = i+2
                break
            
        for i in range(num_line, num_line+len(self.temp)):
            self.D.append(float(lines[i].split()[-1]))    
        self.D = np.array(self.D)
        
    def error_nu(self, nu):
        rand = RandomWalk(self.temp, self.lattice)
        f_cor = self.f0 * np.exp(-self.Ea_f / (self.kb * self.temp))
        D_neb = rand.D_rand(nu=nu) * f_cor
        return np.sum((self.D - D_neb)**2)
    
    def plot_nu(self):
        power = int(np.log10(self.nu_rep))
        nu_lower = 10 ** (power)
        nu_upper = 10 ** (power+1)
        
        D_lower = self.rand.D_rand(nu=nu_lower)
        D_opt = self.rand.D_rand(nu=self.nu_rep)
        D_upper = self.rand.D_rand(nu=nu_upper)
        
        temp_rand = np.linspace(self.temp[0], self.temp[-1], 1000)
        f_cor = self.f0 * np.exp(-self.Ea_f / (self.kb * temp_rand))
        
        D_lower = f_cor * D_lower
        D_opt = f_cor * D_opt
        D_upper = f_cor * D_upper
        
        plt.style.use('default')
        plt.rcParams['figure.figsize'] = (3.8, 3.8)
        plt.rcParams['font.size'] = 11
        fig, ax = plt.subplots()
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.2)
            
        for i, temp in enumerate(self.temp):
            ax.scatter(1/temp, np.log(self.D[i]), color=self.cmap(i), marker='s', s=50)
            
        plt.plot(1/temp_rand, np.log(D_upper), 'k:', linewidth=1, label=f'{nu_upper :.2e}')
        plt.plot(1/temp_rand, np.log(D_opt), 'k-', linewidth=1, label=f'{self.nu_rep :.2e}')
        plt.plot(1/temp_rand, np.log(D_lower), 'k:', linewidth=1, label=f'{nu_lower :.2e}')
        
        if len(self.temp) >= 3:
            x = np.array([self.temp[0], self.temp[int(len(self.temp)/2)], self.temp[-1]])
        else:
            x = self.temp
            
        x_str = [f'1/{temp}' for temp in x]
        x = 1/x
        plt.xticks(x, x_str)
        plt.xlabel('1/T (1/K)', fontsize=14)
        plt.ylabel(r'ln $D$', fontsize=14)
        
        plt.legend(loc='upper right', fancybox=True, framealpha=1, edgecolor='inherit',
                   labelspacing = 0.3, columnspacing=0.5, borderpad=0.2, handlelength=0.6,
                   fontsize=11)

        plt.savefig('nu_rep.svg', transparent=True, dpi=600, bbox_inches="tight")
        print('nu_rep.svg is created.')
        
    def plot_z(self):
        time = self.data.potim * self.data.nsw / 1000
        count = np.sum(self.count, axis=1)
        num_label = np.array([len(label) for label in self.data.label])
        num_label_err = np.array([len(label) for label in self.label_err])
        count = count / (num_label - num_label_err)
        t_res = time / count

        plt.style.use('default')
        plt.rcParams['figure.figsize'] = (3.8, 3.8)
        plt.rcParams['font.size'] = 11
        fig, ax = plt.subplots()
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.2)
            
        for i, temp in enumerate(self.temp):
            ax.bar(temp, t_res[i], width=50, edgecolor='k', color=self.cmap(i))
            ax.scatter(temp, t_res[i], marker='o', edgecolors='k', color='k')    
            
        # Arrhenius fitting
        slop, intercept = np.polyfit(1/self.temp, np.log(t_res), deg=1)

        x = np.linspace(0.99*self.temp[0], 1.01*self.temp[-1], 1000)
        ax.plot(x, np.exp(intercept)*np.exp(slop/x), 'k:')

        plt.xlabel('T (K)', fontsize=14)
        plt.ylabel(r'<$t_{res}$> (ps)', fontsize=14)
        plt.savefig('t_res.svg', transparent=True, dpi=600, bbox_inches="tight")
        print('t_res.svg is created')
        
        t0 = np.exp(intercept) * 1e-12
        self.z_rep = 1/(self.nu_rep * t0)
        
        # use Ea_act from residence time
        self.Ea_hop_vhp = slop * self.kb
        
        # save result
        with open('t_res.txt', 'w', encoding='UTF-8') as f:
            f.write(f'pre-exponential for t_res = {t0 * 1e12 :.6e} ps\n')
            f.write(f'Ea for t_res = {self.kb * slop :.6f} eV\n')
            f.write(f'representative z = {self.z_rep :.6f}\n\n')
            f.write('Raw data :\n')
            f.write(f'     T(K)  t_res(ps)\n')
            for i, temp in enumerate(self.temp):
                f.write(f'{str(temp):>9s}  {t_res[i]:.3f}\n')   
        print('t_res.txt is created')
            
    def save_parameter(self):
        with open('parameter.txt', 'w', encoding='UTF-8') as f:
            f.write('Effective diffusoin parameters related to D0\n')
            f.write(f"  D0     = {self.D0 :.5e} m2/s # pre-exponential for D\n")
            f.write(f"  nu_eff = {self.nu_rep :.5e} Hz   # jump attempt frequency\n")
            f.write(f"  z_eff  = {self.z_rep :<6.5f}          # coordination number\n")
            f.write(f"  a_eff  = {self.a_rep :.5f} Å        # hopping distance\n")
            f.write(f"  f0     = {self.f0 :.5f}          # pre-exponential for correlation factor\n")
            f.write("\n")
            f.write("Effective diffusion parameters related to Ea\n")
            f.write(f"  Ea_D        = {self.Ea_D :.5f} eV      # obtained from Einstein relation\n")
            f.write(f"  Ea          = {self.Ea_hop_rand :.5f} eV      # hopping barrier (random walk)\n")
            f.write(f"  Ea_act      = {self.Ea_hop_vhp :.5f} eV      # hopping barrier (actual diffusion)\n")
            f.write(f"  Ea_f        = {self.Ea_f :.5e} eV  # activation for correlation factor\n")
            f.write(f"  Ea_f1       = {self.Ea_f1 :.5e} eV  # contribution of deviation from Boltzmann distribution\n")
            f.write(f"  Ea_f2       = {self.Ea_f2 :.5e} eV  # contribution of thermal fluctuation\n")
            
        print('parameter.txt is created.')
        
        
    def plot_prob(self):
        path_name = self.cor[0].path_name[:-1] # exclude unkown

        rand = RandomWalk(self.temp, self.lattice)
        rand.get_probability()
        prob_rand = rand.prob_path

        col = [[] for i in range(len(self.lattice.site_names))]
        name = [[] for i in range(len(self.lattice.site_names))]
        for i, n in enumerate(path_name):
            site = self.lattice.path[i]['site_init']
            idx = self.lattice.site_names.index(site)
            col[idx].append(i)
            name[idx].append(n)

        prob = []
        for i in range(len(self.lattice.site_names)):
            count = self.count[:, np.array(col[i])]
            prob.append(count / np.sum(count, axis=1).reshape(-1, 1))

        bar_width=0.35
        cmap = plt.get_cmap("Set2")
        x = np.arange(len(self.temp))

        for i, p in enumerate(prob): 
            plt.style.use('default')
            plt.rcParams['figure.figsize'] = (6, 3.8)
            plt.rcParams['font.size'] = 12
            fig, ax = plt.subplots()
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(1.2)
            
            for j in range(len(x)):
                for k in range(len(name[i])):
                    label = name[i][k] if j==0 else None
                    ax.bar(x[j]-bar_width/2, p[j][k], bottom=np.sum(p[j, :k]),
                        width=bar_width, edgecolor='k', color=cmap(k), label=label)
                    ax.bar(x[j]+bar_width/2, prob_rand[i][j][k], bottom=np.sum(prob_rand[i][j, :k]),
                        width=bar_width, edgecolor='k', color=cmap(k), hatch='//')
                
            plt.xlabel('T (K)', fontsize=16)
            plt.ylabel('Probability', fontsize=16)
            plt.xticks(x, self.temp, fontsize=14)
            plt.ylim([0, 1])
            plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), 
                       fancybox=True, framealpha=1, 
                       edgecolor='inherit', ncol=1,
                       labelspacing = 0.3, columnspacing=0.5, 
                       borderpad=0.3, handlelength=0.6,
                       fontsize=14)
        
            plt.savefig(f'prob_{self.lattice.site_names[i]}.svg', 
                        transparent=True, dpi=600, bbox_inches="tight")
            print(f'prob_{self.lattice.site_names[i]}.svg is created.')
            plt.close()
            
    def mean_Ea_hop_vhp(self):
        path_name = self.cor[0].path_name[:-1]
        col = [[] for i in range(len(self.lattice.site_names))]
        name = [[] for i in range(len(self.lattice.site_names))]
        for i, n in enumerate(path_name):
            site = self.lattice.path[i]['site_init']
            idx = self.lattice.site_names.index(site)
            col[idx].append(i)
            name[idx].append(n)
        prob, Ea_hop_vhp= [], []
        for i in range(len(self.lattice.site_names)):
            count = self.count[:, np.array(col[i])]
            Ea = np.array([self.lattice.path[idx]['Ea'] for idx in col[i]])
            prob_site = count / np.sum(count, axis=1).reshape(-1, 1)
            prob.append(prob_site)
            Ea_hop_vhp.append(np.average(np.dot(prob_site,Ea)))

        rand = RandomWalk(self.temp, self.lattice)
        rand.get_probability()
        rand.prob_path = prob
        rand.D_rand()
        rand.linear_fitting()
        self.Ea_hop_vhp = rand.Ea


    def plot_z_fixed_Ea(self):
        '''
        Ea for residence time is fixed to <Ea_hop>_vhp.
        '''
        # residence time from MD
        time = self.data.potim * self.data.nsw / 1000
        count = np.sum(self.count, axis=1)
        num_label = np.array([len(label) for label in self.data.label])
        num_label_err = np.array([len(label) for label in self.label_err])
        count = count / (num_label - num_label_err)
        self.t_res = time / count 
        
        # z_rep
        result = minimize_scalar(self.error_z)
        t0 = result.x
        self.z_rep = 1/(self.nu_rep * t0) * 1e12
        
        # plot residence time
        plt.style.use('default')
        plt.rcParams['figure.figsize'] = (3.8, 3.8)
        plt.rcParams['font.size'] = 11
        fig, ax = plt.subplots()
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.2)
            
        for i, temp in enumerate(self.temp):
            ax.bar(temp, self.t_res[i], width=50, edgecolor='k', color=self.cmap(i))
            ax.scatter(temp, self.t_res[i], marker='o', edgecolors='k', color='k')    
       
        x = np.linspace(0.99*self.temp[0], 1.01*self.temp[-1], 1000)
        ax.plot(x, t0*np.exp(self.Ea_hop_vhp/(self.kb*x)), 'k:')
        
        plt.xlabel('T (K)', fontsize=14)
        plt.ylabel(r'<$t_{res}$> (ps)', fontsize=14)
        plt.savefig('t_res.svg', transparent=True, dpi=600, bbox_inches="tight")
        print('t_res.svg is created')
         
        # save result
        with open('t_res.txt', 'w', encoding='UTF-8') as f:
            f.write(f'pre-exponential for t_res = {t0 :.6e} ps\n')
            f.write(f'Ea for t_res (fixed to <Ea_hop>_vhp)= {self.Ea_hop_vhp :.6f} eV\n')
            f.write(f'representative z = {self.z_rep :.6f}\n\n')
            f.write('Raw data :\n')
            f.write(f'     T(K)  t_res(ps)\n')
            for i, temp in enumerate(self.temp):
                f.write(f'{str(temp):>9s}  {self.t_res[i]:.5f}\n')
        print('t_res.txt is created')
        
    def error_z(self, t0):
        self.t_res_vhp = t0 * np.exp(self.Ea_hop_vhp / (self.kb * self.temp))     
        return np.sum((self.t_res_vhp - self.t_res)**2)


class Trajectory:
    def __init__(self,
                 xdatcar,
                 force,
                 poscar,
                 symbol,
                 outcar,
                 interval,
                 correction=True,
                 label=False,
                 verbose=False):
        
        self.xdatcar = xdatcar
        self.force = force
        self.poscar = poscar
        self.symbol = symbol
        self.outcar = outcar
        self.correction = correction
        self.label = label
        self.verbose = verbose
        
        # check file
        if not os.path.isfile(self.xdatcar):
            print(f'{self.xdatcar} is not found.')
            sys.exit(0)
        if not os.path.isfile(self.force):
            print(f'{self.force} is not found.')
            sys.exit(0)
        if not os.path.isfile(self.xdatcar):
            print(f'{self.xdatcar} is not found.')
            sys.exit(0)
        
        # read outcar
        self.potim = None
        self.read_outcar()
        self.interval = int(interval * 1000 / self.potim)
        
        # instantiation
        self.lattice = None
        self.traj = None
        self.get_lattice_and_traj()
        
        # corrections
        if self.correction:
            self.do_correction()
            
        # animation
        self.save_animation()
        
    def read_outcar(self):
        with open(self.outcar, 'r') as f:
            for line in f:
                if 'POTIM' in line:
                    self.potim = float(line.split()[2])
                    break
                
    def get_lattice_and_traj(self):
        self.lattice = Lattice(poscar_perf=self.poscar, 
                               symbol=self.symbol)
        
        self.traj = LatticeHopping(xdatcar=self.xdatcar,
                                   lattice=self.lattice,
                                   force=self.force,
                                   interval=self.interval,
                                   verbose=self.verbose)
        
    def do_correction(self):
        check_multivac, check_TS = True, True
        _traj = copy.deepcopy(self.traj)
        # multi-vacancy
        try:
            self.traj.correct_multivacancy(start=1)
            self.traj.check_multivacancy()
            if self.traj.multi_vac == False:
                print('Correction for multi-vacancy : success')
            else:
                print('Correction for multi-vacancy : fail')
                check_multivac = False
        except:
            print('Correction for multi-vacancy : fail')
            check_multivac = False
        
        # TS criteria
        try:
            self.traj.correct_transition_state()
            print('Correction for TS criteria : success')
        except:
            print('Correction for TS criteria : fail')
            check_TS = False
        
        if not(check_multivac and check_TS):
            print('The raw trajectory without corrections will be used.')
            self.traj = _traj
            
    def save_animation(self):
        print('\nInformation on animation')
        print(f'  NSW = {self.traj.nsw} ({self.traj.nsw * self.potim / 1000} ps)')
        print(f'  Interval = {self.interval} step ({self.interval * self.potim / 1000} ps)')
        print(f'  Total step = {self.traj.num_step} (={self.traj.nsw}/{self.traj.interval})')
        print('')
        step = input('Enter init and final steps (int; ex. 0 100 / 0 -1 for all): ')
        try:
            step = list(map(int, step.split()))
        except:
            print('The step number must be integer.')
            sys.exit(0)
            
        if step[-1] > self.traj.num_step:
            print(f'    The final step should be less than {self.traj.num_step}')
            print(f'    The final step is set to {self.traj.num_step}')
            step[-1] = self.traj.num_step
        step = 'all' if step[-1]==-1 else np.arange(step[0], step[-1])
        
        fps = input('Enter fps (int; ex. 20): ')
        try:
            fps = int(fps)
        except:
            print('The fps must be integer.')
        
        print('')
        self.traj.animation(step=step,
                            potim=self.potim,
                            foldername='snapshot',
                            fps=fps,
                            dpi=300,
                            label=self.label)
        
        print('snapshot directory was created.')
        
        
class PathAnalyzer:
    def __init__(self,
                 poscar,
                 symbol,
                 xdatcar,
                 outcar,
                 neb,
                 force,
                 interval,
                 verbose=False):
        
        self.poscar = poscar
        self.symbol = symbol
        self.xdatcar = xdatcar
        self.outcar = outcar
        self.neb = neb
        self.force = force
        self.verbose = verbose

        # check file
        if not os.path.isfile(self.poscar):
            print(f'{self.poscar} is not found.')
            sys.exit(0)
        if not os.path.isfile(self.neb):
            print(f'{self.neb} is not found.')
            sys.exit(0)
        if not os.path.isfile(self.xdatcar):
            print(f'{self.xdatcar} is not found.')
            sys.exit(0)
        if not os.path.isfile(self.outcar):
            print(f'{self.outcar} is not found.')
            sys.exit(0)
        if not os.path.isfile(self.force):
            print(f'{self.force} is not found.')
            sys.exit(0)
        
        # read outcar
        self.potim = None
        self.read_outcar()
        self.interval = int(interval * 1000 / self.potim)
        
        # lattice
        self.lattice = Lattice(self.poscar, symbol=self.symbol)
        self.read_neb()
        
        # traj
        self.traj = LatticeHopping(xdatcar=self.xdatcar,
                                   lattice=self.lattice,
                                   force=self.force,
                                   interval=self.interval,
                                   verbose=self.verbose)
        self.do_correction()
        
        # analyzer
        self.analyzer = Analyzer(traj=self.traj,
                                 lattice=self.lattice,
                                 verbose=self.verbose)
        self.analyzer.get_path_vacancy(verbose=self.verbose)
        try:
            self.analyzer.correct_multipath()
            print('Correction for multi-path : success')
        except:
            print('Correction for multi-vacancy : fail')
        print('')
        
        # print summary
        print(f'{RED}{BOLD}Summary{RESET}')
        self.analyzer.print_summary(disp=False,
                                    save_figure=True,
                                    save_text=True)
        
        if len(self.analyzer.step_unknown) > 0:
            print('unknown steps :', end=' ')
            for step in self.analyzer.step_unknown:
                print(step, end=' ')
            print('')

        print('')
        print('counts.png is created.')
        print('counts.txt is created.')
        
    def read_outcar(self):
        with open(self.outcar, 'r') as f:
            for line in f:
                if 'POTIM' in line:
                    self.potim = float(line.split()[2])
                    break    
        
    def read_neb(self):
        with open(self.neb, 'r', encoding='UTF-8') as f:
            reader = csv.reader(f)
            next(reader) # pass header
            path_neb = [row for row in reader]
            
        for path in path_neb:
            if len(path[0]) > 0:
                path[3] = float(path[3])
                path[4] = float(path[4])
                path[5] = float(path[5])
                path[6] = int(path[6])
                self.lattice.add_path(*path)
                
    def do_correction(self):
        check_multivac, check_TS = True, True
        # multi-vacancy
        try:
            self.traj.correct_multivacancy(start=1)
            self.traj.check_multivacancy()
            if self.traj.multi_vac == False:
                print('Correction for multi-vacancy : success')
            else:
                print('Correction for multi-vacancy : fail')
                check_multivac = False
        except:
            print('Correction for multi-vacancy : fail')
            check_multivac = False
        
        # TS criteria
        try:
            self.traj.correct_transition_state()
            print('Correction for TS criteria : success')
        except:
            print('Correction for TS criteria : fail')
            check_TS = False
            
        if not(check_multivac and check_TS):
            print('Correction Failure : VacHopPy will be terminated.')
            sys.exit(0)
            
            

class GetFingerPrint:
    def __init__(self, 
                 poscar,
                 prefix='fingerprint', 
                 disp=False):
        self.poscar = poscar
        self.prefix = prefix
        self.disp = disp
        
        if not os.path.isfile(self.poscar):
            print(f'{self.poscar} is not found')
            sys.exit(0)
        
        # read poscar
        self.atom = None
        self.read_poscar()
        
        # input parameters
        self.pair = []
        self.get_pair()
        
        # params for fingerprint
        self.Rmax = None
        self.delta = None
        self.sigma = None
        self.get_params()
        
        # fingerprint
        self.fingerprint = []
        self.get_fingerprint()
        self.fingerprint = np.array(self.fingerprint)
        
        # concat fingerprints
        self.fingerprint_concat = self.fingerprint.reshape(1,-1).squeeze()
        
        # save fingerprint
        self.save_fingerprint()
        
        
    def read_poscar(self):
        with open(self.poscar, 'r') as f:
            lines = [line for line in f]
        self.atom = lines[5].split()
    
    def get_pair(self):
        pair = input('input A and B (ex. Hf-O / Hf-Hf,Hf-O / all) : ')
        pair = pair.replace(" ", "")
        
        if pair == 'all':
            self.pair.extend(combinations_with_replacement(self.atom, 2))
        else:
            pair = pair.split(',')
            for p in pair:
                atoms = p.split('-')
                if len(atoms) == 2 and all(atom in self.atom for atom in atoms):
                    self.pair.append(tuple(atoms))
                else:
                    print(f'Invalid pair : {p}')
                    sys.exit(0)
                    
    def get_params(self):
        params = input("input Rmax, delta, and sigma (ex. 15, 0.01, 0.3) : ")
        params = list(map(float, params.replace(" ", "").split(',')))
        self.Rmax, self.delta, self.sigma = params[0], params[1], params[2]
        
    def get_fingerprint(self):
        for (A, B) in self.pair:
            finger = FingerPrint(A, B, 
                                 self.poscar,
                                 self.Rmax, self.delta, self.sigma)
            self.fingerprint.append(finger.fingerprint)
    
    def save_fingerprint(self):
        # save figure
        R = np.linspace(0, self.Rmax, len(self.fingerprint[0]))
        for i in range(len(self.pair)):
            x = R + i*self.Rmax*np.ones_like(R)
            plt.plot(x, self.fingerprint[i], label=f'{self.pair[i][0]}-{self.pair[i][1]}')
        
        plt.axhline(0, 0, 1, color='k', linestyle='--', linewidth=1)
        plt.xlabel("r (Å)", fontsize=13)    
        plt.ylabel('Intensity', fontsize=13)
        plt.legend(fontsize=12)
        plt.savefig(f'{self.prefix}.png', dpi=300)
        print(f'{self.prefix}.png is created.')
        if self.disp:
            plt.show()
        plt.close()
        
        R = np.linspace(0, self.Rmax * len(self.pair), len(self.fingerprint_concat))
        with open(f'{self.prefix}.txt', 'w') as f:
            f.write(f'# Rmax, delta, sigma = {self.Rmax}, {self.delta}, {self.sigma}\n')
            f.write('# pair : ')
            for (A, B) in self.pair:
                f.write(f'{A}-{B}, ')
            f.write('\n')
            for x, y in zip(R, self.fingerprint_concat):
                f.write(f'  {x:2.6f}\t{y:2.6f}\n')
        print(f'{self.prefix}.txt is created.')

    
class GetCosineDistance:
    def __init__(self, fp1, fp2):
        if not os.path.isfile(fp1):
            print(f'{fp1} is not found')
            sys.exit(0)
        if not os.path.isfile(fp2):
            print(f'{fp2} is not found')
            sys.exit(0)
        
        self.fp1 = np.loadtxt(fp1, skiprows=2)[:,1]
        self.fp2 = np.loadtxt(fp2, skiprows=2)[:,1]
        
        if self.fp1.shape != self.fp2.shape:
            print('Size of two fingerprints should be the same')
            sys.exit(0)
        
        self.d_cos = CosineDistance(self.fp1, self.fp2)
        print(f'd_cos = {self.d_cos}')
        
        
def find_last_direct_line(xdatcar):
    try:
        with open(xdatcar, 'rb') as f:
            f.seek(0, os.SEEK_END)
            file_size = f.tell()
            buffer_size = 1024
            buffer = b""
            offset = file_size
            while offset > 0:
                read_size = min(buffer_size, offset)
                offset -= read_size
                f.seek(offset)
                chunk = f.read(read_size)
                buffer = chunk + buffer 
                lines = buffer.split(b'\n')
                buffer = lines[0]
                for line in reversed(lines[1:]):
                    decoded_line = line.decode('utf-8', errors='ignore').strip()
                    if 'Direct' in decoded_line:
                        final_step = int(decoded_line.split('=')[-1])
                        return final_step
            if buffer:
                decoded_line = buffer.decode('utf-8', errors='ignore').strip()
                if 'Direct' in decoded_line:
                    final_step = int(decoded_line.split('=')[-1])
                    return final_step       
        print("No line containing 'Direct' found")
    except FileNotFoundError:
        print(f"{xdatcar} is not found")
    except Exception as e:
        print(f"An error occurred: {e}")


class PhaseTransition:
    def __init__(self,
                 xdatcar,
                 outcar,
                 interval,
                 poscar_mother,
                 prefix1='poscars',
                 prefix2='fingerprints'):
        
        self.xdatcar = xdatcar
        self.outcar = outcar
        self.interval = interval
        self.poscar_mother = poscar_mother
        self.prefix1 = prefix1
        self.prefix2 = prefix2
        
        if not os.path.isfile(self.xdatcar):
            print(f'{self.xdatcar} is not found')
            sys.exit(0)
            
        if not os.path.isfile(self.outcar):
            print(f'{self.outcar} is not found')
            sys.exit(0)
            
        if not os.path.isfile(self.poscar_mother):
            print(f'{self.poscar_mother} is not found')
            sys.exit(0)
        
        if not os.path.isdir(self.prefix1):
            os.makedirs(self.prefix1)
            print(f'{self.prefix1} directory is created.')
        
        if not os.path.isdir(self.prefix2):
            os.makedirs(self.prefix2)
            print(f'{self.prefix2} directory is created.')
            
        # read outcar
        self.potim = None
        self.read_outcar()
        self.interval = int(interval * 1000 / self.potim)
        
        # read xdatcar
        self.nsw = None
        self.position = []
        self.read_xdatcar()
        
        self.pair = []
        self.pair.extend(combinations_with_replacement(self.atom_species, 2))
        
        # input steps
        self.step = None
        self.input_step()
        
        # params
        self.Rmax = None
        self.delta = None
        self.sigma = None
        self.get_params()
        
        # fingerprint of mother phase
        self.finger_mother = self.get_fingerprint(self.poscar_mother)
        self.R = np.linspace(0, self.Rmax*len(self.pair), len(self.finger_mother))
        self.save_fingerprint(self.finger_mother, 'fingerprint_mother.txt')
        
        # fingerprint of xdatcar
        self.d_cos = np.zeros_like(self.step, dtype=float)
        desc = 'progress'
        print('')
        for i, s in enumerate(tqdm(self.step, 
                                   bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'% (Fore.GREEN, Fore.GREEN),
                                   ascii=False,
                                   desc=f'{GREEN}{BOLD}{desc}')):
            # save poscar
            label = format(s, self.digit)
            filename = os.path.join(self.prefix1, f"POSCAR_{label}")
            self.save_poscar(s, filename)
            
            # get fingerpint
            finger = self.get_fingerprint(filename)
            finger_out = f'fingerprint_{label}.txt'
            self.save_fingerprint(finger, finger_out)
            
            # cosine distance
            self.d_cos[i] = CosineDistance(self.finger_mother, finger)
        
        # plot cosine distance
        print('')
        self.save_distance()
            
    def read_outcar(self):
        with open(self.outcar, 'r') as f:
            for line in f:
                if 'POTIM' in line:
                    self.potim = float(line.split()[2])
                    break
                
    def read_xdatcar(self):
        self.nsw = find_last_direct_line(self.xdatcar)
        if self.nsw % self.interval != 0:
            print(f'{self.nsw} step (nsw) is not divided by {self.interval} step (interval)')
            sys.exit(0)
        
        with open(self.xdatcar, 'r') as f:
            lines = np.array([s.strip() for s in f])
        
        self.lattice = np.array([s.split() for s in lines[2:5]], dtype=float)
        self.lattice *= float(lines[1])
        
        self.atom_species = np.array(lines[5].split())
        self.num_species = len(self.atom_species)
        
        self.num_atoms = np.array(lines[6].split(), dtype=int)
        num_atoms_tot = np.sum(self.num_atoms)
        self.num_step = int(self.nsw / self.interval)
        
        digit = int(np.log10(self.num_step)) + 1
        self.digit = f'0{digit}'
        
        for i, spec in enumerate(self.atom_species):           
            atom = {}
            atom['species'] = spec
            atom['num'] = self.num_atoms[i]
            
            traj = np.zeros((atom['num'], self.num_step, 3)) 

            for j in range(atom['num']):
                start = np.sum(self.num_atoms[:i]) + j + 8
                end = lines.shape[0] + 1
                step = num_atoms_tot + 1
                coords = [s.split() for s in lines[start:end:step]]
                coords = np.array(coords, dtype=float)
                
                displacement = np.zeros_like(coords)
                displacement[0,:] = 0
                displacement[1:,:] = np.diff(coords, axis=0)

                # correction for periodic boundary condition
                displacement[displacement>0.5] -= 1.0
                displacement[displacement<-0.5] += 1.0
                displacement = np.cumsum(displacement, axis=0)
                coords = coords[0] + displacement

                # averaged coordination
                coords = coords.reshape(self.num_step, self.interval, 3)
                coords = np.average(coords, axis=1)

                # wrap back into cell
                coords = coords - np.floor(coords)
                traj[j] = coords

            atom['traj'] = traj
            self.position += [atom]
        
    def save_poscar(self, step, filename):
        with open(filename, 'w') as f:
            f.write(f"step_{step}. generated by vachoppy.\n")
            f.write("1.0\n")

            for lat in self.lattice:
                f.write("%.6f %.6f %.6f\n"%(lat[0], lat[1], lat[2]))

            for atom in self.position:
                f.write(f"{atom['species']} ")
            f.write('\n')
            for atom in self.position:
                f.write(f"{atom['num']} ")
            f.write('\n')
            f.write("Direct\n")
            for atom in self.position:
                for traj in atom['traj'][:,step,:]:
                    f.write("%.6f %.6f %.6f\n"%(traj[0], traj[1], traj[2]))
            
    def input_step(self):
        print('')
        step = input(f'Input range of step (ex. 10-100) (available: 0 - {self.num_step}): ')
        step = list(map(int, step.replace(" ", "").split('-')))
        self.step = np.arange(step[0], step[1])
        
    def get_params(self):
        params = input("input Rmax, delta, and sigma (ex. 15, 0.01, 0.3) : ")
        params = list(map(float, params.replace(" ", "").split(',')))
        self.Rmax, self.delta, self.sigma = params[0], params[1], params[2]
        
    def get_fingerprint(self, poscar):
        fingerprint = []
        for (A, B) in self.pair:
            finger = FingerPrint(A, B, poscar, self.Rmax, self.delta, self.sigma)
            fingerprint.append(finger.fingerprint)
        fingerprint = np.array(fingerprint)
        return fingerprint.reshape(1,-1).squeeze()
        
    def save_fingerprint(self, fp, filename):
        file_out = os.path.join(self.prefix2, filename)
        with open(file_out, 'w') as f:
            f.write(f'# Rmax, delta, sigma = {self.Rmax}, {self.delta}, {self.sigma}\n')
            f.write('# pair : ')
            for (A, B) in self.pair:
                f.write(f'{A}-{B}, ')
            f.write('\n')
            for x, y in zip(self.R, fp):
                f.write(f'  {x:2.6f}\t{y:2.6f}\n')
                
    def save_distance(self):
        # plt.figure(figsize=(15, 5))
        plt.scatter(self.step, self.d_cos, s=25)
        plt.xlabel("Step", fontsize=13)
        plt.ylabel('Cosine distnace', fontsize=13)
        # plt.ylim([-0.05, 1.05])
        plt.savefig('cosine_distance.png', dpi=300)
        plt.show()
        plt.close()
        print('cosine_distance.png is created.')
        
        with open('cosine_distance.txt', 'w') as f:
            f.write(f'# Rmax, delta, sigma = {self.Rmax}, {self.delta}, {self.sigma}\n')
            f.write('# pair : ')
            for (A, B) in self.pair:
                f.write(f'{A}-{B}, ')
            f.write('\n')
            for x, y in zip(self.step, self.d_cos):
                f.write(f'  {x}\t{y:.6f}\n')
        print('cosine_distance.txt is created.')