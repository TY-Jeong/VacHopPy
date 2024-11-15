import os
import sys
import csv
import pandas as pd
import numpy as np
from tqdm import tqdm
from colorama import Fore
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

from vachoppy.inout import DataInfo
from vachoppy.einstein import *
from vachoppy.trajectory import *

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
            f.write(f'x_vac = {self.x_vac :.3f}\n')
            if self.x_vac == 1:
                f.write('target : atom\n')
            else:
                f.write('target : vacancy\n')
            f.write('\nParameters for diffusion coefficient : \n')
            f.write(f'  D0 = {np.exp(intercept)} m2/s\n')
            f.write(f'  Ea = {-kb * slop} eV\n')
            f.write('\n')
            f.write('Raw data\n')
            f.write('T(K) \tD(m2/s)\n')
            for temp, D in zip(self.data.temp, self.D):
                f.write(f'{temp} \t{D :.6e}\n')
                   
        print('Einstein.txt is created.')
        print('parameters for diffusion coefficient : ')
        print(f'  D0 = {np.exp(intercept) :.3e} m2/s')
        print(f'  Ea = {-kb * slop :.3f} eV')
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
                 fix_Ea_t_res=False,
                 tolerance_Ea_f=0.1):
        
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
        
        print(f'{CYAN}{BOLD}Representative mass transport parameters.{RESET}')
        
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
            self.plot_z()
        
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
                f.write(f"      Mean correlation factor : {self.cor[i].f_avg:.3f}\n")
                f.write(f"Cumulative correlation factor : {self.cor[i].f_cum:.3f}\n")
                f.write('\n')
                f.write("{:<10} {:<10}\n".format('Label', 'f_cor'))
                for label, f_cor in zip(self.cor[i].label_success, self.cor[i].f_ensemble):
                    f.write("{:<10} {:<10.3f}\n".format(label, f_cor))
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
            if 'Ea' in line:
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
            f.write('Representative parameters related to D0\n')
            f.write(f"  D0     = {self.D0 :.3e} m2/s # obtained from Einstein relation\n")
            f.write(f"  nu_rep = {self.nu_rep :.3e} Hz   # attempt frequency\n")
            f.write(f"  z_rep  = {self.z_rep :<6.3f}         # coordination number\n")
            f.write(f"  a_rep  = {self.a_rep :.3f} Å        # hopping distance\n")
            f.write(f"  f0     = {self.f0 :.3f}          # pre-exponential for correlation factor\n")
            f.write("\n")
            f.write("Representative Parameters related to Ea\n")
            f.write(f"  Ea_D        = {self.Ea_D :.3f} eV  # obtained from Einstein relation\n")
            f.write(f"  Ea_hop_rand = {self.Ea_hop_rand :.3f} eV  # hopping barrier averaged using Boltzmann\n")
            f.write(f"  Ea_hop_vhp = {self.Ea_hop_vhp :.3f} eV   # hopping barrier averaged using probabilities from vachoppy\n")
            f.write(f"  Ea_f        = {self.Ea_f :.3e} eV  # activation for correlation factor\n")
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
                f.write(f'{str(temp):>9s}  {self.t_res[i]:.3f}\n')
        print('t_res.txt is created')
        
    def error_z(self, t0):
        self.t_res_vhp = t0 * np.exp(self.Ea_hop_vhp / (self.kb * self.temp))     
        return np.sum((self.t_res_vhp - self.t_res)**2)
            
