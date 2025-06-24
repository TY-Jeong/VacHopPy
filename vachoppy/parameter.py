import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from tabulate import tabulate
from scipy.optimize import minimize_scalar
from vachoppy.trajectory import *


class Calculator:
    def __init__(self,
                 data,
                 temp,
                 label,
                 lattice,
                 interval,
                 num_vac,
                 tolerance: float = 1e-3,
                 use_incomplete_encounter: bool = True):
        """
        Arguments
        ---------
        data : DataInfo
            DataInfo object
        temp : float
            Temperature (K)
        label : str
            Label
        lattice : Lattice
            Lattice object
        interval : float
            Time interval (ps)
        num_vac : int
            Number of vacancies
        tolerance : flaot, optional
            Tolerance for numerical accuracy
        use_incomplete_encounter : bool, optional
            If true, incomplete encounters are used together in computations.
        """
        
        # arguments
        self.data = data
        self.temp = temp
        self.label = label
        self.lattice = lattice
        self.interval = interval
        self.num_vac = num_vac
        self.tolerance = tolerance
        self.use_incomplete_encounter = use_incomplete_encounter
        
        self.index_temp = None
        self.index_label = None
        self.pos = None
        self.cond = None
        self.force = None
        self.get_address()
        
        # quantities from trajectory
        self.potim = None
        self.num_step = None
        self.transient_vacancy = None
                
        # quantities from analyzer
        self.num_path = len(lattice.path_name)
        self.hopping_history = None
        self.counts = None
        self.path_unknown = None
        self.unknown_name = None
        self.counts_unknown = None
        self.residence_time = None
        self.msd_rand = None
        
        # quantities from encounter
        self.f_cor = None
        self.encounter_num = None
        self.encounter_msd = None
        self.encounter_path_names = None
        self.encounter_path_counts = None
        self.encounter_path_distance = None
        self.extract_AIMD_quantities()
        
        # check success
        self.success = True
        self.fail_reason = None
        
    def get_address(self):
        for i, temp_i in enumerate(self.data.temp):
            if abs(temp_i - self.temp) < self.tolerance:
                self.index_temp = i
                break
        
        for i, label_i in enumerate(self.data.label[self.index_temp]):
            if label_i == self.label:
                self.index_label = i
                break
            
        if None in [self.index_temp, self.index_label]:
            print(f"No data exists (T : {self.temp} K, Label : {self.label})")
            sys.exit(0)
            
        self.pos = self.data.pos[self.index_temp][self.index_label]
        self.cond = self.data.cond[self.index_temp][self.index_label]
        self.force = self.data.force[self.index_temp][self.index_label]
         
    def extract_AIMD_quantities(self):
        # trajectory
        try:
            trajectory = Trajectory(
                interval=self.interval,
                num_vac=self.num_vac,
                lattice=self.lattice,
                pos_file=self.pos,
                force_file=self.force,
                cond_file=self.cond,
                verbose=False
            )
        
        except SystemExit:
            self.success = False
            self.fail_reason = "Error by trajectory.Trajectory"
            return
        
        except BaseException as e:
            self.success = False
            self.fail_reason = "Error by trajectory.Trajectory"
            print(e)
            return
        
        self.potim = trajectory.potim
        self.num_step = trajectory.num_step
        self.transient_vacancy = trajectory.transient_vacancy

        # analyzer
        try:
            analyzer = TrajectoryAnalyzer(
                lattice=self.lattice,
                trajectory=trajectory,
                tolerance=self.tolerance,
                verbose=False
            )
        except SystemExit:
            self.success = False
            self.fail_reason = "Error by trajectory.TrajectoryAnalyzer"
            return
        
        except BaseException as e:
            self.success = False
            self.fail_reason = "Error by trajectory.TrajectoryAnalyzer"
            print(e)
            return
        
        self.hopping_history = analyzer.hopping_history
        self.counts = analyzer.counts
        self.path_unknown = analyzer.path_unknown
        self.unknown_name = analyzer.unknown_name
        self.counts_unknown = analyzer.counts_unknown
        self.residence_time = analyzer.residence_time
        self.msd_rand = analyzer.msd_rand
        
        # encounter
        try:
            encounter = Encounter(
                analyzer=analyzer,
                use_incomplete_encounter=self.use_incomplete_encounter,
                verbose=False
            )
        except SystemExit:
            self.success = False
            self.fail_reason = "Error by trajectory.TrajectoryAnalyzer"
            return
        
        except BaseException as e:
            self.success = False
            self.fail_reason = "Error by trajectory.Encounter"
            print(e)
            return
        
        self.f_cor = encounter.f_cor
        self.encounter_num = encounter.num_encounter
        self.encounter_msd = encounter.msd
        self.encounter_path_names = encounter.path_name
        self.encounter_path_counts = encounter.path_count
        self.encounter_path_distance = encounter.path_distance

        
class Calculator_fail:
    def __init__(self, 
                 data, 
                 temp,
                 label):
        self.data = data
        self.temp = temp
        self.label = label
        self.success = False
        self.fail_reason = 'Unknown reason'


class ParameterExtractor:
    def __init__(self,
                 data,
                 results : list,
                 tolerance: float = 1e-3,
                 verbose: bool = True,
                 figure: bool = True,
                 file_out: str = 'parameter.txt',
                 inset_correlatoin_factor: bool = True):
        """
        Arguments
        ---------
        data : DataInfo
            DataInfo object
        results : list
            List of Calculator objects
        tolerance : float, optional
            Numerial tolerance
        verbose : bool, optional
            Verbosity flag
        figure : bool, optional
            If True, figures will generated
        file_out : str, optional
            File name which the computational summary will be saved in
        inset_correlatoin_factor : bool, optional
            If True, Arrhenius plot for correlation factor is shown together
        """
        # arguments
        self.data = data
        self.results = results
        self.lattice = results[-1].lattice
        self.tolerance = tolerance
        self.verbose = verbose
        self.figure = figure
        
        self.num_vac = results[-1].num_vac
        self.temp = np.array(self.data.temp, dtype=np.float64)
        self.num_temp = len(self.temp)
        
        # constants
        self.kb = 8.61733326e-5
        self.cmap = plt.get_cmap("Set1")
        
        # unknown paths
        self.path_unknown = self.lattice.path_unknown
        
        # labels for successful results
        self.num_label = [
            [result.temp for result in self.results].count(temp) 
            for temp in self.temp
        ]
        self.index = np.concatenate(([0], np.cumsum(self.num_label)))
        
        # correlation factor
        self.f_ind = None
        self.f_avg = None
        self.f_cum = None
        self.correlation_factor()
        
        # diffusion coefficient
        self.D_rand = None
        self.D0_rand = None
        self.Ea = None
        self.random_walk_diffusion_coefficient()
        
        # residence time
        self.tau = None
        self.tau0 = None
        self.residence_time()
        
        # hopping distance
        self.a_eff = np.sqrt(6 * self.D0_rand * self.tau0) * 1e4
        
        # <z>
        self.z_mean = None
        self.mean_number_of_equivalent_paths()
        
        if verbose:
            if file_out is None:
                self.summary()
            else:
                with open(file_out, 'w') as f:
                    original_stdout = sys.stdout
                    sys.stdout = f
                    try:
                        self.summary()
                    finally:
                        sys.stdout = original_stdout
        
        if figure:
            self.save_figure(inset_correlatoin_factor=inset_correlatoin_factor)  
        
    def correlation_factor(self):
        self.f_ind = []
        self.f_avg = []
        self.f_cum = []
        
        # individual correlatoin factor
        self.f_ind = np.array([result.f_cor for result in self.results], dtype=np.float64)
        
        # average at each temperature
        self.f_avg = [
            np.mean([f for f in self.f_ind[start:end] if not np.isnan(f)])
            for start, end in zip(self.index[:-1], self.index[1:])
        ]
        
        # cumulative correlatoin factor
        for i in range(self.num_temp):
            num_encounter = []
            msd_encounter = []
            msd_encounter_rand = []
            
            for j in range(self.index[i], self.index[i+1]):
                
                if self.results[j].f_cor is None:
                    continue
                
                num_encounter.append(self.results[j].encounter_num)
                msd_encounter.append(self.results[j].encounter_msd)
                msd_encounter_rand.append(
                    np.sum(
                        self.results[j].encounter_path_distance**2 * self.results[j].encounter_path_counts
                    )
                )
                
            msd_encounter = np.array(msd_encounter)
            num_encounter = np.array(num_encounter)
            msd_encounter_rand = np.array(msd_encounter_rand)

            msd_encounter = np.sum(msd_encounter * num_encounter) / np.sum(num_encounter)
            msd_encounter_rand = np.sum(msd_encounter_rand) / np.sum(num_encounter)
            
            f_cum_i = msd_encounter / msd_encounter_rand
            self.f_cum.append(f_cum_i)

    def random_walk_diffusion_coefficient(self):
        # random walk diffusion coefficient
        self.D_rand = []
        for start, end in zip(self.index[:-1], self.index[1:]):
            t_i = np.sum(
                np.array(
                    [result.interval * result.num_step for result in self.results[start:end]]
                )
            )
            msd_rand_i = np.sum(
                np.array(
                    [result.msd_rand for result in self.results[start:end]]
                )
            )
            self.D_rand.append(msd_rand_i / (6 * t_i))
        self.D_rand = np.array(self.D_rand) * 1e-8
        
        # Arrhenius fit
        slop, intercept = np.polyfit(1/self.temp, np.log(self.D_rand), deg=1)
        self.Ea = -slop * self.kb
        self.D0_rand = np.exp(intercept)
        
    def residence_time(self):
        self.tau = []
        for start, end in zip(self.index[:-1], self.index[1:]):
            t_i = np.sum(
                np.array(
                    [result.interval * result.num_step for result in self.results[start:end]]
                )
            )
            count_i = np.sum(
                [np.sum(result.counts) / result.num_vac for result in self.results[start:end]]
            )
            self.tau.append(t_i / count_i)
        self.tau = np.array(self.tau)
        
        # Arrhenius fit
        error_tau = lambda tau0: np.linalg.norm(
            self.tau - tau0 * np.exp(self.Ea / (self.kb * self.temp))
        )
        tau0_opt = minimize_scalar(error_tau)
        self.tau0 = tau0_opt.x # ps
            
    def mean_number_of_equivalent_paths(self):
        self.z_mean = []
        z = np.array(
            [path['z'] for path in self.lattice.path], dtype=np.float64
        )
        
        for start, end in zip(self.index[:-1], self.index[1:]):
            counts = np.array(
                [np.sum(result.counts, axis=0) for result in self.results[start:end]]
            )
            counts = np.sum(counts, axis=0)
            
            self.z_mean.append(np.sum(counts) / np.sum(counts / z))
            
    def print_lattice_info(self):
        path_site = [path['site_init'] for path in self.lattice.path]
        num_path_site = [path_site.count(name) for name in self.lattice.site_name]
        
        print('Lattice information :')
        print('  Vacancy type : ', self.lattice.symbol)
        print(f"  Number of vacancies : {self.num_vac}")
        print('  Number of sites : ', len(self.lattice.site_name))
        print('  Number of hopping paths : ', end='')
        for num in num_path_site:
            print(int(num), end=' ')
        print('\n  Number of unknown paths =', len(self.path_unknown))
        print('')
        
        print('Vacancy hopping paths : ')
        header = ['path', 'a (Å)', 'z', 'initial site', 'final site']
        data = [
            [path['name'],  
             f"{path['distance']:.3f}", path['z'], 
             path['site_init'] + f" [{path['coord_init'][0]:.5f} {path['coord_init'][1]:.5f} {path['coord_init'][2]:.5f}]",
             path['site_final'] + f" [{path['coord_final'][0]:.5f} {path['coord_final'][1]:.5f} {path['coord_final'][2]:.5f}]"]
            for path in self.lattice.path
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        print('Unknown hopping paths : ')
        header = ['path', 'a (Å)', 'z', 'initial site', 'final site']
        data = [
            [path['name'],  
             f"{path['distance']:.3f}", '-', 
             path['site_init'] + f" [{path['coord_init'][0]:.5f} {path['coord_init'][1]:.5f} {path['coord_init'][2]:.5f}]",
             path['site_final'] + f" [{path['coord_final'][0]:.5f} {path['coord_final'][1]:.5f} {path['coord_final'][2]:.5f}]"]
            for path in self.path_unknown
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
    def print_simulation_condition(self):
        print(f'Time ineterval used for average (t_interval) : {self.results[0].interval} ps')
        print('')
        
        print('Simulation temperatures (K) : ', end='')
        for temp in self.temp:
            print(temp, end=' ')
        print('\n')   
        
        print('Labels of AIMD data :')
        label_all = list(
            {element for sublabel in self.data.label for element in sublabel}
        )
        label_all.sort()
        
        header = ['label'] + [f"{temp}K" for temp in self.temp]
        data = [
            [label] + ['O' if label in self.data.label[j] else 'X' for j in range(len(self.temp))]
            for label in label_all
        ]
        data.append(['Num'] + self.num_label)
        
        print(tabulate(data, headers=header, tablefmt="simple", stralign='center', numalign='center'))
        print('')
        
    def print_effective_hopping_parameter(self):
        print('Effective hopping parameters : ')
        header = ['parameter', 'value', 'description']
        parameter = [
            'Drand_0 (m2/s)', 'tau0 (ps)', 'Ea (eV)', 'a (Å)', 'f', '<z>'
        ]
        value = [
            f"{self.D0_rand:.5e}", 
            f"{self.tau0:.5f}", 
            f"{self.Ea:.5f}", 
            f"{self.a_eff:.5f}",
            f"{np.average(self.f_cum):.5f}", 
            f"{np.average(self.z_mean):.5f}"
        ]
        description = [
            'pre-exponential for random walk diffusivity',
            'pre-exponential for residence time',
            'hopping barrier',
            'hopping distance',
            'correlation factor',
            'mean number of equivalent paths per path type'
        ]
        data = [[p, v, d] for p, v, d in zip(parameter, value, description)]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')   

    def print_diffusion_coefficient(self):
        print('Random walk diffusion coefficient : ')
        header = ['T (K)', 'D_rand (m2/s)']
        data = [
            [f"{temp}", f"{D:.5e}"] for temp, D in zip(self.temp, self.D_rand)
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
    
    def print_residence_time(self):
        print('Residence time : ')
        header = ['T (K)', 'tau (ps)']
        data = [
            [f"{temp}", f"{tau:.5f}"] for temp, tau in zip(self.temp, self.tau)
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
    
    def print_correlation_factor(self):
        label_all = list(
            {element for sublabel in self.data.label for element in sublabel}
        )
        label_all.sort()
        
        # rearrange f_ind
        f_ind_re = []
        index = [0] + list(np.cumsum(self.num_label))
        for i, f in enumerate(range(len(self.temp))):
            f_ind_re_i = ['-'] * len(label_all)
            for j, cal in enumerate(self.results[index[i]:index[i+1]]):
                f_ind_re_i[label_all.index(cal.label)] = f"{self.f_ind[index[i]+j]:.5f}"
            f_ind_re.append(f_ind_re_i)
        f_ind_re = [list(x) for x in zip(*f_ind_re)] # transpose
        
        # check use_incomplete_encounter
        if self.results[-1].use_incomplete_encounter:
            print("Caution! use_incomplete_encounter = True was set.")
            print("Incompleted encounters are also used for correlatoin factor calculations")
            print("use_incomplete_encounter = True is recommended only when simulation time is not enough to catch sufficient encounters")
            print('')
            
        print('Cumulative correlation factors : ')
        print('(Note: use these values for your work)')
        header = ['T (K)', 'f']
        data = [
            [f"{temp}", f"{f:.5f}"] for temp, f in zip(self.temp, self.f_cum)
        ]
        data.append(['Average', f"{np.average(self.f_cum):.5f}"])
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        print("Individual correlation factors :")
        print('(Note: use these values only for convergence tests)')
        header = ['label'] + [f"f({int(T)}K)" for T in self.temp]
        data = [
            [label] + f_ind_re_i for label, f_ind_re_i in zip(label_all, f_ind_re)
        ]
        data.append(
            ['Average'] + [f"{f:.5f}" for f in self.f_avg]
        )
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
    
    def print_counts(self):
        name = self.lattice.path_name + self.results[-1].unknown_name
        counts = []
        for start, end in zip(self.index[:-1], self.index[1:]):
            num_normal = len(self.lattice.path)
            counts_i = np.zeros((self.num_vac, len(name)))
            
            for result in self.results[start:end]:
                # normal path
                counts_i[:, :num_normal] += result.counts
                
                # unknown path
                num_unknown = result.counts_unknown.shape[-1]
                if num_unknown > 0:
                    counts_i[:, num_normal:num_normal+num_unknown] += result.counts_unknown       
            
            counts_i = np.sum(counts_i, axis=0)
            counts.append(list(counts_i))          
                    
        print("Counts for each hopping path :")
        header = ['T (K)'] + name
        data = [
            [temp] + counts_i for temp, counts_i in zip(self.temp, counts)
        ]
        counts = np.array(counts)
        data.append(['Total'] + np.sum(counts, axis=0).tolist())
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
    
    def print_total_residence_time(self):
        t_simul = []
        t_reside = []
        for start, end in zip(self.index[:-1], self.index[1:]):
            t_simul_i = 0
            t_reside_i = np.zeros_like(self.lattice.site_name, dtype=np.float64)
            
            for result in self.results[start:end]:
                t_simul_i += result.num_step * result.interval * result.num_vac
                t_reside_i += np.sum(result.residence_time, axis=0)
            t_simul.append(t_simul_i)
            t_reside.append(t_reside_i)
        t_simul = np.array(t_simul)
        t_reside = np.array(t_reside)       
        
        print('Time vacancy remained at each site (ps) :')
        header = ['T (K)'] + self.lattice.site_name + ['Total']
        data = [
            [temp] + t_reside_i.tolist() + [t_simul_i]
            for temp, t_reside_i, t_simul_i in zip(self.temp, t_reside, t_simul)
        ]
        data.append(['Total'] + np.sum(t_reside, axis=0).tolist() + [np.sum(t_simul)])
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
    
    def print_lattice_point(self):
        site_num = []
        site_type = []
        site_coord = []
        for i, site in enumerate(self.lattice.lattice_point):
            site_num.append(i)
            site_type.append(site['site'])
            site_coord.append(site['coord'])
        
        print('Lattice point information : ')
        header = ['name', 'site_type', 'frac_coord']
        data = [
            [f"{self.lattice.symbol}{num}", t, f"[{c[0]:.5f} {c[1]:.5f} {c[2]:.5f}]"]
            for num, t, c in zip(site_num, site_type, site_coord)
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
    
    def print_hopping_history(self):
        print("Computational details :")
        for result in self.results:
            print(f"-----------------(T = {result.temp} K Label = {result.label})-----------------")
            print("Vacancy hopping history :")
            for i in range(self.num_vac):
                print(f"# Vacancy{i+1}")
                header = ['num', 'time (ps)', 'path', 'a (Å)', 'initial site', 'final site']
                data = [
                    [
                        f"{j+1}",
                        f"{path['step'] * result.interval:.2f}",
                        f"{path['name']}",
                        f"{path['distance']:.5f}",
                        f"{path['site_init']} [{', '.join(f'{x:.5f}' for x in self.lattice.lattice_point[path['index_init']]['coord'])}]",
                        f"{path['site_final']} [{', '.join(f'{x:.5f}' for x in self.lattice.lattice_point[path['index_final']]['coord'])}]"
                    ] for j, path in enumerate(result.hopping_history[i])
                ]
                print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
                print('')
            
            print("Summary :")
            
            name_all = self.lattice.path_name + result.unknown_name
            counts_all = np.hstack((result.counts, result.counts_unknown))
            counts_all = np.array(counts_all, dtype=np.int32)
            vacancy_name = [f"Vacancy{i+1}" for i in range(result.num_vac)]
            print("# Path counts :")
            header = ['path'] + vacancy_name
            data = np.vstack((name_all, counts_all)).T
            print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
            print('')
            
            print("# Transient vacancies :")
            step_transient = [k for k, d in result.transient_vacancy.items() if len(d) > 0]
            header = ['step', 'number of transient vacancies', 'site indidces']
            data = [
                [str(k), 
                 len(result.transient_vacancy[k]),
                 " ".join(map(str, result.transient_vacancy[k]))
                ] for k in step_transient
            ]
            print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
            print('') 
            
            print("# Vacancy residence time (ps) :")
            header = ['site'] + vacancy_name
            data = np.vstack((self.lattice.site_name, result.residence_time)).T
            print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
            print('')
    
    def save_figure(self,
                    inset_correlatoin_factor=True):
        # D_rand
        plt.style.use('default')
        plt.rcParams['figure.figsize'] = (3.8, 3.8)
        plt.rcParams['font.size'] = 11
        fig, ax = plt.subplots()
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.2)
        
        for i in range(len(self.temp)):
            plt.scatter(1/self.temp[i], np.log(self.D_rand[i]), 
                        color=self.cmap(i), marker='s', s=50, label=str(int(self.temp[i])))
        slop, intercept = np.polyfit(1/self.temp, np.log(self.D_rand), deg=1)
        x = np.linspace(np.min(1/self.temp), np.max(1/self.temp), 100)
        plt.plot(x, slop*x + intercept, 'k:', linewidth=1)
        plt.xlabel('1/T (1/K)', fontsize=14)
        plt.ylabel(r'ln $D_{rand}$ ($m^{2}$/s)', fontsize=14)
        num_data = len(self.D_rand)
        ncol = int(np.ceil(num_data / 5))
        plt.legend(loc='best', fancybox=True, framealpha=1, edgecolor='inherit',
                   ncol=ncol, labelspacing = 0.3, columnspacing=0.5, borderpad=0.2, handlelength=0.6,
                   fontsize=11, title='T (K)', title_fontsize=11)
        if num_data >= 3:
            x = np.array([self.temp[0], self.temp[int(num_data/2)], self.temp[-1]])
        else:
            x = self.temp
        x_str = [f"1/{int(T)}" for T in x]
        x = 1/x
        plt.xticks(x, x_str)
        plt.savefig('D_rand.png', transparent=False, dpi=300, bbox_inches="tight")
        plt.close()
        
        # tau
        plt.style.use('default')
        plt.rcParams['figure.figsize'] = (3.8, 3.8)
        plt.rcParams['font.size'] = 11
        fig, ax = plt.subplots()
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.2)
            
        for i, temp in enumerate(self.temp):
            ax.bar(temp, self.tau[i], width=50, edgecolor='k', color=self.cmap(i))
            ax.scatter(temp, self.tau[i], marker='o', edgecolors='k', color='k')
        x = np.linspace(0.99*self.temp[0], 1.01*self.temp[-1], 1000)
        ax.plot(x, self.tau0 * np.exp(self.Ea/(self.kb*x)), 'k:')
        plt.xlabel('T (K)', fontsize=14)
        plt.ylabel(r'$\tau$ (ps)', fontsize=14)
        plt.savefig('tau.png', transparent=False, dpi=300, bbox_inches="tight")
        plt.close()
        
        # correlation factor
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
        
        if inset_correlatoin_factor:
            slop, intercept = np.polyfit(1/self.temp, np.log(self.f_cum), deg=1)
            Ea_f = -slop * self.kb
            f0 = np.exp(intercept)
            
            # inset graph
            axins = ax.inset_axes([1.125, 0.615, 0.35, 0.35])
            x_ins = np.linspace(1/self.temp[-1], 1/self.temp[0], 100)
            axins.plot(x_ins, slop * x_ins + intercept, 'k:')
            for i in range(len(self.temp)):
                axins.scatter(1/self.temp[i], np.log(self.f_cum[i]), color=self.cmap(i), marker='s')
            axins.set_xlabel('1/T', fontsize=12)
            axins.set_ylabel(r'ln $f$', fontsize=12)
            axins.set_xticks([])
            axins.set_yticks([])
            
            ax.text(1.125, 0.4, 
                    r"$E_{\mathrm{a}} = $" + f"{Ea_f:.3f} eV", 
                    transform=ax.transAxes, fontsize=10)
            ax.text(1.125, 0.3, 
                    r"$f_{\mathrm{0}} = $" + f"{f0:.3f}", 
                    transform=ax.transAxes, fontsize=10)

        plt.savefig('f_cor.png', transparent=False, dpi=300, bbox_inches="tight")
        plt.close()
        print('')    
        
    def summary(self):
        self.print_lattice_info()
        self.print_simulation_condition()
        self.print_effective_hopping_parameter()
        self.print_diffusion_coefficient()
        self.print_residence_time()
        self.print_correlation_factor()
        self.print_counts()
        self.print_total_residence_time()
        self.print_lattice_point()
        self.print_hopping_history()

        
class PostProcess:
    def __init__(self, 
                 file_params='parameter.txt',
                 file_neb = 'neb.csv',
                 verbose=False):
        # check file
        if os.path.isfile(file_params):
            self.file_params = file_params
        else:
            print(f"{file_params} is not found.")
            sys.exit(0)
        if os.path.isfile(file_neb):
            self.file_neb = file_neb
        else:
            print(f"{file_neb} is not found.")
            sys.exit(0)
        self.verbose = verbose
        self.kb = 8.61733326e-5
        
        # read parameter file
        self.num_sites = None
        self.num_paths = None
        self.path_names = []
        self.z = []
        self.temp = None
        self.times = []
        self.counts = []
        self.D0_eff = None
        self.Ea_eff = None
        self.tau0_eff = None
        self.a_eff = None
        self.f = []
        self.f_mean = None
        self.read_parameter()
        
        # read neb file
        self.Ea = None
        self.read_neb()
        
        # P_site
        self.P_site = self.times / np.sum(self.times, axis=1).reshape(-1,1)
        
        # P_esc
        self.P_esc = np.exp(-self.Ea/(self.kb * self.temp[:, np.newaxis]))
        self.P_esc_eff = np.exp(-self.Ea_eff / (self.kb * self.temp))
        
        # P = P_site * P_esc
        self.P = None
        self.get_P()
        
        # z_mean
        self.z_mean = None
        self.z_mean_rep = None # from total counts from all temperatures
        self.get_z_mean()
        
        # z_eff
        self.z_eff = np.sum(self.P * self.z, axis=1) / self.P_esc_eff
        self.z_eff_rep = np.average(self.z_eff)
        
        # <m>
        self.m_mean = self.z_eff / self.z_mean
        self.m_mean_rep = np.average(self.m_mean)
        
        # nu
        self.nu = None
        self.nu_eff = None
        self.nu_eff_rep = None # simple average of nu_eff
        self.get_nu()
        
        if self.verbose:
            self.summary()
        
    def read_parameter(self):
        with open(self.file_params, 'r') as f:
            lines = [line.strip() for line in f]
            
        for i, line in enumerate(lines):
            if "Lattice information :" in line:
                self.num_sites = int(lines[i+3].split()[-1])
                self.num_paths = list(map(int, lines[i+4].split()[-self.num_sites:]))
                self.num_paths = np.array(self.num_paths)
                
            if "Vacancy hopping paths :" in line:
                for j in range(np.sum(self.num_paths)):
                    contents = lines[i+j+3].split()
                    self.path_names.append(contents[0])
                    self.z.append(int(contents[2]))
                self.z = np.array(self.z, dtype=float)
                    
            if "Simulation temperatures (K) :" in line:
                self.temp = np.array(list(map(float, lines[i].split()[4:])), dtype=float)
                
            if "Time vacancy remained at each site (ps) :" in line:
                for j in range(len(self.temp)):
                    self.times.append(
                        list(map(float, lines[i+j+3].split()[1:1+self.num_sites]))
                    )
                self.times = np.array(self.times)
                
            if "Counts for each hopping path :" in line:
                for j in range(len(self.temp)):
                    self.counts.append(
                        list(map(int, lines[i+j+3].split()[1:1+np.sum(self.num_paths)]))
                    )
                self.counts = np.array(self.counts, dtype=float)
                
            if "Effective hopping parameters :" in line:
                self.D0_eff = float(lines[i+3].split()[2])
                self.tau0_eff = float(lines[i+4].split()[2])
                self.Ea_eff = float(lines[i+5].split()[2])
                self.a_eff = float(lines[i+6].split()[2])
                self.f_mean = float(lines[i+7].split()[1]) # no unit
                
            if "Cumulative correlation factors :" in line:
                self.f =[float(lines[i+j+4].split()[1]) for j in range(len(self.temp))]
    
    def read_neb(self):
        neb = pd.read_csv(self.file_neb, header=None).to_numpy()
        self.Ea = np.zeros(len(self.path_names), dtype=float)
        for name_i, Ea_i in neb:
            index = self.path_names.index(name_i)
            self.Ea[index] = float(Ea_i)
    
    def get_P(self):
        P_site_extend = []
        for p_site in self.P_site:
            P_site_i = []
            for p, m in zip(p_site, self.num_paths):
                P_site_i += [float(p)] * m
            P_site_extend.append(P_site_i)
        self.P = np.array(P_site_extend) * self.P_esc
    
    def get_z_mean(self):
        self.z_mean = np.sum(self.counts, axis=1) / np.sum(self.counts / self.z, axis=1)
        self.z_mean_rep = np.sum(self.counts) / np.sum(np.sum(self.counts, axis=0) / self.z)
        
    def get_nu(self):
        times_extend = []
        for time in self.times:
            times_i = []
            for t, m in zip(time, self.num_paths):
                times_i += [float(t)] * m
            times_extend.append(times_i)
        self.nu = self.counts / (self.z * self.P_esc * times_extend) 
        self.nu_eff = np.sum(self.counts, axis=1) / (np.sum(self.times, axis=1) * np.sum(self.P * self.z, axis=1))
        self.nu_eff_rep = np.average(self.nu_eff)
        
    def summary(self):
        # effective parameters
        print("Effective hopping parameters :")
        header = ['parameter', 'value', 'description']
        parameter = ["Drand_0 (m2/s)", 
                     "tau0 (ps)", 
                     "Ea (eV)", 
                     "a (Å)", 
                     "z", 
                     "nu (THz)", 
                     "f", 
                     "<z>", 
                     "<m>"]
        value = [f"{self.D0_eff:.5e}", 
                 f"{self.tau0_eff:.5e}", 
                 f"{self.Ea_eff:.5f}", 
                 f"{self.a_eff:.5f}", 
                 f"{self.z_eff_rep:.5f}", 
                 f"{self.nu_eff_rep:.5f}", 
                 f"{self.f_mean:.5f}", 
                 f"{self.z_mean_rep:.5f}", 
                 f"{self.m_mean_rep:.5f}"]
        description = ['pre-exponential for random walk diffusivity',
                       'pre-exponential for residence time',
                       'hopping barrier',
                       'hopping distance',
                       'number of equivalent paths (coordination number)',
                       'jump attempt frequency',
                       'correlation factor',
                       'mean number of equivalent paths per path type',
                       'mean number of path types (=z/<z>)']
        data = [[p, v, d] for p, v, d in zip(parameter, value, description)]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        # temperature dependence
        print("Effective hopping parameters with respect to temperature :")
        header = ["T (K)", 
                  "z", 
                  "nu (THz)", 
                  "f",
                  "<z>", 
                  "<m>"]
        data = [
            [int(T), 
             f"{z:.5f}", 
             f"{nu:.5f}",
             f"{f:.5f}",
             f"{z_mean:.5f}", 
             f"{m_mean:.5f}"
             ] 
            for T, z, nu, f, z_mean, m_mean in \
                zip(self.temp, self.z_eff, self.nu_eff, self.f, self.z_mean, self.m_mean) 
        ]
        data.append(['Average', 
                     f"{np.average(self.z_eff):.5f}",
                     f"{np.average(self.nu_eff):.5f}",
                     f"{np.average(self.f):.5f}",
                     f"{np.average(self.z_mean):.5f}", 
                     f"{np.average(self.m_mean):.5f}"])
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        # nu with respect to temperature
        print("Jump attempt frequency (THz) with respect to temperature :")
        print("(Note: only paths with sufficient sampling are reliable)")
        header = ["T (K)"] + self.path_names
        data = [
            [int(T)] + list(nu) for T, nu in zip(self.temp, self.nu)
        ]
        data.append(['Average'] + [f"{nu:.5f}" for nu in np.average(self.nu, axis=0)])
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        # P_site with respect to temperature
        print("P_site with respect to temperature :")
        header = ["T (K)"] + [f"site{i+1}" for i in range(self.num_sites)]
        data = [
            [int(T)] + [f"{p:.5e}" for p in p_i] for T, p_i in zip(self.temp, self.P_site)
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        # P_esc with respect to temperature
        print("P_esc with respect to temperature :")
        header = ["T (K)"] + self.path_names
        data = [
            [int(T)] + [f"{p:.5e}" for p in p_i] for T, p_i in zip(self.temp, self.P_esc)
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')