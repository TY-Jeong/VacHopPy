import os
import sys
import time
import copy
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from colorama import Fore
from tabulate import tabulate
from scipy.optimize import minimize_scalar

from vachoppy.inout import *
from vachoppy.trajectory import *

try:
    from mpi4py import MPI
    PARALELL = True
except:
    PARALELL = False

# color map for tqdm
BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m' # Green color
RED = '\033[91m'   # Red color
RESET = '\033[0m'  # Reset to default color



def VacancyHopping_serial(data, 
                          lattice, 
                          interval):
    results = []
    failure = []
    task_size = len(data.datainfo)
    for i in tqdm(range(task_size),
                  bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'%(Fore.GREEN, Fore.RESET),
                  ascii=False,
                  desc=f'{RED}{BOLD}Progress{RESET}'):
        
        try:
            cal = Calculator(
                data=data, index=i, lattice=lattice, interval=interval
            )
        except SystemExit:
            cal = Calculator_fail(data=data, index=i)
        
        if cal.success:
            results.append(cal)
        else:
            failure.append(
                f"  T={cal.temp}K,  Label={cal.label} ({cal.fail_reason})"
            )
    # sort by (temp, label)   
    index = [data.datainfo.index([cal.temp, cal.label]) for cal in results]
    results = [x for _, x in sorted(zip(index, results))]
    # print failed calculations
    if len(failure) > 0:
        print(f"Error reports :")
        for x in failure:
            print(x)
    print('')
    return results


def VacancyHopping_parallel(data, 
                            lattice, 
                            interval):
    time_i = time.time()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    task_size = len(data.datainfo)
    
    if rank==0:
        task_queue = list(range(task_size))
        print(f"Number of AIMD data : {len(task_queue)}")
        results, failure = [], []
        completed_task, terminated_worker, active_workers = 0, 0, size - 1

        while completed_task < task_size or terminated_worker < active_workers:
            status = MPI.Status()
            worker_id, task_result = comm.recv(
                source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status
            )
            
            if status.Get_tag() == 4:
                terminated_worker += 1
                continue
                
            if task_result is not None:
                completed_task += 1
                if task_result.success:
                    results.append(task_result)
                    state = 'success'
                else:
                    failure.append(
                        f"  T={task_result.temp}K,  Label={task_result.label} ({task_result.fail_reason})"
                    )
                    state = 'fail'
                print(f"Progress: {completed_task}/{task_size} finished ({state}), " +
                      f"T={task_result.temp}K Label={task_result.label}, " + 
                      f"remaining workers = {active_workers - terminated_worker}/{active_workers}")
                
            if task_queue:
                new_task = task_queue.pop()
                comm.send(new_task, dest=worker_id, tag=1)
            else:
                comm.send(None, dest=worker_id, tag=0)
                
        while terminated_worker < active_workers:
            worker_id, _ = comm.recv(source=MPI.ANY_SOURCE, tag=4)
            terminated_worker += 1

    else:
        comm.send((rank, None), dest=0, tag=2)
        while True:
            task = comm.recv(source=0, tag=MPI.ANY_TAG)
            if task is None:
                comm.send((rank, None), dest=0, tag=4)
                break
            try:
                cal = Calculator(
                    data=data,
                    index=task,
                    lattice=lattice,
                    interval=interval
                )
            except SystemExit:
                cal = Calculator_fail(data=data, index=task)
            finally:
                comm.send((rank, cal), dest=0, tag=3)
            
    if rank==0:
        index = [data.datainfo.index([cal.temp, cal.label]) for cal in results]
        results = [x for _, x in sorted(zip(index, results))]
        time_f = time.time()
        if failure:
            print(f"\nError reports :")
            for x in failure:
                print(x)
        print('')
        print(f"Total time taken: {time_f - time_i} s")
        return results


class Calculator_fail:
    def __init__(self, data, index):
        self.success = False
        self.fail_reason = 'Unknown reason'
        self.temp = data.datainfo[index][0]
        self.label = data.datainfo[index][1]


class Calculator:
    def __init__(self,
                 data,
                 index,
                 lattice,
                 interval):
        """
        data : vachoppy.inout.DataInfo
        index : index in data.datainfo
        lattice : vachoppy.trajectory.Lattice
        interval : time interval (ps)
        """
        self.data = data
        self.index = index
        self.lattice = lattice
        self.interval = interval
        self.temp, self.label = self.data.datainfo[self.index]
        self.num_path = len(self.lattice.path_names)
        
        # check interval
        potim = self.data.potim[list(self.data.temp).index(self.temp)]
        if (self.interval * 1000) % potim != 0:
            print(f"unvalid interval ({self.temp}K): interval should be a multiple of potim")
            sys.exit(0)
        else:
            self.step_interval = int(self.interval*1000 / potim)
            
        # quantities
        self.path_vac = None
        self.counts = None
        self.unknown = None
        self.t_reside = None
        self.msd_rand = None
        self.encounter_num = None
        self.encounter_msd = None
        self.encounter_path_names = None
        self.encounter_counts = None
        
        # check success
        self.success = True
        self.fail_reason = None
        
        # get quantities
        self.get_quantities()
        
    def get_quantities(self):
        xdatcar = os.path.join(self.data.prefix1,
                               f"{self.data.prefix2}.{self.temp}K",
                               f"XDATCAR_{self.label}")
        force = os.path.join(self.data.prefix1,
                             f"{self.data.prefix2}.{self.temp}K",
                             f"FORCE_{self.label}")
        
        # instantiate VacancyHopping
        try:
            traj = Trajectory(
                xdatcar=xdatcar,
                lattice=self.lattice,
                force=force if self.data.force is not None else None,
                interval=self.step_interval,
                verbose=False
            )
            traj.correct_multivacancy(start=1)
            traj.check_multivacancy()
        except SystemExit:
            self.fail_reason = "Error by trajectory.Trajectory"
            self.success = False
            return
        except BaseException as e:
            self.success = False
            self.fail_reason = "Error by trajectory.Trajectory"
            return
        
        if traj.multi_vac is True:
            self.success = False
            self.fail_reason = "Multi-vacancy issue is not resolved"
            return
        
        if self.data.force is not None:
            try:
                traj.correct_transition_state()
            except SystemExit:
                self.success = False
                self.fail_reason = "Error during TS correction"
                return
            except BaseException as e:
                self.success = False
                self.fail_reason = "Error during TS correction"
                return
        
        # instantiate Analyzer
        try:
            anal = TrajectoryAnalyzer(
                traj=traj,
                lattice=self.lattice,
                verbose=False
            )
        except SystemExit:
            self.success = False
            self.fail_reason = "Error by trajectory.TrajectoryAnalyzer"
            return
        except BaseException as e:
            self.success = False
            self.fail_reason = "Error by trajectory.TrajectoryAnalyzer"
            return
        
            
        self.path_vac = anal.path_vac
        self.counts = anal.counts[:self.num_path]
        self.unknown = anal.path[self.num_path:]
        self.t_reside = anal.total_reside_steps * self.interval
        self.msd_rand = anal.msd_rand
        
        # instantiate Encounter
        try:
            enc = Encounter(
                analyzer=anal,
                verbose=False
            )
        except SystemExit:
            self.success = False
            self.fail_reason = "Error by trajectory.Encounter"
            return
        except BaseException as e:
            self.success = False
            self.fail_reason = "Error by trajectory.Encounter"
            return
        
        self.encounter_num = enc.num_enc
        self.encounter_msd = enc.msd
        self.encounter_path_names = enc.path_names
        self.encounter_path_counts = enc.path_counts
        self.encounter_path_distance = enc.path_dist

class ParameterExtractor:
    def __init__(self,
                 results,
                 data,
                 lattice,
                 tolerance=1e-3,
                 verbose=True,
                 figure=True):
        
        # save arguments
        self.results = copy.deepcopy(results)
        self.data = copy.deepcopy(data)
        self.lattice = copy.deepcopy(lattice)
        self.tolerance = tolerance
        self.verbose = verbose
        self.figure = figure
        
        self.kb = 8.61733326e-5
        self.cmap = plt.get_cmap("Set1")
        self.temp = self.data.temp
        
        # classify unknown paths
        self.unknown_prefix = 'unknown'
        self.unknown_paths = []
        self.merge_unknown_paths()
        self.unknown_paths = sorted(self.unknown_paths, key=lambda x: x['name'])
        
        # labels for successful data
        self.num_label = [0] * len(self.temp)
        for cal in self.results:
            self.num_label[list(self.temp).index(cal.temp)] += 1
        self.label_all = sorted(list(set([cal.label for cal in self.results])))
        
        # correlation factor
        self.f_ind = []
        self.f_avg = []
        self.f_cum = []
        self.calculate_correlation_factor()
        
        # diffusivity
        self.D_rand = []
        self.Ea = None
        self.D0_rand = None
        self.calculate_diffusivity()
        
        # residence time
        self.tau = []
        self.tau0 = None
        self.calculate_residence_time()
        
        # hopping distance
        self.a_eff = np.sqrt(6*self.D0_rand*self.tau0) * 1e4
        
        # <z>
        self.z_mean = []
        self.calculate_z_mean()
        
        # print results
        if self.verbose:
            self.print_lattice_info()
            self.print_simulation_condition()
            self.print_effective_hopping_parameters()
            self.print_diffusivity()
            self.print_residence_time()
            self.print_correlation_factor()
            self.print_counts()
            self.print_total_resideence_time()
            self.print_lattice_point()
            self.print_hopping_history()
        
        # save figures
        if self.figure:
            self.save_figure()
        
    def merge_unknown_paths(self):
        unknown_paths = sum([cal.unknown for cal in self.results],[])
        unknown_distances = np.array([path['distance'] for path in unknown_paths])
        unknown_distances.sort()
                
        unknown_distance_unique = []
        for distance in unknown_distances:
            if not unknown_distance_unique or \
                np.abs(unknown_distance_unique[-1] - distance) > self.tolerance:
                unknown_distance_unique.append(distance)
        unknown_distance_unique = np.array(unknown_distance_unique)
        
        check = [True] * len(unknown_distance_unique)
        for path in unknown_paths:
            index = np.argmin(np.abs(unknown_distance_unique - path['distance']))
            path['name'] = self.unknown_prefix + str(index+1)
            if check[index]:
                path_dic = {}
                path_dic['name'] = path['name']
                path_dic['distance'] = path['distance']
                path_dic['site_init'] = path['site_init']
                path_dic['site_final'] = path['site_final']
                path_dic['coord_init'] = path['coord_init']
                path_dic['coord_final'] = path['coord_final']
                self.unknown_paths.append(path_dic)
                check[index] = False
    
    def calculate_diffusivity(self):
        distance = np.array([path['distance'] for path in self.lattice.path])
        counts = np.array([cal.counts for cal in self.results])
        t = self.data.nsw * self.data.potim * self.num_label
        
        index = [0] + list(np.cumsum(self.num_label))
        self.D_rand = np.array(
            [np.sum(distance**2 * counts[index[i]:index[i+1]]) / (6*t[i]) for i in range(len(self.temp))]
        ) * 1e-5
        
        # Arrhenius fit
        slop, intercept = np.polyfit(1/self.temp, np.log(self.D_rand), deg=1)
        self.Ea = -slop * self.kb
        self.D0_rand = np.exp(intercept)
    
    def calculate_residence_time(self):
        index = [0] + list(np.cumsum(self.num_label))
        counts = np.array([cal.counts for cal in self.results])
        t = self.data.nsw * self.data.potim * self.num_label
        
        self.tau = np.array(
            [t[i] / np.sum(counts[index[i]:index[i+1]]) for i in range(len(self.temp))]
        ) * 1e-3
        
        # fitting
        error_tau = lambda tau0: np.linalg.norm(
            self.tau - tau0 * np.exp(self.Ea / (self.kb * self.temp))
        )
        result = minimize_scalar(error_tau)
        self.tau0 = result.x # ps
        
    def calculate_z_mean(self):
        index = [0] + list(np.cumsum(self.num_label))
        counts = np.array([cal.counts for cal in self.results])
        z = np.array(
            [path['z'] for path in self.lattice.path], dtype=float
        )
        self.z_mean = np.array(
            [np.sum(counts[index[i]:index[i+1]]) for i in range(len(self.temp))]
        )
        self.z_mean /= np.array(
            [np.sum(counts[index[i]:index[i+1]]/z) for i in range(len(self.temp))]
        )
    
    def calculate_correlation_factor(self):
        distance = np.array(
            [path['distance'] for path in self.lattice.path] +
            [path['distance'] for path in self.unknown_paths], dtype=float
        )
        
        # refine encounter_path_counts
        for cal in self.results:
            counts_modified = np.zeros_like(distance, dtype=float)
            for i, d in enumerate(cal.encounter_path_distance):
                index = np.where(np.abs(distance -d) < self.tolerance)[0]
                if index.size > 0:
                    counts_modified[index[0]] = cal.encounter_path_counts[i]
            cal.encounter_path_counts = counts_modified
            
        # correlation factor for individual data    
        self.f_ind = np.array([
            cal.encounter_msd / np.sum(distance**2 * cal.encounter_path_counts)
            for cal in self.results
        ])
        
        # average at each temperature
        index = [0] + list(np.cumsum(self.num_label))
        self.f_avg = np.array([
            np.average(self.f_ind[index[i]:index[i+1]]) for i in range(len(self.temp))
        ])
        
        # cumulative correlation factor
        for i in range(len(self.temp)):
            num_enc = np.array(
                [cal.encounter_num for cal in self.results[index[i]:index[i+1]]]
            )
            msd = np.array(
                [cal.encounter_msd for cal in self.results[index[i]:index[i+1]]]
            )
            msd_cum = np.sum(msd * num_enc) / np.sum(num_enc)
            counts = np.array(
                [cal.encounter_path_counts for cal in self.results[index[i]:index[i+1]]]
            )
            counts_cum = np.sum(counts * num_enc.reshape(-1,1), axis=0) / np.sum(num_enc)
            self.f_cum.append(msd_cum / np.sum(distance**2 * counts_cum))
        self.f_cum = np.array(self.f_cum)

    def print_correlation_factor(self):
        # rearrange f_ind
        f_ind_re = []
        index = [0] + list(np.cumsum(self.num_label))
        for i, f in enumerate(range(len(self.temp))):
            f_ind_re_i = ['-'] * len(self.label_all)
            for j, cal in enumerate(self.results[index[i]:index[i+1]]):
                f_ind_re_i[self.label_all.index(cal.label)] = f"{self.f_ind[index[i]+j]:.5f}"
            f_ind_re.append(f_ind_re_i)
        f_ind_re = [list(x) for x in zip(*f_ind_re)] # transpose
        
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
            [label] + f_ind_re_i for label, f_ind_re_i in zip(self.label_all, f_ind_re)
        ]
        data.append(
            ['Average'] + [f"{f:.5f}" for f in self.f_avg]
        )
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
    def print_diffusivity(self):
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
        
    def print_counts(self):
        name = self.lattice.path_names + [path['name'] for path in self.unknown_paths]
        index = [0] + list(np.cumsum(self.num_label))
        counts = []
        for i in range(len(self.temp)):
            counts_i = np.zeros_like(name, dtype=float)
            for cal in self.results[index[i]:index[i+1]]:
                counts_i[:len(cal.counts)] += cal.counts
                if len(cal.unknown) > 0:
                    for path in cal.unknown:
                        counts_i[name.index(path['name'])] += 1
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
        
    def print_total_resideence_time(self):
        index = [0] + list(np.cumsum(self.num_label))
        t_reside = []
        for i in range(len(self.temp)):
            t_reside_i = np.zeros_like(self.lattice.site_names, dtype=float)
            for cal in self.results[index[i]:index[i+1]]:
                t_reside_i += cal.t_reside 
            t_reside.append(t_reside_i)
        t_reside = np.array(t_reside)
        
        print('Time vacancy remained at each site (ps) :')
        header = ['T (K)'] + self.lattice.site_names
        data = [
            [temp] + t_reside_i.tolist()
            for temp, t_reside_i in zip(self.temp, t_reside)
        ]
        data.append(['Total'] + np.sum(t_reside, axis=0).tolist())
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
    def print_effective_hopping_parameters(self):
        print('Effective hopping parameters : ')
        header = ['parameter', 'value', 'description']
        parameter = [
            'Drand_0 (m2/s)', 'tau0 (ps)', 'Ea (eV)', 'a (Å)', 'f', '<z>'
        ]
        value = [
            f"{self.D0_rand:.5e}", f"{self.tau0:.5f}", f"{self.Ea:.5f}", f"{self.a_eff:.5f}",
            f"{np.average(self.f_cum):.5f}", f"{np.average(self.z_mean):.5f}"
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
        
    def print_lattice_info(self):
        num_path_site = np.zeros(len(self.lattice.site_names))
        for path in self.lattice.path:
            num_path_site[self.lattice.site_names.index(path['site_init'])] += 1
        
        print('Lattice information :')
        print('  Number of sites =', len(self.lattice.site_names))
        print('  Number of hopping paths = ', end='')
        for num in num_path_site:
            print(int(num), end=' ')
        print('')
        print('  Number of unknown paths =', len(self.unknown_paths))
        print('')
        print('Vacancy hopping paths : ')
        header = [
            'path', 'a(Å)', 'z', 'initial site', 'final site'
        ]
        data = [
            [path['name'],  f"{path['distance']:.3f}", path['z'], 
             path['site_init']+f" [{path['coord_init'][0]:.5f} {path['coord_init'][1]:.5f} {path['coord_init'][2]:.5f}]",
             path['site_final']+f" [{path['coord_final'][0]:.5f} {path['coord_final'][1]:.5f} {path['coord_final'][2]:.5f}]"]
            for path in self.lattice.path
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        print('Non-vacancy hopping paths : ')
        header = [
            'path', 'a(Å)', 'z', 'initial site', 'final site'
        ]
        data = [
            [path['name'],  f"{path['distance']:.3f}", '-', 
             path['site_init']+f" [{path['coord_init'][0]:.5f} {path['coord_init'][1]:.5f} {path['coord_init'][2]:.5f}]",
             path['site_final']+f" [{path['coord_final'][0]:.5f} {path['coord_final'][1]:.5f} {path['coord_final'][2]:.5f}]"]
            for path in self.unknown_paths
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
    def print_simulation_condition(self):
        print(f'Time ineterval used for average (t_interval) = {self.results[0].interval} ps')
        print('')
        print('Simulation temperatures (K) : ', end='')
        for temp in self.temp:
            print(temp, end=' ')
        print('\n')
        print('Labels of AIMD data :')
        label_check = [
            ['O' if label in data else 'X' for data in self.data.label] for label in self.label_all
        ]
        num_label = [len(label) for label in self.data.label]
        header = ['label'] + [f"{int(T)}K" for T in self.temp]
        data = [
            [label] + check for label, check in zip(self.label_all, label_check)
        ]
        data.append(['Num'] + num_label)
        print(tabulate(data, headers=header, tablefmt="simple", stralign='center', numalign='center'))
        print('')
        
    def print_hopping_history(self):
        unknown_name = [path['name'] for path in self.unknown_paths]
        unknown_distance = [f"{path['distance']:.3f}" for path in self.unknown_paths]
        unknown_init = [path['site_init'] for path in self.unknown_paths]
        unknown_final = [path['site_final'] for path in self.unknown_paths]
        print("Computational details :")
        for cal in self.results:
            print(f"-----------------(T = {cal.temp} K Label = {cal.label})-----------------")
            print("Vacancy hopping history :")
            header = ['num', 'path', 'a (Å)', 'step', 'init', 'final']
            data=[
                [str(i+1), path['name'], f"{path['distance']:.3f}", f"{path['step']}",
                 f"{path['index_init']} ({path['site_init']})", f"{path['index_final']} ({path['site_final']})"]
                for i, path in enumerate(cal.path_vac)
            ]
            print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
            print('')
            
            print("Path counts :")
            header = ['path', 'count', 'a (Å)', 'init', 'final', 'z']
            data =[
                [name, count, f"{self.lattice.path[i]['distance']:.3f}", 
                 self.lattice.path[i]['site_init'], self.lattice.path[i]['site_final'], self.lattice.path[i]['z']]
                for i, [name, count] in enumerate(zip(self.lattice.path_names, cal.counts))
            ]
            
            unknown_count = [
                [path['name'] for path in cal.unknown].count(name) for name in unknown_name
            ]
            unknown_data = [
                [unknown_name[i], unknown_count[i], unknown_distance[i], unknown_init[i], unknown_final[i], '-']
                for i in range(len(unknown_name))
            ]
            data += unknown_data
            print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
            print('')
            
            print("Total steps a vacancy remained at each site :")
            header = ['site', 'total steps']
            data = [
                [name, f"{int(step/cal.interval)} ({step:.2f} ps)"] 
                for name, step in zip(self.lattice.site_names, cal.t_reside)
            ]
            print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
            print('')
            
            print("Encounter data for correlation factor calculation :")
            print("  MSD for a random walk (Å2) :", f"{cal.msd_rand:.3f}")
            print("  Number of encounters :", cal.encounter_num)
            print("  MSD of encounters (Å2):", f"{cal.encounter_msd:.3f}")
            print("  mean hopping counts per encounter :", f"{np.sum(cal.encounter_path_counts):.3f}")
            print('')
            
    def print_lattice_point(self):
        site_num = []
        site_type = []
        site_coord = []
        for i, site in enumerate(self.lattice.lat_points):
            site_num.append(i+1)
            site_type.append(site['site'])
            site_coord.append(site['coord'])
        
        print('Lattice point information : ')
        print('(Note: name is the same as in VESTA)')
        header = ['name', 'site_type', 'frac_coord']
        data = [
            [f"{self.lattice.symbol}{num}", t, f"[{c[0]:.5f} {c[1]:.5f} {c[2]:.5f}]"]
            for num, t, c in zip(site_num, site_type, site_coord)
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
    def save_figure(self):
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

        # inset graph
        # axins = ax.inset_axes([1.125, 0.615, 0.35, 0.35])
        # x_ins = np.linspace(1/self.temp[-1], 1/self.temp[0], 100)
        # axins.plot(x_ins, -(self.Ea_f/self.kb) * x_ins + np.log(self.f0), 'k:')
        # for i in range(len(self.temp)):
        #     axins.scatter(1/self.temp[i], np.log(self.f_cum[i]), color=self.cmap(i), marker='s')
        # axins.set_xlabel('1/T', fontsize=12)
        # axins.set_ylabel(r'ln $f$', fontsize=12)
        # axins.set_xticks([])
        # axins.set_yticks([])
        plt.savefig('f_cor.png', transparent=False, dpi=300, bbox_inches="tight")
        plt.close()
        print('')


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
                self.num_sites = int(lines[i+1].split()[-1])
                self.num_paths = list(map(int, lines[i+2].split()[-self.num_sites:]))
                self.num_paths = np.array(self.num_paths)
                
            if "Vacancy hopping paths :" in line:
                for j in range(np.sum(self.num_paths)):
                    contents = lines[i+j+3].split()
                    self.path_names.append(contents[0])
                    self.z.append(int(contents[2]))
                self.z = np.array(self.z, dtype=float)
                    
            if "Simulation temperatures (K) :" in line:
                self.temp = np.array(list(map(int, lines[i].split()[4:])), dtype=float)
                
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