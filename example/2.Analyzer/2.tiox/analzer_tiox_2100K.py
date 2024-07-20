import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from vachoppy import trajectory as traj

temp = 2100

bar_width=0.2,
bar_color='c'

dir_fig = f'fig_{temp}'
dir_txt = f'txt_{temp}'

if not os.path.isdir(dir_fig):
    os.mkdir(dir_fig)
    
if not os.path.isdir(dir_txt):
    os.mkdir(dir_txt)

# NEB results 
num_path = 3
site_init = 'VO'
site_final = ['VO', 'VO', 'VO']
path_name = ['OP', 'IP1', 'IP2']
d = [2.80311, 2.56299, 2.96677]
Ea = [1.07, 0.96, 2.16]
    
def analyzer(num):
    poscar_perf = f'../../xdatcar.tiox.{temp}K/POSCAR_SUPERCELL'
    xdatcar = f"../../xdatcar.tiox.{temp}K/XDATCAR_{num}"

    # make trajectory object
    traj_tio2 = traj.LatticeHopping(poscar_perf=poscar_perf,
                                    xdatcar=xdatcar)
    traj_tio2.check_connectivity()
    traj_tio2.check_unique_vac()

    # analyzer
    anal_tio2 = traj.Analyzer(traj=traj_tio2)

    # NEB results
    for i in range(num_path):
        anal_tio2.add_path(path_name=path_name[i],
                        site_init=site_init,
                        site_final=site_final[i],
                        distance=d[i],
                        Ea=Ea[i],
                        dE=0)
        
    # anal_tio2.print_path()

    # lattice information
    for lat_point in anal_tio2.lat_points:
        lat_point['site'] = 'VO'
            
    # search path of vacancy
    anal_tio2.search_path_vac()
    anal_tio2.unwrap_path()
    print('')

    fig = os.path.join(dir_fig, f"counts_{num}.png")
    txt = os.path.join(dir_txt, f"counts_{num}.txt")

    anal_tio2.print_summary(figure=fig,
                            text=txt,
                            bar_width=bar_width,
                            bar_color=bar_color,
                            disp=False)
    p_name = []
    for p_vac in anal_tio2.path_vac:
        p_name += [p_vac['name']]
    
    return p_name

num = [format(i+1, '02') for i in range(10)]
path_counts = np.zeros(len(path_name))

for n in num:
    p_name = analyzer(num=n)
    for i, name in enumerate(path_name):
        path_counts[i] += p_name.count(name)

x = np.arange(len(path_name))
plt.bar(x, path_counts, color='violet', width=bar_width)
plt.xticks(x, path_name)
plt.xlabel('Path', fontsize=13)
plt.ylabel('Counts', fontsize=13)

fig = fig = os.path.join(dir_fig, f"counts_tot.png")
plt.savefig(fig, dpi=300)
plt.show()

txt = os.path.join(dir_txt, f"counts_tot.txt")
with open(txt, 'w') as f:
    f.write(f"total counts = {np.sum(path_counts)}\n\n")
    f.write("path\tcounts\n")
    for name, count in zip(path_name, path_counts):
        f.write(f"{name}\t{int(count)}\n")
        
        