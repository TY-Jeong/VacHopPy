from vachoppy import trajectory as traj

# make trajectory object
traj_hfo2 = traj.LatticeHopping(poscar_perf='POSCAR_SUPERCELL',
                                xdatcar='XDATCAR_2000K_01')
traj_hfo2.check_unique_vac()

# analyzer
anal_hfo2 = traj.Analyzer(traj=traj_hfo2)

# cn3
num_path = 7
site_init = 'cn3'
site_final = ['cn4', 'cn3', 'cn3', 'cn3', 'cn4', 'cn4', 'cn4']
d = [2.54239, 2.57364, 2.78548, 2.83698, 2.93743, 2.96476, 2.98909]
Ea = [0.74, 0.84, 0.85, 1.35, 1.91, 2.07, 2.01]

for i in range(num_path):
    path_name = f"A{i+1}"
    dE = 0 if site_init==site_final[i] else 0.65
    anal_hfo2.add_path(path_name=path_name,
                       site_init=site_init,
                       site_final=site_final[i],
                       distance=d[i],
                       Ea=Ea[i],
                       dE=dE)

## cn4
num_path = 7
site_init = 'cn4'
site_final = ['cn3', 'cn4', 'cn4', 'cn4', 'cn3', 'cn3', 'cn3']
d = [2.54239, 2.57563, 2.6619, 2.72384, 2.93743, 2.96476, 2.98909]
Ea = [0.08, 0.32, 0.86, 0.98, 1.25, 1.42, 1.36]

for i in range(num_path):
    path_name = f"B{i+1}"
    dE = 0 if site_init==site_final[i] else -0.65
    anal_hfo2.add_path(path_name=path_name,
                       site_init=site_init,
                       site_final=site_final[i],
                       distance=d[i],
                       Ea=Ea[i],
                       dE=dE)

# lattice information
anal_hfo2.site_names
for lat_point in anal_hfo2.lat_points:
    x_coord = lat_point['coord'][0]
    if 0.13796 < x_coord < 0.36204 or 0.63796 < x_coord < 0.86204:
        lat_point['site'] = 'cn4'
    else:
        lat_point['site'] = 'cn3'
        
## search path of vacancy
anal_hfo2.search_path_vac()
print('')
# anal_hfo2.print_path_vac()
# print('')
# print(anal_hfo2.path_vac[1])
# for p_vac in anal_hfo2.path_vac:
#     print(p_vac['step'], end=' ')
# print('')