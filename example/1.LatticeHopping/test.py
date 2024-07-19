from vachoppy import trajectory as tj

hfo2 = tj.LatticeHopping(poscar_perf='POSCAR_SUPERCELL',
                         xdatcar='XDATCAR_2000K_01',
                         interval=50,
                         target='O')

# print(hfo2.lat_points.shape)
# print(hfo2.traj_vac_C[0])
# print(hfo2.occ_lat_point.shape)
# print(hfo2.num_atoms[hfo2.idx_target])
# hfo2.animation()
hfo2.save_poscar(step=100, vac=True)
hfo2.save_poscar(step=200, 
                 vac=True,
                 on_lat=True)
