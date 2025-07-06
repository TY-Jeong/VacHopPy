import os
import re
import sys
import json
import shutil
import numpy as np
from lxml import etree
from collections import Counter
from pymatgen.io.vasp.outputs import Vasprun


def parse_in_file(lammps_in, data_override=None):
    with open(lammps_in) as f:
        lines = f.readlines()
    
    potim = nblock = temp = None
    atom_symbols = []
    lammps_dump = None
    lammps_data = None
    use_restart = False
    
    for line in lines:
        line = line.strip()
        # timestep
        if line.startswith("timestep"):
            potim = float(line.split()[1])
        # dump
        elif line.startswith("dump") and "custom" in line:
            parts = line.split()
            try:
                nblock = int(parts[4])
                lammps_dump = parts[5]
            except:
                pass
        # temp
        elif line.startswith("fix") and "nvt" in line and "temp" in line:
            match = re.search(r"temp\s+([0-9.]+)\s+([0-9.]+)", line)
            if match:
                t1, t2 = float(match.group(1)), float(match.group(2))
                if t1 != t2:
                    raise ValueError(f"TEBEG ({t1}) and TEEND ({t2}) are not equal.")
                temp = t1
        # atom symbols
        elif line.startswith("pair_coeff") and "*" in line:
            atom_symbols = line.split()[4:]
        # data file
        elif line.startswith("read_data"):
            lammps_data = line.split()[1]
        elif line.startswith("read_restart"):
            use_restart = True
    
    if use_restart and data_override:
        lammps_data = data_override
    
    if None in [potim, nblock, temp] or not lammps_dump or not atom_symbols or not lammps_data:
        raise ValueError("Failed to parse required variables from input file.")
    
    atom_symbol_map = {i + 1: sym for i, sym in enumerate(atom_symbols)}
    return potim, nblock, temp, atom_symbol_map, lammps_dump, lammps_data


def extract_from_lammps(species: str,
                        in_file: str,
                        data_file_override=None,
                        prefix_pos: str = "pos",
                        prefix_force: str = "force",
                        prefix_cond: str = "cond"):

    potim_ps, nblock, temperature, type_symbol_map, dump_file, data_file = parse_in_file(
        in_file, data_override=data_file_override)
    potim_fs = potim_ps * 1000.0  # ps → fs

    with open(data_file) as f:
        lines = f.readlines()

    num_atoms = None
    atom_data_raw_start = None
    for i, line in enumerate(lines):
        line_strip = line.strip()
        # parse total number of atoms
        if num_atoms is None and line_strip.lower().endswith("atoms") and any(c.isdigit() for c in line_strip):
            num_atoms = int(line_strip.split()[0])
        # start line for atomic coordinates
        elif line_strip.startswith("Atoms"):
            atom_data_raw_start = i
            break

    if num_atoms is None or atom_data_raw_start is None:
        raise ValueError("Failed to find number of atoms or Atoms section in data file.")

    atom_data = []
    line_idx = atom_data_raw_start + 1
    while len(atom_data) < num_atoms and line_idx < len(lines):
        line_stripped = lines[line_idx].strip()
        if line_stripped != "":
            parts = line_stripped.split()
            atom_id = int(parts[0])
            atom_type = int(parts[1])
            atom_data.append((atom_id, atom_type))
        line_idx += 1

    atom_data.sort(key=lambda x: x[0])
    species_list = [type_symbol_map[atom_type] for _, atom_type in atom_data]

    atom_counts = dict(Counter(species_list))
    target_indices = [i for i, s in enumerate(species_list) if s == species]
    if not target_indices:
        raise ValueError(f"No atoms with symbol '{species}' found.")

    frame_count = 0
    box_lengths = None
    lattice_matrix = None

    with open(dump_file) as f:
        for line in f:
            if "ITEM: BOX BOUNDS" in line:
                bounds = []
                tilt_factors = [0.0, 0.0, 0.0]
                for i in range(3):
                    parts = list(map(float, f.readline().strip().split()))
                    bounds.append(parts[:2])
                    if len(parts) == 3:
                        tilt_factors[i] = parts[2]
                lx = bounds[0][1] - bounds[0][0]
                ly = bounds[1][1] - bounds[1][0]
                lz = bounds[2][1] - bounds[2][0]
                xy, xz, yz = tilt_factors

                if lattice_matrix is None:
                    lattice_matrix = [
                        [lx,  0.0, 0.0],
                        [xy,  ly,  0.0],
                        [xz,  yz,  lz]
                    ]
                    box_lengths = np.array([lx, ly, lz])
            elif "ITEM: ATOMS" in line:
                # skip atom lines
                for _ in range(num_atoms):
                    f.readline()
                frame_count += 1

    n_frames = frame_count
    n_target = len(target_indices)

    pos_memmap = np.lib.format.open_memmap(f"{prefix_pos}.npy",
                                           mode='w+',
                                           dtype=np.float64,
                                           shape=(n_frames, n_target, 3))
    force_memmap = np.lib.format.open_memmap(f"{prefix_force}.npy",
                                             mode='w+',
                                             dtype=np.float64,
                                             shape=(n_frames, n_target, 3))

    with open(dump_file) as f:
        t = 0
        prev_pos_frac = None
        while True:
            line = f.readline()
            if not line:
                break

            if "ITEM: BOX BOUNDS" in line:
                for _ in range(3):
                    f.readline()

            elif "ITEM: ATOMS" in line:
                headers = line.strip().split()[2:]
                col_idx = {name: idx for idx, name in enumerate(headers)}
                atoms = np.zeros((num_atoms, 6))
                for _ in range(num_atoms):
                    parts = f.readline().strip().split()
                    atom_id = int(parts[col_idx['id']])
                    index = atom_id - 1
                    atoms[index, :3] = [float(parts[col_idx[c]]) for c in ('x', 'y', 'z')]
                    atoms[index, 3:] = [float(parts[col_idx[c]]) for c in ('fx', 'fy', 'fz')]

                pos_frac_frame = atoms[target_indices, :3] / box_lengths[np.newaxis, :]

                if t == 0:
                    pos_frac_unwrapped = pos_frac_frame
                else:
                    delta = pos_frac_frame - prev_pos_frac
                    delta[delta > 0.5] -= 1
                    delta[delta < -0.5] += 1
                    pos_frac_unwrapped = pos_memmap[t - 1] + delta

                pos_memmap[t] = pos_frac_unwrapped
                force_memmap[t] = atoms[target_indices, 3:]
                prev_pos_frac = pos_frac_frame
                t += 1

    pos_memmap.flush()
    force_memmap.flush()
    print(f"{prefix_pos}.npy is created.")
    print(f"{prefix_force}.npy is created.")

    cond = {
        "symbol": species,
        "nsw": n_frames,
        "potim": potim_fs,
        "nblock": nblock,
        "temperature": temperature,
        "atom_counts": atom_counts,
        "lattice": lattice_matrix
    }
    with open(f"{prefix_cond}.json", "w") as f:
        json.dump(cond, f, indent=2)
    print(f"{prefix_cond}.json is created.")


# def parse_in_file(in_file, data_override=None):
#     with open(in_file) as f:
#         lines = f.readlines()

#     potim = nblock = temperature = None
#     type_symbols = []
#     dump_file = None
#     data_file = None
#     restart_used = False

#     for line in lines:
#         line = line.strip()
#         if line.startswith("timestep"):
#             potim = float(line.split()[1])
#         elif line.startswith("dump") and "custom" in line:
#             parts = line.split()
#             try:
#                 nblock = int(parts[4])
#                 dump_file = parts[5]
#             except:
#                 pass
#         elif line.startswith("fix") and "nvt" in line and "temp" in line:
#             match = re.search(r"temp\s+([0-9.]+)\s+([0-9.]+)", line)
#             if match:
#                 t1, t2 = float(match.group(1)), float(match.group(2))
#                 if t1 != t2:
#                     raise ValueError(f"TEBEG ({t1}) and TEEND ({t2}) are not equal.")
#                 temperature = t1
#         elif line.startswith("pair_coeff") and "*" in line:
#             parts = line.split()
#             type_symbols = parts[4:]
#         elif line.startswith("read_data"):
#             data_file = line.split()[1]
#         elif line.startswith("read_restart"):
#             restart_used = True

#     if restart_used and data_override:
#         data_file = data_override

#     if None in [potim, nblock, temperature] or not dump_file or not data_file or not type_symbols:
#         raise ValueError("Failed to parse required variables from input file.")

#     type_symbol_map = {i + 1: sym for i, sym in enumerate(type_symbols)}
#     return potim, nblock, temperature, type_symbol_map, dump_file, data_file

# def extract_from_lammps(species: str,
#                         in_file: str,
#                         data_file_override=None,
#                         prefix_pos: str = "pos",
#                         prefix_force: str = "force",
#                         prefix_cond: str = "cond"):

#     potim_ps, nblock, temperature, type_symbol_map, dump_file, data_file = parse_in_file(in_file, data_override=data_file_override)
#     potim_fs = potim_ps * 1000.0  # convert ps to fs

#     with open(data_file) as f:
#         lines = f.readlines()

#     box = np.zeros((3, 2))
#     tilt = [0.0, 0.0, 0.0]
#     num_atoms = 0
#     atom_data = []
#     for i, line in enumerate(lines):
#         if "atoms" in line:
#             num_atoms = int(line.strip().split()[0])
#         elif "xlo xhi" in line:
#             box[0] = list(map(float, line.strip().split()[0:2]))
#         elif "ylo yhi" in line:
#             box[1] = list(map(float, line.strip().split()[0:2]))
#         elif "zlo zhi" in line:
#             box[2] = list(map(float, line.strip().split()[0:2]))
#         elif "xy xz yz" in line:
#             tilt = list(map(float, line.strip().split()[0:3]))
#         elif "Atoms" in line:
#             atom_data = lines[i + 2:i + 2 + num_atoms]
#             break

#     species_list = []
#     for line in atom_data:
#         parts = line.strip().split()
#         atom_type = int(parts[1])
#         species_list.append(atom_type)

#     atom_symbols = [type_symbol_map[t] for t in species_list]
#     atom_counts = dict(Counter(atom_symbols))
#     target_indices = [i for i, s in enumerate(atom_symbols) if s == species]
#     if not target_indices:
#         raise ValueError(f"No atoms with symbol '{species}' found.")

#     with open(dump_file) as f:
#         lines = f.readlines()

#     frames = []
#     i = 0
#     while i < len(lines):
#         if "ITEM: TIMESTEP" in lines[i]:
#             i += 2
#         elif "ITEM: NUMBER OF ATOMS" in lines[i]:
#             i += 2
#         elif "ITEM: BOX BOUNDS" in lines[i]:
#             i += 4
#         elif "ITEM: ATOMS" in lines[i]:
#             headers = lines[i].strip().split()[2:]
#             id_index = headers.index("id")
#             x_index = headers.index("x")
#             y_index = headers.index("y")
#             z_index = headers.index("z")
#             fx_index = headers.index("fx")
#             fy_index = headers.index("fy")
#             fz_index = headers.index("fz")

#             atoms = np.zeros((num_atoms, 6))
#             for j in range(num_atoms):
#                 parts = lines[i + 1 + j].strip().split()
#                 idx = int(parts[id_index]) - 1
#                 atoms[idx, 0] = float(parts[x_index])
#                 atoms[idx, 1] = float(parts[y_index])
#                 atoms[idx, 2] = float(parts[z_index])
#                 atoms[idx, 3] = float(parts[fx_index])
#                 atoms[idx, 4] = float(parts[fy_index])
#                 atoms[idx, 5] = float(parts[fz_index])
#             frames.append(atoms)
#             i += num_atoms + 1
#         else:
#             i += 1

#     n_frames = len(frames)
#     n_target = len(target_indices)

#     pos = np.zeros((n_frames, n_target, 3))
#     force = np.zeros((n_frames, n_target, 3))

#     for t in range(n_frames):
#         atoms = frames[t]
#         for j, idx in enumerate(target_indices):
#             pos[t, j, :] = atoms[idx, 0:3]
#             force[t, j, :] = atoms[idx, 3:6]

#     # Box lengths and triclinic matrix
#     lx = box[0, 1] - box[0, 0]
#     ly = box[1, 1] - box[1, 0]
#     lz = box[2, 1] - box[2, 0]
#     xy, xz, yz = tilt
#     box_lengths = np.array([lx, ly, lz])

#     lattice_matrix = [
#         [lx,  0.0, 0.0],
#         [xy,  ly,  0.0],
#         [xz,  yz,  lz]
#     ]

#     # Unwrap position using fractional coordinates
#     pos_frac = pos / box_lengths[np.newaxis, np.newaxis, :]
#     disp = np.zeros_like(pos_frac)
#     disp[1:] = np.diff(pos_frac, axis=0)
#     disp[disp > 0.5] -= 1
#     disp[disp < -0.5] += 1
#     disp = np.cumsum(disp, axis=0)
#     pos_frac_unwrapped = pos_frac[0] + disp

#     # Save fractional coordinates (unwrapped)
#     np.save(f"{prefix_pos}.npy", pos_frac_unwrapped)
#     np.save(f"{prefix_force}.npy", force)
#     print(f"{prefix_pos}.npy is created.")
#     print(f"{prefix_force}.npy is created.")

#     # nsw: 실제 저장된 frame 개수 (position 수)로 갱신
#     nsw = n_frames

#     cond = {
#         "symbol": species,
#         "nsw": nsw,
#         "potim": potim_fs,
#         "nblock": nblock,
#         "temperature": temperature,
#         "atom_counts": atom_counts,
#         "lattice": lattice_matrix
#     }
#     with open(f"{prefix_cond}.json", "w") as f:
#         json.dump(cond, f, indent=2)
#     print(f"{prefix_cond}.json is created.")


def extract_from_vasp(species: str,
                      vasprun: str = "vasprun.xml",
                      prefix_pos: str = "pos",
                      prefix_force: str = "force",
                      prefix_cond: str = "cond"):
    """
    Extracting VacHopPy input files from vasprun.xml.
    Atomic trajectories are unwrapped according to PBC condition.

    Args:
        species (str): atom species (e.g., "O").
        vasprun (str, optional): vasprun.xml file. Defaults to "vasprun.xml".
        prefix_pos (str, optional): prefix for pos file. Defaults to "pos".
        prefix_force (str, optional): prefix for force file. Defaults to "force".
        prefix_cond (str, optional): prefix for cond file. Defaults to "cond".
    """

    if not os.path.isfile(vasprun):
        print(f"{vasprun} is not found.")
        sys.exit(0)

    v = Vasprun(vasprun, 
                parse_dos=False, 
                parse_eigen=False,
                parse_potcar_file=False)
    structure = v.final_structure
    atom_symbols = [str(site.specie) for site in structure.sites]

    # Atom count dictionary
    from collections import Counter
    atom_counts = dict(Counter(atom_symbols))

    # Target indices for selected species
    target_indices = [i for i, sym in enumerate(atom_symbols) if sym == species]
    if not target_indices:
        raise ValueError(f"No atoms with symbol '{species}' found.")

    iterations = v.ionic_steps
    nsw = len(iterations)
    n_atoms = len(target_indices)

    pos = np.zeros((nsw, n_atoms, 3), dtype=np.float64)
    force = np.zeros((nsw, n_atoms, 3), dtype=np.float64)

    for step_idx, step in enumerate(iterations):
        for j, atom_idx in enumerate(target_indices):
            pos[step_idx, j, :] = step["structure"].sites[atom_idx].frac_coords
            force[step_idx, j, :] = step["forces"][atom_idx]

    # PBC refinement
    displacement = np.zeros_like(pos)
    displacement[0:] = 0
    displacement[1:, :] = np.diff(pos, axis=0)
    displacement[displacement>0.5] -= 1.0
    displacement[displacement<-0.5] += 1.0
    displacement = np.cumsum(displacement, axis=0)
    pos = pos[0] + displacement
    
    # save positions and forces
    np.save(f"{prefix_pos}.npy", pos)
    np.save(f"{prefix_force}.npy", force)
    print(f"{prefix_pos}.npy is created.")
    print(f"{prefix_force}.npy is created.")

    # Extract metadata
    incar = v.incar
    potim = float(incar.get("POTIM", -1))
    nblock = int(incar.get("NBLOCK", -1))
    tebeg = float(incar.get("TEBEG", -1))
    teend = float(incar.get("TEEND", -1))

    if tebeg != teend:
        raise ValueError(f"TEBEG ({tebeg}) and TEEND ({teend}) are not equal.")

    lattice = structure.lattice.matrix.tolist()  # 3x3 list

    cond = {
        "symbol": species,
        "nsw": nsw,
        "potim": potim,
        "nblock": 1, # for vasprun.xml, all steps are stored
        "temperature": tebeg,
        "atom_counts": atom_counts,
        "lattice": lattice
    }
    
    cond_file = f"{prefix_cond}.json"
    with open(cond_file, "w") as f:
        json.dump(cond, f, indent=2)

    print(f"{cond_file} is created.")


def combine_vasprun(vasprun1, 
                    vasprun2,
                    vasprun_out = "vasprun_combined.xml"):
    # Load and merge files
    v1 = etree.parse(vasprun1)
    v2 = etree.parse(vasprun2)

    r1 = v1.getroot()
    r2 = v2.getroot()

    # remove duplicated <parameters> from v2
    p2 = r2.find("parameters")
    if p2 is not None:
        r2.remove(p2)

    # append calculations from v2
    calcs2 = r2.findall("calculation")
    for c in calcs2:
        r1.append(c)

    # recalculate total steps
    nsw_total = len(r1.findall("calculation"))

    # fix ALL <i name="NSW"> across entire tree
    for elem in r1.findall(".//i[@name='NSW']"):
        elem.text = str(nsw_total)

    # save result
    etree.ElementTree(r1).write(vasprun_out, pretty_print=True)
    

def crop_vasprun(vasprun,
                 nsw_crop,
                 vasprun_out="vasprun_cropped.xml"):
    print(f"Cropping first {nsw_crop} iterations from {vasprun}...")
    
    tree = etree.parse(vasprun)
    root = tree.getroot()

    calculations = root.findall(".//calculation")
    total = len(calculations)
    for calc in calculations[nsw_crop:]:
        calc.getparent().remove(calc)

    print(f"  Original steps: {total}, Cropped to: {nsw_crop}")

    # fix ALL <i name="NSW"> across entire tree
    for elem in root.findall(".//i[@name='NSW']"):
        elem.text = str(nsw_crop)

    tree.write(vasprun_out, pretty_print=True, encoding='utf-8', xml_declaration=True)
    print(f"Saved cropped file to {vasprun_out}")
 

def CosineDistance(fp1, fp2):
    """
    fp1 : fingerprint of structure 1
    fp2 : fingerprint of structure 2
    """
    dot = np.dot(fp1, fp2)
    norm1 = np.linalg.norm(fp1, ord=2)
    norm2 = np.linalg.norm(fp2, ord=2)

    return 0.5 * (1 - dot/(norm1*norm2))
