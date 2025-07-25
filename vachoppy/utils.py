import os
import re
import sys
import json
import shutil
import numpy as np
from lxml import etree
from collections import Counter
from pymatgen.io.vasp.outputs import Vasprun

import re
import numpy as np
import json
import os
from tqdm import tqdm
from colorama import Fore

from vachoppy.trajectory import *

# color map for tqdm
BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m' # Green color
RED = '\033[91m'   # Red color
RESET = '\033[0m'  # Reset to default color


def get_path_single(interval,
                    poscar_lattice="POSCAR_LATTICE",
                    pos_file="pos.npy",
                    force_file="force.npy",
                    cond_file="cond.json",
                    use_incomplete_encounter=False,
                    rmax=3.25):
    
    with open('cond.json', 'r') as f:
        data = json.load(f)
    counts_data = data.get("atom_counts", {})
    
    structure = Structure.from_file(poscar_lattice)
    counts_perf = structure.composition.get_el_amt_dict()
    for sym in counts_perf.keys():
        if abs(counts_data[sym] - counts_perf[sym]) > 1e-3:
            symbol = sym
            num_vac = round(counts_perf[sym] - counts_data[sym])
            break
    
    lattice = Lattice(
        poscar_lattice=poscar_lattice,
        symbol=symbol,
        rmax=rmax,
        verbose=True
    )
    
    traj = Trajectory(
        interval=interval,
        num_vac=num_vac,
        lattice=lattice,
        pos_file=pos_file,
        force_file=force_file,
        cond_file=cond_file,
        verbose=True
    )
    
    anal = TrajectoryAnalyzer(
        lattice=lattice,
        trajectory=traj,
        verbose=True
    )
    
    enc = Encounter(
        analyzer=anal,
        use_incomplete_encounter=use_incomplete_encounter
    )


def get_lattice_from_bounds(bounds: list[str]) -> list[list[float]]:
    tokens = [list(map(float, line.split())) for line in bounds]
    if len(tokens[0]) == 2:  # Orthogonal
        (xlo, xhi), (ylo, yhi), (zlo, zhi) = tokens
        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo
        lattice = [
            [lx, 0.0, 0.0],
            [0.0, ly, 0.0],
            [0.0, 0.0, lz]
        ]
    elif len(tokens[0]) == 3:  # Triclinic
        (xlo, xhi, xy), (ylo, yhi, xz), (zlo, zhi, yz) = tokens
        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo
        lattice = [
            [lx, 0.0, 0.0],
            [xy, ly, 0.0],
            [xz, yz, lz]
        ]
    else:
        raise ValueError("Unknown BOX BOUNDS format.")
    return lattice

def get_nsw_from_dump(dump_file, nblock, num_atoms):
    frame_size = num_atoms + 9
    read_bytes = frame_size * 2 * 100
    with open(dump_file, "rb") as f:
        f.seek(0, os.SEEK_END)
        file_size = f.tell()
        f.seek(max(file_size - read_bytes, 0), os.SEEK_SET)
        tail = f.read().decode(errors="ignore").splitlines()

    last_step = None
    for i in reversed(range(len(tail))):
        if tail[i].startswith("ITEM: TIMESTEP"):
            last_step = int(tail[i + 1].strip())
            break

    with open(dump_file, "r") as f:
        for line in f:
            if line.startswith("ITEM: TIMESTEP"):
                first_step = int(next(f).strip())
                break

    if first_step is None or last_step is None:
        raise ValueError("Could not find timestep info")

    return 1 + (last_step - first_step) // nblock

def extract_from_lammps(species, 
                        log_path, 
                        pos_prefix='pos', 
                        force_prefix='force', 
                        cond_prefix='cond'):
    timestep = None
    dump_file = None
    nblock = None
    temperature = None
    symbol_map = {}

    with open(log_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.strip().startswith("timestep"):
            timestep = float(line.strip().split()[1]) * 1000  # ps -> fs
        elif line.startswith("dump "):
            match = re.search(r'dump\s+\S+\s+all\s+custom\s+(\d+)\s+(\S+)', line)
            if match:
                nblock = int(match.group(1))
                dump_file = match.group(2)
        elif line.startswith("fix "):
            match = re.search(r'temp\s+([0-9.]+)\s+[0-9.]+\s+[0-9.]+', line)
            if match:
                temperature = float(match.group(1))
        elif line.startswith("pair_coeff "):
            parts = line.strip().split()
            i = 1
            while parts[i] == '*':
                i += 1
            symbols = parts[i + 1:]
            for i, s in enumerate(symbols):
                symbol_map[i + 1] = s

    if timestep is None:
        raise RuntimeError("timestep not found in log file")
    if dump_file is None:
        raise RuntimeError("dump file not found in log file")
    if nblock is None:
        raise RuntimeError("nblock not found in log file")
    if temperature is None:
        raise RuntimeError("temperature not found in log file")
    if not symbol_map:
        raise RuntimeError("pair_coeff symbols not found in log file")

    selected_types = [t for t, sym in symbol_map.items() if sym == species]
    if not selected_types:
        raise ValueError(f"species '{species}' not found in symbol_map")

    num_atoms = None
    lattice = None
    check_num_atoms, check_lattice, check_idx = False, False, False
    element_counts = {sym: 0 for sym in symbol_map.values()}
    
    with open(dump_file, 'r') as f:
        while True:
            line = f.readline()
            if (not check_num_atoms) and line.startswith("ITEM: NUMBER OF ATOMS"):
                num_atoms = int(f.readline().strip())
                check_num_atoms = True
            elif (not check_lattice) and line.startswith("ITEM: BOX BOUNDS"):
                bounds_lines = [f.readline().strip() for _ in range(3)]
                lattice = get_lattice_from_bounds(bounds_lines)
                check_lattice = True
            elif (not check_idx) and line.startswith("ITEM: ATOMS"):
                headers = line.split()[2:]
                type_idx = headers.index("type")
                xs_idx = headers.index("xs")
                ys_idx = headers.index("ys")
                zs_idx = headers.index("zs")
                ix_idx = headers.index("ix")
                iy_idx = headers.index("iy")
                iz_idx = headers.index("iz")
                fx_idx = headers.index("fx")
                fy_idx = headers.index("fy")
                fz_idx = headers.index("fz")
                check_idx = True
                
                for _ in range(num_atoms):
                    atom_line = f.readline().split()
                    atom_type = int(atom_line[type_idx])
                    symbol = symbol_map[atom_type]
                    element_counts[symbol] += 1
            
            if check_num_atoms and check_lattice and check_idx:
                break

    nsw = get_nsw_from_dump(dump_file, nblock, num_atoms)
    
    pos_memmap = np.lib.format.open_memmap(f"{pos_prefix}.npy", 
                                           mode='w+', 
                                           dtype=np.float64, 
                                           shape=(nsw, element_counts[species], 3))
    
    force_memmap = np.lib.format.open_memmap(f"{force_prefix}.npy", 
                                             mode='w+', 
                                             dtype=np.float64, 
                                             shape=(nsw, element_counts[species], 3))

    frame_idx = 0
    with open(dump_file, 'r') as f:
        for _ in tqdm(range(nsw),
                      bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'%(Fore.GREEN, Fore.RESET),
                      ascii=False,
                      desc=f"{RED}{BOLD}Reading {dump_file}{RESET}"):
            while True:
                line = f.readline()
                if not line:
                    break
                if line.startswith("ITEM: TIMESTEP"):
                    _ = f.readline()
                elif line.startswith("ITEM: NUMBER OF ATOMS"):
                    _ = f.readline()
                elif line.startswith("ITEM: BOX BOUNDS"):
                    _ = [f.readline() for _ in range(3)]
                elif line.startswith("ITEM: ATOMS"):
                    # _ = f.readline()
                    count = 0
                    for _ in range(num_atoms):
                        atom_line = f.readline().split()
                        atom_type = int(atom_line[type_idx])
                        if atom_type in selected_types:
                            frac_coord = [
                                float(atom_line[xs_idx]) + int(atom_line[ix_idx]),
                                float(atom_line[ys_idx]) + int(atom_line[iy_idx]),
                                float(atom_line[zs_idx]) + int(atom_line[iz_idx])
                            ]
                            force = [
                                float(atom_line[fx_idx]),
                                float(atom_line[fy_idx]),
                                float(atom_line[fz_idx])
                            ]
                            pos_memmap[frame_idx, count] = frac_coord
                            force_memmap[frame_idx, count] = force
                            count += 1
                    frame_idx += 1
                    break

    pos_memmap.flush()
    force_memmap.flush()

    atom_counts = {}
    for i in sorted(symbol_map.keys()):
        sym = symbol_map[i]
        atom_counts[sym] = element_counts.get(sym, 0)

    cond = {
        "symbol": species,
        "nsw": (nsw - 1) * nblock + 1,
        "potim": timestep,
        "nblock": nblock,
        "temperature": temperature,
        "atom_counts": atom_counts,
        "lattice": lattice
    }

    with open(cond_prefix + ".json", "w") as f:
        json.dump(cond, f, indent=2)
        
    print(f"{pos_prefix}.npy is created.")
    print(f"{force_prefix}.npy is created.")
    print(f"{cond_prefix}.json is created.")


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
                    vasprun_out="vasprun_combined.xml"):
    # Use recover=True to handle incomplete XML files
    parser = etree.XMLParser(recover=True)

    # Load and merge files
    v1 = etree.parse(vasprun1, parser)
    v2 = etree.parse(vasprun2, parser)

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
    
    # tree = etree.parse(vasprun)
    parser = etree.XMLParser(recover=True)
    tree = etree.parse(vasprun, parser)
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
