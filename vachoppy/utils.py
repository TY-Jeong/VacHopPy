import os
import sys
import json
import shutil
import numpy as np
from pymatgen.io.vasp.outputs import Vasprun


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

    v = Vasprun(vasprun, parse_dos=False, parse_eigen=False)
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


def concat_xdatcar(xdatcar1, xdatcar2, xdatcar_out):
    if not os.path.isfile(xdatcar1):
        print(f'{xdatcar1} is not found')
        sys.exit(0)

    if not os.path.isfile(xdatcar2):
        print(f'{xdatcar2} is not found')
        sys.exit(0)
        
    nsw1 = find_last_direct_line(xdatcar1)
    nsw2 = find_last_direct_line(xdatcar2)
    
    # copy xdatcar1
    shutil.copyfile(xdatcar1, xdatcar_out)

    # concatnate xdatcars
    with open(xdatcar2, 'r') as f:
        lines = [line.strip() for line in f]

    num_atoms = np.array(lines[6].split(), dtype=int)
    num_atoms = np.sum(num_atoms)

    coords = np.zeros((nsw2, num_atoms, 3))
    for i in range(nsw2):
        idx_i = 8 + (num_atoms+1)*i
        idx_f = idx_i + num_atoms

        coord = lines[idx_i:idx_f]
        coord = [list(map(float, s.split())) for s in coord]
        coord = np.array(coord)
        coords[i] = coord

    with open(xdatcar_out, 'a') as f:
        for idx, coord in enumerate(coords, start=1):
            step = str(nsw1 + idx).rjust(6)
            f.write(f'Direct configuration={step}\n')
            for x, y, z in coord:
                f.write("   %.8f  %.8f  %.8f\n"%(x,y,z))


def concat_force(force1, force2, force_out):
    if not os.path.isfile(force1):
        print(f'{force1} is not found')
        sys.exit(0)

    if not os.path.isfile(force2):
        print(f'{force2} is not found')
        sys.exit(0)

    with open(force1, 'r') as f:
        lines = [s.split()[0].isalpha() for s in f]
        lines = np.array(lines)
    nsw1 = len(np.where(lines==True)[0])

    with open(force2, 'r') as f:
        lines = [s for s in f]

    shutil.copyfile(force1, force_out)
    with open(force_out, 'a') as f:
        for s in lines:
            if s.split()[0].isalpha():
                step = int(s.strip().split()[1])
                step += nsw1
                f.write(f'Iteration {step}\n')
            else:
                f.write(s)


def update_outcar(outcar1, outcar2, outcar_out):
    if not os.path.isfile(outcar1):
        print(f'{outcar1} is not found')
        sys.exit(0)
        
    if not os.path.isfile(outcar2):
        print(f'{outcar2} is not found')
        sys.exit(0)
        
    with open(outcar1, 'r') as f:
        for line in f:
            if 'NSW' in line:
                nsw1 = int(line.split()[2])
                break
            
    with open(outcar2, 'r') as f:
        for line in f:
            if 'NSW' in line:
                nsw2 = int(line.split()[2])
                break
            
    # write outcar_out
    if not os.path.isfile(outcar_out):
        with open(outcar1, 'r') as f1:
            with open(outcar_out, 'w') as f_out:
                for line in f1:
                    if 'NSW' in line:
                        line = f"   NSW    =  {nsw1+nsw2}    number of steps for IOM\n"
                    if 'Iteration' in line:
                        break
                    f_out.write(line)
 

def CosineDistance(fp1, fp2):
    """
    fp1 : fingerprint of structure 1
    fp2 : fingerprint of structure 2
    """
    dot = np.dot(fp1, fp2)
    norm1 = np.linalg.norm(fp1, ord=2)
    norm2 = np.linalg.norm(fp2, ord=2)

    return 0.5 * (1 - dot/(norm1*norm2))

