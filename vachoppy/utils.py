import os
import sys
import shutil
import numpy as np


def concatXDATCAR(xdatcar1, 
                  xdatcar2, 
                  outcar1, 
                  outcar2, 
                  xdatcar_out, 
                  outcar_out):

    if not os.path.isfile(xdatcar1):
        print(f"{xdatcar1} does not exist.")
        sys.exit(0)

    if not os.path.isfile(xdatcar2):
        print(f"{xdatcar2} does not exist.")
        sys.exit(0)

    if not os.path.isfile(outcar1):
        print(f"{outcar1} does not exist.")
        sys.exit(0)

    if not os.path.isfile(outcar2):
        print(f"{outcar2} does not exist.")
        sys.exit(0)

    # read outcars
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
                    f_out.write(line)
                
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


def concatForce(force1, 
                force2, 
                force_out):
    if not os.path.isfile(force1):
        print(f'there is no {force1}.')
        sys.exit(0)

    if not os.path.isfile(force2):
        print(f'there is no {force2}.')
        sys.exit(0)
        
    with open(force1, 'r') as f:
        lines = [s.split()[0].isalpha() for s in f]
        lines = np.array(lines)
    nsw1 = len(np.where(lines==True)[0])

    with open(force2, 'r') as f:
        lines = [s for s in f]
        
    shutil.copyfile(force1, force_new)
    with open(force_new, 'a') as f:
        for s in lines:
            if s.split()[0].isalpha():
                step = int(s.strip().split()[1])
                step += nsw1
                f.write(f'Iteration {step}\n')
            else:
                f.write(s)
                

def extractForce(file_in, 
                 file_out='force.dat'):
    """
    extract force profile from vasprun.xml
    """
    # read vasprun.xml
    with open(file_in, 'r') as f:
        lines = [s.strip() for s in f]
    
    # system info
    nsw, num_atoms = None, None
    for line in lines:
        if 'NSW' in line:
            nsw = int(line.split()[-1].split('<')[0])
        if '<atoms>' in line:
            num_atoms = int(line.split()[1])
        if nsw and num_atoms:
            break
    
    # save forces
    step = 0
    forces = np.zeros((nsw, num_atoms, 3))
    for i, line in enumerate(lines):
        if 'forces' in line:
            force = [list(map(float, s.split()[1:4])) for s in lines[i+1:i+1+num_atoms]]
            force = np.array(force)
            forces[step] = force
            step += 1

    # write out_file
    with open(file_out, 'w') as f:
        for i, force in enumerate(forces, start=1):
            f.write(f"Iteration {i}\n")
            for fx, fy, fz in force:
                fx = str(fx).rjust(12)
                fy = str(fy).rjust(12)
                fz = str(fz).rjust(12)
                f.write(f"{fx} {fy} {fz}\n")