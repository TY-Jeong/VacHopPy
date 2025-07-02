# VacHopPy 

---
**VacHopPy** is a Python package for analyzing vacancy hopping mechanisms based on *Ab initio* molecular dynamics (AIMD) simulations. A detailed explanation on **VacHopPy** framwork is available [**here**](https://arxiv.org/abs/2503.23467).

<div align="center">
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/main/imgs/logo.png" width="550"/>
</p>
</div>


## Features

* Tracking of **vacancy trajectories** in AIMD simulations
* Extraction of **effective hopping parameter** set
* Assessment of lattice stability or **phase transitions**

<br>

**Effective hopping parameter** set, a key improvement of **VacHopPy**, is a single, representative set of hopping parameters, which is determined by integrating all possible hopping paths in a given system considering energetic and geometric properties. Hence, the effective hopping parameters are suitable for multiscaling modeling, bridging the *ab initio* calculations and device-scale simulations (e.g., continuum models).

The list of effective hopping parameters, which can be obtained using **VacHopPy** is summarized below:




<div align="center">

| Symbol | Description                     |
|--------|---------------------------------|
| D      | Diffusion coefficient (m²/s)    |
| f      | Correlation factor              |
| τ      | Residence time (ps)             |
| a      | Hopping distance (Å)            |
| Eₐ     | Hopping barrier (eV)            |
| z      | Number of equivalent paths      |
| ν      | Attempt frequency (THz)         |
</div>


Within the **VacHopPy** framework, the temperature dependencies of overall diffusion behavior, represented by *D* and *τ*, are simplified to an Arrehnius equation composed of the effective hopping parameters. Please see the original paper for a detailed description.



## Contents

* Installation
* Requirement
* List of commands
* How to implement
  1. Preparation
  2. Vacancy trajectory visualization
  3. Extraction of effective hopping parameters
  4. Assessment of lattice stability
* Reference
* License



## Installation

This package can be easily installed via pip.

```bash
pip install vachoppy
```

If you need parallelized calculations, download the `mpi4py` package:

```bash
pip install mpi4py
```

# Requirements

**VacHopPy** requires **Python 3.10 or higher** and the following Python packages:
* numpy
* lxml
* tqdm
* colorama
* matplotlib >= 3.10.0
* scipy
* tabulate
* pymatgen >= 2024.6.10

Current version of **VacHopPy** was developed and tested using **VASP 5.4.4** and **LAMMPS (12 Jun 2025)**.

## Available commands

**VacHopPy** provides a command-line interface (CLI). Belows are available CLI commands:

<div align=center>
<table>
    <tr>
        <th scope="col">Option 1</td>
        <th scope="col">Option 2</td>
        <th scope="col">Use</td>
    </tr>
    <tr>
        <td rowspan="4">-m<br>(main)</td>
        <td>t</td>
        <td>Make an animation for vacancy trajectories</td>
    </tr>
    <tr>
        <!-- <td>2</td> -->
        <td>p</td>
        <td>Calculate effective hopping parameters (excluding z and ν)</td>
    </tr>
    <tr>
        <!-- <td>4</td> -->
        <td>pp</td>
        <td>Calculate z and ν (post-processing for -m p mode)</td>
    </tr>
    <tr>
        <!-- <td>5</td> -->
        <td>f</td>
        <td>Perform fingerprint analyses</td>
    <tr>
        <td rowspan="6">-u<br>(utility)</td>
        <td>extract_data</td>
        <td>Extract AIMD data (cond.json, pos.npy, force.npy)</td>
    </tr>
    <tr>
        <!-- <td>2</td> -->
        <td>combine_vasprun</td>
        <td>Combine two successive vasprun.xml</td>
    </tr>
    <tr>
        <!-- <td>2</td> -->
        <td>crop_vasprun</td>
        <td>Crop vasprun.xml</td>
    </tr>
    <tr>
        <!-- <td>3</td> -->
        <td>fingerprint</td>
        <td>Extract fingerprint</td>
    </tr>
    <tr>
        <!-- <td>4</td> -->
        <td>cosine_distance</td>
        <td>Calculate cosine distance</td>
    </tr>
</table>
</div>

For detailed descriptions, please use `-h` flag:

```bash
vachoppy -h # list of available commands
vachoppy -m p -h # explanation for '-m p' mode
```

For time-consuming commands, `vachoppy -m p` and `vachoppy -m f`, parallelization is supported by **mpirun**. For parallelization, please specify `--parallel` flag:

```bash
vachoppy -m p 0.1 # serial
mpirun -np {num_nodes} vachoppy -m p 0.1 --parallel # parallel
```


Below is a summary of the main commands (only main modules are shown for clarity):


<p align="center">
  <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/main/imgs/Flowchart.JPG" width="800"/>
</p>


## How to implement

Example files can be downloaded from:

* **Example1** : AIMD data extraction from VASP / LAMMPS outputs. [download (100 MB)](https://drive.google.com/file/d/1HGPto-X8RyVGjUuvjAmdd-5Yb4CHxyE6/view?usp=sharing)
* **Example2** : Vacancy hopping in rutile TiO<SUB>2</SUB> [download (30 GB)](https://drive.google.com/file/d/1F92sJWvM_5gIvE76sLAbCjxv0U0L9ngG/view?usp=sharing)
* **Example3** : Phase transition of monoclinic HfO<SUB>2</SUB> at 2200 K  [download (300 MB)](https://drive.google.com/file/d/1F53tN4NGgw5jcU24qBHQf6XGHNI_i5TR/view?usp=sharing)

## 1. Preparation

>Download **Example1** directory linked above.

### Input data

To run **VacHopPy**, the user needs three types of input data: **AIMD data**, **POSCAR_LATTICE**, and **neb.csv** (*optional*). The current version of **VacHopPy** supports AIMD simulations performed under the NVT ensemble only.

#### (1) AIMD data

AIMD data can be extracted from the **VASP** or **LAMMPS** output files using the `vachoppy -u extract_data` command.

- **Extracting AIMD data from VASP**

Navigate to the `Example1/vasp` directory and run:

```bash
vachoppy -u extract_data O vasprun.xml 
```
Here, the arguments are:
* atom symbol = O
* MD result = vasprun.xml

The `vasprun.xml` file contains the result of an AIMD simulation of a rutile TiO₂ system with two oxygen vacancies at 2100 K.
The `atom symbol` specifies the type of atom to track, as well as the type of vacancy.
When the atom symbol is set to `O`, this command extracts the positions and forces of oxygen atoms from the `vasprun.xml` file.

This command generates the following three output files:

* **cond.json** <br> cond.json contains simulation conditions, such as temperature, time, atom numbers, and lattice parameters.

* **pos.npy** (numpy binary) <br> pos.npy contains evolution in atomic positions during the simulation. The raw trajectories are refined with consideration of periodic boundary condition (PBC).

* **force.npy** (numpy binary) <br> force.npy contains force vectors acting on atoms.

To reduce the size of the input data, only atoms of the specified type (as defined by the atom symbol) are included in the output files.


- **Extract AIMD data from LAMMPS**

Starting from version 2.0.0, **VacHopPy** supports integration with **LAMMPS**.
To extract AIMD data from LAMMPS output files, use the `-l` flag:

Navigate to the `Example1/lammps` directory and run:

```bash
vachoppy -u extract_data O lammps.in -l
```
Here, the arguments are:
* atom symbol = O
* MD result = lammps.in

This command reads both the **LAMMPS data file** (e.g., coo.lammps) and the **LAMMPS dump file** (e.g., lammps.dump) specified in the `lammps.in`, and produces the same three output files: cond.json, pos.npy, and force.npy.

Using LAMMPS with **machine learning potentials (MLPs)** can greatly reduce the computational cost required to sample sufficient hopping events in MD simulations. Although LAMMPS is not an AIMD engine, it can still provide reliable MD trajectories when used with a well-trained MLP.


#### (2) POSCAR_LATTICE
This file contains the perfect crystal structure without vacnacies. Its lattice parameters must match those of input structure (e.g., POSCAR) of the MD simulations. This file is used to define the lattice points for vacancy identification.


#### (3) neb.csv (*optional*)
This file contains **hopping barriers (Eₐ)** for all vacancy hopping paths in the system. Below is an example of `neb.csv` (example system: rutile TiO<SUB>2</SUB>):

```ruby
# neb.csv
A1,0.8698
A2,1.058
A3,1.766
```
Here, the **first column** corresponds to the **path names**, and the **second column** contains the **Eₐ values**. The user can obtain a list of possible vacancy hopping paths by running the `vachoppy -m t` or `vachoppy -m p` command. For the `vachoppy -m t` command (when used with the `-v` flag), hopping path information is saved to the **trajectory.txt** file. For the `vachoppy -m p` command, the information is saved to the **parameter.txt file**.

>**Note1**: the **neb.csv** file is only required to extract the effective values for **coordination number (z)** and **attempt frequency (ν)** (by running the `vachoppy -m pp` command).

>**Note2**: It is highly recommended to **perform NEB calculations using a larger supercell** than that used in MD simulations. In MD simulations, thermal fluctuations attenuate interactions with periodic images and provide a broader sampling of atomic configurations, which helps approximate the effects of a larger supercell.


#### (4) File organization
Since AIMD simulations typically cover timescales shorter than a nanosecond, a single AIMD simulation may contain a few hopping events. However, since **VacHopPy** computes the effective hopping parameters in static manner, sufficient sampling of hopping events is necessary to ensure reliablilty. To address this, **VacHopPy** processes multiple AIMD datasets simultaneously. Each AIMD dataset is distinguished by a label appended after an underscore in the file names (e.g., cond_{label}.json, pos_{label}.npy, force_{label}.npy). Below is an example of the recommended file structure:

```bash
{where VacHopPy is executed}
 ┣ traj # name (traj) is specified by -p1 flag
 ┃ ┣ traj.1900K # prefix (traj) is specifed by -p2 flag
 ┃ ┃ ┣ cond_01.json, pos_01.npy, force_03.npy  
 ┃ ┃ ┣ cond_02.json, pos_02.npy, force_03.npy  
 ┃ ┃ ┗ cond_02.json, pos_02.npy, force_04.npy  
 ┃ ┣ traj.2000K
 ┃ ┃ ┣ cond_01.json, pos_01.npy, force_03.npy  
 ┃ ┃ ┣ cond_02.json, pos_02.npy, force_03.npy  
 ┃ ┃ ┗ cond_02.json, pos_02.npy, force_04.npy  
 ┃ ┗ traj.2100K
 ┃   ┣ cond_01.json, pos_01.npy, force_03.npy  
 ┃   ┣ cond_02.json, pos_02.npy, force_03.npy  
 ┃   ┗ cond_02.json, pos_02.npy, force_04.npy  
 ┃
 ┣ POSCAR_LATTICE
 ┗ neb.csv
```

The name of the outer directory is specified by the `-p1` flag (default: traj), and the prefix of the inner directories is specified by the `-p2` flag (default: traj). Each inner directory must contain AIMD datasets generated at the same temperature.


### Hyperparameter: t<SUB>interval</SUB>
To run **VacHopPy**, the user needs to determine one hyperparameter, **t<SUB>interval</SUB>**, in advance. This parameter defines the time interval for averaging atomic positions and forces. Thermal fluctuations in AIMD simulations can make it difficult to precisely determine atomic occupancies. However, since these fluctuations are random, they can be effectively averaged out over time. **VacHopPy** processes AIMD data by dividing it into segments of length of t<SUB>interval</SUB>. Each segment corresponds to one analysis step, representing the averaged structure over that time interval. The total number of steps is given by t<SUB>simulation</SUB>/t<SUB>interval</SUB>, where t<SUB>simulation</SUB> is the total AIMD simulation time.



<div align=center>
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/main/imgs/t_interval.jpg" width="800"/>
</p>
</div>

Choosing an appropriate t<SUB>interval</SUB> is crucial for reliable analysis. The t<SUB>interval</SUB> should be large enough to mitigate thermal fluctuations but short enough to prevent multiple hopping events from being included in a single step. A typical value is around 0.05-0.1 ps, through it may vary depending on the system. 

One recommended approach for determining the optimal t<SUB>interval</SUB> is through convergence tests using the correlation factor (f). Below is an example of a convergence test (example system: rutile TiO<SUB>2</SUB>):

<div align=center>
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/refs/heads/main/imgs/converge.png" width="550"/>
</p>
</div>

The left and right figures show the convergences of f with respect to the number of AIMD datasets (N<SUB>cell</SUB>) and t<SUB>interval</SUB>, respectively, at each temperature. The results confirm that convergence is achieved at **N<SUB>cell</SUB>=20** and **t<SUB>interval</SUB>=0.07 ps**. 


## 2. Vacancy trajectory visualization

>Download **Example2** directory linked above.

Navigate to the `Example2` directory and run:
```bash
 vachoppy -m t 0.07 2100 07 -v 
 ```

Here, the arguments are:

* t<SUB>interval</SUB> = 0.07 ps
* temperature = 2100 K
* label = 07

The `-v` flag (verbosity flag) enables the generation of the **trajectory.txt** file, which contains the identified vacancy hopping paths and hopping history. You can adjust the resolution of the animation using the `--dpi` flag (default: 300). The number and type of vacancies are automatically determined by comparing the POSCAR_LATTICE file with the AIMD datasets.


**Output:**

The trajectory animation is saved as **traj.gif**, with individual snapshots stored in the **snapshot** directory. Below is an example of traj.gif:

<div align=center>
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/main/imgs/traj.gif" width="550"/>
</p>
</div>

In this animation, the solid box represents the lattice (here, rutile TiO<SUB>2</SUB> 2×2×3 supercell), and the color-coded circles indicate the lattice points corresponding to the selected atom type (here, oxygen). The **yellow-colored circles** mark the vacancy positions, while other colors denote occupied lattice points. Atomic movements are depicted with arrows matching the color of the moving atoms.

Occasionally, spurious vacancies may appear; hence, the detected number of vacancies exceeds the initial count in the system. This scenario arise when two or more atoms are assigned to the same lattice site and are commonly associated with non-vacancy-hopping mechanism (e.g., kick-out mechanism). Because such spurious vacancies typically vanish within a few hundred femtoseconds, they are referred to as **transient vacancies**.


In this animation, the transient vacancies are represented as **orange-colored circles**. Below three successive snapshots shows the kick-out mechanism:

<div align=center>
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/refs/heads/main/imgs/kick-out.png" width="800"/>
</p>
</div>



## 3. Extraction of effective hopping parameters

Navigate to the `Example2` directory and run:

For serial computation:

```bash
vachoppy -m p 0.07
```

For parallel computation:
```bash
mpirun -np {num_nodes} vachoppy -m p 0.07 --parallel
```

Here, the arguments are:

* t<SUB>interval</SUB> = 0.07 ps

For serial computation, the process is displayed via a progress bar. For parallel computation, process is recorded in the **VACHOPPY_PROGRESS** file in real time. The number and type of vacancies are automatically determined by comparing the POSCAR_LATTICE file with the AIMD datasets.

**Output:**

All results are stored in the **parameter.txt** file, which includes:

1.	A list of vacancy hopping paths in the system

2.	Effective hopping parameters (excluding z and ν)

3.	Vacancy hopping history for each AIMD dataset


To print the **vacancy hopping paths**, use:

```bash
awk '/Vacancy hopping paths :/ {f=1} f; /^$/ {f=0}' parameter.txt
```

To print the **effective hopping parameters**, use:
```bash
awk '/Effective/ {f=1} f; /^$/ {f=0}' parameter.txt
```

Below is the expected output:
<div align=center>
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/refs/heads/main/imgs/parameter.png" width="650"/>
</p>
</div>

The names of hopping paths are automatically assigned in ascending order of hopping distance. Symmetrically distinct sites are labeled with different alphabets (e.g., A1, A2, … for site 1 and B1, B2, … for site 2, where site 1 and site 2 are symmetrically distinct).

The correlation factor (f) is inherently temperature-dependent, but an average value across all simulated temperatures is displayed. To print f values for each temperature, use:

```bash
awk '/Cumulative/ {f=1} f; /^$/ {f=0}' parameter.txt
```
Below is the expected output:
<div align=center>
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/refs/heads/main/imgs/correlation_factor.png" width="650"/>
</p>
</div>


## 3-1. Extract effective values for z and ν (*optional*)

To obtain the effective values of z and ν, users must first perform NEB calculations for all vacancy hopping paths. The results of the NEB calculations should be stored in a file named **neb.csv** (see above).


Navigate to the `Example2` directory and run:

```bash
vachoppy -m pp
```

This command reads `parameter.txt` and `neb.csv` files and outputs `postprocess.txt` which contains the complete set of the effective hopping parameters.

To print the **effective hopping parameters**, use:

```bash
awk '/hopping parameter/ {f=1} f; /^$/ {f=0}' postprocess.txt
```
Below is the expected output:
<div align=center>
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/refs/heads/main/imgs/postprocess.png" width="650"/>
</p>
</div>
In the upper table, the effective hopping parameters averaged over all simulated temperatures are shown, while the temperature-dependent parameters are presented in the lower table.

<br>
<br>

Additionally, **VacHopPy** provides **individual attempt frequencies** for each hopping path. To print them, use
```bash
awk '/Jump attempt frequency/ {f=1} f; /^$/ {f=0}' postprocess.txt
```
Attempt frequencies are estimated based on statistical approach. Therefore, only the values for hopping paths with a sufficient number of hopping events can be considered reliable.


## 4. Assessment of lattice stability

>Download and unzip the **Example3** file linked above.

**VacHopPy** employs the fingerprint analysis proposed by Oganov *et al.* to assess lattice stability. The key quantities used in this analysis are the **fingerprint vector (*ψ*)** and the **cosine distance (d<SUB>cos</SUB>)**. Detailed descriptions can be found in [**this paper**](https://www.sciencedirect.com/science/article/pii/S0010465510001840).

### Fingerprint vector (*ψ*)

To construct *ψ*, three parameters are required: 

1. Threshold radius (**R<SUB>max</SUB>**) 
2. Bin size (**Δ**)
3. Standard deviation for Gaussian-smeared delta function (**σ**) 

A well-defined *ψ* satisfies *ψ*(r=0) = -1 and converges to 0 as r → ∞. Therefore, the user needs to set these parameters appropriately to ensure these conditions are met. The *ψ* can be generated by using `vachoppy -u fingerprint` command: 


Navigate to the `Example3` directory and run:
```bash
vachoppy -u fingerprint POSCAR_MONO 20.0 0.04 0.04 -d 
```

Here, the arguments are:

* Atomic structure = POSCAR_MONO (monoclinic HfO<SUB>2</SUB>)
* R<SUB>max</SUB> = 20 Å
* Δ =0.04 Å
* σ = 0.04 Å 


Using `-d` flag displays the resulting *ψ* in a pop-up window. Below is an example output (*ψ* for monoclinic HfO<SUB>2</SUB>):

<div align=center>
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/main/imgs/fingerprint_mono.png" width="550"/>
</p>
</div>

To enhance robustness of *ψ*, **VacHopPy** considers all possible atom pairs (e.g., Hf-Hf, Hf-O, and O-O) and concatenates them to construct a single well-defined *ψ*.

---

### Cosine distance (d<SUB>cos</SUB>)

Cosine distance (**d<SUB>cos</SUB>(x)**) quantifies structural similarity to a reference phase x, where a lower d<SUB>cos</SUB>(x) indicates a greater similarity. By analyzing variations in d<SUB>cos</SUB>(x) over time, users can **assess lattice stability** or **explore phase transitions** occurred in the AIMD simulations.


#### (1) Assessment of lattice stability

Navigate to the `Example3` directory and run:

For serial computation:
```bash
vachoppy -m f 0.07 20 0.04 0.04 -v vasprun_1600K.xml -p POSCAR_MONO
```

For parallel computation:
```bash
mpirun -np {num_nodes} vachoppy -m f 0.07 20 0.04 0.04 -v vasprun_1600K.xml -p POSCAR_MONO --parallel
```

Here, the arguments are:

* t<SUB>interval</SUB> = 0.07 ps
* R<SUB>max</SUB> = 20 Å
* Δ = 0.04 Å
* σ = 0.04 Å 

The `-v` flag specifies **vasprun.xml** file (default: vasprun.xml), where **vasprun_1600K.xml** contains the AIMD trajectory at 1600 K. The `-p` flag specifies the **reference phase** (default: POSCAR_REF), where **POSCAR_MONO** contains monoclinic HfO<SUB>2</SUB> lattice.

Results are stored in **cosine_distance.txt** and **cosine_distance.png**. To avoid overwriting, rename **cosine_distance.txt** to **dcos_1600K_mono.txt**.

----
Next, run:

For serial computation:
```bash
vachoppy -m f 0.07 20 0.04 0.04 -v vasprun_2200K.xml -p POSCAR_MONO
```
For parallel computation:
```bash
mpirun -np {num_nodes} vachoppy -m f 0.07 20 0.04 0.04 -v vasprun_2200K.xml -p POSCAR_MONO --parallel 
```

Here, **vasprun_2200K.xml** contains the AIMD data at 2200 K. Rename **cosine_distance.txt** to **dcos_2200K_mono.txt**.

----

For comparison, plot **dcos_1600K_mono.txt** and **dcos_2200K_mono.txt** together using `plot.py`:

```bash
# plot.py
import sys
import numpy as np
import matplotlib.pyplot as plt

data, num_data = [], len(sys.argv)-1
for i in range(num_data):
    data.append(np.loadtxt(sys.argv[i+1], skiprows=2))

plt.rcParams['figure.figsize'] = (6, 2.5)
plt.rcParams['font.size'] = 10

space = 0.09
for i, data_i in enumerate(data):
    data_i[:,1] -= np.average(data_i[:,1]) - space * i
    plt.scatter(data_i[:,0], data_i[:,1], s=10)
    
plt.yticks([])
plt.xlabel('Time (ps)', fontsize=12)
plt.ylabel(r'$d_{cos}$($x$)', fontsize=12)
plt.legend(loc='center right')
plt.show()
```

```bash
python plot.py dcos_1600K_mono.txt dcos_2200K_mono.txt
```

<div align=center>
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/refs/heads/main/imgs/dcos_1.png" width="550"/>
</p>
</div>

In this figure, the d<SUB>cos</SUB> data at each temperature is arranged vertically; hence, the absolute y-values are meaningless. Instead, the focus is on the relative change in d<SUB>cos</SUB> over time. 

* At 1600 K, d<SUB>cos</SUB> remains nearly constant, indicating structural stability.
* At 2200 K, d<SUB>cos</SUB> exhibits substantial fluctuations near 20 ps, suggesting that the monoclinic lattice becomes unstable at high temperatures.

It is important to note that the lattice parameters were constrained to those of monoclinic lattice since the AIMD simulations were performed under **NVT ensemble**. As a result, any lattice distortion is not sustained but instead revert to the original lattice, producing peaks in the d<SUB>cos</SUB> trace.

In unstable lattices, such as monoclinic HfO<SUB>2</SUB> at 2200 K, vacancies are poorly defined since atomic vibration centers may shift away from the original lattice point. Consequently, vacancy trajectory determination (`vachoppy -m t`) and effective hopping parameter extraction (`vachoppy -m p`) may lack accuracy.

---

#### (2) Exploring phase transition

By varying the reference phase, users can explore phase transitions occuring in AIMD simulations.

Navigate to the `Example3` directory and run:

For serial computation:
```bash
vachoppy -m f 0.07 20 0.04 0.04 -v vasprun_2200K.xml -p POSCAR_TET
```

For parallel computation:
```bash
mpirun -np {num_nodes} vachoppy -m f 0.07 20 0.04 0.04 -v vasprun_2200K.xml -p POSCAR_TET --parallel 
```
Here, **POSCAR_TET** contains the atomic structure of **tetragonal HfO<SUB>2</SUB>**. To avoid overwriting, rename **cosine_distance.txt** to **dcos_2200K_tet.txt**.

---

Next, run:

For serial computation:
```bash
vachoppy -m f 0.07 20 0.04 0.04 -v vasprun_2200K.xml -p POSCAR_AO
```

For parallel computation:
```bash
mpirun -np {num_nodes} vachoppy -m f 0.07 20 0.04 0.04 -v vasprun_2200K.xml -p POSCAR_AO --parallel 
```

Here, **POSCAR_AO** contains the atomic structure of **antipolar orthorhombic HfO<SUB>2</SUB>**. Rename **cosine_distance.txt** to **dcos_2200K_ao.txt**.

---

To compare the results, run `plot.py`:

```bash
python plot.py dcos_2200K_mono.txt dcos_2200K_tet.txt dcos_2200K_ao.txt
```

<div align=center>
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/refs/heads/main/imgs/dcos_2.png" width="550"/>
</p>
</div>

As before, the d<SUB>cos</SUB> data is arranged vertically, so the absolute y-values are not meaningful. Instead, the focus is on the relative change in d<SUB>cos</SUB> over time. 

* As d<SUB>cos</SUB>(*mono*) increases, 
* d<SUB>cos</SUB>(*tet*) decreases, 
* while d<SUB>cos</SUB>(*ao*) remain nearly constant. 

This result clearly suggests that the phase transition is directed toward the **tetragonal phase**.


## References
If you used **VacHopPy** package, please cite [**this paper**](https://arxiv.org/abs/2503.23467)

If you used `vachoppy -m f` or `vachoppy -u fingerprint` commands, also cite [**this paper**](https://www.sciencedirect.com/science/article/pii/S0010465510001840).
