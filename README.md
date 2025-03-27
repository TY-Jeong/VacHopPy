# VacHopPy 

---
**VacHopPy** is a Python package for analyzing vacancy hopping mechanisms based on *Ab initio* molecular dynamics (AIMD) simulations. A detailed explanation on **VacHopPy** framwork is available in **here(*will be updated*)**.


<div align=center>
<p>
    <img src="./imgs/logo.png" width="550"/>
</p>
</div>


## Features

* Tracking of **vacancy trajectories** in AIMD simulations
* Extraction of **effective hopping parameter** set
* Assessment of lattice stability or **phase transitions**

<br>

**Effective hopping parameter** set, a key improvement of **VacHopPy**, is a single, representative set of hopping parameters, which is determined by integrating all possible hopping paths in a given system considering energetic and geometric properties. It includes the following parameters (the bar expression ($`\bar{x}`$) is used to emphasize the parameters are effective values):

<div align="center">

|<center>Symbol</center>|<center>Description</center>|
|:---:|---|
| $`\bar{a}`$ |Hopping distance (Å)|
| $`\bar{E}_{a}`$ |Hopping barrier (eV)|
| $`\bar{z}`$ |Coordination number|
| $`\bar{ν}`$ |Jump attempt frequency (THz)|
| $`f`$ |Correlation factor|

</div>


Effective hopping parameter set is beneficial for **multiscale modeling**, which bridges the *ab intio* calculations and larger-scale simulations, such as TCAD, continuum models, and KMC methods, since required mass transport quantities are expressed in a simple Arrhenius equation. For example, diffusion coefficient ($D$) is given by:

```math
D = \frac{1}{6}\bar{z}\bar{a}^{2}\bar{ν}f \cdot \exp(-\frac{\bar{E}_{a}}{k_{B}T})
```

<br>

where $k_B$ represents the Boltzmann constant. Note that the exact expression of $D$ consist multiple exponential terms, each corresponds to a distinct vacancy hopping path.


## Contents

* Installation
* List of commands
* How to implement
  1. Preparation
  2. Vacancy trajectory visualization
  3. Extraction of effective hopping parameters
  4. Assessment of lattice stability
* Reference
* License



## Installation

This package can be easily installed via pip. The current version of  **VacHopPy** was developed based on VASP 5.4.4.

```bash
pip intall vachoppy
```

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
        <td>Calculate effective hopping parameters (excluding &#772;z and &#772;ν)</td>
    </tr>
    <tr>
        <!-- <td>4</td> -->
        <td>pp</td>
        <td>Calculate &#772;z and &#772;ν (post-processing for -m p option)</td>
    </tr>
    <tr>
        <!-- <td>5</td> -->
        <td>f</td>
        <td>Perform fingerprint analyses</td>
    <tr>
        <td rowspan="6">-u<br>(utility)</td>
        <td>extract_force</td>
        <td>Extract FORCE file from vasprun.xml</td>
    </tr>
    <tr>
        <!-- <td>2</td> -->
        <td>concat_xdatcar</td>
        <td>Concatenate two XDATCAR files</td>
    </tr>
    <tr>
        <!-- <td>3</td> -->
        <td>concat_force</td>
        <td>Concatenate two FORCE files</td>
    </tr>
    <tr>
        <!-- <td>4</td> -->
        <td>update_outcar</td>
        <td>Combine two OUTCAR files</td>
    </tr>
    <tr>
        <!-- <td>5</td> -->
        <td>fingerprint</td>
        <td>Extract fingerprint</td>
    </tr>
    <tr>
        <!-- <td>6</td> -->
        <td>cosine_distance</td>
        <td>Calculate cosine distance</td>
    </tr>
</table>
</div>

For detailed descriptions, please use `-h` options:

```bash
vachoppy -h # list of available commands
vachoppy -m p -h # explanation for '-m p' option
```

For time-consuming commaands, `vachoppy -m p` and `vachoppy -m f`, parallelization is supported by **mpirun**. For parallelization, please specify `--parallel` option:

```bash
vachoppy -m p O 0.1 # serial computation
mpirun -np 10 vachoppy -m p O 0.1 --parallel # parallel computation with 10 cpu nodes.
```


Belows is summary of the main commands (only main modules are shown for clarity):
<div align=center>
<p>
    <img src="./imgs/flowchart2.svg" width="800"/>
</p>
</div>


## How to implement

Example files can be downloaded from:

* **Example1** : Vacancy hopping in rutile TiO<SUB>2</SUB> [download (29 GB)](https://drive.google.com/file/d/1SudMlQk40cJnVgrkklK6b4nhiF3YWPOY/view?usp=sharing)
* **Example2** : Phase transition of monoclinic HfO<SUB>2</SUB> at 2200 K  [download (102 MB)](https://drive.google.com/file/d/1SuxEHmGdVNkk-mogdWWDOOUPZqX74QG5/view?usp=sharing)

## 1. Preparation
### Input data
To run **VacHopPy**, the user needs four types of input data: **XDATCAR**, **OUTCAR**, **FORCE**, and **POSCAR_LATTICE**. In current version, **VacHopPy** supports only single-vacancy simulation, with multi-vacancy support planned for a future update.

#### (1) XDATCAR and OUTCAR 
XDATCAR and OUTCAR are standard VASP output files containing atomic trajectory from AIMD simulation and simulation conditions, respectively. The input atomic structure (*i.e.*, POSCAR) must include a **single vacancy**, and the AIMD simulation should be performed under the **NVT ensemble** (set **MDALGO = 2**, **NBLOCK = 1**).

#### (2) FORCE (*optinal*)
FORCE file contains force vectors acting on atoms and can be extracted from the **vasprun.xml** file (a standard VASP output) using the following command:
```bash
vachoppy -u extract_force -in vasprun.xml -out FORCE
```
The FORCE file helps assign atoms to corresponding lattice points based on their relatice positions to transition states. If the FORCE file is not provided, atomic occupancies will be determined solely based on proximity. Once atoms are assigned, the vacancy position is identified as an unoccupied lattice point.

#### (3) POSCAR_LATTICE
POSCAR_LATTICE contains the perfect crystal structure without a vacancy. Its lattice parameters must match those of input structure (POSCAR) of the AIMD simulations. This file is used to define the lattice points for vacancy identification.

#### (4) File organization
Since AIMD simulations commonly cover timescales shorter than nanoseconds, a single AIMD simulation may contain a few hopping events. However, since **VacHopPy** computes the effective hopping parameters in static manner, sufficient sampling of hopping events is necessary to ensure reliablilty. To address this, **VacHopPy** processes multiple AIMD datasets simultaneously. Each AIMD dataset is distinguished by a number appended after an underscore in the XDATCAR and FORCE file names (e.g., XDATCAR_01, FORCE_01). Below is an example of the recommended file structure:

```bash
Example1
 ┣ traj
 ┃ ┣ traj.1900K # AIMD simulations conducted at 1900 K
 ┃ ┃ ┣ XDATCAR_01, FORCE_01 # Simiulations in the same directory should be 
 ┃ ┃ ┣ XDATCAR_02, FORCE_02 # conducted under the same conditions
 ┃ ┃ ┣ XDATCAR_03, FORCE_03
 ┃ ┃ ┗ OUTCAR
 ┃ ┣ traj.2000K
 ┃ ┃ ┣ XDATCAR_01, FORCE_01
 ┃ ┃ ┣ XDATCAR_02, FORCE_02
 ┃ ┃ ┣ XDATCAR_03, FORCE_03
 ┃ ┃ ┗ OUTCAR
 ┃ ┗ traj.2100K
 ┃ ┃ ┣ XDATCAR_01, FORCE_01
 ┃ ┃ ┣ XDATCAR_02, FORCE_02
 ┃ ┃ ┣ XDATCAR_03, FORCE_03
 ┃ ┃ ┗ OUTCAR
 ┗ POSCAR_LATTICE # POSCAR of the perfect crystal
```

Simulations at the same temperature should be conducted under identical conditions. In other words, **NSW** and **POTIM** tags in INCAR should be the same. Therefore, only one OUTCAR file is needed per temperature directory.


### Hyperparameter: t<SUB>interval</SUB>
To run **VacHopPy**, the user needs to determine one hyperparameter, **t<SUB>interval</SUB>**, in advance. This parameter defines the time interval for averaging atomic positions and forces. Thermal fluctuations in AIMD simulations can obscure precise atomic occupancy determination. However, since these fluctuations are random, they can be effectively averaged out over time. **VacHopPy** processes AIMD data by dividing it into segments of length of t<SUB>interval</SUB>. Each segment represents a single step in the analysis. The total number of steps is given by t<SUB>simulation</SUB>/t<SUB>interval</SUB>, where t<SUB>simulation</SUB> is the total AIMD simulation time.


<div align=center>
<p>
    <img src="./imgs/t_interval.jpg" width="800"/>
</p>
</div>

Choosing an appropriate t<SUB>interval</SUB> is crucial for reliable analysis. The t<SUB>interval</SUB> should be large enough to mitigate thermal fluctuations but short enough to prevent multiple hopping events from being included in a single step. A typical value is around 0.1 ps, through it may vary depending on the system. 

One recommended approach for determining the optimal t<SUB>interval</SUB> is through convergence tests using the correlation factor ($f$). Below is an example of a convergence test (example system: rutile TiO<SUB>2</SUB>):

<div align=center>
<p>
    <img src="./imgs/convergence_test.jpg" width="550"/>
</p>
</div>

The left and rigut figures show the convergences of $f$ with respect to the number of AIMD datasets (N<SUB>AIMD</SUB>) and t<SUB>interval</SUB>, respectively, at each temperature. The results confirm that convergence is achieved at **N<SUB>AIMD</SUB>=20** and **t<SUB>interval</SUB>=0.1 ps**. (You can obtain the same reulsts using the data in **Example** 1 above)


## 2. Vacancy trajectory visualization

>Download and unzip the **Example1** file linked above.

Navigate to the `Example1` directory and run:
```bash
 vachoppy -m t O 0.1 2100 03 -v # vacancy type, t_interval, temperature, label
 ```

Here, the arguments are:

* vacancy type : O (oxygen vacancy)
* t<SUB>interval</SUB> = 0.1 ps
* temperature = 2100 K
* label = 03

Therefore, this command visualizes the oxygen vacancy trajectory of **XDATCAR_03** in **traj.2100K** directory. Using the `-v` option (verbosity tag) prints the computational details, including vacancy hopping paths and vacancy hopping histories.


**Output:**

The result is saved as `traj.gif`, with individual snapshots stored in the `snapshot` directory. Below is an example of traj.gif:

<div align=center>
<p>
    <img src="./imgs/traj.gif" width="550"/>
</p>
</div>

In this animation, the solid box represents the lattice (here, rutile TiO<SUB>2</SUB>), and the color-coded circles indicate the lattice points corresponding to the selected atom type (here, oxygen). The **yellow-colored circle** marks the vacancy position (*i.e.*, the unoccupied lattice point), while other colors denote occupied lattice points. Atomic movements are depicted with arrows matching the color of the moving atoms. The user can adjust the resolution of the animation using the ``--dpi`` option (**default: 300**).


## 3. Extraction of effective hopping parameters

Navigate to the `Example1` directory and run:

```bash
# For serial computation
vachoppy -m p O 0.1 # symbol, t_interval

# For parallel computation
mpirun -np 10 vachoppy -m p O 0.1 --parallel # 10: number of cpu nodes
```
Here, the arguments are:

* vacancy type = O (oxygen vacancy)
* t<SUB>interval</SUB> = 0.1 ps


For serial computation, the process is displayed via a progress bar. For parallel computation, process is recorded in the `VACHOPPY_PROGRESS` file in real time.

**Output:**

All results are stored in `parameter.txt` file, which includes: 

1. A list of **vacancy hopping paths** in the system 
2. **Effective hopping parameters** (except for $\bar{z}$ and $\bar{\nu}$)
3. **Vacancy hopping history** for each AIMD dataset.

To find the **effective hopping parameters**, search for **Effective hopping parameters** in `parameter.txt` file:

<div align=center>
<p>
    <img src="./imgs/parameters.png" width="650"/>
</p>
</div>

---

The `vachoppy -m p` command extracts effective hopping parameters except for $\bar{z}$ and $\bar{\nu}$. To calculate $\bar{z}$ and $\bar{\nu}$, the user needs an additional input data, `neb.csv`. This file contains **hopping barriers ($E_{a}$)** for all vacancy hopping paths in the system. Below is an example of `neb.csv` (example system: rutile TiO<SUB>2</SUB>):

```bash
# neb.csv
A1,0.8698
A2,1.058
A3,1.766
```
Here, the **first column** corresponds to the **path names**, and the **second column** contains the **$E_{a}$ values**. The user can find the hopping path information in the `parameter.txt` file under **Vacancy hopping paths**:

<div align=center>
<p>
    <img src="./imgs/hopping_paths.png" width="650"/>
</p>
</div>

**Recommendation:** It is highly recommended to **perform NEB calculations using a larger supercell** than that used in AIMD simulations. In AIMD, thermal fluctuations attenuate interactions with periodic images and provide a broader sampling of atomic configurations, which helps approximate the effects of a larger supercell.


---

Navigate to the `Example1` directory and run:

```bash
vachoppy -m pp
```

This command reads `parameter.txt` and `neb.csv` files and outputs `postprocess.txt` which contains the complete set of the effective hopping parameters. To find the final values, search for **Effective hopping parameters** in the `postprocess.txt`:

<div align=center>
<p>
    <img src="./imgs/postprocess.png" width="650"/>
</p>
</div>

Additionally, **VacHopPy** provides **individual jump atempt frequencies** for each hopping paths. Find **Jump attempt frequency (THz)** in `postprocess.txt`:

<div align=center>
<p>
    <img src="./imgs/nu.png" width="650"/>
</p>
</div>


## 4. Assessment of lattice stability

>Download and unzip the **Example2** file linked above.

**VacHopPy** employs the fingerprint analysis proposed by Oganov *et al.* to assess lattice stability. The key quantities used in this analysis are the **fingerprint vector (*ψ*)** and the **cosine distance (d<SUB>cos</SUB>)**. Detailed descriptions can be found in [**this paper**](https://www.sciencedirect.com/science/article/pii/S0010465510001840).

### Fingerprint vector (*ψ*)

To construct *ψ*, three parameters are required: 

1. Threshold radius (**R<SUB>max</SUB>**) 
2. Bin size (**Δ**)
3. Standard deviation for Gaussian-smeared delta function (**σ**) 

A well-defined *ψ* satisfies *ψ*(r=0) = -1 and converges to 0 as r → ∞. Therefore, the user needs to set these parameters appropriately to ensure these conditions are met. The *ψ* can be generated by using `vachoppy -u fingerprint` command: 


Navigate to the `Example2` directory and run:
```bash
vachoppy -u fingerprint POSCAR_MONO 20.0 0.04 0.04 -d # POSCAR, R_max, Δ, σ
```

Here, the arguments are:

* Atomic structure = POSCAR_MONO (monoclinic HfO<SUB>2</SUB>)
* R<SUB>max</SUB> = 20 Å
* Δ =0.04 Å
* σ = 0.04 Å 


Using `-d` option displays the resulting *ψ* in a pop-up window. Below is an example output (*ψ* for monoclinic HfO<SUB>2</SUB>):

<div align=center>
<p>
    <img src="./imgs/fingerprint_mono.png" width="550"/>
</p>
</div>

To enhance robustness of *ψ*, **VacHopPy** considers all possible atom pairs (e.g., Hf-Hf, Hf-O, and O-O) and concatenates them to construct a single well-defined *ψ*.

---

### Cosine distance (d<SUB>cos</SUB>)

Cosine distance (**d<SUB>cos</SUB>($x$)**) quantifies structural similarity to a reference phase $x$, where a lower d<SUB>cos</SUB>($x$) indicates a greater similarity. By analyzing variations in d<SUB>cos</SUB>($x$) over time, users can **assess lattice stability** or **explore phase transitions** occurred in the AIMD simulations.


#### (1) Assessment of lattice stability

Navigate to the `Example2` directory and run:

```bash
# For serial computation
vachoppy -m f 0.05 20 0.04 0.04 -x XDATCAR_1600K -p POSCAR_MONO -o OUTCAR_1600K

# For parallel computation
mpirun -np 10 vachoppy -m f 0.05 20 0.04 0.04 -x XDATCAR_1600K -p POSCAR_MONO -o OUTCAR_1600K --parallel
```

Here, the arguments are:

* t<SUB>interval</SUB> = 0.05 ps
* R<SUB>max</SUB> = 20 Å
* Δ = 0.04 Å
* σ = 0.04 Å 

The `-x` option specifies **XDATCAR** file (**default: XDATCAR**), where `XDATCAR_1600K` contains the AIMD trajectory at 1600 K. The `-p` option specifies the **reference phase** (**default: POSCAR_REF**), where `POSCAR_MONO` contains **monoclinic HfO<SUB>2</SUB>** lattice. The `-o` option specifies the **OUTCAR** file (**default: OUTCAR**). 

Results are stored in `cosine_distance.txt` and `cosine_distance.png`. To prevent overwriting, rename `cosine_distance.txt` to `dcos_1600K_mono.txt`.

----
Next, run:

```bash
# For serial computation
vachoppy -m f 0.05 20 0.04 0.04 -x XDATCAR_2200K -p POSCAR_MONO -o OUTCAR_2200K

# For parallel computation
mpirun -np 10 vachoppy -m f 0.05 20 0.04 0.04 -x XDATCAR_2200K -p POSCAR_MONO -o OUTCAR_2200K --parallel 
```

Here, `XDATCAR_2200K` and `OUTCAR_2200K` contain the AIMD trajecoty and simulation conditions at 2200 K, respectively. Rename `cosine_distance.txt` to `dcos_2200K_mono.txt`.

----

For comparison, plot `dcos_1600K_mono.txt` and `dcos_2200K_mono.txt` together using `plot.py`:

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
    <img src="./imgs/dcos_1.png" width="550"/>
</p>
</div>

In this figure, the d<SUB>cos</SUB> data at each temperature is arranged vertically; hence, the absolute y-values are meaningless. Instead, the focus is on the relative change in d<SUB>cos</SUB> over time. 

* At 1600 K, d<SUB>cos</SUB> remains nearly constant, indicating structural stability.
* At 2200 K, d<SUB>cos</SUB> exhibits substantial fluctuations near 20 ps, suggesting that the monoclinic lattice becomes unstable at high temperatures.

It is important to note that the lattice parameters were contrained to those of monoclinic lattice since the AIMD simulations were performed under **NVT ensmeble**. As a result, any lattice distortion is not sustained but instead revert to the original lattice, producing peaks in the d<SUB>cos</SUB> trace.

In unstable lattices, such as monoclinic HfO<SUB>2</SUB> at 2200 K, vacancies are poorly defined since atomic vibrtaion centers may shift away from the original lattice point. Consequently, vacancy trajectory determination (`vachoppy -m t`) and effective hopping parameter extraction (`vachoppy -m p`) may lack accuracy.

---

#### (2) Exploring phase transition

By varying the reference phase, users can explore phase transitions occuring in AIMD simulatoins.

Navigate to the `Example2` directory and run:
```bash
# For serial computation
vachoppy -m f 0.05 20 0.04 0.04 -x XDATCAR_2200K -p POSCAR_TET -o OUTCAR_2200K

# For parallel computation
mpirun -np 10 vachoppy -m f 0.05 20 0.04 0.04 -x XDATCAR_2200K -p POSCAR_TET -o OUTCAR_2200K --parallel 
```
Here, `POSCAR_TET` contains the atomic structure of **tetragonal HfO<SUB>2</SUB>**. To prevent overwriting, rename `cosine_distance.txt` to `dcos_2200K_tet.txt`.

---

Next, run:
```bash
# For serial computation
vachoppy -m f 0.05 20 0.04 0.04 -x XDATCAR_2200K -p POSCAR_AO -o OUTCAR_2200K

# For parallel computation
mpirun -np 10 vachoppy -m f 0.05 20 0.04 0.04 -x XDATCAR_2200K -p POSCAR_AO -o OUTCAR_2200K --parallel 
```
Here, `POSCAR_AO` contains the atomic structure of **antipolar orthorhombic HfO<SUB>2</SUB>**. Rename `cosine_distance.txt` to `dcos_2200K_ao.txt`.

---

To compare the results, run `plot.py`:

```bash
python plot.py dcos_2200K_mono.txt dcos_2200K_tet.txt dcos_2200K_ao.txt
```

<div align=center>
<p>
    <img src="./imgs/dcos_2.png" width="550"/>
</p>
</div>

As before, the d<SUB>cos</SUB> data is arranged vertically, so the absolute y-values are not meaningful. Instead, the focus is on the relative change in d<SUB>cos</SUB> over time. 

* As d<SUB>cos</SUB>(*mono*) increases, 
* d<SUB>cos</SUB>(*tet*) decreases, 
* while d<SUB>cos</SUB>(*ao*) remain nearly constant. 

This result clearly suggets that the phase transition is directed toward the **tetragonal phase**.


## Reference
If you used **VacHopPy** package, plese cite **this paper(*will be updated*: please contact helianthus312@gmail.com)**

If you used `vachoppy -m f` or `vachoppy -u fingerprint` commands, also cite [**this paper**](https://www.sciencedirect.com/science/article/pii/S0010465510001840).
