# Vacancy Trajectory Analysis

To follow this tutorial, please download and unzip `Example2/` from [this link](https://drive.google.com/file/d/1xBj3iP4eUInB2OKxCHstel4CTyXUDU9t/view?usp=sharing):

----

`VacHopPy` provides the `trajectory` command to analyze and visualize the vacancy trajectories from a single HDF5 file. 

First, navigate into the `Example2/` directory you downloaded. You will find the two essential files:

```bash
cd path/to/Example2/
ls
# >> POSCAR_TiO2  TRAJ_O_01.h5
```

  * **`POSCAR_TiO2`** 

    This file contains the crystal structure of the **perfect, vacancy-free** material. While this example uses the VASP POSCAR format, any format supported by the **Atomic Simulation Environment (ASE)** is compatible. `VacHopPy` uses this file to determine the **lattice sites** and potential **vacancy hopping paths**.

  * **`TRAJ_O_01.h5`**: 
  
    This is the HDF5 file containing the MD trajectory data. The example system is **rutile TiO₂** containing **two oxygen vacancies**, simulated at **2100 K**.

You can inspect the contents of the HDF5 file using the `show` command.

```bash
vachoppy show TRAJ_O_01.h5
```

This command prints a summary of the stored simulation data:

```{code-block} bash
:class: scrollable-output

==================================================
  Trajectory File: TRAJ_O_01.h5
==================================================

[Simulation Parameters]
  - Atomic Symbol:      O
  - Number of Frames:   100000
  - Temperature:        2100.0 K
  - Time Step:          1.0 fs

[Composition]
  - Counts:             O: 46, Ti: 24
  - Total Atoms:        70

[Lattice Vectors (Ang)]
  [  9.29130,   0.00000,   0.00000]
  [  0.00000,   9.29130,   0.00000]
  [  0.00000,   0.00000,   8.90104]

[Stored Datasets]
  - positions:          Shape = (100000, 46, 3)
  - forces:             Shape = (100000, 46, 3)
==================================================
```

---

## How to Obtain Vacancy Trajectory

To run the analysis, execute the `trajectory` command with the paths to the **HDF5** and **structure file**, along with the **symbol of the mobile species**. Using the `--unwrap` flag is recommended to generate continuous trajectories that are not folded back into the unit cell.

```bash
vachoppy trajectory TRAJ_O_01.h5 POSCAR_TiO2 O --unwrap
```

---

## Understanding the Output

The command produces two types of output: a detailed summary printed to the terminal and two generated files.

### Terminal Output

The command prints a multi-step summary directly to the terminal, including:

* **Step 1**: Automatic estimation of the analysis time interval (`t_interval`).

* **Step 2**: Progress of the main trajectory analysis.

* **Step 3**: A summary of all unique hopping paths that were observed.

* **Step 4**: A detailed, time-ordered history of every hop for each vacancy.

```{code-block} bash
:class: scrollable-output

[STEP1] Automatic t_interval Estimation:
============================================================
               Automatic t_interval Estimation
============================================================
  Estimating from TRAJ_O_01.h5
             -> t_interval : 0.050 ps
============================================================
    Adjusting t_interval to the nearest multiple of dt
============================================================
    - dt                  : 0.0010 ps
    - Original t_interval : 0.0501 ps
    - Adjusted t_interval : 0.0500 ps (50 frames)
============================================================


[STEP2] Identifying Vacancy Trajectory:
Analyze Trajectory: 100%|##############################| 1/1 [00:00<00:00, 1952.66it/s]


[STEP3] Summary of Hopping Paths:
====================================================================================================================
                                              Hopping Path Information
====================================================================================================================
Name    a (Ang)    z    Count    Initial Site (Fractional Coordinate)    Final Site (Fractional Coordinate)
------  ---------  ---  -------  --------------------------------------  ------------------------------------
A1      2.56255    1    44       site1 [0.09751, 0.40249, 0.16667]       site1 [-0.09751, 0.59751, 0.16667]
A2      2.80308    8    50       site1 [0.09751, 0.40249, 0.16667]       site1 [0.34751, 0.34751, 0.33333]
A3      2.96701    2    2        site1 [0.09751, 0.40249, 0.16667]       site1 [0.09751, 0.40249, -0.16667]
====================================================================================================================


[STEP4] Summary of Hopping Histories:
====================================================================================================================
                                           Hopping Sequence of Vacancy 0
====================================================================================================================
Num    Time (ps)    Path    a (Ang)    Initial Site (Fractional Coordinate)    Final Site (Fractional Coordinate)
-----  -----------  ------  ---------  --------------------------------------  ------------------------------------
1      9.25         A1      2.56255    site1 [0.09751, 0.40249, 0.16667]       site1 [0.90249, 0.59751, 0.16667]
2      15.75        A2      2.80308    site1 [0.90249, 0.59751, 0.16667]       site1 [0.84751, 0.34751, 0.33333]
3      16.45        A1      2.56255    site1 [0.84751, 0.34751, 0.33333]       site1 [0.65249, 0.15249, 0.33333]
4      20           A1      2.56255    site1 [0.65249, 0.15249, 0.33333]       site1 [0.84751, 0.34751, 0.33333]
5      20.25        A1      2.56255    site1 [0.84751, 0.34751, 0.33333]       site1 [0.65249, 0.15249, 0.33333]
6      24.55        A1      2.56255    site1 [0.65249, 0.15249, 0.33333]       site1 [0.84751, 0.34751, 0.33333]
7      27.55        A3      2.96701    site1 [0.84751, 0.34751, 0.33333]       site1 [0.84751, 0.34751, 0.66667]
8      27.6         A2      2.80308    site1 [0.84751, 0.34751, 0.66667]       site1 [0.90249, 0.59751, 0.83333]
9      29           A2      2.80308    site1 [0.90249, 0.59751, 0.83333]       site1 [0.65249, 0.65249, 0.00000]
10     29.15        A2      2.80308    site1 [0.65249, 0.65249, 0.00000]       site1 [0.90249, 0.59751, 0.83333]
11     39.65        A2      2.80308    site1 [0.90249, 0.59751, 0.83333]       site1 [0.65249, 0.65249, 0.66667]
12     40.8         A2      2.80308    site1 [0.65249, 0.65249, 0.66667]       site1 [0.90249, 0.59751, 0.50000]
13     41.25        A2      2.80308    site1 [0.90249, 0.59751, 0.50000]       site1 [0.84751, 0.84751, 0.33333]
14     41.8         A2      2.80308    site1 [0.84751, 0.84751, 0.33333]       site1 [0.90249, 0.59751, 0.16667]
15     41.9         A2      2.80308    site1 [0.90249, 0.59751, 0.16667]       site1 [0.84751, 0.84751, 0.33333]
16     42.25        A2      2.80308    site1 [0.84751, 0.84751, 0.33333]       site1 [0.59751, 0.90249, 0.16667]
17     43.35        A2      2.80308    site1 [0.59751, 0.90249, 0.16667]       site1 [0.65249, 0.65249, 0.33333]
18     43.75        A2      2.80308    site1 [0.65249, 0.65249, 0.33333]       site1 [0.90249, 0.59751, 0.50000]
19     46.85        A2      2.80308    site1 [0.90249, 0.59751, 0.50000]       site1 [0.65249, 0.65249, 0.33333]
20     47.4         A2      2.80308    site1 [0.65249, 0.65249, 0.33333]       site1 [0.90249, 0.59751, 0.50000]
21     50.85        A2      2.80308    site1 [0.90249, 0.59751, 0.50000]       site1 [0.65249, 0.65249, 0.33333]
22     51.8         A1      2.56255    site1 [0.65249, 0.65249, 0.33333]       site1 [0.84751, 0.84751, 0.33333]
23     56.05        A2      2.80308    site1 [0.84751, 0.84751, 0.33333]       site1 [0.90249, 0.59751, 0.50000]
24     57.95        A2      2.80308    site1 [0.90249, 0.59751, 0.50000]       site1 [0.15249, 0.65249, 0.33333]
25     58.35        A2      2.80308    site1 [0.15249, 0.65249, 0.33333]       site1 [0.09751, 0.40249, 0.16667]
26     61.65        A1      2.56255    site1 [0.09751, 0.40249, 0.16667]       site1 [0.90249, 0.59751, 0.16667]
27     82.05        A2      2.80308    site1 [0.90249, 0.59751, 0.16667]       site1 [0.15249, 0.65249, 0.33333]
28     83.55        A2      2.80308    site1 [0.15249, 0.65249, 0.33333]       site1 [0.40249, 0.59751, 0.50000]
29     83.75        A1      2.56255    site1 [0.40249, 0.59751, 0.50000]       site1 [0.59751, 0.40249, 0.50000]
30     84           A1      2.56255    site1 [0.59751, 0.40249, 0.50000]       site1 [0.40249, 0.59751, 0.50000]
31     84.6         A1      2.56255    site1 [0.40249, 0.59751, 0.50000]       site1 [0.59751, 0.40249, 0.50000]
32     85.15        A2      2.80308    site1 [0.59751, 0.40249, 0.50000]       site1 [0.34751, 0.34751, 0.33333]
33     85.35        A1      2.56255    site1 [0.34751, 0.34751, 0.33333]       site1 [0.15249, 0.15249, 0.33333]
34     92.55        A2      2.80308    site1 [0.15249, 0.15249, 0.33333]       site1 [0.09751, 0.90249, 0.16667]
35     92.85        A2      2.80308    site1 [0.09751, 0.90249, 0.16667]       site1 [0.15249, 0.65249, 0.33333]
36     94.15        A1      2.56255    site1 [0.15249, 0.65249, 0.33333]       site1 [0.34751, 0.84751, 0.33333]
37     94.35        A2      2.80308    site1 [0.34751, 0.84751, 0.33333]       site1 [0.40249, 0.59751, 0.50000]
38     95.35        A2      2.80308    site1 [0.40249, 0.59751, 0.50000]       site1 [0.15249, 0.65249, 0.33333]
39     95.4         A2      2.80308    site1 [0.15249, 0.65249, 0.33333]       site1 [0.09751, 0.40249, 0.50000]
====================================================================================================================

====================================================================================================================
                                           Hopping Sequence of Vacancy 1
====================================================================================================================
Num    Time (ps)    Path    a (Ang)    Initial Site (Fractional Coordinate)    Final Site (Fractional Coordinate)
-----  -----------  ------  ---------  --------------------------------------  ------------------------------------
1      5.4          A2      2.80308    site1 [0.90249, 0.59751, 0.83333]       site1 [0.65249, 0.65249, 0.00000]
2      9.25         A1      2.56255    site1 [0.65249, 0.65249, 0.00000]       site1 [0.84751, 0.84751, 0.00000]
3      9.9          A1      2.56255    site1 [0.84751, 0.84751, 0.00000]       site1 [0.65249, 0.65249, 0.00000]
4      10.1         A1      2.56255    site1 [0.65249, 0.65249, 0.00000]       site1 [0.84751, 0.84751, 0.00000]
5      10.3         A1      2.56255    site1 [0.84751, 0.84751, 0.00000]       site1 [0.65249, 0.65249, 0.00000]
6      10.75        A1      2.56255    site1 [0.65249, 0.65249, 0.00000]       site1 [0.84751, 0.84751, 0.00000]
7      12.05        A1      2.56255    site1 [0.84751, 0.84751, 0.00000]       site1 [0.65249, 0.65249, 0.00000]
8      12.4         A2      2.80308    site1 [0.65249, 0.65249, 0.00000]       site1 [0.59751, 0.90249, 0.83333]
9      14.3         A2      2.80308    site1 [0.59751, 0.90249, 0.83333]       site1 [0.65249, 0.15249, 0.66667]
10     14.4         A1      2.56255    site1 [0.65249, 0.15249, 0.66667]       site1 [0.84751, 0.34751, 0.66667]
11     16.25        A2      2.80308    site1 [0.84751, 0.34751, 0.66667]       site1 [0.09751, 0.40249, 0.50000]
12     17.05        A2      2.80308    site1 [0.09751, 0.40249, 0.50000]       site1 [0.84751, 0.34751, 0.66667]
13     20.05        A2      2.80308    site1 [0.84751, 0.34751, 0.66667]       site1 [0.90249, 0.09751, 0.50000]
14     20.35        A1      2.56255    site1 [0.90249, 0.09751, 0.50000]       site1 [0.09751, 0.90249, 0.50000]
15     21.9         A2      2.80308    site1 [0.09751, 0.90249, 0.50000]       site1 [0.84751, 0.84751, 0.66667]
16     24.65        A1      2.56255    site1 [0.84751, 0.84751, 0.66667]       site1 [0.65249, 0.65249, 0.66667]
17     26.2         A2      2.80308    site1 [0.65249, 0.65249, 0.66667]       site1 [0.90249, 0.59751, 0.50000]
18     26.65        A1      2.56255    site1 [0.90249, 0.59751, 0.50000]       site1 [0.09751, 0.40249, 0.50000]
19     26.85        A1      2.56255    site1 [0.09751, 0.40249, 0.50000]       site1 [0.90249, 0.59751, 0.50000]
20     28.45        A2      2.80308    site1 [0.90249, 0.59751, 0.50000]       site1 [0.84751, 0.84751, 0.66667]
21     28.75        A1      2.56255    site1 [0.84751, 0.84751, 0.66667]       site1 [0.65249, 0.65249, 0.66667]
22     29.15        A3      2.96701    site1 [0.65249, 0.65249, 0.66667]       site1 [0.65249, 0.65249, 0.33333]
23     31.6         A2      2.80308    site1 [0.65249, 0.65249, 0.33333]       site1 [0.40249, 0.59751, 0.16667]
24     34.7         A1      2.56255    site1 [0.40249, 0.59751, 0.16667]       site1 [0.59751, 0.40249, 0.16667]
25     34.85        A2      2.80308    site1 [0.59751, 0.40249, 0.16667]       site1 [0.65249, 0.65249, 0.33333]
26     42.3         A1      2.56255    site1 [0.65249, 0.65249, 0.33333]       site1 [0.84751, 0.84751, 0.33333]
27     42.7         A1      2.56255    site1 [0.84751, 0.84751, 0.33333]       site1 [0.65249, 0.65249, 0.33333]
28     43           A1      2.56255    site1 [0.65249, 0.65249, 0.33333]       site1 [0.84751, 0.84751, 0.33333]
29     44.2         A1      2.56255    site1 [0.84751, 0.84751, 0.33333]       site1 [0.65249, 0.65249, 0.33333]
30     44.85        A1      2.56255    site1 [0.65249, 0.65249, 0.33333]       site1 [0.84751, 0.84751, 0.33333]
31     45           A1      2.56255    site1 [0.84751, 0.84751, 0.33333]       site1 [0.65249, 0.65249, 0.33333]
32     45.1         A1      2.56255    site1 [0.65249, 0.65249, 0.33333]       site1 [0.84751, 0.84751, 0.33333]
33     45.7         A1      2.56255    site1 [0.84751, 0.84751, 0.33333]       site1 [0.65249, 0.65249, 0.33333]
34     45.8         A1      2.56255    site1 [0.65249, 0.65249, 0.33333]       site1 [0.84751, 0.84751, 0.33333]
35     50.3         A1      2.56255    site1 [0.84751, 0.84751, 0.33333]       site1 [0.65249, 0.65249, 0.33333]
36     50.8         A1      2.56255    site1 [0.65249, 0.65249, 0.33333]       site1 [0.84751, 0.84751, 0.33333]
37     51           A2      2.80308    site1 [0.84751, 0.84751, 0.33333]       site1 [0.59751, 0.90249, 0.16667]
38     52.15        A2      2.80308    site1 [0.59751, 0.90249, 0.16667]       site1 [0.65249, 0.65249, 0.33333]
39     55.25        A2      2.80308    site1 [0.65249, 0.65249, 0.33333]       site1 [0.59751, 0.40249, 0.16667]
40     56.55        A1      2.56255    site1 [0.59751, 0.40249, 0.16667]       site1 [0.40249, 0.59751, 0.16667]
41     60.6         A2      2.80308    site1 [0.40249, 0.59751, 0.16667]       site1 [0.65249, 0.65249, 0.00000]
42     60.7         A2      2.80308    site1 [0.65249, 0.65249, 0.00000]       site1 [0.59751, 0.40249, 0.16667]
43     61.2         A2      2.80308    site1 [0.59751, 0.40249, 0.16667]       site1 [0.84751, 0.34751, 0.00000]
44     62.2         A2      2.80308    site1 [0.84751, 0.34751, 0.00000]       site1 [0.59751, 0.40249, 0.16667]
45     62.85        A2      2.80308    site1 [0.59751, 0.40249, 0.16667]       site1 [0.34751, 0.34751, 0.00000]
46     64.35        A2      2.80308    site1 [0.34751, 0.34751, 0.00000]       site1 [0.40249, 0.09751, 0.16667]
47     69.85        A2      2.80308    site1 [0.40249, 0.09751, 0.16667]       site1 [0.34751, 0.34751, 0.33333]
48     72.5         A1      2.56255    site1 [0.34751, 0.34751, 0.33333]       site1 [0.15249, 0.15249, 0.33333]
49     75.35        A1      2.56255    site1 [0.15249, 0.15249, 0.33333]       site1 [0.34751, 0.34751, 0.33333]
50     85.1         A2      2.80308    site1 [0.34751, 0.34751, 0.33333]       site1 [0.40249, 0.59751, 0.50000]
51     85.3         A1      2.56255    site1 [0.40249, 0.59751, 0.50000]       site1 [0.59751, 0.40249, 0.50000]
52     87.2         A1      2.56255    site1 [0.59751, 0.40249, 0.50000]       site1 [0.40249, 0.59751, 0.50000]
53     92.2         A1      2.56255    site1 [0.40249, 0.59751, 0.50000]       site1 [0.59751, 0.40249, 0.50000]
54     94.05        A1      2.56255    site1 [0.59751, 0.40249, 0.50000]       site1 [0.40249, 0.59751, 0.50000]
55     94.3         A2      2.80308    site1 [0.40249, 0.59751, 0.50000]       site1 [0.15249, 0.65249, 0.33333]
56     94.55        A1      2.56255    site1 [0.15249, 0.65249, 0.33333]       site1 [0.34751, 0.84751, 0.33333]
57     95.45        A2      2.80308    site1 [0.34751, 0.84751, 0.33333]       site1 [0.40249, 0.59751, 0.50000]
====================================================================================================================

Results are saved in 'trajectory.json'.
Trajectory is saved in 'trajectory.html'.

Execution Time: 8.251 seconds
Peak RAM Usage: 0.056 GB
```

If you manually specify the `t_interval` using the `--t_interval` flag, the automatic estimation process (Step 1) will be skipped.


### Generated Files

In addition to the terminal output, the command generates two files in your current directory:

* `trajectory.json` 

    Contains the raw, time-series **PBC-unwrapped Cartesian coordinates** for each vacancy, along with simulation metadata. This file is ideal for custom plotting or further analysis.

* `trajectory.html` 

    An **interactive 3D plot** of the vacancy trajectories, generated using Plotly.

### Interactive 3D Visualization

You can open `trajectory.html` in any web browser to rotate, zoom, and inspect the vacancy paths. Below is an example of the interactive plot embedded directly in this documentation.

<div class="plotly-iframe-container">
  <iframe src="../../_static/trajectory.html" width="100%" height="400px" frameborder="0"></iframe>
</div>

-----

## Working with trajectory.json File

The `trajectory.json` file provides a straightforward way to access the raw trajectory data for your own custom analysis. You can load and process it using a simple Python script.

```python
import json
import numpy as np

with open('trajectory.json', 'r') as f:
    contents = json.load(f)

# Vacancy trajectories
trajectory_vac0 = np.array(contents['Vacancy0'], dtype=np.float64)
trajectory_vac1 = np.array(contents['Vacancy1'], dtype=np.float64)

# Summary
print("\nSimulation Parameters")
print(f"  - Atomic Symbol     : {contents['symbol']}")
print(f"  - Temperature       : {contents['temperature']} K")
print(f"  - Num. of Vacancies : {contents['num_vacancies']}")
print(f"  - t_interval        : {contents['t_interval']} ps")

print("\nLattice Vectors (Å)")
for vector in contents['lattice']:
    print(f"  [{vector[0]:>9.5f}, {vector[1]:>9.5f}, {vector[2]:>9.5f}]")
    
print("\nVacancy Trajectories")
print(f"  - Shape of the 1st Trajectory : {trajectory_vac0.shape}")
print(f"  - Shape of the 2nd Trajectory : {trajectory_vac1.shape}")
```

Example output:

```bash
Simulation Parameters
  - Atomic Symbol     : O
  - Temperature       : 2100.0 K
  - Num. of Vacancies : 2
  - t_interval        : 0.05 ps

Lattice Vectors (Å)
  [  9.29130,   0.00000,   0.00000]
  [  0.00000,   9.29130,   0.00000]
  [  0.00000,   0.00000,   8.90104]

Vacancy Trajectories
  - Shape of the 1st Trajectory : (2000, 3)
  - Shape of the 2nd Trajectory : (2000, 3)
```

Each vacancy trajectory is stored as a List with the shape (n_steps, 3), representing the PBC-unwrapped Cartesian coordinates in Ångstroms at each analysis step.

To convert the step index back to simulation time in picoseconds, simply multiply it by the `t_interval` value stored in the file. For example, the coordinates at `trajectory_vac0[10]` correspond to the time `10 * contents['t_interval']` ps.
