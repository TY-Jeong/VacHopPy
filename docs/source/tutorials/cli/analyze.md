# Hopping Parameter Extraction

To follow this tutorial, please download and unzip `Example3/` from [this link](https://drive.google.com/file/d/1xBj3iP4eUInB2OKxCHstel4CTyXUDU9t/view?usp=sharing):

----

## How to Extract Hopping Parameters

The core feature of `VacHopPy`—calculating **effective hopping parameters**—is automated by the `analyze` command. The basic syntax is as follows:

```bash
vachoppy analyze [PATH_TRAJ] [PATH_STRUCTURE] [SYMBOL] --neb [NEB_CSV]
```

### Command Arguments

This command takes the following primary arguments:

* **`PATH_TRAJ`**

    Path to either a single HDF5 trajectory file or a root directory containing multiple HDF5 files. If files from multiple temperatures are provided, an Arrhenius-type analysis is automatically performed.

* **`PATH_STRUCTURE`**

    Path to a structure file of the perfect, vacancy-free material. While this example uses VASP POSCAR, any format supported by the **Atomic Simulation Environment (ASE)** is compatible. `VacHopPy` uses this file to determine lattice sites and potential hopping paths.

* **`SYMBOL`**

    The chemical symbol of the diffusing species (e.g., O for oxygen).

* **`--neb [NEB_CSV]`** (Optional)

    Path to a CSV file containing the activation energy barrier for each hopping path, as pre-calculated from a method like Nudged Elastic Band (NEB).


### The NEB Data File

This optional CSV file provides essential energy information for calculating certain parameters like the attempt frequency ($\nu$) and coordination number ($z$). The format is a simple two-column list of path names and their corresponding energy barriers in eV.

Example `neb_TiO2.csv`:

```bash
A1,0.8698
A2,1.058
A3,1.766
```
The required path names (e.g., `A1`, `A2`) can be identified using the `Site` class.


```python
from vachoppy.core import Site

site = Site([PATH_STRUCTURE], [SYMBOL]) 
site.summary()
```

Running the `summary()` method provides the necessary structural information about crystal sites and potential hopping paths.

Example output:

```bash
==============================================================================================================
                                             Site Analysis Summary
==============================================================================================================
  - Structure File       : POSCAR_TiO2
  - Diffusing Symbol     : O
  - Inequivalent Sites   : 1 found
  - Inequivalent Paths   : 3 found (with Rmax = 3.25 Å)

-- Hopping Path Information --

Name    Init Site    Final Site    a (Å)    z    Initial Coord (Frac)      Final Coord (Frac)
------  -----------  ------------  -------  ---  ------------------------  -------------------------
A1      site1        site1         2.5626   1    [0.0975, 0.4025, 0.1667]  [-0.0975, 0.5975, 0.1667]
A2      site1        site1         2.8031   8    [0.0975, 0.4025, 0.1667]  [0.3475, 0.3475, 0.3333]
A3      site1        site1         2.967    2    [0.0975, 0.4025, 0.1667]  [0.0975, 0.4025, -0.1667]
==============================================================================================================
```


### List of Extractable Parameters

The specific effective hopping parameters that `VacHopPy` can calculate depend entirely on the data you provide. There are two factors:

* **Number of Temperatures**

    Do your HDF5 files represent a single temperature or a range of temperatures?

* **NEB Data** 

    Have you provided a CSV file with pre-calculated hopping barriers?

The following table summarizes what can be calculated under different conditions. (**O** = Calculable, **X** = Not Calculable).

<div align="center">

| Parameter | Symbol | Single T | Multiple T | NEB Data Required? |
| :--- | :---: | :---: | :---: | :---: |
| Diffusivity | $D$ | O | O | No |
| Residence Time | $\tau$ | O | O | No |
| Correlation Factor | $f$ | O | O | No |
| Hopping Barrier | $E_a$ | X | O | No |
| Hopping Distance | $a$ | O | O | No |
| Coordination Number | $z$ | X | O | **Yes** |
| Attempt Frequency | $\nu$ | O | O | **Yes** |

</div>


-----

## Understanding the Output

Navigate into the `Example2/` directory you downloaded. You will find the three files (or directories):

