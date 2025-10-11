# Getting Started

## Installation
You can install `VacHopPy` using pip:
```bash
pip install vachoppy
```

Alternatively, you can install the latest version in development from [VacHopPy GitHub](https://github.com/TY-Jeong/VacHopPy):

```bash
git clone git@github.com:TY-Jeong/VacHopPy.git
cd VacHopPy
pip install -e .
```

---
## Installation test

Once the `VacHopPy` package is successfully installed, you can import `VacHopPy` within Python interpreter:

```python
import vachoppy
```

Or the following command can be executed:
```bash
vachoppy -h
```

---

## Input File Preparation

`VacHopPy` does not directly read raw MD trajectory files like `vasprun.xml` for its main analysis. Instead, it uses the **HDF5 format** as its primary input.

HDF5 allows for highly efficient, streaming-based data access. As a result, `VacHopPy` can process massive trajectory datasets quickly while consuming minimal RAM (typically only a few gigabytes), even for simulations that are hundreds of gigabytes in size.

Before running an analysis, you must first convert your MD trajectory into this format.

### Converting Your Trajectory to HDF5

You can convert your files using either a simple command-line tool or the Python script.

#### Via the Command-Line Interface (CLI)

The most straightforward method is the built-in `convert` command.

```bash
vachoppy convert vasprun_TiO2.xml 2000.0 1.0 --label 2000K
```

The arguments are as follows:

* `vasprun_TiO2.xml`: The path to your source MD trajectory file.

* `2000.0`: The simulation temperature in Kelvin.

* `1.0`: The time step (dt) between frames in femtoseconds.

* `2000K`: An optional suffix to append to the output filenames.

For efficient storage, the converter automatically splits the trajectory by chemical species. The command above would generate two separate files: `TRAJ_Ti_2000K.h5` and `TRAJ_O_2000K.h5`. This conversion process is powered by Atomic Simulation Environment (**ASE**), enabling compatibility with a wide range of file [**formats supported by the ASE**](https://ase-lib.org/ase/io/io.html).

#### Via the Python script

The CLI command is a convenient wrapper for the parse_md function. You can achieve the same result within a Python script:

```python
from vachoppy.core import parse_md

parse_md(
  filename='vasprun_TiO2.xml',
  temperature=2000.0,
  dt=1.0,
  label='2000K'
)
```

````{note}
For **LAMMPS dump** format, `VacHopPy` automatically uses [**MDAnalysis**](https://www.mdanalysis.org) as a backend due to the limitaion of ASE. For more details and advanced options, consult the help message with `vachoppy convert -h` or refer to `parse_lammps` module in the API documentation.

