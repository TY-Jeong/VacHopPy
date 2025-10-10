# VacHopPy (v 3.0.0)

---
**VacHopPy** is a Python package for analyzing vacancy hopping mechanisms based on molecular dynamics (MD) simulations. A detailed explanation on **VacHopPy** framwork is available [**here**](https://arxiv.org/abs/2503.23467).

<div align="center">
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/main/imgs/logo.png" width="550"/>
</p>
</div>


## Features

* Tracking of **vacancy trajectories** in MD simulations
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
| ν<SUP>*</SUP>| Atomic vibration frequency (THz)         |
</div>



## Installation

This package can be easily installed via pip.

```bash
pip install vachoppy
```

**VacHopPy** requires **Python 3.10 or higher** and the following Python packages:
* ase
* h5py
* tqdm
* numpy
* scipy
* joblib
* pandas
* Pillow
* plotly
* tabulate
* mdanalysis
* matplotlib >= 3.10.0
* pymatgen >= 2024.6.10




## References
If you used **VacHopPy** package, please cite [**this paper**](https://arxiv.org/abs/2503.23467)
