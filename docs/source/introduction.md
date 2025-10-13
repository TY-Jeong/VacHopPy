# Introduction

## List of Extractable Parameters

The key parameters extractable with `VacHopPy` are summarized below.

<div align="center">

| Symbol | Parameter                  | Unit |
| :----: | :------------------------- | :--: |
|  $D$   | Vacancy diffusivity        | m²/s |
| $\tau$ | Vacancy residence time     |  ps  |
|  $f$   | Correlation factor         |  -   |
| $E_a$  | Hopping barrier            |  eV  |
|  $a$   | Hopping distance           |  Å   |
|  $z$   | Coordination number        |  -   |
| $\nu$  | Attempt frequency          | THz  |
| $\nu^*$| Atomic vibration frequency | THz  |

</div>

All parameters listed above are calculated as **effective values**. For a detailed explanation of the concept of effective hopping parameters, please refer to the following section.


---

## What are effective hopping parameters?

In *ab initio* calculations, such as the nudged elastic band (NEB) method, hopping parameters are typically determined for a specific, predefined migration path. Consequently, in a system with multiple distinct migration pathways, one obtains a collection of different hopping parameter sets, each corresponding to a unique path.

However, to compare these computational results with experimental data or to use them in continuum-scale models, a single, representative set of parameters that encapsulates the overall diffusion behavior is required. We term this representative set the **effective hopping parameters**, which serve as ready-to-use inputs for continuum-level analysis.

The effective hopping parameters self-consistently satisfy the following fundamental diffusion equations:

$$
D = \frac{1}{6} f z a^2 \nu \cdot \exp(-E_a / k_B T)
$$

$$
\Gamma = \frac{1}{\tau} = z \nu \cdot \exp(-E_a / k_B T)
$$

For a more detailed theoretical background, please refer to [this paper](https://arxiv.org/abs/2503.23467).

---


## Understanding Data Requirements

`VacHopPy` is designed to be flexible, capable of processing either a single MD trajectory or an ensemble of trajectories simultaneously. Furthermore, it can analyze data from simulations performed at multiple temperatures to extract temperature-dependent properties.

However, not all parameters can be calculated from the same set of data. The table below summarizes the data requirements for each parameter, where 'O' indicates that the parameter can be calculated, and 'X' indicates that it cannot.


<div align="center">

| Symbol | Single Temperature | Multiple Temperatures | Note              |
| :----: | :----------------: | :-------------------: | :---------------: |
|  $D$   | O                  | O                     | -                 |
| $\tau$ | O                  | O                     | -                 |
|  $f$   | O                  | O                     | -                 |
| $E_a$  | X                  | O                     | -                 |
|  $a$   | O                  | O                     | -                 |
|  $z$   | X                  | O                     | NEB data required |
| $\nu$  | O                  | O                     | NEB data required |
| $\nu^*$| O                  | O                     | -                 |

</div>

As shown, parameters like the effective hopping barrier ($E_a$) require sampling MD trajectories across a range of temperatures to perform an Arrhenius-type analysis.

The "NEB data required" note signifies that for certain parameters ($z$ and $\nu$), additional data from NEB calculations is necessary. `VacHopPy` needs the hopping barrier for each distinct migration path in the system. However, since these energy barriers can be highly sensitive to the local atomic environment, this approach is considered reliable only for **monovacancy systems**.


### Identifying Hopping Paths for NEB

To assist in setting up these required NEB calculations, `VacHopPy` provides the `Site` class to identify potential hopping paths in your crystal structure. 

```python
from vachoppy.core import Site

site = Site('POSCAR_TiO2', 'O') 
site.summary()
```
The `Site` class is initialized with two arguments:

1. The **path to a structure file** of the perfect, vacancy-free crystal (e.g., 'POSCAR_TiO2'). It supports any format compatible with the Atomic Simulation Environment ([ASE](https://ase-lib.org/ase/io/io.html)).

2. The **chemical symbol** of the species whose vacancy is being analyzed (e.g., 'O' for an oxygen vacancy).

Running the `summary()` method provides the necessary structural information about crystal sites and potential migration paths to serve as a starting point for your NEB calculations.

---

## Additional Tools

In addition to the primary parameter extraction, `VacHopPy` also provides several secondary tools for more in--depth analysis:

* **Vacancy Trajectory Visualization**

  Visualize the movement and pathways of vacancies over time.

* **MSD-based Diffusivity Calculation**

  Calculate atomic diffusivity using mean squared displacement (MSD) analysis.

* **Structural Evolution Tracking**

  Identify and track the evolution of the material's structure throughout the simulation.

