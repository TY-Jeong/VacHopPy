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

## Additional Tools

In addition to the primary parameter extraction, `VacHopPy` also provides several secondary tools for more in--depth analysis:

* **Vacancy Trajectory Visualization**

  Visualize the movement and pathways of vacancies over time.

* **MSD-based Diffusivity Calculation**

  Calculate atomic diffusivity using mean squared displacement (MSD) analysis.

* **Structural Evolution Tracking**

  Identify and track the evolution of the material's structure throughout the simulation.

