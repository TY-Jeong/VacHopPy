
<div align="center">
<p>
    <img src="https://raw.githubusercontent.com/TY-Jeong/VacHopPy/main/imgs/logo.png" width="550"/>
</p>
</div>

**`VacHopPy`** is a Python package for analyzing **vacancy-mediated diffusion** from molecular dynamics (MD) simulations. The package facilitates a reliable **multiscale modeling** workflow by calculating key **effectuve hopping parameters** that can be directly incorporated into continuum-scale models.

---
# Key Features

* **Effective Parameter Extraction** 

  Derives a single set of effective hopping parameters from the MD simulations, 
  providing ready-to-use inputs for continuum-scale models.

* **Ensemble Analysis**

  Simultaneously processes an ensemble of multiple MD trajectories. 
  This approach efficiently samples rare hopping events in high-barrier systems, 
  avoiding the need for a single, prohibitively long simulation.

* **Memory-Efficient Processing**

  Reads and interprets large-scale trajectories via a streaming approach, 
  enabling the analysis of massive datasets within just a few gigabytes of RAM.

---

# Effective Hopping Parameters

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

`VacHopPy` was developed with **multiscale modeling** in mind. It calculates the various parameters required for continuum-scale models, providing each as an effective value. These effective parameters are delivered in simple, ready-to-use forms—such as a single value or an Arrhenius form—making them directly suitable for continuum-level analysis.

The key effective hoppping parameters extractable with `VacHopPy` are summarized below.

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

---
# Getting Started

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

# Command-Line Interface (CLI)

In addition to its Python library, `VacHopPy` provides a powerful Command-Line Interface (CLI) to perform key tasks directly from your terminal. This allows you to run analyses quickly without writing custom scripts.

The general command structure is:

```bash
vachoppy <command> [options]
```

-----
# Available Commands

To see a list of all available commands and their brief descriptions, use the `-h` or `--help` option:

```bash
vachoppy -h
```

Here is a summary of the main commands available:

<div align="center">

| Command | Description |
|---|---|
| `trajectory` | Identify and visualize vacancy trajectories from a single trajectory file. |
| `analyze` | Extract hopping parameters from an ensemble of trajectories. |
| `vibration` | Extract atomic vibration frequency from a single trajectory. |
| `distance` | Trace change in cosine distance from a reference structure over time. |
| `fingerprint` | Calculate and plot the fingerprint for a single static structure. |
| `msd` | Calculate diffusivity from mean squared displacement (Einstein relation). |
| `convert` | Convert various MD trajectory formats to the standard HDF5 format. |
| `concat` | Concatenate two successive HDF5 trajectory files into a new one. |
| `show` | Display a metadata summary of a HDF5 trajectory file. |

</div>

----
# Documentation

You can find detailed description for the CLI commands and Python API on this [web site](https://vachoppy.readthedocs.io/en/latest/)

-----
# References
If you used `VacHopPy` package, please cite [**this paper**](https://arxiv.org/abs/2503.23467)
