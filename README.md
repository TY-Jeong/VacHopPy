# VacHopPy 

---
**VacHopPy** is a Python package to analyze trajectory of vacancy hopping mechanisms, based on *Ab initio* molecular dynamics (AIMD) simulations.


<div align=center>
<p>
    <img src="./imgs/logo.png" width="550"/>
</p>
</div>


A Key improvement in **VacHopPy** is introduction of an **effective hopping parameter** set. The hopping parameters inherently depend on the selction of hopping paths, resulting in multiple sets of hopping parameters within a given lattice. However, for macroscopic simulations (e.g., TCAD, continuum models, KMC methods), a well-defined, single effective hopping parameter set is required. This is because macroscopic equations typically do not account for multiple hopping paths. To sum up, the **effective hopping parameter** set is a single, consolidated parameter set that represents all hopping paths within a given lattice, making it suitable for multiscale modeling.


<div align=center>
<p>
    <img src="./imgs/fig2.jpg" width="550"/>
</p>
</div>


## Features

* Simultaneous calculation of **effective hopping parameters**
* Tracking of **vacancy trajectories** in AIMD simulations
* Assessment of lattice stability or **phase transitions**

 Belows are a **list of effective hopping parameters** which can be obtained from **VacHopPy**:

<div align="center">

|<center>Symbol</center>|<center>Description</center>|
|:---:|---|
|D<SUB>0</SUB>|Pre-exponential of diffusion coefficient (m<SUP>2</SUP>/s)|
|E<SUB>a</SUB>|Hopping barrier (eV)|
|a|Hopping distance (Å)|
|z|Coordination number|
|ν|Jump attempt frequency (THz)|
|f|Correlation factor|

</div>



## Contents

* Installation
* List of commands
* How to implement
  * Vacancy trajectory determination
    * Making animation
    * Distribution of hopping path
  * Effective hopping parameter calculation
    * Diffusion coefficient 
    * Atomic vibration coefficient
  * Assessment of lattice stability
  

## Installation

This package can be easily installed via pip

```ruby
pip intall vachoppy
```

The latest **VacHopPy** was developed based on VASP 5.4.4 and Python 3.12.4

## List of commands

**VacHopPy** provides a command-line interface (CLI). Belows are available CLI commands:

<div align=center>
<table>
    <tr>
        <th scope="col">Option 1</td>
        <th scope="col">Option 2</td>
        <th scope="col">Use</td>
    </tr>
    <tr>
        <td rowspan="5">-m<br>(main)</td>
        <td>p</td>
        <td>Calculate effective hopping parameters (excluding z and ν)</td>
    </tr>
    <tr>
        <!-- <td>2</td> -->
        <td>pp</td>
        <td>Calculate z and ν (post-processing for `-m p` option)</td>
    </tr>
    <tr>
        <!-- <td>3</td> -->
        <td>e</td>
        <td>Calculate diffusion coefficient using Einstein relation</td>
    </tr>
    <tr>
        <!-- <td>4</td> -->
        <td>t</td>
        <td>Make an animation for vacancy trajectories</td>
    </tr>
    <tr>
        <!-- <td>5</td> -->
        <td>f</td>
        <td>Perform fingerprint analyses</td>
    </tr>
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
```ruby
vachoppy -h # list of available commands
vachoppy -m p -h # explanation for '-m p' option
```

## How to implement

Example files can be downloaded from:

* **Example 1** : Vacancy trajectory determination & Effective hopping parameter calculations (28 GB) click
* **Example 2** : Assessment of lattice stability or phase transitions (10 GB) click

## 0. Preparation
**VacHopPy** reads AIMD simulation data in VASP format (XDATCAR, OUTCAR, and FORCE). **XDATCAR** and **OUTCAR** are the typical VASP output files, contain information on atomic positions and simulation conditions, respectively. **FORCE** (optinal) includes force vectors and can be extracted from **vasprun.xml** file using `vachoppy -u extract_force` command. If FORCE files are included in the input dataset, the trajectory is determined based on transition state (TS) distribution; otherwise, the trajectory is determined based simply on proximity.

> In current version, **VacHopPy** supports only AIMD simulations conducted using the **NVT ensmeble**. Each ensemble cell should contains a single vacancy. (Support for multi vacancies will be added in a future update) 

Since AIMD simulations are commonly conducted on time scales shorter than nanoseconds, a single AIMD simulation includes a limited number of hopping events. To overcome this limitation, **VacHopPy** can simultaneously process multiple bundles of AIMD simulation results, each belonging to the same NVT ensemble group. Each bundle is distinguished by a number appended after an underscore in the XDATCAR and FORCE file names (e.g., XDATCAR_01, FORCE_01). Below is an example of file tree:


```ruby
Example1
 ┣ traj
 ┃ ┣ traj.1900K # AIMD simulations conducted at 1900 K
 ┃ ┃ ┣ XDATCAR_01, FORCE_01 # Simiulations should be 
 ┃ ┃ ┣ XDATCAR_02, FORCE_02 # conducted in the same condition
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
 ┗ POSCAR_LATTICE # POSCAR of perfect cell
```

In this example, AIMD simulations were performed at three temperatures (1900 K, 2000 K, and 2100 K), and three cells are emplyed for each temperature. The simulations in the same temperature should be conducted with the same conditions. Hence, only one OUTCAR file exist in each subdirectory.

## 1. Vacancy trajectory determination
User can obtain vacancy trajectory using:
```ruby
 vachoppy -m t O 0.1 2100 03 # symbol, t_width, temperature, label
 ```

 Output:
<div align=center>
<p>
    <img src="./imgs/traj.gif" width="550"/>
</p>
</div>


## 2. Effective hopping parameter calculation
Use:
```ruby
vachoppy -m p O 0.1 # symbol, t_width
```
This command will provides effective hopping parameters of an oxygen vacancy.

## 3. Diffusion coefficient using Einstein relation
Use:
```ruby
vachoppy -m e O 50 --skip 1
```
This command will provides diffusion coefficient of an oxygen vacancy at each temperature.

## 4. Assessment of lattice stability
