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

To ease usage, **VacHopPy** provides a command-line interface (CLI). Belows are available CLI commands:

<div align=center>
<table style="border:1px solid black;margin-left:auto;margin-right:auto;">
    <tr>
        <th scope="col">Option 1</td>
        <th scope="col">Option 2</td>
        <th scope="col">Use</td>
        <th scope="col">Note</td>
    </tr>
    <tr>
        <td rowspan="5">-m<br>(main)</td>
        <td>p</td>
        <td>Calculate effective hopping parameters (excluding z and ν)</td>
        <td> </td>
    </tr>
    <tr>
        <!-- <td>2</td> -->
        <td>pp</td>
        <td>Calculate z and ν (post-processing for `-m p` option)</td>
        <td>neb.csv is needed</td>
    </tr>
    <tr>
        <!-- <td>3</td> -->
        <td>e</td>
        <td>Calculate diffusion coefficient using Einstein relation</td>
        <td> </td>
    </tr>
    <tr>
        <!-- <td>4</td> -->
        <td>t</td>
        <td>Make an animation for vacancy trajectories</td>
        <td> </td>
    </tr>
    <tr>
        <!-- <td>5</td> -->
        <td>f</td>
        <td>Perform fingerprint analyses</td>
        <td> </td>
    </tr>
    <tr>
        <td rowspan="6">-u<br>(util)</td>
        <td>extract_force</td>
        <td>Extract FORCE file from vasprun.xml</td>
        <td> </td>
    </tr>
    <tr>
        <!-- <td>2</td> -->
        <td>concat_xdatcar</td>
        <td>Concatenate two XDATCAR files</td>
        <td> </td>
    </tr>
    <tr>
        <!-- <td>3</td> -->
        <td>concat_force</td>
        <td>Concatenate two FORCE files</td>
        <td> </td>
    </tr>
    <tr>
        <!-- <td>4</td> -->
        <td>update_outcar</td>
        <td>Combine two OUTCAR files</td>
        <td> </td>
    </tr>
    <tr>
        <!-- <td>5</td> -->
        <td>fingerprint</td>
        <td>Extract fingerprint</td>
        <td> </td>
    </tr>
    <tr>
        <!-- <td>6</td> -->
        <td>cosine_distance</td>
        <td>Calculate cosine distance</td>
        <td> </td>
    </tr>
</table>
</div>


The user can see a brief explanation for each command using `-h` option:
```ruby
vachoppy -h # overall explanation
vachoppy -m p -h # specific explanation for -m p option
```



