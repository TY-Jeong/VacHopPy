# VacHopPy 
---
**VacHopPy** is a Python package to analyze trajectory of vacancy hopping in materials. This module was developed based on VASP 5.4.4 and Python 3.11.5.

## Features
* Tracking vacancy in *ab-inito* molecular dynamics (AIMD).
* Calculating mass transport parameters of vacancies.
* Analyzing structure changes during AIMD.

## Contents
* Getting started
* Step by step examples
  * Tracking vacancy
    * Making animation
    * Distribution of hopping path
  * Mass transsport parameter calculation
    * Diffusion coefficient 
    * Atomic vibration coefficient
  * Fingerprint of structure
  
## Getting started
*will be updated*

## Step by step examples
---
## Tracking vacancy
>  For the following example, the user can use files in `VacHopPy/example/1.LatticeHopping/data.zip`, which contains MD trajectory of monoclinic HfO<SUB>2</SUB> with one V<SUB>O</SUB><SUP>2+</SUP> at the temperature of 2200 K. 

The core module for vacancy trajectory analysis is `trajectory.LatticeHopping`.
```ruby
from vachoppy import trajectory

traj = trajectory.LatticeHopping(poscar_perf='data/POSCAR_novcac'
                                 xdatcar='data/XDATCAR_01'
                                 interval=50, 
                                 target='O')
```
The `poscar_perf` is the path to POSCAR  of perfect crystalline and is used to construct lattice points. The `xdatcar` is XDATCAR of defective cell containing one vacancy, obtaind from AIMD. To alleviate thermal fluctuation, the positions of each atoms are averaged at specified `interval`. Therefore, if the user ran AIMD for 15000 iterations and used an interval of 50, VacHopPy considers the trajectory of AIMD to be composed of 300 (=15000/50) steps and generate 300 averaged structures from each step. Note that the **position of each atoms are determined to the nearest lattice point** from the averaged position at each step. The `target` specifies which atom will be traced. For example, `target='O'` will trace oxgen ions and an oxygen vacancy.

The method of allocating the nearest lattice point might cause unexpected problem, such as multiple vacancies existing simultaneously (**multi-vacancy problem**) or multiple sequential hopping occuring within one step (**multi-hopping problem**). VacHopPy provides correction functions for the problems, will be addressed later.

> The `trajectory` module assumes that lattice points are maintained during MD simulation. It is well known that monoclinic lattice of HfO<SUB>2</SUB> becomes unstable, so that the MD trajectory at 2200 K is inappropriate. Nevertheless, the following examples were written using the MD trajetory at 2200 K in order to explain how to correct for the unexpected problems.


### Making animation
Using `trajectory.LatticeHopping.animation` method, the user can easily generate the animation of vacancy movements. By default arguments, the animation is saved in `traj.gif` and the snapshots are saved in the `./gif` directory.
```ruby
traj.animation(potim=2, # time interval in AIMD
               fps=30)
```
The user can customize the animation by specifying arguments desribed in ***link***. Below is an example animation of oxygen vacancy hopping in monoclinic HfO<SUB>2</SUB> at the temperature of 2200K.
<div align=center>
<p>
    <img src="./imgs/traj.gif" width="550" height="412" />
</p>
</div>

The **oxygen vacancy** is represented by a **yellow color**, while other colored points correspond to occupied lattice points by oxgen ions. The movements of oxygen ions are displayed with arrows of the same color.

Note that not only hoppings through the vacancy, but also site-exchanges between oxygein ions are also observed, indicating that the lattice hopping is not the only mechanism of oxygen transport.

For example, below is a snapshot at step 119. The oxygen ions labelled by 34(magenta), 35(green), 58(light green), 59(grey) exchanged their site, with hourglass-shape. (The labels can be displayed by setting `label=True`.)

<div align=center>
<p>
    <img src="./imgs/snapshot_119.png" width="550" height="412" />
</p>
</div>

with `trajectory.LatticeHopping.save_traj_on_lat` module, the user can focus on the interested atoms and steps. For example, the movements of ions within the hourglass will be shown by the following code.

```ruby
traj.save_traj_on_lat(lat_point=[34, 35, 58, 59],
                      step=[118, 119],
                      vac=True,
                      label=False)
```

This code will save the snapshots in `./traj_on_lat` directory. In the snapshots, the initial sites are represented by background colors, as shown below.
<div align=center>
<p>
    <img src="./imgs/traj_118.png" width="550" height="412" /> 
    <img src="./imgs/traj_119.png" width="550" height="412" />
</p>
</div>




### Correction for multi-vacancy problem
The method of allocating the oxygen (or other) ions to the nearest lattice points has a limitation that can result in an unphysical situation where two or more vacancies exist simultaneously, eventhough there was only one vacancy in the cell. This situation is usually observed at high temperature, especially when the existing crystalline structure becomes unstable. 

In this situation, a plausible vacancy site should be carefully identified. VacHopPy provide a function that find the plausible vacancy site based on the connectivity between the previous vacancy site. For example, 

### Trajectory of each atom


