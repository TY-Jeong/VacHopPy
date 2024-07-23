# VacHopPy 
---
**VacHopPy** is a Python package to analyze trajectory of vacancy hopping in crystalline materials. This module was developed based on VASP 5.4.4. 

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
>  Note: The user can use files in `VacHopPy/example/1.LatticeHopping/data.zip` for the following example.

The core module to get vacancy trajectory is `trajectory.LatticeHopping`.
```ruby
from vachoppy import trajectory

traj = trajectory.LatticeHopping(poscar_perf, # usee data/POSCAR_novac file
                                 xdatcar,     # use data/XDATCAR_01 file
                                 interval=50, 
                                 target='O' )
```
The `poscar_perf` is the path to POSCAR  of perfect crystalline and is used to construct lattice points. The `xdatcar` is XDATCAR of defective cell containing one vacancy, obtaind from AIMD. To alleviate thermal fluctuation, the positions of each atoms are averaged at specified `interval`. Therefore, if the user ran AIMD for 15000 iterations and used an interval of 50, VacHopPy considers the trajectory of AIMD to be composed of 300 (=15000/50) steps and generate 300 averaged structures from each step. Note that the **position of each atoms are determined to the nearest lattice point** from the averaged position at each step. The `target` specifies which atom will be traced. For example, `target = 'O'` will trace oxgen ions and an oxygen vacancy.

The determination based on averaged position might cause unintended situatoins, such as multiple vacancies existing simultaneously or multiple sequential hopping occuring within one step. To address the problems, VacHopPy provides several correction functions, which will be described latter.

### Making animation
Using `trajectory.LatticeHopping.animation` method, the user can easily generate the animation of vacancy movements. By default arguments, the animation is saved in `traj.gif` and the snapshots are saved in the `./gif` directory.
```ruby
traj.animation(potim=2, # time interval in AIMD
               fps=30)
```
The user can adjust the animation by specifying arguments desribed in ***link***. Below is an example animation of oxygen vacancy hopping in monoclinic HfO<SUB>2</SUB> at the temperature of 2200K.
<div align=center>
<p>
    <img src="./imgs/traj.gif" width="550" height="412" />
</p>
</div>

The oxygen vacancy is represented by a yellow color, while other colored points correspond to occupied lattice points by oxgen ions. The movements of oxygen ions are displayed with arrows of the same color.

It is well known that monoclinc HfO<SUB>2</SUB> becomes unstable at the higher temperatures over 2000K. Hence the movements of oxygen ions irrelevant to the vacancy are also observed 

### Trajectory of each atom


