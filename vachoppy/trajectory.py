"""
vachoppy.trajectory
===================

Contains the core classes that form the computational engine for vacancy
diffusion analysis in `VacHopPy`.

This module provides a hierarchy of classes designed to process trajectory
data in progressive stages, from the lowest level of identifying individual
atomic hops to the highest level of calculating ensemble-averaged physical
properties and performing Arrhenius fits.

While these classes can be used individually for advanced or custom workflows,
they are typically orchestrated by the high-level setup function in the
`vachoppy.core` module.

Core Classes
------------
- **Trajectory**: The foundational class that reads a single HDF5 trajectory,
  determines the site occupation of each atom at each time step, and
  reconstructs the raw vacancy movement paths.
- **TrajectoryAnalyzer**: Analyzes the output of a `Trajectory` object to
  identify and quantify discrete hopping events, calculating statistics like
  hop counts and site residence times.
- **Encounter**: Processes results from `TrajectoryAnalyzer` to identify
  atom-vacancy encounters and calculate the correlation factor (f).
- **CalculatorSingle**: A high-level orchestrator that combines the three
  classes above to run the full pipeline on a *single* trajectory file. It
  serves as the worker process for ensemble calculations.
- **CalculatorEnsemble**: The top-level analysis class that manages an ensemble
  of trajectories. It uses `TrajectoryBundle` for file discovery and runs
  multiple `CalculatorSingle` instances in parallel to compute
  temperature-dependent properties.
- **TrajectoryBundle**: A helper class used by `CalculatorEnsemble` to
  discover, validate, and group a collection of trajectory files for
  multi-temperature analysis.

Typical Usage
-------------
Users should typically **not** import classes from this module directly.
Instead, the intended entry point is the `vachoppy.core.Calculator`
function, which handles the instantiation and execution of the appropriate
classes from this module behind the scenes.
"""

from __future__ import annotations

__all__ =['Trajectory', 'TrajectoryAnalyzer', 'Encounter', 'TrajectoryBundle', 'CalculatorSingle', 'CalculatorEnsemble']

import os
import h5py
import json
import math
import copy   
import itertools
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

from PIL import Image
from typing import Tuple, Optional
from pathlib import Path
from tqdm.auto import tqdm
from tabulate import tabulate
from joblib import Parallel, delayed
from scipy.optimize import minimize_scalar

from collections import defaultdict
from itertools import permutations

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

from vachoppy.frequency import AttemptFrequency
from vachoppy.utils import monitor_performance


class Arrow3D(FancyArrowPatch):
    """A custom 3D arrow patch for matplotlib."""
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)
    
    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)


class Trajectory:
    """Analyzes a single MD trajectory to trace atomic and vacancy hops.

    This class reads a trajectory from an HDF5 file, identifies which lattice
    site each atom occupies at each time step, validates hops using a robust
    transition state criterion, and reconstructs the movement paths of vacancies.
    It serves as the primary engine for analyzing vacancy diffusion events and
    provides methods for both interactive (Plotly) and static (matplotlib)
    visualization of the results.

    Args:
        path_traj (str):
            Path to the HDF5 trajectory file.
        site (Site):
            An initialized `Site` object containing pre-analyzed information
            about the crystal lattice sites.
        t_interval (float):
            The time interval in picoseconds (ps) between each analysis step.
            This determines the coarse-graining level of the analysis.
        force_margin (float, optional):
            The minimum force magnitude (eV/Å) required for the directional
            transition state check. Defaults to 0.05.
        cos_margin (float, optional):
            The required cosine difference to reject a hop during the transition
            state check. A higher value makes the check stricter. Defaults to 0.1.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Attributes:
        occupation (numpy.ndarray):
            A 2D array of shape (num_atoms, num_steps) storing the lattice site
            index occupied by each atom at each step.
        hopping_history (list[dict]):
            A chronological list of all validated hop events. Each entry contains
            details such as the time step, hopping atom index, and the initial
            and final site indices.
        vacancy_trajectory_index (dict):
            A dictionary mapping each time step to the indices of occupied
            vacancy sites.
        vacancy_trajectory_coord_cart (dict):
            A dictionary mapping each time step to the Cartesian coordinates (Å)
            of occupied vacancy sites within the unit cell.
        unwrapped_vacancy_trajectory_coord_cart (dict):
            A dictionary mapping each time step to the continuous, unwrapped
            Cartesian coordinates of vacancies, representing the true diffusion path.
        hopping_sequence (dict):
            A dictionary mapping each time step to the reconstructed vacancy hop
            paths between sites.

    Raises:
        FileNotFoundError: If the `path_traj` file is not found.
        ValueError: If the trajectory file has an invalid format, missing data,
            or if parameters (`t_interval`, `symbol`, etc.) are inconsistent.
        IOError: If the HDF5 file cannot be read.
    """

    def __init__(self,
                 path_traj: str,
                 site: Site,
                 t_interval: float,
                 force_margin: float = 0.05,
                 cos_margin: float = 0.1,
                 verbose: bool = True):
        
        self._validate_traj(path_traj)
        self.cmap = self._custrom_cmap()
        
        self.path_traj = path_traj
        self.site = site
        self.force_margin = force_margin
        self.cos_margin = cos_margin
        self.t_interval = t_interval
        self.verbose = verbose
        
        # Read cond
        self.total_frames = None    # nsw
        self.temperature = None     
        self.dt = None         
        self.symbol = None
        self.num_frames = None  
        self.num_steps = None   
        self.num_atoms = None  
        self.frame_interval = None 
        self.lattice_parameter = None
        self._read_cond()
        
        # Read site
        self.lattice_sites = None   
        self.lattice_sites_cart = None  
        self.num_lattice_sites = None   
        self.num_vacancies = None      
        self._read_site()
        
        self.occupation = None
        self._get_occupation(self.force_margin, self.cos_margin) 
        
        self.trace_arrows = None
        self._get_trace_arrows() 
        
        # Vacancy trajectory
        self.hopping_sequence = {}
        self.vacancy_trajectory_index = {}
        self.vacancy_trajectory_coord_cart = {}
        self.transient_vacancy = {}
        self._get_vacancy_trajectory()
        
        # Unwrapped trajectory
        self.unwrapped_vacancy_trajectory_coord_cart = None
        self._get_unwrapped_vacancy_trajectory()
        
    @monitor_performance
    def animate_vacancy_trajectory(self,
                                   vacancy_indices: int | list,
                                   step_init: int = 0,
                                   step_final: int | None = None,
                                   unwrap: bool = False,
                                   max_trace_length: int = 10,
                                   update_alpha: float = 0.8,
                                   fps: int = 5,
                                   save: bool = True,
                                   filename: str = "trajectory_video.html",
                                   verbose: bool = True) -> None:
        """Generates an interactive 3D animation of vacancy trajectories using Plotly.

        The animation includes play/pause controls and a slider to navigate through
        time. The real-time position of each vacancy is shown, and its recent path
        fades over time to indicate the direction of movement.

        Args:
            vacancy_indices (int | list):
                Index or list of indices for the vacancies to be animated.
            step_init (int, optional):
                The starting step for the animation. Defaults to 0.
            step_final (int | None, optional):
                The ending step for the animation. If None, animates to the end.
                Defaults to None.
            unwrap (bool, optional):
                If True, plots the continuous, unwrapped path. Defaults to False.
            max_trace_length (int, optional):
                The number of past path segments to display as a fading tail.
                Defaults to 10.
            update_alpha (float, optional):
                Fading rate for the tail. Closer to 0 makes paths fade faster.
                Defaults to 0.8.
            fps (int, optional):
                Frames per second for animation playback. Defaults to 5.
            save (bool, optional):
                If True, saves the animation as a standalone HTML file.
                Defaults to True.
            filename (str, optional):
                The name of the output HTML file if `save` is True.
                Defaults to "trajectory_video.html".
            verbose (bool, optional):
                Verbosity flag. Defaults to True.

        Returns:
            None: This method saves an HTML file or displays a plot.

        Raises:
            ValueError: If the specified `step_init` or `step_final` is out of bounds.
        """
        if unwrap:
            coord_source = self.unwrapped_vacancy_trajectory_coord_cart
            title_prefix = "Unwrapped Vacancy Animation"
        else:
            coord_source = self.vacancy_trajectory_coord_cart
            title_prefix = "Vacancy Animation"

        if not coord_source:
            print("Vacancy trajectory data is not available."); return

        if step_final is None: step_final = self.num_steps - 1
        if not (0 <= step_init <= step_final < self.num_steps):
            raise ValueError(f"Invalid step range [{step_init}, {step_final}].")
        
        if isinstance(vacancy_indices, int):
            indices_to_plot = [vacancy_indices]
        else:
            indices_to_plot = vacancy_indices
        
        available_steps = sorted([s for s in coord_source.keys() 
                                  if s is not None and step_init <= s <= step_final])

        all_paths_coords = []
        for vac_idx in indices_to_plot:
            path = [
                coord_source[s][vac_idx] 
                for s in available_steps if vac_idx < len(coord_source.get(s, []))
            ]
            if path: all_paths_coords.append(np.array(path, dtype=np.float64))

        static_traces = []
        if unwrap and all_paths_coords:
            inv_lattice = np.linalg.inv(self.lattice_parameter)
            frac_coords = np.dot(np.vstack(all_paths_coords), inv_lattice)
            cell_indices = np.floor(frac_coords).astype(int)
            unique_cells = np.unique(cell_indices, axis=0)
            a, b, c = self.lattice_parameter
            for cell_vec in unique_cells:
                i, j, k = cell_vec
                origin = i * a + j * b + k * c
                is_center_cell = (i == 0 and j == 0 and k == 0)
                if not is_center_cell:
                    v = np.array([
                        origin, origin+a, origin+b, origin+c, 
                        origin+a+b, origin+b+c, origin+c+a, origin+a+b+c
                    ])
                    edges = [(0,1), (0,2), (0,3), (1,4), (1,6), (2,4), 
                             (2,5), (3,5), (3,6), (4,7), (5,7), (6,7)]
                    x_e, y_e, z_e = [], [], []; [ (x_e.extend([v[s][0], v[e][0], None]), 
                                                   y_e.extend([v[s][1], v[e][1], None]), 
                                                   z_e.extend([v[s][2], v[e][2], None])) 
                                                 for s, e in edges ]
                    static_traces.append(go.Scatter3d(x=x_e, 
                                                      y=y_e, 
                                                      z=z_e, 
                                                      mode='lines', 
                                                      line=dict(color='lightgrey', width=1), 
                                                      showlegend=False)
                                         )
                    supercell_sites = self.lattice_sites_cart + origin
                    static_traces.append(go.Scatter3d(x=supercell_sites[:, 0], 
                                                      y=supercell_sites[:, 1], 
                                                      z=supercell_sites[:, 2], 
                                                      mode='markers', 
                                                      marker=dict(color='grey', size=3, opacity=0.2), 
                                                      showlegend=False)
                                         )
            v = np.array(
                [np.zeros(3)] + list(self.lattice_parameter) + 
                [self.lattice_parameter[0] + self.lattice_parameter[1]] + 
                [self.lattice_parameter[1] + self.lattice_parameter[2]] + 
                [self.lattice_parameter[2] + self.lattice_parameter[0]] + 
                [np.sum(self.lattice_parameter, axis=0)]
            )
            edges = [(0,1), (0,2), (0,3), (1,4), (1,6), (2,4), 
                     (2,5), (3,5), (3,6), (4,7), (5,7), (6,7)]
            x_e, y_e, z_e = [], [], []; [ (x_e.extend([v[s][0], v[e][0], None]), 
                                           y_e.extend([v[s][1], v[e][1], None]), 
                                           z_e.extend([v[s][2], v[e][2], None])) 
                                         for s, e in edges ]
            static_traces.append(go.Scatter3d(x=x_e, 
                                              y=y_e, 
                                              z=z_e, 
                                              mode='lines', 
                                              line=dict(color='black', width=2.5), 
                                              showlegend=False)
                                 )
            static_traces.append(go.Scatter3d(x=self.lattice_sites_cart[:, 0], 
                                              y=self.lattice_sites_cart[:, 1], 
                                              z=self.lattice_sites_cart[:, 2], 
                                              mode='markers', 
                                              marker=dict(color='dimgrey', size=4, opacity=0.8), 
                                              showlegend=False)
                                 )
        else:
            v = np.array(
                [np.zeros(3)] + list(self.lattice_parameter) + 
                [self.lattice_parameter[0] + self.lattice_parameter[1]] + 
                [self.lattice_parameter[1] + self.lattice_parameter[2]] + 
                [self.lattice_parameter[2] + self.lattice_parameter[0]] + 
                [np.sum(self.lattice_parameter, axis=0)]
            )
            edges = [(0,1), (0,2), (0,3), (1,4), (1,6), (2,4), 
                     (2,5), (3,5), (3,6), (4,7), (5,7), (6,7)]
            x_e, y_e, z_e = [], [], []; [ (x_e.extend([v[s][0], v[e][0], None]), 
                                           y_e.extend([v[s][1], v[e][1], None]), 
                                           z_e.extend([v[s][2], v[e][2], None])) 
                                         for s, e in edges ]
            static_traces.append(go.Scatter3d(x=x_e, 
                                              y=y_e, 
                                              z=z_e, 
                                              mode='lines', 
                                              line=dict(color='black', width=2), 
                                              showlegend=False)
                                 )
            static_traces.append(go.Scatter3d(x=self.lattice_sites_cart[:, 0], 
                                              y=self.lattice_sites_cart[:, 1], 
                                              z=self.lattice_sites_cart[:, 2], 
                                              mode='markers', 
                                              marker=dict(color='dimgrey', size=4, opacity=0.8), 
                                              showlegend=False)
                                 )

        frames = []
        colors = ['#636EFA', '#EF553B', '#00CC96', 
                  '#AB63FA', '#FFA15A', '#19D3F3']
        
        iterable = enumerate(available_steps)
        iterable = tqdm(
            iterable,
            total=len(available_steps),
            desc="Make Animation",
            ascii=True,
            bar_format='{l_bar}{bar:30}{r_bar}')

        for s_idx, step in iterable:
            dynamic_traces = []
            for i, path_coords in enumerate(all_paths_coords):
                vac_idx = indices_to_plot[i]
                solid_color = colors[i % len(colors)]
                
                dynamic_traces.append(go.Scatter3d(
                    x=[path_coords[s_idx, 0]], y=[path_coords[s_idx, 1]], z=[path_coords[s_idx, 2]],
                    mode='markers',
                    marker=dict(symbol='circle', 
                                size=5, # Adjust this value to change vacancy size
                                color=solid_color, 
                                line=dict(width=2, color='black')),
                    name=f'Vacancy {vac_idx}',
                    legendgroup=f'vacancy_{vac_idx}'
                ))
                
                current_alpha = 1.0
                trace_length = min(s_idx, max_trace_length)
                for j in range(max_trace_length):
                    if j < trace_length:
                        p_start, p_end = path_coords[s_idx - 1 - j], path_coords[s_idx - j]
                        h = solid_color.lstrip('#'); r, g, b = tuple(int(h[k:k+2], 16) for k in (0, 2, 4))
                        rgba_color = f'rgba({r}, {g}, {b}, {current_alpha})'
                        dynamic_traces.append(go.Scatter3d(
                            x=[p_start[0], p_end[0]], 
                            y=[p_start[1], p_end[1]], 
                            z=[p_start[2], p_end[2]],
                            mode='lines', 
                            line=dict(color=rgba_color, width=8), 
                            showlegend=False
                        ))
                        current_alpha *= update_alpha
                    else:
                        dynamic_traces.append(go.Scatter3d(x=[None], 
                                                           y=[None], 
                                                           z=[None], 
                                                           mode='lines', 
                                                           line=dict(color='rgba(0,0,0,0)'), 
                                                           showlegend=False)
                                              )
            
            frames.append(go.Frame(data=static_traces + dynamic_traces, name=str(step)))

        fig = go.Figure(data=frames[0].data if frames else static_traces)
        fig.update(frames=frames)
        
        def frame_args(duration):
            return {"frame": {"duration": duration}, "mode": "immediate", "transition": {"duration": 0}}
        fig.update_layout(
            updatemenus=[{
                "buttons": [
                    {"args": [None, frame_args(1000 / fps)], 
                     "label": "▶ Play", 
                     "method": "animate"},
                    {"args": [[None], frame_args(0)], 
                     "label": "❚❚ Pause", 
                     "method": "animate"},
                ], "direction": "left", 
                "pad": {"r": 10, "t": 70}, 
                "type": "buttons", 
                "x": 0.1, 
                "xanchor": "right", 
                "y": 0, 
                "yanchor": "top"
            }],
            sliders=[{
                "active": 0, 
                "yanchor": "top", 
                "xanchor": "left",
                "currentvalue": {"font": {"size": 16}, 
                                 "prefix": "Step: ", 
                                 "visible": True, 
                                 "xanchor": "right"},
                "transition": {"duration": 0}, 
                "pad": {"b": 10, "t": 50}, 
                "len": 0.9, 
                "x": 0.1, 
                "y": 0,
                "steps": [
                    {"args": [[f.name], frame_args(0)], 
                     "label": f.name, "method": "animate"} for f in fig.frames
                    ]
            }]
        )

        fig.update_layout(
            title_text=f'{title_prefix} (Steps {step_init}-{step_final})',
            scene=dict(xaxis_title='x (Å)', 
                       yaxis_title='y (Å)', 
                       zaxis_title='z (Å)', 
                       aspectmode='data'),
            showlegend=True, margin=dict(l=0, r=0, b=0, t=40)
        )

        if save:
            fig.write_html(filename)
            if verbose: print(f"\n'{filename}' created.\n")
        else:
            fig.show()

    def plot_vacancy_trajectory(self, 
                                vacancy_indices: int | list, 
                                step_init: int = 0, 
                                step_final: int | None = None,
                                unwrap: bool = False,
                                alpha: float = 0.6,
                                disp: bool = True,
                                save: bool = False,
                                filename: str = "trajectory.html") -> None:
        """Generates an interactive 3D plot of vacancy trajectories using Plotly.

        This method creates a static but interactive (rotatable, zoomable) 3D plot
        of the complete vacancy paths over a specified time range.

        Args:
            vacancy_indices (int | list):
                Index or list of indices for the vacancies to be plotted.
            step_init (int, optional):
                The starting step for the trajectory plot. Defaults to 0.
            step_final (int | None, optional):
                The ending step for the plot. If None, plots to the end.
                Defaults to None.
            unwrap (bool, optional):
                If True, plots the continuous, unwrapped path. Defaults to False.
            alpha (float, optional):
                Opacity for trajectory lines. Defaults to 0.6.
            disp (bool, optional):
                If True, displays the plot interactively. Defaults to True.
            save (bool, optional):
                If True, saves the plot as a standalone HTML file. Defaults to False.
            filename (str, optional):
                The name of the output HTML file if `save` is True.
                Defaults to "trajectory.html".
        
        Returns:
            None: This method saves an HTML file or displays a plot.

        Raises:
            ValueError: If the specified `step_init` or `step_final` is out of bounds.
        """
        if unwrap:
            if self.unwrapped_vacancy_trajectory_coord_cart is None:
                self._get_unwrapped_vacancy_trajectory()
            coord_source = self.unwrapped_vacancy_trajectory_coord_cart
            title_prefix = "Unwrapped Vacancy Trajectory"
        else:
            coord_source = self.vacancy_trajectory_coord_cart
            title_prefix = "Vacancy Trajectory"

        if not coord_source:
            print("Vacancy trajectory data is not available.")
            return

        if step_final is None:
            step_final = self.num_steps - 1
        
        if not (0 <= step_init <= step_final < self.num_steps):
            raise ValueError(f"Invalid step range [{step_init}, {step_final}].")
            
        fig = go.Figure()

        if isinstance(vacancy_indices, int):
            indices_to_plot = [vacancy_indices]
        else:
            indices_to_plot = vacancy_indices
        
        available_steps = sorted([s for s in coord_source.keys() 
                                  if s is not None and step_init <= s <= step_final])

        all_paths_coords = []
        for vac_idx in indices_to_plot:
            path_coords = [coord_source[s][vac_idx] 
                           for s in available_steps if vac_idx < len(coord_source[s])]
            if path_coords:
                all_paths_coords.append(np.array(path_coords, dtype=np.float64))

        if unwrap and all_paths_coords:
            self._plot_supercell(fig, np.vstack(all_paths_coords))
            v = np.array(
                [np.zeros(3)] + list(self.lattice_parameter) + 
                [self.lattice_parameter[0] + self.lattice_parameter[1]] + 
                [self.lattice_parameter[1] + self.lattice_parameter[2]] + 
                [self.lattice_parameter[2] + self.lattice_parameter[0]] + 
                [np.sum(self.lattice_parameter, axis=0)]
            )
            edges = [(0,1), (0,2), (0,3), (1,4), (1,6), (2,4), 
                     (2,5), (3,5), (3,6), (4,7), (5,7), (6,7)]
            x_e, y_e, z_e = [], [], []; [ (x_e.extend([v[s][0], v[e][0], None]), 
                                           y_e.extend([v[s][1], v[e][1], None]), 
                                           z_e.extend([v[s][2], v[e][2], None])) 
                                         for s, e in edges ]
            fig.add_trace(go.Scatter3d(x=x_e, 
                                       y=y_e, 
                                       z=z_e, 
                                       mode='lines', 
                                       line=dict(color='black', width=2.5), 
                                       name='Center Cell')
                          )
        else:
            v = np.array(
                [np.zeros(3)] + list(self.lattice_parameter) + 
                [self.lattice_parameter[0] + self.lattice_parameter[1]] + 
                [self.lattice_parameter[1] + self.lattice_parameter[2]] + 
                [self.lattice_parameter[2] + self.lattice_parameter[0]] + 
                [np.sum(self.lattice_parameter, axis=0)]
            )
            edges = [(0,1), (0,2), (0,3), (1,4), (1,6), (2,4), 
                     (2,5), (3,5), (3,6), (4,7), (5,7), (6,7)]
            x_e, y_e, z_e = [], [], []; [ (x_e.extend([v[s][0], v[e][0], None]), 
                                           y_e.extend([v[s][1], v[e][1], None]), 
                                           z_e.extend([v[s][2], v[e][2], None])) 
                                         for s, e in edges ]
            fig.add_trace(go.Scatter3d(x=x_e, 
                                       y=y_e, 
                                       z=z_e, 
                                       mode='lines', 
                                       line=dict(color='black', width=2), 
                                       name='Unit Cell')
                          )
        
        fig.add_trace(go.Scatter3d(x=self.lattice_sites_cart[:, 0], 
                                   y=self.lattice_sites_cart[:, 1], 
                                   z=self.lattice_sites_cart[:, 2], 
                                   mode='markers', 
                                   marker=dict(color='dimgrey', size=4, opacity=0.8), 
                                   name='Lattice Sites')
                      )

        colors = ['#636EFA', '#EF553B', '#00CC96', 
                  '#AB63FA', '#FFA15A', '#19D3F3']
        start_marker, end_marker = 'diamond', 'x'

        for i, path_coords in enumerate(all_paths_coords):
            vac_idx = indices_to_plot[i]
            
            solid_color = colors[i % len(colors)]
            
            h = solid_color.lstrip('#')
            r, g, b = tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
            rgba_color = f'rgba({r}, {g}, {b}, {alpha})'
            
            fig.add_trace(go.Scatter3d(
                x=path_coords[:, 0], y=path_coords[:, 1], z=path_coords[:, 2],
                mode='lines',
                line=dict(color=rgba_color, width=8),
                name=f'Vacancy {vac_idx} Path'
            ))
            
            fig.add_trace(go.Scatter3d(
                x=[path_coords[0, 0]], y=[path_coords[0, 1]], z=[path_coords[0, 2]],
                mode='markers',
                marker=dict(symbol=start_marker, 
                            size=5, 
                            color=solid_color, 
                            line=dict(width=2, color='black')),
                name=f'Vac {vac_idx} Start'
            ))
            fig.add_trace(go.Scatter3d(
                x=[path_coords[-1, 0]], y=[path_coords[-1, 1]], z=[path_coords[-1, 2]],
                mode='markers',
                marker=dict(symbol=end_marker, 
                            size=5, 
                            color=solid_color, 
                            line=dict(width=2, color='black')),
                name=f'Vac {vac_idx} End'
            ))

        fig.update_layout(
            title_text=f'{title_prefix} (Steps {step_init}-{step_final})',
            scene=dict(xaxis_title='x (Å)', 
                       yaxis_title='y (Å)', 
                       zaxis_title='z (Å)', 
                       aspectmode='data'),
            showlegend=True,
            margin=dict(l=0, r=0, b=0, t=40)
        )
        
        if save:
            fig.write_html(filename)
            if self.verbose: print(f"'{filename}' created.")
                
        if disp: fig.show()
    
        # TrajAnalyzer 클래스 내부에 추가될 헬퍼(worker) 메서드
    
    def _create_snapshot(self, 
                         step: int, 
                         foldername: str, 
                         dpi: int, 
                         label: bool, 
                         vac: bool, 
                         update_alpha: float, 
                         legend: bool, 
                         atom_indices: np.ndarray) -> str:
        """Creates and saves a single snapshot image for a given time step."""
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

        self._plot_lattice(ax, label=label)

        for i, atom_idx in enumerate(atom_indices):
            site_idx = self.occupation[atom_idx, step]
            coords = self.lattice_sites_cart[site_idx]
            ax.scatter(*coords, s=100, facecolor=self.cmap[atom_idx % len(self.cmap)],
                    edgecolor='k', alpha=0.8, label=f"Atom {atom_idx}")
        
        alpha = 1.0
        for i in reversed(range(step + 1)):
            if alpha < 0.1: break
            for arrow in self.trace_arrows.get(i, []):
                arrow_prop = dict(mutation_scale=15, arrowstyle='->',
                                color=arrow['c'], alpha=alpha,
                                shrinkA=0, shrinkB=0, lw=1.5)
                disp_arrow = Arrow3D(*arrow['p'].T, **arrow_prop)
                ax.add_artist(disp_arrow)
            alpha *= update_alpha

        if vac:
            vac_coords = self.vacancy_trajectory_coord_cart.get(step, [])
            if len(vac_coords) > 0:
                ax.scatter(*vac_coords.T, s=120, facecolor='yellow', edgecolor='k',
                        marker='o', alpha=0.8, zorder=10, label='Vacancy')
            trans_vac_indices = self.transient_vacancy.get(step, [])
            if len(trans_vac_indices) > 0:
                trans_vac_coords = self.lattice_sites_cart[trans_vac_indices]
                ax.scatter(*trans_vac_coords.T, s=120, facecolor='orange', edgecolor='k',
                        marker='s', alpha=0.8, zorder=10, label='Transient Vacancy')

        time = step * self.t_interval
        time_tot = (self.num_steps - 1) * self.t_interval
        ax.set_title(f"Time: {time:.2f} ps / {time_tot:.2f} ps (Step: {step}/{self.num_steps-1})")

        if legend:
            ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1.0))

        snapshot_path = os.path.join(foldername, f"snapshot_{step:04d}.png")
        plt.savefig(snapshot_path, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        
        return snapshot_path

    @ monitor_performance
    def animate_occupation(self,
                           index: list | str = 'all',
                           step_init: int = 0,
                           step_final: int | None = None,
                           vac: bool = True,
                           gif: bool = True,
                           filename: str = 'occupation_video.gif',
                           foldername: str = 'snapshots',
                           update_alpha: float = 0.9,
                           fps: int = 10,
                           loop: int = 0,
                           dpi: int = 150,
                           n_jobs: int = -1,
                           legend: bool = False,
                           label: bool = False,
                           verbose: bool = True) -> None:
        """Generates a GIF animation visualizing the time evolution of site occupations.

        This method creates a series of 3D snapshot images for each time step and
        compiles them into a GIF. Each frame displays the crystal lattice and shows
        the position of color-coded atoms based on which site they occupy at that
        moment in time.

        This provides a dynamic, visual representation of how individual atoms hop
        between lattice sites and how vacancies consequently move through the
        crystal structure.

        Args:
            index (list | str, optional):
                Indices of atoms to display. Defaults to 'all'.
            step_init (int, optional):
                The starting step for the animation. Defaults to 0.
            step_final (int | None, optional):
                The ending step for the animation. If None, animates to the end.
                Defaults to None.
            vac (bool, optional):
                If True, highlights vacancy and transient vacancy sites.
                Defaults to True.
            gif (bool, optional):
                If True, creates a GIF file from the snapshots. Defaults to True.
            filename (str, optional):
                Output GIF file name. Defaults to 'occupation_video.gif'.
            foldername (str, optional):
                Directory to save temporary snapshot images. Defaults to 'snapshots'.
            update_alpha (float, optional):
                Fading rate for atomic hop trace arrows. Defaults to 0.9.
            fps (int, optional):
                Frames per second for the output GIF. Defaults to 10.
            loop (int, optional):
                Number of loops for the GIF (0 for infinite). Defaults to 0.
            dpi (int, optional):
                Resolution in dots per inch for snapshots. Defaults to 150.
            n_jobs (int, optional):
                The number of CPU cores for parallel processing. -1 uses all
                available cores. Defaults to -1.
            legend (bool, optional):
                If True, displays a legend for atoms. Defaults to False.
            label (bool, optional):
                If True, displays numeric labels for lattice sites. Defaults to False.

        Returns:
            None: This method saves image files and a GIF to disk.
        """
        if not os.path.isdir(foldername):
            os.mkdir(foldername)

        if step_final is None: step_final = self.num_steps - 1
        if not (0 <= step_init <= step_final < self.num_steps):
            raise ValueError(f"Invalid step range [{step_init}, {step_final}].")
   
        step_indices = np.arange(step_init, step_final)
        atom_indices = np.arange(self.num_atoms) if index == 'all' else np.array(index)
        
        files = Parallel(n_jobs=n_jobs)(
            delayed(self._create_snapshot)(
                step, foldername, dpi, label, vac, update_alpha, legend, atom_indices
            ) for step in tqdm(step_indices, 
                               desc="Make Animation",
                               bar_format='{l_bar}{bar:30}{r_bar}',
                               ascii=True))
        files.sort()

        if gif:
            print(f"\nMerging {len(files)} snapshots into a GIF...")
            imgs = [Image.open(file) for file in files]
            if imgs:
                imgs[0].save(fp=filename, format='GIF', append_images=imgs[1:],
                            save_all=True, duration=int(1000 / fps), loop=loop)
                print(f"Successfully created '{filename}'.\n")
    
    def distance_PBC(self, 
                     coord1: np.ndarray, 
                     coord2: np.ndarray) -> float | np.ndarray:
        """Calculates the PBC-aware distance between fractional coordinates.

        This method computes the shortest distance between one or more initial
        points (`coord1`) and a final point (`coord2`), respecting the periodic
        boundary conditions of the lattice.

        Args:
            coord1 (np.ndarray):
                Initial coordinate(s) in fractional form. Can be a single point
                (1D array) or multiple points (2D array).
            coord2 (np.ndarray):
                Final coordinate in fractional form (1D array).

        Returns:
            float | np.ndarray:
                The calculated distance(s) in Cartesian units. Returns a float
                if `coord1` is 1D, or a 1D array of distances if `coord1` is 2D.
        """
        dist_frac = coord1 - coord2
        dist_frac -= np.round(dist_frac) 
        dist_cart = np.dot(dist_frac, self.lattice_parameter)
        return np.linalg.norm(dist_cart, axis=-1)     

    def _custrom_cmap(self):
        """Color map for visualization"""
        cmap = [
            'blue',
            'red',
            'teal',
            'indigo',
            'lime',
            'darkgoldenrod',
            'cyan',
            'hotpink',
            'dodgerblue',
            'dimgray',
            'forestgreen',
            'slateblue'
        ]
        return cmap      
    
    def _plot_lattice(self, ax, label=False) -> None:
        """Helper method to plot the unit cell and lattice sites on a 3D axis."""
        coord_origin = np.zeros(3)
        
        def plot_edge(start, end):
            edge = np.vstack([start, end]).T
            ax.plot(edge[0], edge[1], edge[2], 'k-', marker='none', lw=0.5)

        a, b, c = self.lattice_parameter
        # Unit cell vertices
        v = [coord_origin, a, b, c, a+b, b+c, c+a, a+b+c]
        
        # Edges of the unit cell
        edges = [
            (v[0], v[1]), (v[0], v[2]), (v[0], v[3]), (v[1], v[4]),
            (v[1], v[6]), (v[2], v[4]), (v[2], v[5]), (v[3], v[5]),
            (v[3], v[6]), (v[4], v[7]), (v[5], v[7]), (v[6], v[7])
        ]

        for start, end in edges:
            plot_edge(start, end)

        ax.scatter(*self.lattice_sites_cart.T, s=20, facecolor='none', edgecolors='k', alpha=0.5)
        
        if label:
            for i, coord in enumerate(self.lattice_sites_cart):
                ax.text(*coord, s=f"{i}", fontsize='xx-small', ha='center', va='center')

        ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)'); ax.set_zlabel('z (Å)')
        ax.set_aspect('equal', adjustable='box')

    def _plot_supercell(self, fig, all_path_coords):
        """Plots only the visited supercells and their sites for context."""
        if all_path_coords.size == 0: return

        inv_lattice = np.linalg.inv(self.lattice_parameter)
        
        frac_coords = np.dot(all_path_coords, inv_lattice)
        cell_indices = np.floor(frac_coords).astype(int)
        
        unique_cells = np.unique(cell_indices, axis=0)
        
        a, b, c = self.lattice_parameter
        
        for cell_vec in unique_cells:
            i, j, k = cell_vec
            origin = i * a + j * b + k * c

            if i == 0 and j == 0 and k == 0:
                continue
            
            v = np.array([origin, origin+a, origin+b, origin+c, origin+a+b, 
                          origin+b+c, origin+c+a, origin+a+b+c])
            edges = [(0,1), (0,2), (0,3), (1,4), (1,6), (2,4), (2,5), 
                     (3,5), (3,6), (4,7), (5,7), (6,7)]
            
            x_edges, y_edges, z_edges = [], [], []
            for s, e in edges:
                x_edges.extend([v[s][0], v[e][0], None])
                y_edges.extend([v[s][1], v[e][1], None])
                z_edges.extend([v[s][2], v[e][2], None])
            
            fig.add_trace(go.Scatter3d(
                x=x_edges, y=y_edges, z=z_edges, mode='lines',
                line=dict(color='lightgrey', width=1), showlegend=False
            ))
            
            supercell_sites = self.lattice_sites_cart + origin
            fig.add_trace(go.Scatter3d(
                x=supercell_sites[:, 0], y=supercell_sites[:, 1], z=supercell_sites[:, 2],
                mode='markers', marker=dict(color='grey', size=3, opacity=0.2),
                showlegend=False
            ))
    
    def _validate_traj(self, path_traj: str) -> None:
        """
        Validates the structure and content of the HDF5 trajectory file.

        This method checks for the correct file extension, file existence, and
        the presence of required datasets ('positions', 'forces') and metadata
        attributes ('symbol', 'nsw', 'dt', 'temperature', 'atom_counts', 'lattice').

        Args:
            path_traj (str): The file path to validate.

        Raises:
            ValueError: If the file extension is not '.h5' or if required
                metadata or datasets are missing.
            FileNotFoundError: If the trajectory file does not exist.
            IOError: If the file cannot be read as an HDF5 file.
        """
        if not path_traj.endswith('.h5'):
            raise ValueError(f"Error: Trajectory file must have a .h5 extension, but got '{path_traj}'.")

        if not os.path.isfile(path_traj):
            raise FileNotFoundError(f"Error: Input file '{path_traj}' not found.")
        
        try:
            with h5py.File(path_traj, "r") as f:
                required_datasets = ["positions", "forces"]
                for dataset in required_datasets:
                    if dataset not in f:
                        raise ValueError(f"Error: Required dataset '{dataset}' not found in '{path_traj}'.")

                metadata_str = f.attrs.get("metadata")
                if not metadata_str:
                    raise ValueError(f"Error: Required attribute 'metadata' not found in '{path_traj}'.")
                
                cond = json.loads(metadata_str)
                required_keys = ["symbol", "nsw", "dt", "temperature", "atom_counts", "lattice"]
                for key in required_keys:
                    if key not in cond:
                        raise ValueError(f"Error: Required key '{key}' not found in metadata of '{path_traj}'.")

        except (IOError, OSError) as e:
            raise IOError(f"Error: Failed to read '{path_traj}' as an HDF5 file. Reason: {e}")
          
    def _read_cond(self) -> None:
        """Reads simulation conditions from the HDF5 file's metadata.

        This internal method extracts simulation parameters from the trajectory
        file, validates them against the provided `site` object, calculates
        derived values like the frame interval, and populates the instance
        attributes (e.g., self.dt, self.temperature).

        Raises:
            ValueError: If the symbol or lattice parameters in the trajectory
                file do not match the `site` object, or if `t_interval` is not
                a valid multiple of the simulation timestep `dt`.
            ZeroDivisionError: If the timestep `dt` is zero.
        """
        eps = self.site.eps
        with h5py.File(self.path_traj, "r") as f:
            cond = json.loads(f.attrs.get("metadata"))

            self.total_frames = cond.get("nsw")
            self.temperature = cond.get("temperature")
            self.dt = cond.get("dt")
            self.symbol = cond.get("symbol")
            self.num_atoms = cond.get("atom_counts")[self.symbol]
            self.lattice_parameter = np.array(cond.get("lattice"), dtype=np.float64)
            
            if self.symbol != self.site.symbol:
                raise ValueError(
                    f"Symbol mismatch: Expected '{self.site.symbol}' from site object, "
                    f"but found '{self.symbol}' in '{self.path_traj}'."
                )

            if not np.all(np.abs(self.lattice_parameter - self.site.lattice_parameter) <= eps):
                raise ValueError(
                    f"Lattice parameter mismatch between site object and trajectory file '{self.path_traj}'."
                )
            
            if self.dt <= 0:
                raise ZeroDivisionError(f"Timestep 'dt' must be positive, but got {self.dt}.")
            
            val = (self.t_interval * 1000) / self.dt
            
            if math.isclose(val, round(val), rel_tol=0, abs_tol=eps):
                self.frame_interval = int(round(val))
            else:
                raise ValueError(
                    f"'t_interval' ({self.t_interval} ps) must be a multiple of "+
                    f"'dt' ({self.dt} ps)."
                )
            
            if self.frame_interval == 0:
                raise ValueError(
                    f"'t_interval' ({self.t_interval} ps) is too small "+
                    "compared to 'dt' ({self.dt} ps), resulting in a zero frame interval."
                    )

            self.num_steps = self.total_frames // self.frame_interval
            self.num_frames = self.num_steps * self.frame_interval

    def _read_site(self) -> None:
        self.lattice_sites = np.array(
            [p['coord'] for p in self.site.lattice_sites], dtype=np.float64
        )
        self.lattice_sites_cart = np.array(
            [p['coord_cart'] for p in self.site.lattice_sites], dtype=np.float64
        )
        self.num_lattice_sites = len(self.lattice_sites)
        self.num_vacancies = self.num_lattice_sites - self.num_atoms
    
    def _displacement_PBC(self, 
                          coord1: np.ndarray, 
                          coord2: np.ndarray) -> np.ndarray:
        """Calculates the Cartesian displacement vector between two fractional
        coordinates under Periodic Boundary Conditions (PBC)."""
        disp_frac = coord2 - coord1
        disp_frac -= np.round(disp_frac)
        return np.dot(disp_frac, self.lattice_parameter)
    
    def _nearest_lattice_points(self, 
                                coords: np.ndarray) -> int:
        """Finds the index of the nearest lattice site for each atom."""
        delta = coords[:, np.newaxis, :] - self.lattice_sites[np.newaxis, :, :]
        delta -= np.round(delta)
        dist_sq = np.sum(np.dot(delta, self.lattice_parameter)**2, axis=2)
        return np.argmin(dist_sq, axis=1)
    
    def _get_occupation(self,
                        force_margin : float,
                        cos_margin : float) -> None:
        """Analyzes atomic trajectory to determine site occupations over time.

        This method iterates through the trajectory and assigns each atom to its
        nearest lattice site. For atoms that appear to have hopped between
        timesteps, a robust transition state (TS) criterion is applied to
        validate the hop and prevent miscounting oscillations as jumps.

        The TS criterion has two conditions:
        1. Force Magnitude Check: The force on the atom must be greater than
           `force_margin` for its direction to be considered reliable.
        2. Cosine Difference Check: The force vector must be pointing back
           towards the initial site more significantly than towards the final
           site (by a margin of `cos_margin`) for the hop to be rejected.

        The final occupation history is stored in the `self.occupation` attribute.

        Args:
            force_margin (float, optional): The minimum force magnitude (in eV/Å)
                required to perform the directional TS check.
            cos_margin (float, optional): The required difference between the cosine
                of the angle to the initial site and the final site to reject a hop.
        """
        with h5py.File(self.path_traj, 'r') as f:
            pos_data = f['positions']
            force_data = f['forces']
            
            check_init = False
            occupation = np.zeros((self.num_steps, self.num_atoms), dtype=np.int16)

            for i in range(self.num_steps):
                start = i * self.frame_interval
                end = start + self.frame_interval
                
                pos_chunk_raw = pos_data[start:end]
                force_chunk_raw = force_data[start:end]
                pos_chunk = np.average(pos_chunk_raw.astype(np.float64), axis=0)
                force_chunk = np.average(force_chunk_raw.astype(np.float64), axis=0)
                
                occupation_i = self._nearest_lattice_points(pos_chunk)
                
                if not check_init:
                    if len(set(occupation_i)) == self.num_atoms:
                        for j in range(i + 1):
                            occupation[j] = occupation_i
                        check_init = True
                    continue
                
                indices_move_atom = np.where(occupation_i != occupation[i-1])[0]
                
                for index in indices_move_atom:
                    site_init = occupation[i-1][index]
                    site_final = occupation_i[index]
                    
                    force_atom = force_chunk[index]
                    p_atom = pos_chunk[index]
                    p_init = self.lattice_sites[site_init]
                    p_final = self.lattice_sites[site_final]
                    
                    r_init = self._displacement_PBC(p_atom, p_init)
                    r_final = self._displacement_PBC(p_atom, p_final)
                    
                    norm_f = np.linalg.norm(force_atom)
                    norm_init = np.linalg.norm(r_init)
                    norm_final = np.linalg.norm(r_final)
                    
                    eps = 1e-12
                    if norm_f < eps or norm_init < eps:
                        cos_init = np.nan
                    else:
                        cos_init = np.dot(force_atom, r_init) / (norm_f * norm_init)
                        
                    if norm_f < eps or norm_final < eps:
                        cos_final = np.nan
                    else:
                        cos_final = np.dot(force_atom, r_final) / (norm_f * norm_final)
                        
                    accept_hop = True
                    
                    if norm_f > force_margin and not (np.isnan(cos_init) or np.isnan(cos_final)):
                        if (cos_init - cos_final) > cos_margin:
                            accept_hop = False
                            
                    if not accept_hop:
                        occupation_i[index] = site_init
                        
                occupation[i] = occupation_i
        self.occupation = occupation.T
    
    def _get_trace_arrows(self) -> None:
        """
        Identifies all atomic hops from the occupation data and formats them
        for visualization.
        """
        if self.occupation is None:
            raise RuntimeError("Occupation data is not available. Run _get_occupation() first.")
        
        change_in_occ = np.diff(self.occupation, axis=1)
        move_atom_indices, move_step_indices = np.where(change_in_occ != 0)
        move_step_indices += 1
        
        trace_arrows = {}
        for step, atom_idx in zip(move_step_indices, move_atom_indices):
            site_idx_init = self.occupation[atom_idx][step-1]
            site_idx_final = self.occupation[atom_idx][step]
            
            arrow = {
                'c': self.cmap[atom_idx % len(self.cmap)],
                'lattice_point': [site_idx_init, site_idx_final],
                'p': np.vstack((
                    self.lattice_sites_cart[site_idx_init],
                    self.lattice_sites_cart[site_idx_final]
                ))
            }

            if step in trace_arrows:
                trace_arrows[step].append(arrow)
            else:
                trace_arrows[step] = [arrow]
        
        for step in range(self.num_steps):
            if step not in trace_arrows:
                trace_arrows[step] = []
                
        self.trace_arrows = trace_arrows
    
    def _trace_vacancy_paths(self, site_init, site_final, paths):
        """Helper method to find connected paths for vacancies."""
        path_map = defaultdict(list)
        for to_site, from_site in paths:
            path_map[from_site].append(to_site)

        site_final_set = set(site_final)
        candidate_routes = {s: [] for s in site_init}

        # Step 1: Collect all possible routes from each start site using DFS
        for s_init in site_init:
            stack = [(s_init, [s_init])]
            while stack:
                current, route = stack.pop()
                if current in site_final_set and not (len(route) == 1 and current == s_init):
                    candidate_routes[s_init].append(route)
                for next_site in path_map.get(current, []):
                    if next_site not in route:
                        stack.append((next_site, route + [next_site]))

        # Step 2: Find a valid combination of routes that uses unique end sites
        for ordering in permutations(site_init):
            used_finals = set()
            results, used_paths = [], set()
            for s in ordering:
                found = False
                for route in candidate_routes[s]:
                    if route[-1] not in used_finals:
                        results.append(route)
                        used_finals.add(route[-1])
                        for i in range(len(route) - 1):
                            used_paths.add((route[i+1], route[i]))
                        found = True
                        break
                if not found: results.append(None)

            if len(results) == len(site_init) and None not in results:
                reordered = [None] * len(site_init)
                for i, s in enumerate(ordering):
                    reordered[site_init.index(s)] = results[i]
                unused_paths = list(set(map(tuple, paths)) - used_paths)
                return reordered, unused_paths
        
        # Fallback if no perfect permutation is found
        fallback, used_paths, used_finals = [], set(), set()
        for s in site_init:
            found = False
            for route in candidate_routes.get(s, []):
                if route[-1] not in used_finals:
                    fallback.append(route); used_finals.add(route[-1])
                    for i in range(len(route) - 1): used_paths.add((route[i+1], route[i]))
                    found = True; break
            if not found: fallback.append(None)
        unused_paths = list(set(map(tuple, paths)) - used_paths)
        return fallback, unused_paths

    def _get_vacancy_trajectory(self) -> None:
        """
        Calculates the trajectory of vacancies based on atomic occupation data.
        This method tracks how empty lattice sites move over time, handling
        simple hops and complex, multi-atom movements.
        """
        self.transient_vacancy = {0: np.array([], dtype=np.int16)}
        step_transient = {0: False}
        all_lattice_sites = np.arange(self.num_lattice_sites)
        
        # Find the first stable step with the correct number of vacancies
        step_init = 0
        while step_init < self.num_steps:
            site_vac = np.setdiff1d(all_lattice_sites, self.occupation[:, step_init])
            if len(site_vac) == self.num_vacancies:
                break
            step_init += 1
        
        # Back-fill the initial steps with the first stable vacancy configuration
        for step in range(step_init + 1):
            self.vacancy_trajectory_index[step] = copy.deepcopy(site_vac)
            self.vacancy_trajectory_coord_cart[step] = self.lattice_sites_cart[site_vac]
            
        # Main loop to trace vacancy movement
        for step in range(step_init + 1, self.num_steps):
            site_vac_new = np.setdiff1d(all_lattice_sites, self.occupation[:, step])
            step_transient[step] = len(site_vac_new) > self.num_vacancies
            
            site_init = np.setdiff1d(site_vac, site_vac_new)   # Vacancies that disappeared
            site_final = np.setdiff1d(site_vac_new, site_vac) # Vacancies that appeared

            if len(site_init) == 0: # No hop
                self.vacancy_trajectory_index[step] = copy.deepcopy(site_vac)
                self.vacancy_trajectory_coord_cart[step] = self.lattice_sites_cart[site_vac]
            else: # Hops occurred
                loop = 1
                site_transient = []
                paths = [arrow['lattice_point'] for arrow in self.trace_arrows.get(step, [])]
                while step_transient.get(step - loop, False):
                    paths += [arrow['lattice_point'] for arrow in self.trace_arrows.get(step - loop, [])]
                    site_transient.extend(self.transient_vacancy.get(step - loop, []))
                    loop += 1
                
                effective_site_final = np.array(list(set(list(site_final) + site_transient)))
                path_connect, unused_path = self._trace_vacancy_paths(list(site_init), effective_site_final, paths)
                
                for i, site in enumerate(site_init):
                    path_route = path_connect[i]
                    if path_route is None:
                        raise RuntimeError(
                            f"Failed to find a vacancy trajectory path at step {step}.\n"
                            f"  - Start sites: {site_init}\n"
                            f"  - End sites: {effective_site_final}\n"
                            f"  - Available path segments: {paths}"
                        )
                    
                    site_vac[list(site_vac).index(site)] = path_route[-1]
                    site_final = site_final[site_final != path_route[-1]]
                    effective_site_final = effective_site_final[effective_site_final != path_route[-1]]
                
                site_remain = np.setdiff1d(site_vac, site_vac_new)
                if len(site_remain) > 0:
                    site_unexpect = np.setdiff1d(site_vac_new, site_vac)
                    path_unexpect, _ = self._trace_vacancy_paths(list(site_remain), site_unexpect, unused_path)
                    
                    for i, site in enumerate(site_remain):
                        path_route = path_unexpect[i]
                        if path_route is None: continue # Skip if no path found
                        site_vac[list(site_vac).index(site)] = path_route[-1]
                        site_final = site_final[site_final != path_route[-1]]
                        
                        for path in path_connect:
                            if path and path[-1] == site:
                                path.append(path_route[-1]); break
                
                self.hopping_sequence[step] = copy.deepcopy(path_connect)
                self.vacancy_trajectory_index[step] = copy.deepcopy(site_vac)
                self.vacancy_trajectory_coord_cart[step] = self.lattice_sites_cart[site_vac]
                
            self.transient_vacancy[step] = site_final

    def _get_unwrapped_vacancy_trajectory(self):
        """
        Calculates the continuous, unwrapped Cartesian trajectory for each vacancy.
        It uses the hopping sequence to track displacements across periodic boundaries.
        """
        if not self.vacancy_trajectory_index:
            return

        unwrapped_coords = {step: None for step in range(self.num_steps)}
        
        step_init = 0
        while step_init < self.num_steps:
            if self.vacancy_trajectory_index.get(step_init) is not None:
                break
            step_init += 1
        
        if step_init == self.num_steps: return

        initial_indices = self.vacancy_trajectory_index[step_init]
        unwrapped_coords[step_init] = self.lattice_sites_cart[initial_indices]
        
        prev_unwrapped = unwrapped_coords[step_init]
        prev_indices = initial_indices

        for step in range(step_init + 1, self.num_steps):
            current_indices = self.vacancy_trajectory_index[step]
            current_unwrapped = prev_unwrapped.copy()

            hops = self.hopping_sequence.get(step, [])
            for route in hops:
                if route is None: continue
                
                start_idx, end_idx = route[0], route[-1]

                vac_k_list = np.where(prev_indices == start_idx)[0]
                if len(vac_k_list) == 0: continue
                vac_k = vac_k_list[0]

                start_cart_unwrapped = prev_unwrapped[vac_k]

                disp_cart = self._displacement_PBC(
                    self.lattice_sites[start_idx], 
                    self.lattice_sites[end_idx]
                )

                current_unwrapped[vac_k] = start_cart_unwrapped + disp_cart

            unwrapped_coords[step] = current_unwrapped.copy()
            prev_unwrapped = current_unwrapped
            prev_indices = current_indices
        
        self.unwrapped_vacancy_trajectory_coord_cart = unwrapped_coords


class TrajectoryAnalyzer:
    """Analyzes and quantifies hopping statistics from Trajectory and Site objects.

    This class serves as a post-processor for the `Trajectory` analysis. It takes
    pre-computed lattice information (from a `Site` object) and a detailed hop
    analysis (from a `Trajectory` object) to calculate key diffusion statistics,
    such as hop counts per path, site residence times, and the theoretical
    mean squared displacement from a random walk perspective.

    The main workflow is to initialize the class with a fully processed `Trajectory`
    object. All calculations are performed automatically upon initialization, and
    results can be viewed using the `.summary()` method.

    Args:
        trajectory (Trajectory):
            An initialized and fully processed `Trajectory` object.
        site (Site | None, optional):
            An initialized `Site` object. If None, the `site` object
            associated with the `trajectory` object will be used.
            Defaults to None.
        eps (float, optional):
            A tolerance value for floating-point comparisons when categorizing
            hopping paths. Defaults to 1e-3.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Attributes:
        hopping_history (list[list[dict]]):
            A nested list containing a time-ordered sequence of hop dictionaries
            for each vacancy.
        counts (numpy.ndarray):
            A 2D array of shape (num_vacancies, num_paths) counting the
            occurrences of each pre-defined hop path.
        counts_unknown (numpy.ndarray):
            A 2D array counting occurrences of dynamically discovered "unknown"
            hop paths that were not in the initial `Site` definition.
        residence_time (numpy.ndarray):
            A 2D array of shape (num_vacancies, num_sites) storing the total
            time (in ps) each vacancy spent at each inequivalent site type.
        msd_rand (float):
            The Mean Squared Displacement (Å²) calculated from all observed hops,
            based on random walk theory.

    Raises:
        AttributeError: If the input `trajectory` or `site` objects are missing
            required attributes for the analysis.

    Examples:
        >>> # First, create Site and Trajectory objects
        >>> site_info = Site(path_structure='POSCAR', symbol='O')
        >>> traj_analysis = Trajectory(
        ...     path_traj='TRAJ_O.h5',
        ...     site=site_info,
        ...     t_interval=0.1
        ... )
        >>> # Now, create the analyzer with the trajectory results
        >>> analyzer = TrajectoryAnalyzer(trajectory=traj_analysis, verbose=True)
        # A full summary of hopping statistics is printed upon initialization.
    """
    def __init__(self,
                 trajectory: Trajectory,
                 site: Site | None = None,
                 eps: float = 1e-3,
                 verbose: bool = True):
        
        self.trajectory = trajectory
        
        if site is None:
            self.site = trajectory.site
        else:
            self.site = site
            
        self.num_vacancies = trajectory.num_vacancies
        self.eps = eps
        self.verbose = verbose
        
        # site(lattice)
        self.path = self.site.path
        self.path_name = self.site.path_name
        self.site_name = self.site.site_name
        self.path_distance = np.array([p['distance'] for p in self.site.path])
        self.lattice_sites_info = self.site.lattice_sites

        self.path_unknown = []
        
        # hopping history
        self.hopping_history = [[] for _ in range(self.num_vacancies)]
        self.counts = np.zeros((self.num_vacancies, len(self.path_name)))
        self._hopping_statistics()
        
        # unknown path
        self.unknown_name = [p['name'] for p in self.path_unknown]
        self.counts_unknown = np.zeros((self.num_vacancies, len(self.path_unknown)))
        self._counts_unknown_path()
        
        # random walk msd
        self.msd_rand = None
        self._random_walk_msd()
        
        # vacancy residence time
        self.residence_time = np.zeros((self.num_vacancies, len(self.site_name)))
        self._get_residence_time()
        
        if verbose:
            self.summary()
            
    def _hopping_statistics(self) -> None:
        """Categorizes each hop from the trajectory into known or unknown paths."""
        for step, sequence in self.trajectory.hopping_sequence.items():
            if not sequence: continue
            for path_route in sequence:
                if path_route is None: continue
                
                try:
                    index_vac = list(self.trajectory.vacancy_trajectory_index[step]).index(path_route[-1])
                except ValueError:
                    continue
                
                for i in range(len(path_route) - 1):
                    index_init = path_route[i]
                    index_final = path_route[i+1]
                    
                    distance = self.trajectory.distance_PBC(
                        self.trajectory.lattice_sites[index_init],
                        self.trajectory.lattice_sites[index_final]
                    )
                    
                    # hop 분류
                    check_normal = False
                    index_path = -1
                    site_init_name = self.lattice_sites_info[index_init]['site']
                    site_final_name = self.lattice_sites_info[index_final]['site']

                    for j, p in enumerate(self.path):
                        check1 = abs(p['distance'] - distance) < self.eps
                        check2 = p['site_init'] == site_init_name
                        check3 = p['site_final'] == site_final_name
                        
                        if check1 and check2 and check3:
                            check_normal = True
                            index_path = j
                            break
                    
                    if check_normal:
                        path_info = copy.deepcopy(self.path[index_path])
                        path_info.update({'step': step, 'index_init': index_init, 'index_final': index_final})
                        self.hopping_history[index_vac].append(path_info)
                        self.counts[index_vac, index_path] += 1
                    else:
                        check_unknown = False
                        for j, p in enumerate(self.path_unknown):
                            check1 = abs(p['distance'] - distance) < self.eps
                            check2 = p['site_init'] == site_init_name
                            check3 = p['site_final'] == site_final_name
                            if check1 and check2 and check3:
                                check_unknown = True
                                index_path = j
                                break
                        
                        if check_unknown:
                            path_info = copy.deepcopy(self.path_unknown[index_path])
                            path_info.update({'step': step, 'index_init': index_init, 'index_final': index_final})
                            self.hopping_history[index_vac].append(path_info)
                        else:
                            unknown_new = {
                                'site_init': site_init_name,
                                'site_final': site_final_name,
                                'distance': distance,
                                'coord_init': self.lattice_sites_info[index_init]['coord'],
                                'coord_final': self.lattice_sites_info[index_final]['coord'],
                                'name': f"unknown{len(self.path_unknown)+1}"
                            }
                            self.path_unknown.append(copy.deepcopy(unknown_new))
                            unknown_new.update({'step': step, 'index_init': index_init, 'index_final': index_final})
                            self.hopping_history[index_vac].append(unknown_new)
    
    def _get_residence_time(self) -> None:
        """Calculates the residence time of each vacancy at each inequivalent site."""
        for indices in self.trajectory.vacancy_trajectory_index.values():
            for i, index in enumerate(indices):
                index_site = self.site_name.index(self.site.lattice_sites[index]['site'])
                self.residence_time[i, index_site] += 1
        self.residence_time *= self.trajectory.t_interval

    def _counts_unknown_path(self) -> None:
        """Counts the occurrences of newly found unknown paths."""
        self.counts_unknown = np.zeros((self.num_vacancies, len(self.path_unknown)))
        self.unknown_name = [p['name'] for p in self.path_unknown]
        for index_vac in range(self.num_vacancies):
            for path in self.hopping_history[index_vac]:
                if 'unknown' in path['name']:
                    try:
                        index_path = self.unknown_name.index(path['name'])
                        self.counts_unknown[index_vac, index_path] += 1
                    except ValueError:
                        continue

    def _random_walk_msd(self) -> None:
        """Calculates the Mean Squared Displacement based on random walk theory."""
        if not self.path_unknown:
            distance_all = self.path_distance
            counts_all = self.counts
        else:
            distance_all = np.array(
                list(self.path_distance) + [p['distance'] for p in self.path_unknown]
            )
            counts_all = np.hstack((self.counts, self.counts_unknown))
        
        self.msd_rand = np.average(
            np.sum(distance_all**2 * counts_all, axis=1)
        )

    def summary(self) -> None:
        """Prints a comprehensive, multi-part summary of the hopping analysis.

        The summary includes tables for:
        - Path Counts: The number of times each defined (and unknown) path was traversed by each vacancy.
        - Residence Time: The total time each vacancy spent at each inequivalent site type.
        - Hopping Sequence: A detailed, time-ordered log of every hop for each vacancy.
        """
        # update unknown path information (for bundle application)
        self._counts_unknown_path()
        
        # Path counts
        name_all = self.path_name + self.unknown_name
        counts_all = np.hstack((self.counts, self.counts_unknown))
        counts_all = np.array(counts_all, dtype=np.int32)
        vacancy_name = [f"Vacancy{i+1}" for i in range(self.num_vacancies)]
        print("# Path counts :")
        header = ['path'] + vacancy_name
        data = np.vstack((name_all, counts_all)).T
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        # Vacancy residence time
        print("# Vacancy residence time (ps) :")
        header = ['site'] + vacancy_name
        data = np.vstack((self.site_name, self.residence_time)).T
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')
        
        # Hopping sequence
        print("# Hopping sequence :")
        for i in range(self.num_vacancies):
            print(f"# Vacancy{i+1}")
            header = ['num', 'time (ps)', 'path', 'a (Ang)', 'initial site', 'final site']
            data = [
                [
                    f"{j+1}",
                    f"{path['step'] * self.trajectory.t_interval:.2f}",
                    f"{path['name']}",
                    f"{path['distance']:.5f}",
                    f"{path['site_init']} [{', '.join(f'{x:.5f}' for x in self.site.lattice_sites[path['index_init']]['coord'])}]",
                    f"{path['site_final']} [{', '.join(f'{x:.5f}' for x in self.site.lattice_sites[path['index_final']]['coord'])}]"
                ] for j, path in enumerate(self.hopping_history[i])
            ]
            print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
            print('')
            

class Encounter:
    """Analyzes atom-vacancy encounters to calculate diffusion correlation factors.

    An "encounter" is a sequence of correlated hops involving a **specific atom**
    and a **specific vacancy**. This sequence continues as long as the atom
    interacts with the same vacancy, even if it hops away and returns.

    The encounter is considered **terminated** when the atom performs a hop by
    exchanging with a **different** vacancy. At that moment, a new encounter
    begins with the new vacancy. This class identifies these events to compute
    the correlation factor (f).

    Args:
        analyzer (TrajectoryAnalyzer):
            An initialized `TrajectoryAnalyzer` object containing the full
            hopping statistics from a trajectory.
        use_incomplete_encounter (bool, optional):
            If True, encounters that were still in progress at the end of the
            simulation are included in the statistical analysis. This can improve
            statistics for short simulations but may slightly skew results.
            Defaults to True.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Attributes:
        f_cor (float):
            The calculated correlation factor (f), a measure of how
            correlated successive jumps of an atom are.
        msd (float):
            The mean squared displacement (Å²) of atoms, averaged over all
            analyzed encounters.
        path_count (numpy.ndarray):
            An array counting the total number of times each pre-defined path
            was traversed across all encounters.
        num_encounter (int):
            The total number of encounters used in the analysis.
        encounter_complete (list[dict]):
            A list of all encounter events that were fully completed within the
            simulation time.
        encounter_in_process (list[dict]):
            A list of encounter events that were still ongoing when the
            simulation ended.

    Raises:
        AttributeError:
            If the input `analyzer` object is missing required, pre-computed
            attributes (e.g., `hopping_sequence`).

    Examples:
        >>> # Assuming site_info and traj_analysis objects have been created
        >>> analyzer = TrajectoryAnalyzer(trajectory=traj_analysis)
        >>> encounter_stats = Encounter(analyzer=analyzer, verbose=True)
        # A summary is printed upon initialization if encounters are found.
        # Access results directly:
        # print(f"Correlation Factor (f): {encounter_stats.f_cor:.4f}")
    """
    def __init__(self,
                 analyzer: TrajectoryAnalyzer,
                 use_incomplete_encounter: bool = True,
                 verbose: bool = True):
        self.analyzer = analyzer
        self.trajectory = analyzer.trajectory
        self.site = analyzer.site
        self.eps = analyzer.eps
        self.use_incomplete_encounter = use_incomplete_encounter
        self.verbose = verbose
        
        # path information
        self.path = self.site.path + self.analyzer.path_unknown
        self.path_name = [p['name'] for p in self.path]
        
        # unwrapped vacancy trajectory
        self.vacancy_coord_unwrap = None
        self._get_unwrapped_vacancy_trajectory()
        
        # encounter
        self.encounter_complete = []
        self.encounter_in_process = []
        self._find_encounters()
        self.num_encounter_complete = len(self.encounter_complete)
        self.num_encounter_incomplete = len(self.encounter_in_process)
        
        # encounter information
        self.msd = None
        self.path_count = None
        self.f_cor = None
        self.path_distance = None
        
        self.encounter_all = list(self.encounter_complete)
        if use_incomplete_encounter:
            self.encounter_all += self.encounter_in_process
        self.num_encounter = len(self.encounter_all)
        
        if self.num_encounter == 0:
            if verbose:
                print("Warning: No complete encounters were found to analyze.")
        else:
            self._analyze_encounter()
            self._calculate_correlation_factor()
            if verbose:
                self.summary()
                
    def _find_path_name(self, index_init: int, index_final: int) -> str | None:
        """Finds the pre-defined path name for a hop between two site indices."""
        site_init_name = self.site.lattice_sites[index_init]['site']
        site_final_name = self.site.lattice_sites[index_final]['site']
        
        coord_init = self.site.lattice_sites[index_init]['coord']
        coord_final = self.site.lattice_sites[index_final]['coord']
        distance = self.trajectory.distance_PBC(coord_init, coord_final)
        
        for path in self.path:
            if (path['site_init'] == site_init_name and
                path['site_final'] == site_final_name and
                abs(path['distance'] - distance) < self.eps):
                return path['name']
        return None

    def _get_unwrapped_vacancy_trajectory(self) -> None:
        """Generates a simplified unwrapped trajectory for all vacancies combined."""
        source = self.trajectory.unwrapped_vacancy_trajectory_coord_cart
        num_steps = self.trajectory.num_steps
        num_vac = self.trajectory.num_vacancies
        
        self.vacancy_coord_unwrap = np.zeros((num_steps, num_vac, 3))
        for step, coords in source.items():
            if coords is not None and len(coords) == num_vac:
                self.vacancy_coord_unwrap[step] = coords

    def _find_encounters(self) -> None:
        """
        Identifies and categorizes atom-vacancy encounters from the hopping sequence.
        """
        for step, sequence in self.trajectory.hopping_sequence.items():
            for path_connect in sequence:
                if path_connect is None: continue
                # vacancy index
                try:
                    index_vac = list(
                        self.trajectory.vacancy_trajectory_index[step]
                    ).index(path_connect[-1])
                except (ValueError, IndexError): continue
                
                # decompose path
                trace_arrow = [
                    path_connect[i-1:i+1][::-1] for i in range(len(path_connect)-1, 0, -1)
                ]
                
                # encounter analysis
                coord_init = self.vacancy_coord_unwrap[step][index_vac]
                for path in trace_arrow:
                    # atom index
                    try:
                        index_atom = list(self.trajectory.occupation[:, step]).index(path[-1])
                    except ValueError:
                        loop = 1
                        found_atom = False
                        while step - loop >= 0:
                            match = next(
                                (arrow for arrow in self.trajectory.trace_arrows.get(step-loop, []) 
                                 if path == arrow['lattice_point']), 
                                None
                            )
                            if match:
                                try:
                                    index_atom = list(
                                        self.trajectory.occupation[:, step-loop]
                                    ).index(path[-1])
                                    found_atom = True
                                    break
                                except ValueError: pass
                            loop += 1
                        if not found_atom: continue
                        
                    # path name
                    path_name = self._find_path_name(*path)
                    if path_name is None: continue
                    
                    # unwrapped atomic coord
                    disp_frac = self.trajectory.lattice_sites[path[1]] - self.trajectory.lattice_sites[path[0]]
                    disp_frac -= np.round(disp_frac)
                    coord_final = coord_init + np.dot(disp_frac, self.trajectory.lattice_parameter)
                    
                    # comparison with existing encounters
                    index_encounter = None
                    for i, enc in enumerate(self.encounter_in_process):
                        if enc['index_atom'] == index_atom:
                            index_encounter = i
                            break
                        
                    # case 1. no matching encounter
                    if index_encounter is None:
                        encounter = {'index_atom': index_atom, 'index_vac': index_vac,
                                     'coord_init': coord_init, 'coord_final': coord_final,
                                     'hopping_history': [path_name]}
                        self.encounter_in_process.append(encounter)
                        coord_init = coord_final
                        continue
                    
                    # matching encounter
                    encounter_match = self.encounter_in_process[index_encounter]
                    coord_encounter = encounter_match['coord_final']
                    distance = np.linalg.norm(
                        np.dot(coord_encounter - coord_init, self.trajectory.lattice_parameter)
                    )
                    
                    # case 2. exactly matching encounter
                    if distance < self.eps:
                        # exchange with the associated vacancy : update encoutner
                        if encounter_match['index_vac'] == index_vac:
                            encounter_match['coord_final'] = coord_final
                            encounter_match['hopping_history'].append(path_name)
                            
                        # exchange with a new vacancy : terminate encounter
                        else:
                            # terminate the existing encounter
                            # self.encounter_complete.append(encounter_match.copy())
                            self.encounter_complete.append(copy.deepcopy(encounter_match))
                            del self.encounter_in_process[index_encounter]
                            
                            # initiate a new encounter
                            encounter = {
                                'index_atom': index_atom, 
                                'index_vac': index_vac,
                                'coord_init': coord_init, 
                                'coord_final': coord_final,
                                'hopping_history': [path_name]
                            }
                            self.encounter_in_process.append(encounter)
                            
                    # case 3. PBC matching encounter:
                    else:
                        # terminate the existing encounter
                        # self.encounter_complete.append(encounter_match.copy())
                        self.encounter_complete.append(copy.deepcopy(encounter_match))
                        del self.encounter_in_process[index_encounter]
                        
                        # initiate a new encounter
                        encounter = {
                            'index_atom': index_atom, 
                            'index_vac': index_vac,
                            'coord_init': coord_init, 
                            'coord_final': coord_final,
                            'hopping_history': [path_name]
                        }
                        self.encounter_in_process.append(encounter)
                    
                    coord_init = coord_final

    def _analyze_encounter(self):
        """
        Calculates MSD and path counts from the encounters, considering only
        pre-defined paths from the Site object.
        """
        displacement = []
        self.path_count = np.zeros(len(self.path_name), dtype=np.float64)
        
        for encounter in self.encounter_all:
            for name in encounter['hopping_history']:
                try:
                    index_path = self.path_name.index(name)
                    self.path_count[index_path] += 1
                except ValueError:
                    continue

            disp = encounter['coord_final'] - encounter['coord_init']
            displacement.append(disp)
            
        displacement = np.array(displacement)
        if displacement.size == 0:
            self.msd = 0.0
        else:
            self.msd = np.average(np.sum(displacement**2, axis=1))
            
    def _calculate_correlation_factor(self):
        """Calculates the tracer correlation factor f."""
        self.path_distance = np.array([path['distance'] for path in self.path])
        denominator = np.sum(self.path_distance**2 * (self.path_count / self.num_encounter))
        
        if denominator > 1e-9:
            self.f_cor = self.msd / denominator
        else:
            self.f_cor = np.nan

    def summary(self) -> None:
        """Prints a comprehensive, formatted summary of the encounter analysis.

        The summary includes key calculated values like the correlation factor
        and MSD, statistics on the number and type of encounters, and a
        detailed table of path-wise hop counts within the encounters.
        """
        print("# Encounter Analysis")
        print(f"  Use incomplete encounters : {self.use_incomplete_encounter}")
        print(f"  Correlation factor (f)    : {self.f_cor:.5f}")
        print(f"  Mean Squared Disp. (MSD)  : {self.msd:.5f} Ang^2")
        print(f"  Num. complete encounters  : {self.num_encounter_complete}")
        print(f"  Num. incomplete encounters: {self.num_encounter_incomplete}")
        print(f"  Num. encounters in use    : {self.num_encounter}")
        print(f"  Total hopping events      : {int(np.sum(self.path_count))}")
        if self.num_encounter > 0:
            print(f"  Mean hops per encounter   : {np.sum(self.path_count) / self.num_encounter:.5f}")
        print('') 
        
        print(f"# Pathwise Counts in Encounters")
        header = ['Path', 'a (Ang)', 'Count', 'Count/Encounter']
        data = [
            [
                name,
                f"{dist:.5f}",
                f"{int(count)}",
                f"{count / self.num_encounter:.5f}" if self.num_encounter > 0 else "N/A"
            ]
            for name, dist, count in zip(self.path_name, self.path_distance, self.path_count)
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print('')


class TrajectoryBundle:
    """Manages and validates one or more HDF5 trajectory files.

    This class handles two main scenarios:

    1.  **Single File Mode**: Processes a single, specified HDF5 trajectory file.

    2.  **Directory Search Mode**: Recursively searches a directory for HDF5 files,
        groups them by temperature, and validates their consistency.

    It ensures that all processed files are from consistent simulations (e.g.,
    same timestep, atom counts, lattice). The primary workflow is to initialize
    the class, which triggers file discovery and validation. Results can then be
    saved to a summary file using the `.summary()` method.

    Args:
        path_traj (str):
            Path to a single HDF5 file or a root directory to search for files.
        symbol (str):
            The chemical symbol of the target element, used to filter files.
        prefix (str, optional):
            The prefix of trajectory files to search for (e.g., "TRAJ").
            Used only when searching a directory. Defaults to "TRAJ".
        depth (int, optional):
            The maximum directory depth to search for files. Used only when
            searching a directory. Defaults to 2.
        eps (float, optional):
            The tolerance for comparing floating-point values in metadata
            (e.g., temperature, dt, lattice vectors). Defaults to 1.0e-3.
        verbose (bool, optional):
            Verbosity flag. Defaults to False.

    Attributes:
        temperatures (list[float]):
            A sorted list of unique temperatures found across all valid files.
        traj (list[list[str]]):
            A nested list where each sublist contains the full file paths for a
            corresponding temperature in `self.temperatures`.
        atom_count (int):
            The number of atoms for the specified symbol, confirmed to be
            consistent across all files.
        lattice (numpy.ndarray):
            The 3x3 lattice vectors as a NumPy array, confirmed to be consistent.

    Raises:
        FileNotFoundError:
            If the provided `path_traj` does not exist, or if no valid trajectory
            files matching the criteria are found during a directory search.
        ValueError:
            If `depth` is less than 1, or if inconsistent simulation parameters
            (dt, atom_counts, lattice) are found among the files.
        IOError:
            If a trajectory file or its metadata cannot be read.
    """
    def __init__(self,
                 path_traj: str,
                 symbol: str,
                 prefix: str = "TRAJ",
                 depth: int = 2,
                 eps: float = 1.0e-3,
                 verbose: bool = False):
        
        self.eps = eps
        self.path_traj = Path(path_traj)
        self.depth = depth
        self.symbol = symbol
        self.prefix = prefix
        self.verbose = verbose
        
        if not self.path_traj.exists():
            raise FileNotFoundError(f"Error: '{self.path_traj}' not found.")
        
        if self.depth < 1:
            raise ValueError("Error: depth must be 1 or greater.")

        # List of TRAJ_*_.h5 files 
        self.temperatures = []
        self.traj = []
        self._search_traj()
        
        self.atom_count = None
        self.lattice = None
        self._validate_consistency()
        
        if self.verbose:
            self.summary()
        
    def _validate_consistency(self) -> None:
        """
        Validates that all found trajectory files share consistent simulation parameters
        and sets the class attributes for lattice and atom_count.
        """
        all_paths = list(itertools.chain.from_iterable(self.traj))
        
        if not all_paths:
            return

        ref_path = all_paths[0]
        try:
            with h5py.File(ref_path, "r") as f:
                cond = json.loads(f.attrs["metadata"])
                ref_dt = cond.get("dt")
                ref_atom_counts = cond.get("atom_counts")
                ref_lattice = np.array(cond.get("lattice"), dtype=np.float64)
                
                self.lattice = ref_lattice
                self.atom_count = ref_atom_counts.get(self.symbol)

        except Exception as e:
            raise IOError(f"Could not read reference metadata from '{ref_path}'. Reason: {e}")
        
        if len(all_paths) > 1:
            for path in all_paths[1:]:
                try:
                    with h5py.File(path, "r") as f:
                        cond = json.loads(f.attrs['metadata'])
                        
                        current_dt = cond.get("dt")
                        if abs(current_dt - ref_dt) > self.eps:
                            raise ValueError(
                                f"Inconsistent 'dt' parameter found.\n"
                                f"  - Reference '{ref_path}': {ref_dt}\n"
                                f"  - Conflicting '{path}': {current_dt}"
                            )
                        
                        current_atom_counts = cond.get('atom_counts')
                        if current_atom_counts != ref_atom_counts:
                            raise ValueError(
                                f"Inconsistent 'atom_counts' parameter found.\n"
                                f"  - Reference '{ref_path}': {ref_atom_counts}\n"
                                f"  - Conflicting '{path}': {current_atom_counts}"
                            )
                            
                        current_lattice = np.array(cond.get('lattice'))
                        if not np.all(np.abs(current_lattice - ref_lattice) <= self.eps):
                            raise ValueError(
                                f"Inconsistent 'lattice' parameter found.\n"
                                f"  - Reference '{ref_path}':\n{ref_lattice}\n"
                                f"  - Conflicting '{path}':\n{current_lattice}"
                            )
                except Exception as e:
                    if isinstance(e, ValueError):
                        raise
                    raise IOError(f"Could not read or validate metadata from '{path}'. Reason: {e}")
        
    def _search_traj(self) -> None:
        """
        Searches for HDF5 trajectory files, groups them by temperature,
        and populates self.temperatures and self.traj lists.
        """
        if self.path_traj.is_file():
            candidate_files = [self.path_traj]
                
        elif self.path_traj.is_dir():
            glob_iterators = []
            for i in range(self.depth):
                dir_prefix = '*/' * i
                pattern = f"{dir_prefix}{self.prefix}*.h5"
                glob_iterators.append(self.path_traj.glob(pattern))    
            candidate_files = itertools.chain.from_iterable(glob_iterators)
            
        else:
            raise ValueError(f"Error: The path '{self.path_traj}' is not a regular file or directory.")
        
        found_paths = []
        for file_path in sorted(candidate_files):
            try:
                with h5py.File(file_path, "r") as f:
                    metadata_str = f.attrs.get("metadata")
                    if not metadata_str:
                        if self.verbose:
                            print(f"Warning: File '{str(file_path.resolve())}' is " +
                                  "missing 'metadata' attribute. Skipping.")
                        continue
                    
                    cond = json.loads(metadata_str)
                    file_symbol = cond.get("symbol")
                    file_temp = cond.get("temperature")
                    
                    if file_symbol == self.symbol:
                        if file_temp is None:
                            if self.verbose:
                                print(f"Warning: File '{str(file_path.resolve())}' is missing " +
                                      "'temperature' in metadata. Skipping.")
                            continue
                        
                        if "positions" in f and "forces" in f:
                            found_paths.append((str(file_path.resolve()), float(file_temp)))
                        else:
                            if self.verbose:
                                print(f"Warning: File '{str(file_path.resolve())}' is missing " +
                                      "required datasets ('positions', 'forces'). Skipping.")
                                
            except Exception as e:
                if self.verbose:
                    print(f"Warning: Could not read or verify '{str(file_path.resolve())}'. " +
                          f"Skipping file. Reason: {e}")
        
        if not found_paths:
            raise FileNotFoundError(
                f"Error: No valid trajectory files found for symbol '{self.symbol}' " +
                f"in path '{self.path_traj}' with depth {self.depth}."
            )
            
        sorted_files = sorted(found_paths, key=lambda item: item[1])
        temp_groups, traj_groups = [], []
        for path, temp in sorted_files:
            if not traj_groups or abs(temp - temp_groups[-1]) > self.eps:
                temp_groups.append(temp)
                traj_groups.append([path])
            else:
                traj_groups[-1].append(path)
        
        self.temperatures = temp_groups
        self.traj = traj_groups
    
    def summary(self, filename: str = 'Bundle.txt') -> None:
        """Creates a text file summarizing the found trajectories.

        This method generates a detailed text summary that includes:
        - Consistent system information (symbol, atom count, lattice).
        - A list of all found trajectory files, grouped by temperature.
        - The simulation time for each file and the total time per temperature.

        The summary is saved to the specified file.

        Args:
            filename (str, optional):
                The path and name of the output summary file.
                Defaults to 'bundle.txt'.

        Returns:
            None: This method does not return a value; it writes a file to disk.

        Examples:
            >>> bundle = TrajectoryBundle(...)
            >>> # Save summary to the default 'bundle.txt'
            >>> bundle.summary()

            >>> # Save summary to a custom file
            >>> bundle.summary(filename='simulation_summary.log')
        """
        with open(filename, 'w', encoding='utf-8') as f:
            f.write("--- System Information ---\n")
            f.write(f"Symbol        : {self.symbol}\n")
            f.write(f"Atom Count    : {self.atom_count}\n")
            f.write(f"Lattice (Ang) :\n{self.lattice}\n\n")
            f.write("--- Trajectory Summary ---\n\n")
            
            for temp, paths in zip(self.temperatures, self.traj):
                f.write(f"--- Temperature: {temp:.1f} K (Found {len(paths)} files) ---\n")
                
                total_sim_time = 0.0
                
                for path in paths:
                    try:
                        with h5py.File(path, 'r') as h5f:
                            cond = json.loads(h5f.attrs['metadata'])
                            nsw = cond.get('nsw')
                            dt = cond.get('dt')
                            
                            if nsw is not None and dt is not None:
                                sim_time = (nsw * dt) / 1000.0  # in ps
                                total_sim_time += sim_time
                                f.write(f"  - {path}: {sim_time:.2f} ps\n")
                            else:
                                f.write(f"  - {path}: [nsw/dt not found in metadata]\n")
                    except Exception as e:
                        f.write(f"  - {path}: [Error reading metadata: {e}]\n")
                        
                f.write(f"\n  Total simulation time: {total_sim_time:.2f} ps\n\n")
        print(f"Summary written to '{filename}'")
        

class CalculatorSingle(Trajectory):
    """Orchestrates the full vacancy diffusion analysis for a single trajectory.

    This class provides a high-level interface for analyzing a single simulation run.
    Upon initialization, it immediately performs the foundational analysis by
    processing the trajectory to identify all hopping events and atom-vacancy
    encounters.

    After the object is created, the `.calculate()` method must be called to
    compute the final, derived physical properties (e.g., diffusivity, correlation
    factor, residence time) from the initial results.

    Args:
        path_traj (str):
            Path to the HDF5 trajectory file.
        site (Site):
            An initialized `Site` object containing lattice structure and path info.
        t_interval (float):
            The time interval in picoseconds (ps) for coarse-graining the analysis.
        eps (float, optional):
            A tolerance for floating-point comparisons. Defaults to 1.0e-3.
        use_incomplete_encounter (bool, optional):
            If True, incomplete encounters are included in the final analysis.
            Defaults to True.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Attributes:
        f (float):
            The calculated correlation factor. Populated after `.calculate()` is called.
        D (float):
            The final vacancy diffusivity (m²/s). Populated after `.calculate()` is called.
        D_rand (float):
            The random-walk diffusivity (m²/s). Populated after `.calculate()` is called.
        tau (float):
            The average vacancy residence time (ps). Populated after `.calculate()` is called.
        a (float):
            Effective hopping distance (Å). Populated after `.calculate()` is called.
        hopping_history (list):
            A detailed log of every hop for each vacancy, available after initialization.
        counts (numpy.ndarray):
            An array containing the counts for each predefined hopping path.

    Raises:
        FileNotFoundError:
            If the specified `path_traj` does not exist.
        ValueError:
            If critical data is missing from the trajectory or site information.

    Examples:
        >>> # 1. Initialize the Site and CalculatorSingle objects
        >>> site_info = Site("POSCAR", symbol="O")
        >>> calc = CalculatorSingle(
        ...     path_traj="TRAJ_O_1000K.h5",
        ...     site=site_info,
        ...     t_interval=0.1
        ... )
        >>>
        >>> # 2. Compute the final physical properties
        >>> calc.calculate()
        >>>
        >>> # 3. View the results
        >>> calc.summary()
        >>> print(f"Correlation Factor: {calc.f:.4f}")
    """
    @monitor_performance
    def __init__(self,
                 path_traj: str,
                 site: Site,
                 t_interval: float,
                 eps: float = 1.0e-3,
                 use_incomplete_encounter: bool = True,
                 verbose: bool = True):
        
        self.analyzer = None
        self.encounter = None
        
        self.f = None
        self.D_rand = None
        self.D = None
        self.tau = None
        self.a = None
        
        # Trajectory
        super().__init__(path_traj=path_traj, 
                         site=site, 
                         t_interval=t_interval, 
                         verbose=False)
        
        # TrajectoryAnalyzer
        self.analyzer = TrajectoryAnalyzer(
            trajectory=self, 
            site=self.site, 
            eps=eps, 
            verbose=False
        )
        self.hopping_history = self.analyzer.hopping_history
        self.counts = self.analyzer.counts
        self.path_unknown = self.analyzer.path_unknown
        self.unknown_name = self.analyzer.unknown_name
        self.counts_unknown = self.analyzer.counts_unknown
        self.residence_time = self.analyzer.residence_time
        self.msd_rand = self.analyzer.msd_rand
        
        # Encounter
        self.encounter = Encounter(
            analyzer=self.analyzer,
            use_incomplete_encounter=use_incomplete_encounter,
            verbose=False
        )
        
        self.f = self.encounter.f_cor
        self.verbose = verbose
        self.calculate_is_done = False
        
        # extra properties
        self.a_path = None
        self.z_path = None
        self.temperatures = None
        self.P_site = None
        self.times_site = None
        self.counts_hop = None
     
    def _get_correlation_factor(self):
        """Extracts the correlation factor from the encounter analysis."""
        if self.f is None:
            self.f = self.encounter.f_cor
    
    def _get_random_walk_diffusivity(self):
        """Calculates the random walk diffusivity."""
        total_time = self.t_interval * self.num_steps
        self.D_rand = self.msd_rand / (6 * total_time) * 1e-8 # m2/s
    
    def _get_diffusivity(self):
        """Calculates the final tracer diffusivity."""
        self.D = self.D_rand * self.f
        
    def _get_residence_time(self):
        """Calculates the mean residence time."""
        total_time = self.t_interval * self.num_steps
        total_jumps = np.sum(self.counts) / self.num_vacancies
        self.tau = total_time / total_jumps
        
    def _get_mean_number_of_equivalent_paths(self):
        """Calculates the weighted harmonic mean of equivalent paths (z)."""
        z = np.array([path['z'] for path in self.site.path], dtype=np.float64)
        counts = np.sum(self.counts, axis=0)
        self.z_mean = np.sum(counts) / np.sum(counts / z)
    
    def _get_hopping_distance(self):
        """Calculates hopping distance"""
        self.a = np.sqrt(6 * self.D_rand * self.tau) * 1e4
        
    def _get_extra_properties(self):
        """Calculates extra properties: a_path, z_path, temperatures, counts_hop, times"""
        self.a_path, self.z_path = [], []
        for path in self.site.path:
            self.a_path.append(path['distance'])
            self.z_path.append(path['z'])
        self.temperatures = np.array([self.temperature])
        self.counts_hop = [np.sum(self.counts, axis=0).tolist()]
        
        name_to_type_index = {name: i for i, name in enumerate(self.site.site_name)}
        site_index_to_type_map = np.array(
            [name_to_type_index[site['site']] for site in self.site.lattice_sites]
        )
        all_indices = np.concatenate(list(self.vacancy_trajectory_index.values()))
        all_type_indices = site_index_to_type_map[all_indices]
        count_site = np.bincount(all_type_indices, minlength=len(self.site.site_name))
        self.times_site = count_site * self.t_interval
        self.P_site = self.times_site / np.sum(self.times_site)
        
        self.times_site = [self.times_site.tolist()]
        self.P_site = [self.P_site.tolist()]
    
    def calculate(self, **kwargs) -> None:
        """Calculates the final physical properties from the initial analysis.

        This method should be called after the `CalculatorSingle` object has been
        initialized. It computes derived quantities like diffusivity, correlation
        factor, and residence time from the already-processed hop and encounter data.
        """
        self._get_correlation_factor()
        self._get_random_walk_diffusivity()
        self._get_diffusivity()
        self._get_residence_time()
        self._get_hopping_distance()
        self._get_mean_number_of_equivalent_paths()
        self._get_extra_properties()
        self.calculate_is_done = True
        
    def plot_counts(self,
                    title: str | None = None,
                    save: bool = True, 
                    filename: str = 'counts.png', 
                    dpi: int = 300) -> None:
        """Generates a bar plot showing the total counts for each migration path.

        This plot visualizes the frequency of both predefined (known) and
        dynamically discovered (unknown) hopping events.

        Args:
            title (str | None, optional):
                A custom title for the plot. Defaults to None.
            save (bool, optional):
                If True, saves the figure to a file. Defaults to True.
            filename (str, optional):
                Filename for the saved plot. Defaults to 'counts.png'.
            dpi (int, optional):
                Resolution in dots per inch for the saved figure. Defaults to 300.

        Returns:
            None: This method displays a plot and optionally saves it to a file.
        """
        if self.counts is None or self.counts_unknown is None:
            raise RuntimeError("Path counts have not been calculated. "
                               "Please run the .calculate() method first.")

        name_all = self.site.path_name + [p['name'] for p in self.path_unknown]
        counts_all = np.append(np.sum(self.counts, axis=0), 
                               np.sum(self.counts_unknown, axis=0))
       
        if len(name_all) == 0:
            print("Warning: No path count data available to plot.")
            return

        fig, ax = plt.subplots(figsize=(10, 6))
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.2)
            
        x_pos = np.arange(len(name_all))
        bars = ax.bar(x_pos, counts_all, color='steelblue', edgecolor='k', alpha=0.8)

        ax.set_ylabel('Total Counts', fontsize=13)
        ax.set_xlabel('Path Name', fontsize=13)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(name_all, rotation=45, ha="right", rotation_mode="anchor")

        if title is not None:
            ax.set_title(title, fontsize=12, pad=10)
            
        ax.bar_label(bars, fmt='%d', padding=3, fontsize=10)
        ax.yaxis.grid(True, linestyle='--', alpha=0.7)
        ax.set_axisbelow(True)
        ax.set_ylim(top=ax.get_ylim()[1] * 1.1)
        
        fig.tight_layout()
        if save:
            fig.savefig(filename, dpi=dpi, bbox_inches="tight")
        
        plt.show()
        plt.close(fig)
              
    def summary(self) -> None:
        """Prints a comprehensive, formatted summary of the analysis results.

        The summary includes input parameters and all key calculated physical
        properties like diffusivity, correlation factor, and residence time.
        It automatically calls `.calculate()` if it has not been run yet.

        Returns:
            None: This method prints a summary to the console.
        """
        if not self.analyzer or not self.encounter:
            raise RuntimeError("Primary analysis objects (analyzer, encounter) are missing.")
        if not self.calculate_is_done:
            self.calculate()
            
        print("=" * 60)
        print("Summary for Trajectory dataset")
        print(f"  - Path to TRAJ file   : {self.path_traj}")
        print(f"  - Lattice structure   : {self.site.path_structure}")
        print(f"  - t_interval          : {self.t_interval:.3f} ps ({self.frame_interval} frames)")
        print(f"  - Temperatures (K)    : {self.temperatures.tolist()}")
        print(f"  - Num. of TRAJ files  : [1]")
        print("=" * 60)
        
        print("\n" + "="*16 + " Temperature-Dependent Data " + "="*16)
        headers = ["Temp (K)", "D (m2/s)", "D_rand (m2/s)", "f", "tau (ps)"]
        table_data = zip([self.temperature], [self.D], [self.D_rand], [self.f], [self.tau])
        formats = [".1f", ".3e", ".3e", ".4f", ".4f"]
        table_data = []
        for row in zip([self.temperature], [self.D], [self.D_rand], [self.f], [self.tau]):
            formatted_row = [f"{value:{fmt}}" for value, fmt in zip(row, formats)]
            table_data.append(formatted_row)
        table = tabulate(table_data, headers=headers, 
                            tablefmt="simple", stralign='left', numalign='left')
        print(table)
        print("=" * 60)
        
        print("\n" + "="*17 + " Final Fitted Parameters " + "="*18)
        print(f"Diffusivity (D):")
        print(f"  - Ea          : - ")
        print(f"  - D0          : - ")
        print(f"  - R-squared   : -")
        print(f"Random Walk Diffusivity (D_rand):")
        print(f"  - Ea          : -")
        print(f"  - D0          : -")
        print(f"  - R-squared   : -")
        print(f"Correlation Factor (f):")
        print(f"  - Ea          : -")
        print(f"  - f0          : -")
        print(f"  - R-squared   : -")
        print(f"Residence Time (tau):")
        print(f"  - Ea (fixed)  : -")
        print(f"  - tau0        : -")
        print(f"  - R-squared   : -")
        print(f"Effective Hopping Distance (a) : -")
        print("=" * 60)
        
    def show_hopping_history(self) -> None:
        """Prints a detailed, time-ordered table of the hopping sequence for each vacancy.

        Returns:
            None: This method prints a table to the console.
        """
        for i in range(self.num_vacancies):
            print("=" * 116)
            print(" "*43 + f"Hopping Sequence of Vacancy {i}")
            print("=" * 116)
            
            header = ['Num', 'Time (ps)', 'Path', 'a (Ang)', 
                      'Initial Site (Fractional Coordinate)', 
                      'Final Site (Fractional Coordinate)']
            data = [
                [
                    f"{j+1}",
                    f"{path['step'] * self.t_interval:.2f}",
                    f"{path['name']}",
                    f"{path['distance']:.5f}",
                    f"{path['site_init']} [{', '.join(f'{x:.5f}' for x in self.site.lattice_sites[path['index_init']]['coord'])}]",
                    f"{path['site_final']} [{', '.join(f'{x:.5f}' for x in self.site.lattice_sites[path['index_final']]['coord'])}]"
                ] for j, path in enumerate(self.hopping_history[i])
            ]
            print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
            print("=" * 116 + "\n")
            
    def show_hopping_paths(self) -> None:
        """Prints a summary table of all observed hopping paths and their counts.

        This includes both the paths pre-defined in the `Site` object and any
        new "unknown" paths discovered during the trajectory analysis.

        Returns:
            None: This method prints a table to the console.
        """
        
        path_info = []
        for i, path in enumerate(self.site.path):
            p = {
                'name': path['name'],
                'a': path['distance'],
                'z': path['z'],
                'count': np.sum(self.counts, axis=0)[i],
                'site_init': f"{path['site_init']} [{', '.join(f'{x:.5f}' for x in path['coord_init'])}]",
                'site_final': f"{path['site_final']} [{', '.join(f'{x:.5f}' for x in path['coord_final'])}]"
            }
            path_info.append(p)
        
        for i, path in enumerate(self.path_unknown):
            p = {
                'name': path['name'],
                'a': path['distance'],
                'z': '-',
                'count': int(np.sum(self.counts_unknown, axis=0)[i]),
                'site_init': f"{path['site_init']} [{', '.join(f'{x:.5f}' for x in path['coord_init'])}]",
                'site_final': f"{path['site_final']} [{', '.join(f'{x:.5f}' for x in path['coord_final'])}]"
            }
            path_info.append(p)
            
        print("=" * 116)
        print(" " * 46 + "Hopping Path Information")
        print("=" * 116)
        header = ['Name',
                  'a (Ang)',
                  'z',
                  'Count',
                  'Initial Site (Fractional Coordinate)',
                  'Final Site (Fractional Coordinate)']
        data = [
            [
                p['name'],
                p['a'],
                p['z'],
                p['count'],
                p['site_init'],
                p['site_final']
            ] for p in path_info
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print("=" * 116 + "\n")
       
    def save_trajectory(self, filename: str = "trajectory.json") -> None:
        """Saves the unwrapped Cartesian trajectories of vacancies to a JSON file.

        The output JSON contains simulation metadata and the time-series coordinates
        for each vacancy, keyed as 'Vacancy0', 'Vacancy1', etc.

        Args:
            filename (str, optional):
                The name of the output JSON file. Defaults to 'trajectory.json'.

        Returns:
            None: This method writes a file to disk.
        """
        trajectories = {i:[] for i in range(self.num_vacancies)}
        for coords in self.unwrapped_vacancy_trajectory_coord_cart.values():
            for i, coord in enumerate(coords):
                trajectories[i].append(coord.tolist())
                
        contents = {
            'traj': str(Path(self.path_traj).resolve()),
            'symbol': self.symbol,
            't_interval': self.t_interval,
            'temperature': self.temperature,
            'num_vacancies': self.num_vacancies,
            'lattice': self.lattice_parameter.tolist()
        }
        for i in range(self.num_vacancies):
            contents[f'Vacancy{i}'] = trajectories[i]
            
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(contents, f, indent=2)

        if self.verbose:
            print(f"Trajectoris saved to '{filename}'.")

class CalculatorEnsemble(TrajectoryBundle):
    """Orchestrates the full analysis pipeline for an ensemble of trajectories.

    This class is the main user interface for multi-temperature diffusion analysis.
    It inherits from `TrajectoryBundle`, using it to automatically discover, validate,
    and group HDF5 trajectory files by temperature from a given path.

    The primary workflow involves two steps:

    1.  **Initialization**: Create an instance of the class, which prepares the
        dataset.

    2.  **Calculation**: Call the `.calculate()` method to run the entire
        computational pipeline. This process analyzes each trajectory in parallel,
        aggregates the results, and performs Arrhenius fits to determine activation
        energies and pre-exponential factors.

    Once calculated, results can be viewed with the `.summary()` method, visualized
    with a suite of plotting methods (e.g., `.plot_D()`), or saved to a file.

    Args:
        path_traj (str):
            Path to a single HDF5 file or a root directory to search for files.
        site (Site):
            An initialized `Site` object with lattice and path information.
        t_interval (float):
            The time interval in picoseconds (ps) for coarse-graining the analysis.
        prefix (str, optional):
            Prefix of trajectory files to search for (used in directory scans).
            Defaults to "TRAJ".
        depth (int, optional):
            Maximum directory depth for file searches. Defaults to 2.
        use_incomplete_encounter (bool, optional):
            If True, incomplete encounters are included in the analysis.
            Defaults to True.
        eps (float, optional):
            A tolerance for floating-point comparisons. Defaults to 1.0e-3.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Attributes:
        temperatures (numpy.ndarray):
            Array of unique temperatures (K) found in the analysis.

        The following attributes are populated after calling the `.calculate()` method:

        D (numpy.ndarray):
            Temperature-dependent vacancy diffusivity (m²/s).
        Ea_D (float):
            Activation energy (eV) for diffusivity from Arrhenius fit.
            This value corresponds to the diffusion barrier.
        D0 (float):
            Pre-exponential factor (m²/s) for diffusivity.
        D_rand (numpy.ndarray):
            Temperature-dependent random-walk diffusivity (m²/s).
        Ea_D_rand (float):
            Activation energy (eV) for random-walk diffusivity.
            This value corresponds to the hopping barrier.
        D_rand0 (float):
            Pre-exponential factor (m²/s) for random-walk diffusivity.
        f (numpy.ndarray):
            Temperature-dependent correlation factor.
        Ea_f (float):
            Activation energy (eV) for the correlation factor.
        f0 (float):
            Pre-exponential factor (ps) for the correlation factor
        tau (numpy.ndarray):
            Temperature-dependent average vacancy residence time (ps).
        Ea_tau (float):
            Activation energy (eV) for residence time.
            This value is the same with the Ea_D_rand.
        tau0 (float):
            Pre-exponential factor (ps) for the vacancy residence time.
        a (numpy.ndarray):
            Effective hopping distance (Å).
        calculators (list[CalculatorSingle]):
            A list of the successfully completed `CalculatorSingle` instances,
            ordered corresponding to the `all_traj_paths` attribute.

        The following attributes are populated after calling the `.calculate_attempt_frequency()` method:

        z (numpy.ndarray):
            Temperature-dependent effective coordination number.
        nu (numpy.ndarray):
            Temperature-dependent effective attempt frequency (THz).
        nu_path (numpy.ndarray):
            Path-wise attempt frequency (THz), with shape
            (n_temperatures, n_paths).
            
        The following attributes are populated after calling the `.decompose_diffusivity()` method:
        
        Dx (numpy.ndarray):
            Temperature-dependent vacancy diffusivity in the x-direction (m²/s).
        D0_x (float):
            Pre-exponential factor for Dx (m²/s) from Arrhenius fit.
        Ea_x (float):
            Activation energy for Dx (eV) from Arrhenius fit.
        Dy (numpy.ndarray):
            Temperature-dependent vacancy diffusivity in the y-direction (m²/s).
        D0_y (float):
            Pre-exponential factor for Dy (m²/s) from Arrhenius fit.
        Ea_y (float):
            Activation energy for Dy (eV) from Arrhenius fit.
        Dz (numpy.ndarray):
            Temperature-dependent vacancy diffusivity in the z-direction (m²/s).
        D0_z (float):
            Pre-exponential factor for Dz (m²/s) from Arrhenius fit.
        Ea_z (float):
            Activation energy for Dz (eV) from Arrhenius fit.

    Raises:
        FileNotFoundError:
            If `path_traj` does not exist or no valid files are found.
        ValueError:
            If inconsistent simulation parameters are found among trajectories.

    Examples:
        >>> # 1. Initialize the Site and CalculatorEnsemble objects
        >>> site_info = Site("POSCAR", symbol="O")
        >>> ensemble = CalculatorEnsemble(
        ...     path_traj="path/to/trajectories/",
        ...     site=site_info,
        ...     t_interval=0.1
        ... )
        >>>
        >>> # 2. Run the entire parallel analysis pipeline
        >>> ensemble.calculate(n_jobs=-1)
        >>>
        >>> # 3. View and visualize the main results
        >>> ensemble.summary()
        >>> ensemble.plot_D(filename="diffusivity_plot.png")
        >>>
        >>> # 4. (Optional) Calculate attempt frequencies with NEB data
        >>> ensemble.calculate_attempt_frequency(neb_csv="path/to/neb_data.csv")
        >>> ensemble.plot_nu()
        >>>
        >>> # 5. (Optional) Decompose diffusivity and plot results
        >>> ensemble.decompose_diffusivity()
        >>> ensemble.plot_D_xyz()
    """
    def __init__(self,
                 path_traj: str,
                 site: Site,
                 t_interval: float,
                 prefix: str = "TRAJ",
                 depth: int = 2,
                 use_incomplete_encounter: bool = True,
                 eps: float = 1.0e-3,
                 verbose: bool = True):
        
        self.kb = 8.61733326e-5
        # self.cmap = plt.get_cmap("Set1")
        
        self.site = site
        self.t_interval = t_interval
        self.use_incomplete_encounter = use_incomplete_encounter
        self.eps = eps
        self.symbol = self.site.symbol
        
        super().__init__(path_traj, 
                         self.symbol, 
                         prefix=prefix, 
                         depth=depth, 
                         eps=self.eps, 
                         verbose=False)
        
        self.temperatures = np.array(self.temperatures, dtype=np.float64)
        self.all_traj_paths = [path for temp_paths in self.traj for path in temp_paths]
        self.verbose = verbose
        
        self.calculators = None
        self.num_vacancies = None
        self.counts = None  # (num_vacancies, num_paths)
        
        self.path_unknown = None
        self.counts_unknown = None
        
        # correlation factor
        self.f = None       # (n_temp,)
        self.f0 = None
        self.Ea_f = None
        self.f_R2 = None
        
        # random walk diffusivity
        self.D_rand = None  # (n_temp,)
        self.D_rand0 = None
        self.Ea_D_rand = None
        self.D_rand_R2 = None
        
        # diffusivity
        self.D = None       # (n_temp,)
        self.D0 = None
        self.Ea_D = None
        self.D_R2 = None
        
        # residence time
        self.tau = None     # (n_temp,)
        self.tau0 = None
        self.Ea_tau = None
        self.tau_R2 = None
        
        # effective hopping distance
        self.a = None
       
        # <z>
        self.z_mean = None
        
        # extra properties for attempt frequency calc.
        self.a_path = None      # (n_path,)
        self.z_path = None      # (n_path,)
        self.P_site = None      # (n_temp, n_site)
        self.times_site = None  # (n_temp, n_site)
        self.counts_hop = None  # (n_temp, n_path)
        
        # attirbutes populated after `.calculate_attempt_frequency()` is called.
        self.attempt_frequency = None
        self.z = None       # (n_temp,)
        self.nu = None      # (n_temp,)
        self.nu_path = None # (n_temp, n_path)
        
        # Decomposed diffusivity and its fits
        self.msd = None
        self.Dx, self.Dy, self.Dz = None, None, None
        self.D0_x, self.Ea_x, self.Dx_R2 = None, None, None
        self.D0_y, self.Ea_y, self.Dy_R2 = None, None, None
        self.D0_z, self.Ea_z, self.Dz_R2 = None, None, None
        
    @monitor_performance
    def calculate(self, n_jobs: int = -1, verbose=True) -> None:
        """Runs the full analysis pipeline on all trajectories in parallel.

        This method performs the following steps:
        1. Runs the `CalculatorSingle` analysis for each trajectory file in parallel.
        2. Gathers and aggregates the results by temperature.
        3. Calculates temperature-averaged physical properties (D, f, tau, etc.).
        4. Performs Arrhenius fits for these properties if more than one temperature
        is available.
        5. Calculates the effective hopping distance from the fit results.

        Args:
            n_jobs (int, optional):
                Number of CPU cores for parallel analysis. -1 uses all available
                cores. Defaults to -1.
            verbose (bool, optional):
                Verbosity flag for the performance monitor decorator.
                Defaults to True.
        
        Returns:
            None: This method populates the instance's attributes with the
                analysis results and does not return any value.
        """
        
        results = Parallel(n_jobs=n_jobs)(
            delayed(self._run_single_task)(path_traj) 
            for path_traj in tqdm(self.all_traj_paths,
                                  desc=f'Analyze Trajectory',
                                  bar_format='{l_bar}{bar:30}{r_bar}',
                                  ascii=True) 
        )
        
        # Successfully terminated calculators
        successful_results = [res for res in results if res[1] is not None]

        # Sort the results in order of traj
        path_order_map = {path: i for i, path in enumerate(self.all_traj_paths)}
        successful_results = sorted(successful_results, key=lambda res: path_order_map[res[0]])
        self.calculators = [calc for _, calc in successful_results]
        if verbose:
            print(f"\nAnalysis complete: {len(successful_results)} successful, " +
                f"{len(results) - len(successful_results)} failed.")
        
        # Count the number of calc at each temperature
        calc_temps = np.array([calc.temperature for calc in self.calculators])
        self.num_calc_temp = np.bincount(
            np.argmin(np.abs(calc_temps[:, np.newaxis] - self.temperatures), axis=1), 
            minlength=len(self.temperatures)
        )
        self.index_calc_temp = np.concatenate(([0], np.cumsum(self.num_calc_temp)))
        
        self.num_vacancies = self.calculators[0].num_vacancies
        self.t_interval = self.calculators[0].t_interval
        self.frame_interval = self.calculators[0].frame_interval
        self.counts = np.sum(np.array([c.counts for c in self.calculators]), axis=0) # (num_vacancies, num_paths)

        self._gather_unknown_paths()
        self._get_correlation_factor()
        self._get_random_walk_diffusivity()
        self._get_diffusivity()
        self._get_residence_time()
        self._get_effective_hopping_distance()
        self._get_mean_number_of_equivalent_paths()
        self._get_extra_properties()
        
        if len(self.temperatures) >= 2:
            self._fit_correlation_factor()
            self._fit_random_walk_diffusivity()
            self._fit_diffusivity()
            self._fit_residence_time()
            
    
    # ===================================================================
    # Internal Helper Methods
    # ===================================================================
         
    def _run_single_task(self, path_traj) -> Tuple[str, Optional['CalculatorSingle']]:
        """
        Worker function to run analysis on a single trajectory file.
        
        Args:
            traj_path (str): The file path of the trajectory to analyze.

        Returns:
            A tuple containing the input path and the resulting CalculatorSingle
            object, or None if the analysis failed.
        """
        try:
            calc = CalculatorSingle(
                path_traj=path_traj,
                site=self.site,
                t_interval=self.t_interval,
                eps=self.eps,
                use_incomplete_encounter=self.use_incomplete_encounter,
                verbose=False
            )
            return (path_traj, calc)
        
        except Exception as e:
            print(f"\nWarning: Failed to analyze '{path_traj}'. Reason: {e}. Skipped.")
            return (path_traj, None)

    def _get_path_key(self, path):
        """Generates a unique, hashable key for a given migration path.

        This helper method creates a tuple that can be used as a dictionary key to
        uniquely identify a path. It combines the initial and final sites with a
        discretized version of the path distance. Discretizing the distance by
        rounding it relative to a small epsilon (`self.eps`) ensures that paths
        with floating-point distances that are very close are treated as identical.

        Args:
            path (dict): 
                A dictionary representing a migration path, requiring
                'site_init', 'site_final', and 'distance' keys.

        Returns:
            tuple: A hashable tuple `(site_init, site_final, discretized_distance)`
                used for unique path identification.
        """
        discretized_distance = round(path['distance'] / self.eps)
        return (path['site_init'], path['site_final'], discretized_distance)

    def _gather_unknown_paths(self):
        """Consolidates and re-labels all 'unknown' paths from multiple calculators.

        This method performs three main tasks to create a unified list of previously
        unidentified migration paths:

        1.  **Path Aggregation**: It iterates through all `path_unknown` lists from
            each calculator. Using a dictionary-based registry for efficient O(1)
            lookups, it identifies unique unknown paths (based on sites and
            distance) and counts their total occurrences across all calculators.

        2.  **Final List Population**: It processes the registry to populate the
            instance's final attributes: `self.path_unknown` (a list of unique
            path dictionaries) and `self.count_unknown` (a corresponding list of
            occurrence counts).

        3.  **History Re-labeling**: It performs a final pass through the
            `hopping_history` of each calculator. Any history event originally
            marked as 'unknown' is updated with the new canonical name
            (e.g., 'unknown1', 'unknown2') assigned during aggregation.

        This method modifies `self.path_unknown`, `self.count_unknown`, and the
        `hopping_history` within each calculator object in place.

        Returns:
            None
        """
        path_registry = {}
        
        for c in self.calculators:
            if not hasattr(c, 'path_unknown') or not c.path_unknown:
                continue

            for path in c.path_unknown:
                key = self._get_path_key(path)
                
                if key in path_registry:
                    path_registry[key]['count'] += 1.0
                else:
                    new_name = f"unknown{len(path_registry) + 1}"
                    path['name'] = new_name
                    path_registry[key] = {
                        'path': path,
                        'count': 1.0
                    }
        
        self.path_unknown = []
        self.counts_unknown = []
        for reg_item in path_registry.values():
            self.path_unknown.append(reg_item['path'])
            self.counts_unknown.append(reg_item['count'])

        for c in self.calculators:
            if not hasattr(c, 'hopping_history'):
                continue

            for history_vacancy in c.hopping_history:
                for history in history_vacancy:
                    if history['name'].startswith('unknown'):
                        key = self._get_path_key(history)
                        if key in path_registry:
                            history['name'] = path_registry[key]['path']['name']

    def _get_correlation_factor(self):
        """Calculates the temperature-dependent correlation factor."""
        
        if self.calculators is None:
            raise RuntimeError("Please call the .calculate() method first.") 

        self.f_ind = np.array([c.f for c in self.calculators], dtype=np.float64)
        num_enc_all = np.array([c.encounter.num_encounter for c in self.calculators])
        msd_enc_all = np.array([c.encounter.msd for c in self.calculators])
        
        msd_enc_rand_all = np.array([
            np.sum(c.encounter.path_distance**2 * c.encounter.path_count)
            if c.f is not None else np.nan
            for c in self.calculators
        ])
        
        self.f_avg = [np.nanmean(self.f_ind[start:end]) 
                      for start, end in zip(self.index_calc_temp[:-1], self.index_calc_temp[1:])]

        self.f = []
        for start, end in zip(self.index_calc_temp[:-1], self.index_calc_temp[1:]):
            num_enc_slice = num_enc_all[start:end]
            msd_enc_slice = msd_enc_all[start:end]
            msd_enc_rand_slice = msd_enc_rand_all[start:end]
            
            valid_mask = ~np.isnan(msd_enc_rand_slice)
            if not np.any(valid_mask):
                self.f.append(np.nan)
                continue

            num_enc_valid = num_enc_slice[valid_mask]
            total_num_enc = np.sum(num_enc_valid)
            
            if total_num_enc == 0:
                self.f.append(np.nan)
                continue
                
            msd_total = np.sum(msd_enc_slice[valid_mask] * num_enc_valid) / total_num_enc
            msd_rand_total = np.sum(msd_enc_rand_slice[valid_mask]) / total_num_enc
            
            self.f.append(msd_total / msd_rand_total)
        self.f = np.array(self.f, dtype=np.float64)
    
    def _get_random_walk_diffusivity(self):
        """Calculates the temperature-dependent random walk diffusivity."""
        if self.calculators is None:
            raise RuntimeError("Please call the .calculate() method first.")
            
        total_time_all = np.array([c.t_interval * c.num_steps for c in self.calculators])
        msd_rand_all = np.array([c.msd_rand for c in self.calculators])
        
        self.D_rand = []
        for start, end in zip(self.index_calc_temp[:-1], self.index_calc_temp[1:]):
            t_i = np.sum(total_time_all[start:end])
            msd_rand_i = np.sum(msd_rand_all[start:end])
            if t_i == 0:
                self.D_rand.append(np.nan)
            else:
                self.D_rand.append(msd_rand_i / (6 * t_i))
                
        self.D_rand = np.array(self.D_rand) * 1e-8
    
    def _get_diffusivity(self):
        """Calculates the temperature-dependent tracer diffusivity."""
        if self.D_rand is None:
            raise RuntimeError(
                "D_rand is not calculated. Please call '_get_random_walk_diffusivity()' first."
            )
        if self.f is None:
            raise RuntimeError(
                "Correlation factor (f) is not calculated. Please call '_get_correlation_factor()' first."
            )
            
        d_rand_arr = np.asarray(self.D_rand)
        f_arr = np.asarray(self.f)
        
        if d_rand_arr.shape != f_arr.shape:
            raise ValueError(
                f"Shape mismatch: D_rand shape {d_rand_arr.shape} does not match f shape {f_arr.shape}."
            )

        self.D = d_rand_arr * f_arr
    
    def _get_residence_time(self):
        """Calculates the temperature-dependent residence time."""
        if self.calculators is None:
            raise RuntimeError("Please call the .calculate() method first.")
        
        total_time_all = np.array([c.t_interval * c.num_steps for c in self.calculators])
        total_jumps_all = np.array(
            [(np.sum(c.counts) + np.sum(c.counts_unknown)) / c.num_vacancies 
             for c in self.calculators]
        )

        self.tau = []
        for start, end in zip(self.index_calc_temp[:-1], self.index_calc_temp[1:]):
            t_i = np.sum(total_time_all[start:end])
            count_i = np.sum(total_jumps_all[start:end])
            
            if count_i == 0:
                self.tau.append(np.nan)
            else:
                self.tau.append(t_i / count_i)
                
        self.tau = np.array(self.tau, dtype=np.float64)

    def _get_effective_hopping_distance(self):
        """Calculates the effective hopping distance from fit parameters."""
        if self.D_rand is None: self._get_random_walk_diffusivity()
        if self.tau is None: self._get_residence_time()
        
        term = 6 * np.asarray(self.D_rand) * np.asarray(self.tau)
        self.a = np.full_like(term, np.nan)
        
        valid_mask = term >= 0

        if not np.all(valid_mask):
            invalid_temps = self.temperatures[~valid_mask]
            print(
                f"Warning: The term (6 * D_rand * tau) was negative for temperatures "
                f"{invalid_temps.tolist()}. 'a' is set to NaN for these points."
            )
            
        self.a[valid_mask] = np.sqrt(term[valid_mask]) * 1e4 # convert to Å
        
    def _get_mean_number_of_equivalent_paths(self):
        """Calculates the temperature-dependent mean number of equivalent paths."""
        if self.calculators is None:
            raise RuntimeError("Please call the .calculate() method first.")

        z = np.array([path['z'] for path in self.site.path], dtype=np.float64)
        all_counts = np.array([np.sum(c.counts, axis=0) for c in self.calculators])

        self.z_mean = []
        for start, end in zip(self.index_calc_temp[:-1], self.index_calc_temp[1:]):
            counts = np.sum(all_counts[start:end], axis=0)

            total_jumps = np.sum(counts)
            if total_jumps == 0:
                self.z_mean.append(np.nan)
                continue
            denominator = np.sum(np.divide(counts, z, where=z!=0))
            
            if denominator == 0:
                self.z_mean.append(np.nan)
            else:
                self.z_mean.append(total_jumps / denominator)
        self.z_mean = np.array(self.z_mean, dtype=np.float64) 
        
    def _get_extra_properties(self):
        """Calculates P_site, times, counts_hop, a_path, z_path"""
        name_to_type_index = {name: i for i, name in enumerate(self.site.site_name)}
        site_index_to_type_map = np.array(
            [name_to_type_index[site['site']] for site in self.site.lattice_sites]
        )
        
        self.counts_hop = []
        self.times_site = []
        self.P_site = []

        for start, end in zip(self.index_calc_temp[:-1], self.index_calc_temp[1:]):
            
            counts_list_temp = []
            time_temp = 0.0
            indices_list_temp = []

            for calc in self.calculators[start:end]:
                counts_list_temp.append(np.sum(calc.counts, axis=0))
                indices_list_temp.append(np.concatenate(list(calc.vacancy_trajectory_index.values())))

            counts_hop_temp = np.sum(np.array(counts_list_temp, dtype=np.float64), axis=0)
            self.counts_hop.append(counts_hop_temp.tolist())
            
            if indices_list_temp:
                all_indices_in_temp = np.concatenate(indices_list_temp)
                all_type_indices = site_index_to_type_map[all_indices_in_temp]
                count_site = np.bincount(all_type_indices, minlength=len(self.site.site_name))
                time_temp = count_site * self.t_interval
                self.times_site.append(time_temp.tolist())
                
                total_counts = np.sum(count_site)
                if total_counts > 0:
                    P_site_temp = count_site.astype(np.float64) / total_counts
                else:
                    P_site_temp = np.zeros_like(self.site.site_name, dtype=np.float64)
                self.P_site.append(P_site_temp.tolist())
            else:
                self.P_site.append(np.zeros_like(self.site.site_name, dtype=np.float64).tolist())
                
        self.a_path, self.z_path = [], []
        for path in self.site.path:
            self.a_path.append(path['distance'])
            self.z_path.append(path['z'])
        
    def _fit_correlation_factor(self):
        """Performs an Arrhenius fit on the correlation factor data."""
        self.Ea_f, self.f0, self.f_R2 = self._Arrhenius_fit(
            self.temperatures, self.f, self.kb
        )
    
    def _fit_random_walk_diffusivity(self):
        """Performs an Arrhenius fit on the random walk diffusivity data."""
        self.Ea_D_rand, self.D_rand0, self.D_rand_R2 = self._Arrhenius_fit(
            self.temperatures, self.D_rand, self.kb
        )
    
    def _fit_diffusivity(self):
        """Performs an Arrhenius fit on the diffusivity data."""
        self.Ea_D, self.D0, self.D_R2 = self._Arrhenius_fit(
            self.temperatures, self.D, self.kb
        )
    
    def _fit_residence_time(self):
        """Fits the residence time data to find the pre-exponential factor tau0."""
        if self.tau is None:
            self._get_residence_time()
            
        if self.Ea_D_rand is None:
            self._get_random_walk_diffusivity()
            self._fit_random_walk_diffusivity()
        self.Ea_tau = self.Ea_D_rand
        
        valid_mask = ~np.isnan(self.tau)
        if np.sum(valid_mask) < 2:
            self.tau0 = np.nan
            self.tau_R2 = np.nan
            return

        tau_valid = self.tau[valid_mask]
        temps_valid = self.temperatures[valid_mask]

        error_tau = lambda tau0: np.linalg.norm(
            tau_valid - tau0 * np.exp(self.Ea_tau / (self.kb * temps_valid))
        )

        tau0_opt = minimize_scalar(error_tau)
        self.tau0 = tau0_opt.x # ps
    
        ss_residual = tau0_opt.fun ** 2
        ss_total = np.sum((tau_valid - np.mean(tau_valid))**2)

        if ss_total == 0:
            self.tau_R2 = 1.0 if ss_residual == 0 else 0.0
        else:
            self.tau_R2 = 1 - (ss_residual / ss_total)

    @staticmethod
    def _Arrhenius_fit(temperatures: np.ndarray, 
                       data: np.ndarray, 
                       kb: float) -> Tuple[float, float, float]:
        """
        Fits data to the Arrhenius equation and returns fitting parameters.

        Args:
            temperatures (np.ndarray): Array of temperatures (in K).
            data (np.ndarray): Array of data to be fitted (e.g., diffusivity).
            kb (float): Boltzmann constant in the desired units (e.g., eV/K).

        Returns:
            Tuple[float, float, float]: A tuple containing:
                - activation_energy (float)
                - pre_exponential_factor (float)
                - r_squared (float)
        """
        if len(temperatures) <= 1:
            raise ValueError("At least two temperature points are required.")
        
        valid_mask = data > 0
        if np.sum(valid_mask) < 2:
            raise ValueError("At least two valid data points (> 0) are required.")
        
        temps_valid = temperatures[valid_mask]
        data_valid = data[valid_mask]
        
        x = 1 / temps_valid
        y = np.log(data_valid)
        slope, intercept = np.polyfit(x, y, deg=1)
        
        activation_energy = -slope * kb
        pre_exponential_factor = np.exp(intercept)
        
        # Calculate R-squared for goodness of fit
        y_predicted = slope * x + intercept
        ss_residual = np.sum((y - y_predicted)**2)
        ss_total = np.sum((y - np.mean(y))**2)
        
        if ss_total == 0:
            r_squared = 1.0 if ss_residual == 0 else 0.0
        else:
            r_squared = 1 - (ss_residual / ss_total)
            
        return activation_energy, pre_exponential_factor, r_squared

    # ===================================================================
    # Plotting Methods
    # ===================================================================
     
    def _create_arrhenius_plot(self, 
                               data: np.ndarray, 
                               temperatures: np.ndarray, 
                               slope: float, 
                               intercept: float,
                               ylabel: str, 
                               title: Optional[str] = None, 
                               figsize: Tuple[float, float] = (7, 6),
                               ax: Optional[plt.Axes] = None) -> Tuple[plt.Figure, plt.Axes]:
        """
        Generic helper function to create a styled Arrhenius plot.
        Can draw on a provided Axes object or create a new figure.
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()
            
        valid_mask = ~np.isnan(data) & (data > 0)
        temps_valid = temperatures[valid_mask]
        data_valid = data[valid_mask]

        # fig, ax = plt.subplots(figsize=figsize)  <-- 이 라인이 문제였습니다. 삭제합니다!

        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.2)

        x_points = 1000 / temps_valid
        y_points = np.log(data_valid)
        
        cmap = plt.get_cmap("viridis", len(self.temperatures))
        temp_color_map = {temp: cmap(i) for i, temp in enumerate(self.temperatures)}

        for i in range(len(temps_valid)):
            temp = temps_valid[i]
            ax.scatter(x_points[i], y_points[i], color=temp_color_map[temp],
                       marker="o", s=100, label=f"{temp:.0f} K", edgecolor='k', alpha=0.8)

        if len(x_points) > 1 and not np.isnan(slope):
            x_fit = np.linspace(np.min(x_points), np.max(x_points), 100)
            ax.plot(x_fit, slope * x_fit + intercept, 'k--', linewidth=1.5)
        
        # 범례가 중복되지 않도록 기존 범례를 지우고 다시 생성
        if ax.get_legend() is not None:
            ax.get_legend().remove()
        ax.legend(title='Temperature', fontsize=11, title_fontsize=12)

        ax.set_xlabel('1000 / T (K⁻¹)', fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        if title:
            ax.set_title(title, fontsize=13, pad=10)
        
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # fig.tight_layout() # 이 함수를 호출하는 상위 함수에서 레이아웃을 조절하므로 주석 처리
        
        return fig, ax
             
    def plot_D_rand(self, 
                    title: str | None = "Random Walk Diffusivity (Arrhenius Plot)",
                    disp: bool = True,
                    save: bool = True, 
                    filename: str = "D_rand.png", 
                    dpi: int = 300) -> None:
        """Generates an Arrhenius plot for the random-walk diffusivity (D_rand).

        Args:
            title (str | None, optional): A custom title for the plot.
            disp (bool, optional): If True, displays the plot. Defaults to True.
            save (bool, optional): If True, saves the plot to a file. Defaults to True.
            filename (str, optional): Filename for the saved plot.
            dpi (int, optional): Resolution for the saved figure.
        """
        if len(self.temperatures) < 2:
            raise ValueError("At least two temperature points are required for an Arrhenius plot.")
        if self.D_rand is None or self.Ea_D_rand is None:
             raise RuntimeError("Please run analysis and fitting for 'D_rand' first.")
         
        slope = -self.Ea_D_rand / (self.kb * 1000) 
        intercept = np.log(self.D_rand0)

        fig, ax = self._create_arrhenius_plot(
            data=self.D_rand, temperatures=self.temperatures,
            slope=slope, intercept=intercept,
            ylabel=r'ln($D_{rand}$) (m$^2$/s)', title=title
        )
        
        text_str = (f'$E_a = {self.Ea_D_rand:.2f}$ eV\n'
                    f'$D_0 = {self.D_rand0:.2e}$ m$^2$/s\n'
                    f'$R^2 = {self.D_rand_R2:.3f}$')
        
        ax.text(0.05, 0.05, text_str, transform=ax.transAxes, fontsize=10,
                verticalalignment='bottom', horizontalalignment='left',
                bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.7))
        
        if save: fig.savefig(filename, dpi=dpi, bbox_inches="tight")
        if disp: plt.show()
        plt.close(fig)
        
    def plot_f(self, 
               title: str | None = "Correlation Factor vs. Temperature",
               disp: bool = True,
               save: bool = True, 
               filename: str = 'f.png', 
               dpi: int = 300) -> None:
        """Plots the correlation factor (f) as a function of temperature.

        Args:
            title (str | None, optional): A custom title for the plot.
            disp (bool, optional): If True, displays the plot. Defaults to True.
            save (bool, optional): If True, saves the plot to a file. Defaults to True.
            filename (str, optional): Filename for the saved plot.
            dpi (int, optional): Resolution for the saved figure.
        """
        if self.f is None:
            raise RuntimeError("Correlation factor has not been calculated. Please run .calculate() first.")

        valid_mask = ~np.isnan(self.f)
        temps_valid = self.temperatures[valid_mask]
        f_valid = self.f[valid_mask]

        fig, ax = plt.subplots(figsize=(7, 6))
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.2)
        
        cmap = plt.get_cmap("viridis", len(self.temperatures))
        temp_color_map = {temp: cmap(i) for i, temp in enumerate(self.temperatures)}

        plt.plot(temps_valid, f_valid, linestyle='--', color='k')    
        
        for i in range(len(temps_valid)):
            temp = temps_valid[i]
            ax.scatter(temp, f_valid[i], color=temp_color_map[temp],
                       marker='s', s=100, label=f"{temp:.0f} K", edgecolor='k', alpha=0.8)
        
        ax.set_ylim([0, 1])
        ax.set_xlabel('Temperature (K)', fontsize=12)
        ax.set_ylabel('Correlation Factor, $f$', fontsize=12)
        if title: ax.set_title(title, fontsize=14, pad=10)
        ax.legend(title='Temperature', fontsize=11, title_fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.7)
        
        if self.Ea_f is not None:
             text_str = (f'$E_a = {self.Ea_f:.3f}$ eV\n'
                         f'$f_0 = {self.f0:.3f}$\n'
                         f'$R^2 = {self.f_R2:.4f}$')
             ax.text(0.05, 0.05, text_str, transform=ax.transAxes, fontsize=10,
                     verticalalignment='bottom', horizontalalignment='left',
                     bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.7))
        
        fig.tight_layout()
        if save: fig.savefig(filename, dpi=dpi, bbox_inches="tight")
        if disp: plt.show()
        plt.close(fig)
        
    def plot_D(self, 
               title: str | None = "Diffusivity (Arrhenius Plot)",
               disp: bool = True,
               save: bool = True, 
               filename: str = "D.png", 
               dpi: int = 300) -> None:
        """Generates an Arrhenius plot for the final diffusivity (D).

        Args:
            title (str | None, optional): A custom title for the plot.
            disp (bool, optional): If True, displays the plot. Defaults to True.
            save (bool, optional): If True, saves the plot to a file. Defaults to True.
            filename (str, optional): Filename for the saved plot.
            dpi (int, optional): Resolution for the saved figure.
        """
        if len(self.temperatures) < 2:
            raise ValueError("At least two temperature points are required for an Arrhenius plot.")
        if self.D is None or self.Ea_D is None:
             raise RuntimeError("Please run analysis and fitting for 'D' first.")
         
        slope = -self.Ea_D / (self.kb * 1000)
        intercept = np.log(self.D0)
        
        fig, ax = self._create_arrhenius_plot(
            data=self.D, temperatures=self.temperatures,
            slope=slope, intercept=intercept,
            ylabel=r'ln(D) (m$^2$/s)', title=title
        )
        
        text_str = (f'$E_a = {self.Ea_D:.2f}$ eV\n'
                    f'$D_0 = {self.D0:.2e}$ m$^2$/s\n'
                    f'$R^2 = {self.D_R2:.3f}$')
        
        ax.text(0.05, 0.05, text_str, transform=ax.transAxes, fontsize=10,
                verticalalignment='bottom', horizontalalignment='left',
                bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.7))
        
        if save: fig.savefig(filename, dpi=dpi, bbox_inches="tight")
        if disp: plt.show()
        plt.close(fig)
        
    def plot_tau(self, 
                 title: str | None = "Residence Time vs. Temperature",
                 disp: bool = True,
                 save: bool = True, 
                 filename: str = 'tau.png', 
                 dpi: int = 300) -> None:
        """Plots the average residence time (tau) vs. temperature.

        Args:
            title (str | None, optional): A custom title for the plot.
            disp (bool, optional): If True, displays the plot. Defaults to True.
            save (bool, optional): If True, saves the plot to a file. Defaults to True.
            filename (str, optional): Filename for the saved plot.
            dpi (int, optional): Resolution for the saved figure.
        """
        if self.tau is None:
            raise RuntimeError("Residence time has not been calculated. Please run .calculate() first.")
        
        valid_mask = ~np.isnan(self.tau)
        temps_valid = self.temperatures[valid_mask]
        tau_valid = self.tau[valid_mask]

        if len(tau_valid) == 0:
            print("Warning: No valid residence time data to plot.")
            return
        
        fig, ax = plt.subplots(figsize=(7, 6))
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.2)
        
        cmap = plt.get_cmap("viridis", len(self.temperatures))
        temp_color_map = {temp: cmap(i) for i, temp in enumerate(self.temperatures)}
        colors = [temp_color_map[t] for t in temps_valid]

        if len(temps_valid) > 1:
            bar_width = np.mean(np.diff(temps_valid)) * 0.7
        else:
            bar_width = temps_valid[0] * 0.1

        ax.bar(temps_valid, tau_valid, width=bar_width, color=colors,
               edgecolor='k', alpha=0.6, zorder=2)
               
        ax.scatter(temps_valid, tau_valid, c=colors,
                   marker='o', s=100, edgecolor='k', alpha=0.9, zorder=3)
        
        if self.tau0 is not None and not np.isnan(self.tau0) and len(temps_valid) > 1:
            x_fit = np.linspace(np.min(temps_valid) * 0.95, np.max(temps_valid) * 1.05, 200)     
            y_fit = self.tau0 * np.exp(self.Ea_tau / (self.kb * x_fit))
            ax.plot(x_fit, y_fit, 'k--', linewidth=1.5, label='Arrhenius Fit', zorder=4)
            
            text_str = (f'$E_a = {self.Ea_tau:.2f}$ eV\n'
                        f'$\\tau_0 = {self.tau0:.2e}$ ps\n'
                        f'$R^2 = {self.tau_R2:.3f}$')
            
            ax.text(0.05, 0.05, text_str, transform=ax.transAxes, fontsize=10,
                    verticalalignment='bottom', horizontalalignment='left',
                    bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.7))

        ax.set_xlabel('Temperature (K)', fontsize=12)
        ax.set_ylabel(r'Residence Time, $\tau$ (ps)', fontsize=12)
        if title: ax.set_title(title, fontsize=13, pad=10)
        ax.grid(True, which='both', linestyle='--', alpha=0.7)
        ax.legend()
        
        fig.tight_layout()
        if save: fig.savefig(filename, dpi=dpi, bbox_inches="tight")
        if disp: plt.show()
        plt.close(fig)

    def plot_a(self,
               title: str | None = "Hopping Distance vs. Temperature",
               disp: bool = True,
               save: bool = True,
               filename: str = 'a.png', 
               dpi: int = 300) -> None:
        """Plots the hopping distance (a) as a function of temperature.

        Args:
            title (str | None, optional): A custom title for the plot.
            disp (bool, optional): If True, displays the plot. Defaults to True.
            save (bool, optional): If True, saves the plot to a file. Defaults to True.
            filename (str, optional): Filename for the saved plot.
            dpi (int, optional): Resolution for the saved figure.
        """
        if self.a is None:
            raise RuntimeError("Hopping distance has not been calculated. Please run .calculate() first.")

        valid_mask = ~np.isnan(self.a)
        temps_valid = self.temperatures[valid_mask]
        a_valid = self.a[valid_mask]

        fig, ax = plt.subplots(figsize=(7, 6))
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.2)
            
        cmap = plt.get_cmap("viridis", len(self.temperatures))
        temp_color_map = {temp: cmap(i) for i, temp in enumerate(self.temperatures)}
        
        plt.plot(self.temperatures, self.a, linestyle='--', color='k')
        
        for i in range(len(temps_valid)):
            temp = temps_valid[i]
            ax.scatter(temp, a_valid[i], color=temp_color_map[temp],
                       marker='s', s=100, label=f"{temp:.0f} K", edgecolor='k', alpha=0.8)
        
        ax.set_xlabel('Temperature (K)', fontsize=12)
        ax.set_ylabel(r'Effective Hopping distance, $a_{eff}$ (Å)', fontsize=12)
        if title: ax.set_title(title, fontsize=14, pad=10)
        ax.legend(title='Temperature', fontsize=11, title_fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.7)
        
        bottom, top = ax.get_ylim()
        
        if (top - bottom) < 0.5:
            center = (top + bottom) / 2
            new_bottom = center - 0.25
            new_top = center + 0.25
            ax.set_ylim(new_bottom, new_top)
        
        fig.tight_layout()
        if save: fig.savefig(filename, dpi=dpi, bbox_inches="tight")
        if disp: plt.show()
        plt.close(fig)
        
    def plot_counts(self,
                    title: str | None = "Total Hopping Counts per Path",
                    disp: bool = True,
                    save: bool = True,
                    filename: str = 'counts.png', 
                    dpi: int = 300) -> None:
        """Generates a bar plot of total hop counts for each path type.

        This plot visualizes the frequency of all observed hopping events,
        summed across all trajectories in the ensemble.

        Args:
            title (str | None, optional): A custom title for the plot.
            disp (bool, optional): If True, displays the plot. Defaults to True.
            save (bool, optional): If True, saves the plot to a file. Defaults to True.
            filename (str, optional): Filename for the saved plot.
            dpi (int, optional): Resolution for the saved figure.
        """
        if self.counts is None or self.counts_unknown is None:
            raise RuntimeError("Path counts have not been calculated.")

        name_all = self.site.path_name + [p['name'] for p in self.path_unknown]
        counts_all = np.append(np.sum(self.counts, axis=0), np.array(self.counts_unknown, dtype=int))
       
        if len(name_all) == 0:
            print("Warning: No path count data available to plot.")
            return

        fig, ax = plt.subplots(figsize=(max(6, len(name_all) * 0.5), 6))
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.2)
            
        x_pos = np.arange(len(name_all))
        bars = ax.bar(x_pos, counts_all, color='steelblue', edgecolor='k', alpha=0.8)

        ax.set_ylabel('Total Counts', fontsize=12)
        ax.set_xlabel('Path Name', fontsize=12)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(name_all, rotation=45, ha="right")
        if title: ax.set_title(title, fontsize=13, pad=10)
            
        ax.bar_label(bars, fmt='%d', padding=3, fontsize=9)
        ax.yaxis.grid(True, linestyle='--', alpha=0.7)
        ax.set_axisbelow(True)
        ax.set_ylim(top=ax.get_ylim()[1] * 1.15)
        
        fig.tight_layout()
        if save: fig.savefig(filename, dpi=dpi, bbox_inches="tight")
        if disp: plt.show()
        plt.close(fig)

    def summary(self):
        """Prints a comprehensive summary of the ensemble analysis results.

        The summary includes a table of temperature-dependent data (D, f, tau)
        and a detailed breakdown of the final fitted Arrhenius parameters for
        each physical quantity.
        """
        if self.calculators is None:
            print("Calculator has not been run. Running .calculate() now...")
            self.calculate()
        
        print("=" * 68)
        print("Summary for Trajectory dataset")
        print(f"  - Path to TRAJ bundle : {self.path_traj} (depth={self.depth})")
        print(f"  - Lattice structure   : {self.site.path_structure}")
        print(f"  - t_interval          : {self.t_interval:.3f} ps ({self.frame_interval} frames)")
        print(f"  - Temperatures (K)    : {self.temperatures.tolist()}")
        print(f"  - Num. of TRAJ files  : {self.num_calc_temp.tolist()}")
        print("=" * 68)
        
        if self.calculators and len(self.temperatures) > 0:
            print("\n" + "="*20 + " Temperature-Dependent Data " + "="*20)
            
            headers = ["T (K)", "D (m2/s)", "D_rand (m2/s)", "f", "tau (ps)", "a (Ang)"]
            
            formats = [".1f", ".3e", ".3e", ".4f", ".4f", ".4f"]
            table_data = []
            for row in zip(self.temperatures, self.D, self.D_rand, self.f, self.tau, self.a):
                formatted_row = [f"{value:{fmt}}" for value, fmt in zip(row, formats)]
                table_data.append(formatted_row)
    
            table = tabulate(table_data, headers=headers, 
                             tablefmt="simple", stralign='left', numalign='left')
            print(table)
            print("=" * 68)
        
        if self.Ea_D is not None:
            print("\n" + "="*21 + " Final Fitted Parameters " + "="*22)
            print(f"Diffusivity (D):")
            print(f"  - Ea          : {self.Ea_D:.3f} eV")
            print(f"  - D0          : {self.D0:.3e} m^2/s")
            print(f"  - R-squared   : {self.D_R2:.4f}")
            print(f"Random Walk Diffusivity (D_rand):")
            print(f"  - Ea          : {self.Ea_D_rand:.3f} eV")
            print(f"  - D0          : {self.D_rand0:.3e} m^2/s")
            print(f"  - R-squared   : {self.D_rand_R2:.4f}")
            print(f"Correlation Factor (f):")
            print(f"  - Ea          : {self.Ea_f:.3f} eV")
            print(f"  - f0          : {self.f0:.3f}")
            print(f"  - R-squared   : {self.f_R2:.4f}")
            print(f"Residence Time (tau):")
            print(f"  - Ea (fixed)  : {self.Ea_D_rand:.3f} eV")
            print(f"  - tau0        : {self.tau0:.3e} ps")
            print(f"  - R-squared   : {self.tau_R2:.4f}")
        else:
            print("\n" + "="*21 + " Final Fitted Parameters " + "="*22)
            print(f"Diffusivity (D):")
            print(f"  - Ea          : - ")
            print(f"  - D0          : - ")
            print(f"  - R-squared   : -")
            print(f"Random Walk Diffusivity (D_rand):")
            print(f"  - Ea          : -")
            print(f"  - D0          : -")
            print(f"  - R-squared   : -")
            print(f"Correlation Factor (f):")
            print(f"  - Ea          : -")
            print(f"  - f0          : -")
            print(f"  - R-squared   : -")
            print(f"Residence Time (tau):")
            print(f"  - Ea (fixed)  : -")
            print(f"  - tau0        : -")
            print(f"  - R-squared   : -")
        print("=" * 68)
            
    def show_hopping_paths(self) -> None:
        """ Prints a summary table of all observed hopping path types and their counts."""
        
        path_info = []
        for i, path in enumerate(self.site.path):
            p = {
                'name': path['name'],
                'a': path['distance'],
                'z': path['z'],
                'count': np.sum(self.counts, axis=0)[i],
                'site_init': f"{path['site_init']} [{', '.join(f'{x:.5f}' for x in path['coord_init'])}]",
                'site_final': f"{path['site_final']} [{', '.join(f'{x:.5f}' for x in path['coord_final'])}]"
            }
            path_info.append(p)
        
        for i, path in enumerate(self.path_unknown):
            p = {
                'name': path['name'],
                'a': path['distance'],
                'z': '-',
                'count': int(np.sum(self.counts_unknown, axis=0)[i]),
                'site_init': f"{path['site_init']} [{', '.join(f'{x:.5f}' for x in path['coord_init'])}]",
                'site_final': f"{path['site_final']} [{', '.join(f'{x:.5f}' for x in path['coord_final'])}]"
            }
            path_info.append(p)
            
        print("=" * 110)
        print(" " * 43 + "Hopping Path Information")
        print("=" * 110)
        header = ['Name',
                  'a (Ang)',
                  'z',
                  'Count',
                  'Initial Site (Fractional Coordinate)',
                  'Final Site (Fractional Coordinate)']
        data = [
            [
                p['name'],
                p['a'],
                p['z'],
                p['count'],
                p['site_init'],
                p['site_final']
            ] for p in path_info
        ]
        print(tabulate(data, headers=header, tablefmt="simple", stralign='left', numalign='left'))
        print("=" * 110 + "\n")
    
    def show_hopping_history(self) -> None:
        """Prints instructions on how to view the history for a single trajectory.

        This method does not print a history table directly, as this is an
        ensemble analysis. Instead, it provides a user guide explaining how to
        access the `show_hopping_history` method of an individual
        `CalculatorSingle` object within the `.calculators` list.
        """
        
        print("\n" + "="*87)
        print(" "*5 + "NOTE: `show_hopping_history()` is a method for a single trajectory analysis!")
        print("="*87)
        print("This method displays the detailed hop-by-hop history for a single simulation.\n"
              "The current object is a 'CalculatorEnsemble' that manages a collection of simulations.\n\n"
              "To view the history for a specific trajectory, you need to select one of the individual\n"
              "calculator objects stored in the `.calculators` list using its index.\n")
        print("You can find the index for a specific file by checking the `.all_traj_paths` attribute.")
        print("The order of files in this list matches the order of the `.calculators` list.\n")
        print("Example command to list all files and their indices:")
        print("  >>> for i, path in enumerate(your_bundle_object.all_traj_paths):")
        print("  ...     print(f'Index {i}: {path_traj}')")
        print("\nOnce you know the index, you can run:")
        print("  # Show history for the first trajectory in the bundle (index 0)")
        print("  >>> your_bundle_object.calculators[0].show_hopping_history()\n")
        print("  # Show history for the third trajectory in the bundle (index 2)")
        print("  >>> your_bundle_object.calculators[2].show_hopping_history()\n")
        print("="*87)
    
    def save_trajectory(self,
                        path_dir: str = 'trajectories',
                        prefix: str = 'trajectory') -> None:
        """Saves the unwrapped vacancy trajectories from all simulations.

        This method creates a directory structure organized by temperature
        (e.g., 'trajectories/300K/'), and saves the unwrapped vacancy trajectory
        from each simulation as a separate JSON file within the appropriate
        subdirectory.

        Args:
            path_dir (str, optional): The root directory for saving trajectories.
            prefix (str, optional): The prefix for the output JSON filenames.
        """
        if not hasattr(self, 'calculators') or not self.calculators:
            raise RuntimeError("Analysis has not been run. Please call the .calculate() method first.")

        for i, (start, end) in enumerate(zip(self.index_calc_temp[:-1], self.index_calc_temp[1:])):
            temp_str = f"{self.temperatures[i]:.0f}K"
            path_dir_temp = os.path.join(path_dir, temp_str)
            os.makedirs(path_dir_temp, exist_ok=True)

            label = 0
            for index in range(start, end):
                label += 1
                filename = f"{prefix}_{label:02d}.json"
                output_path = os.path.join(path_dir_temp, filename)
                self.calculators[index].save_trajectory(filename=output_path)

        if self.verbose:
            print(f"All vacancy trajectories have been saved in the '{path_dir}' directory.")
            
    def save_parameters(self,
                        filename: str = "parameters.json") -> None:
        """Saves all calculated and fitted physical parameters to a JSON file.

        The output JSON includes a 'description' key that explains each parameter,
        making the file self-documenting.

        Args:
            filename (str, optional): The name of the output JSON file.
        """
            
        description = {
            'D'         : 'Diffusivity (m2/s): (n_temperatures,)',
            'D0'        : 'Pre-exponential factor for diffusivity (m2/s)',
            'Ea_D'      : 'Activation barrier for diffusivity (eV)',
            'D_rand'    : 'Random walk diffusivity (m2/s): (n_temperatures,)',
            'D_rand0'   : 'Pre-exponential factor for random walk diffusivity (m2/s)',
            'Ea_D_rand' : 'Activation barrier for random walk diffusivity (eV)',
            'f'         : 'Correlation factor: (n_temperatures,)',
            'f0'        : 'Pre-exponential factor for correlation factor',
            'Ea_f'      : 'Activation barrier for correlation factor (eV)',
            'tau'       : 'Residence time (ps): (n_temperatures,)',
            'tau0'      : 'Pre-exponential factor for residence time (ps)',
            'Ea_tau'    : 'Activation barrier for residence time (eV)',
            'a'         : 'Effective hopping distance (Ang) (n_temperatures,)',
            'a_path'    : 'Path-wise hopping distance (Ang): (n_paths,)',
            'nu'        : 'Effective attempt frequency (THz): (n_temperatures,)',
            'nu_path'   : 'Path-wise attempt frequency (THz): (n_temperatures, n_paths)',
            'Ea_path'   : 'Path-wise hopping barrier (eV): (n_paths)',
            'z'         : 'Effective number of the equivalent paths: (n_temperatures,)',
            'z_path'    : 'Number of the equivalent paths of each path: (n_paths,)',
            'z_mean'    : 'Mean number of the equivalent paths per path type: (n_temperatures,)',
            'm_mean'    : 'Mean number of path types: (n_temperatures,)',
            'P_site'    : 'Site occupation probability: (n_temperatures, n_sites)',
            'P_esc'     : 'Escape probability: (n_temperature, n_paths)',
            'times_site': 'Total residence times at each site (ps): (n_temperature, n_sites)',
            'counts_hop': 'Counts of hops at each temperature: (n_temperature, n_paths)'
        }
        is_list = {
            'D'         : True,
            'D0'        : False,
            'Ea_D'      : False,
            'D_rand'    : True,
            'D_rand0'   : False,
            'Ea_D_rand' : False,
            'f'         : True,
            'f0'        : False,
            'Ea_f'      : False,
            'tau'       : True,
            'tau0'      : False,
            'Ea_tau'    : False,
            'a'         : True,
            'a_path'    : True,
            'nu'        : True,
            'nu_path'   : True,
            'Ea_path'   : True,
            'z'         : False,
            'z_path'    : True,
            'z_mean'    : True,
            'm_mean'    : True,
            'P_site'    : True,
            'P_esc'     : True,
            'times_site': True,
            'counts_hop': True
        }
        contents = {
            'symbol': self.symbol,
            'path': self.site.path_name,
            'site': self.site.site_name,
            'temperatures': self.temperatures.tolist(),
            'num_vacancies': self.num_vacancies,
        }
        
        for param in description.keys():
            if hasattr(self, param):
                value = getattr(self, param)
                if isinstance(value, (np.ndarray, np.generic)):
                    value = value.tolist()
                if is_list[param] and not isinstance(value, list):
                    contents[param] = [value]
                else:
                    contents[param] = value
            else:
                contents[param] = None
        contents['description'] = description
        
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(contents, f, indent=2)
        
        if self.verbose:
            print(f"Parameters saved to '{filename}'")

    # ===================================================================
    # Methods for attempt frequency
    # ===================================================================
            
    def calculate_attempt_frequency(self,
                                    neb_csv: str,
                                    filename: str = "parameters.json") -> None:
        """Calculates attempt frequencies and updates the parameter file.

        This method first saves the current analysis results, then runs the
        `AttemptFrequency` analysis using the provided NEB data, and finally
        updates the parameter file with the newly calculated frequency results.

        Args:
            neb_csv (str):
                Path to the CSV file containing NEB-calculated activation barriers.
            filename (str, optional):
                The name of the parameter JSON file to create and update.
                Defaults to "parameters.json".
        """
        
        self.save_parameters(filename=filename)
        self.attempt_frequency = AttemptFrequency(filename, neb_csv, verbose=False)
        self.attempt_frequency.update_json()
        self.z = self.attempt_frequency.z
        self.nu = self.attempt_frequency.nu
        self.nu_path = self.attempt_frequency.nu_path
        
        if self.verbose:
            self.attempt_frequency.summary()
    
    def plot_nu(self,
                title: str | None = "Attempt Frequency vs. Temperature",
                disp: bool = True,
                save: bool = True,
                filename: str = "attempt_frequency.png",
                dpi: int = 300) -> None:
        """Plots the effective attempt frequency (nu) vs. temperature.

        This method delegates the plotting task to the internal `AttemptFrequency` object.
        `.calculate_attempt_frequency()` must be called first.

        Args:
            title (str | None, optional): A custom title for the plot.
            disp (bool, optional): If True, displays the plot. Defaults to True.
            save (bool, optional): If True, saves the plot. Defaults to True.
            filename (str, optional): Filename for the saved plot.
            dpi (int, optional): Resolution for the saved figure.
        """
        if self.attempt_frequency is None:
            raise RuntimeError("Attempt frequencies have not been calculated. "
                               "Please run the .calculate_attempt_frequency() method first.")
        
        self.attempt_frequency.plot_nu(title=title,
                                       disp=disp,
                                       save=save,
                                       filename=filename,
                                       dpi=dpi)
        
    def plot_z(self,
               title: str | None = "Coordination Number vs. Temperature",
               disp: bool = True,
               save: bool = True,
               filename: str = "coordination_number.png",
               dpi: int = 300) -> None:
        """Plots the effective coordination number (z) vs. temperature.

        This method delegates the plotting task to the internal `AttemptFrequency` object.
        `.calculate_attempt_frequency()` must be called first.
        
        Args:
            title (str | None, optional): A custom title for the plot.
            disp (bool, optional): If True, displays the plot. Defaults to True.
            save (bool, optional): If True, saves the plot. Defaults to True.
            filename (str, optional): Filename for the saved plot.
            dpi (int, optional): Resolution for the saved figure.
        """
        if self.attempt_frequency is None:
            raise RuntimeError("Attempt frequencies have not been calculated. "
                               "Please run the .calculate_attempt_frequency() method first.")
            
        self.attempt_frequency.plot_z(title=title,
                                      disp=disp,
                                      save=save,
                                      filename=filename,
                                      dpi=dpi)
    
    # ===================================================================
    # Methods for decomposing diffusivity
    # ===================================================================    
        
    def decompose_diffusivity(self, verbose: bool = True) -> None:
        """Decomposes the final diffusivity (self.D) into directional components
        (Dx, Dy, Dz) and performs Arrhenius fits on them.

        This method uses a statistically robust approach:

        1. The already-calculated total diffusivity (self.D) is used as the basis.
        2. For each temperature, the MSD component ratios (x, y, z) are calculated
           at each available time step and then time-averaged to get a stable
           mean contribution for each direction.
        3. The total D is apportioned into Dx, Dy, and Dz using these mean ratios,
           ensuring physical consistency.
        4. Finally, an Arrhenius fit is performed on the decomposed diffusivities.

        The results are stored in instance attributes like `self.Dx`, `self.Ea_x`,
        `self.D0_x`, etc.

        Args:
            verbose (bool, optional):
                If True, prints a summary of the fitting results. Defaults to True.
        """
        
        if self.D is None:
            raise RuntimeError("Total diffusivity (self.D) has not been calculated. "
                               "Please run the main .calculate() method first.")
        
        # Temperature-wise MSD calculation
        if not hasattr(self, 'msd') or not self.msd:
            self.msd = []
            for start, end in zip(self.index_calc_temp[:-1], self.index_calc_temp[1:]):
                num_steps = [c.num_steps for c in self.calculators[start:end]]
                if not num_steps:
                    self.msd.append(np.empty((0, 3)))
                    continue

                max_steps_in_segment = max(num_steps)
                squared_disp_sum = np.zeros((max_steps_in_segment, 3))
                step_counts = np.zeros(max_steps_in_segment)

                for c in self.calculators[start:end]:
                    n_vac, n_step = c.num_vacancies, c.num_steps
                    traj_dict = c.unwrapped_vacancy_trajectory_coord_cart
                    if not traj_dict: continue
                    
                    sorted_keys = sorted(traj_dict.keys())
                    positions = np.array([traj_dict[key] for key in sorted_keys])
                    disp = positions - positions[0, :, :]
                    squared_disp = disp ** 2
                    squared_disp_sum[:n_step] += np.sum(squared_disp, axis=1)
                    step_counts[:n_step] += n_vac
                
                counts_reshaped = step_counts.reshape(-1, 1)
                msd_temp = np.divide(squared_disp_sum, counts_reshaped,
                                     out=np.zeros_like(squared_disp_sum),
                                     where=(counts_reshaped != 0))
                self.msd.append(msd_temp)

        # Decompose self.D into xyz
        self.Dx, self.Dy, self.Dz = [], [], []

        for i, msd_T in enumerate(self.msd):
            if msd_T.shape[0] < 2:
                self.Dx.append(np.nan); 
                self.Dy.append(np.nan); 
                self.Dz.append(np.nan)
                continue
            
            msd_total_per_step = np.sum(msd_T, axis=1, keepdims=True)
            ratios_over_time = np.divide(msd_T, msd_total_per_step,
                                         out=np.full_like(msd_T, 1/3),
                                         where=(msd_total_per_step > self.eps))
            
            # Time-averaged ration_xyz
            mean_ratios = np.mean(ratios_over_time, axis=0)

            # D_total = (Dx + Dy + Dz) / 3
            self.Dx.append(3 * self.D[i] * mean_ratios[0])
            self.Dy.append(3 * self.D[i] * mean_ratios[1])
            self.Dz.append(3 * self.D[i] * mean_ratios[2])

        self.Dx = np.array(self.Dx)
        self.Dy = np.array(self.Dy)
        self.Dz = np.array(self.Dz)
            
        # Arrhenius fitting for D_xyz
        if len(self.temperatures) >= 2:
            try:
                self.Ea_x, self.D0_x, self.Dx_R2 = self._Arrhenius_fit(self.temperatures, self.Dx, self.kb)
            except (ValueError, np.linalg.LinAlgError):
                self.Ea_x, self.D0_x, self.Dx_R2 = np.nan, np.nan, np.nan

            try:
                self.Ea_y, self.D0_y, self.Dy_R2 = self._Arrhenius_fit(self.temperatures, self.Dy, self.kb)
            except (ValueError, np.linalg.LinAlgError):
                self.Ea_y, self.D0_y, self.Dy_R2 = np.nan, np.nan, np.nan

            try:
                self.Ea_z, self.D0_z, self.Dz_R2 = self._Arrhenius_fit(self.temperatures, self.Dz, self.kb)
            except (ValueError, np.linalg.LinAlgError):
                self.Ea_z, self.D0_z, self.Dz_R2 = np.nan, np.nan, np.nan

        if verbose:
            print("="*11 + " Temperature-Dependent Decomposed Diffusivity " + "="*11)
            headers = ["T (K)", "Dx (m2/s)", "Dy (m2/s)", "Dz (m2/s)"]
            formats = [".1f", ".3e", ".3e", ".3e"]
            
            table_data = []
            for row in zip(self.temperatures, self.Dx, self.Dy, self.Dz):
                formatted_row = [f"{value:{fmt}}" for value, fmt in zip(row, formats)]
                table_data.append(formatted_row)

            table = tabulate(table_data, headers=headers, 
                            tablefmt="simple", stralign='left', numalign='left')
            print(table)
            print("=" * 68)

            if len(self.temperatures) >= 2:
                print("\n" + "="*19 + " Decomposed Diffusivity Fits " + "="*20)
                print(f"Directional Diffusivity (Dx):")
                print(f"  - Ea        : {self.Ea_x:.3f} eV")
                print(f"  - D0        : {self.D0_x:.3e} m^2/s")
                print(f"  - R-squared : {self.Dx_R2:.4f}")
                print(f"Directional Diffusivity (Dy):")
                print(f"  - Ea        : {self.Ea_y:.3f} eV")
                print(f"  - D0        : {self.D0_y:.3e} m^2/s")
                print(f"  - R-squared : {self.Dy_R2:.4f}")
                print(f"Directional Diffusivity (Dz):")
                print(f"  - Ea        : {self.Ea_z:.3f} eV")
                print(f"  - D0        : {self.D0_z:.3e} m^2/s")
                print(f"  - R-squared : {self.Dz_R2:.4f}")
                print("=" * 68)

    def plot_msd_xyz(self,
                     title: str | None = "MSD Components vs. Time",
                     disp: bool = True,
                     save: bool = True,
                     filename: str = "msd_xyz.png",
                     dpi: int = 300) -> None:
        """
        Plots the x, y, and z components of the MSD vs. time for all temperatures.

        Each component is plotted on a separate subplot in a horizontal layout.
        This plot shows the raw MSD data without any fitting lines.

        Args:
            title (str | None, optional): A custom title for the plot figure.
            disp (bool, optional): If True, displays the plot. Defaults to True.
            save (bool, optional): If True, saves the plot to a file. Defaults to True.
            filename (str, optional): Filename for the saved plot.
            dpi (int, optional): Resolution for the saved figure.
        """
        if not hasattr(self, 'msd') or not self.msd:
            raise RuntimeError("MSD has not been calculated. Please run .decompose_diffusivity() first.")

        fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharey=True)
        components = ['x', 'y', 'z']
        
        cmap = plt.get_cmap("viridis", len(self.temperatures))
        temp_color_map = {temp: cmap(i) for i, temp in enumerate(self.temperatures)}

        for i, ax in enumerate(axes):
            for temp_idx, msd_T in enumerate(self.msd):
                if msd_T.shape[0] == 0: continue
                
                temperature = self.temperatures[temp_idx]
                color = temp_color_map[temperature]
                
                n_steps = msd_T.shape[0]
                time_ps = np.arange(n_steps) * self.t_interval

                ax.plot(time_ps, msd_T[:, i], color=color, alpha=0.7,
                        label=f"{temperature:.0f} K")
            
            ax.set_xlabel('Time (ps)', fontsize=12)
            ax.set_ylabel(f'MSD$_{{{components[i]}}}$ (Å²)', fontsize=12)
            ax.grid(True, linestyle='--', alpha=0.7)
            ax.legend(title='Temperature', fontsize=10)

        if title: fig.suptitle(title, fontsize=16)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        
        if save: fig.savefig(filename, dpi=dpi, bbox_inches="tight")
        if disp: plt.show()
        plt.close(fig)
        
    def plot_D_xyz(self,
                   title: str | None = "Decomposed Diffusivity (Arrhenius Plot)",
                   disp: bool = True,
                   save: bool = True,
                   filename: str = "D_xyz.png",
                   dpi: int = 300) -> None:
        """
        Generates Arrhenius plots for the decomposed directional diffusivities (Dx, Dy, Dz).

        This method creates a figure with three horizontal subplots, one for each
        spatial component, to visualize their temperature dependence and Arrhenius fits.

        Args:
            title (str | None, optional): A custom title for the overall figure.
            disp (bool, optional): If True, displays the plot. Defaults to True.
            save (bool, optional): If True, saves the plot to a file. Defaults to True.
            filename (str, optional): Filename for the saved plot.
            dpi (int, optional): Resolution for the saved figure.
        """
        if len(self.temperatures) < 2:
            raise ValueError("At least two temperature points are required for an Arrhenius plot.")
        if not hasattr(self, 'Ea_x') or self.Ea_x is None:
            raise RuntimeError("Decomposed diffusivities have not been calculated and fitted. "
                               "Please run .decompose_diffusivity() first.")

        components = [
            {'label': 'x', 'data': self.Dx, 'Ea': self.Ea_x, 'D0': self.D0_x, 'R2': self.Dx_R2},
            {'label': 'y', 'data': self.Dy, 'Ea': self.Ea_y, 'D0': self.D0_y, 'R2': self.Dy_R2},
            {'label': 'z', 'data': self.Dz, 'Ea': self.Ea_z, 'D0': self.D0_z, 'R2': self.Dz_R2}
        ]

        fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharey=True)

        for i, comp in enumerate(components):
            ax = axes[i]
            is_fit_valid = not (np.isnan(comp['Ea']) or np.isnan(comp['D0']))
            
            if is_fit_valid:
                slope = -comp['Ea'] / (self.kb * 1000)
                intercept = np.log(comp['D0'])
            else:
                slope, intercept = np.nan, np.nan
            
            self._create_arrhenius_plot(
                data=comp['data'],
                temperatures=self.temperatures,
                slope=slope,
                intercept=intercept,
                ylabel=fr'ln($D_{{{comp["label"]}}}$) (m$^2$/s)',
                title=f'Directional Diffusivity ($D_{{{comp["label"]}}}$)',
                ax=ax
            )

            if is_fit_valid:
                text_str = (f'$E_a = {comp["Ea"]:.2f}$ eV\n'
                            f'$D_0 = {comp["D0"]:.2e}$ m$^2$/s\n'
                            f'$R^2 = {comp["R2"]:.3f}$')
            else:
                text_str = "Fit not available"
                
            ax.text(0.05, 0.05, text_str, transform=ax.transAxes, fontsize=10,
                    verticalalignment='bottom', horizontalalignment='left',
                    bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.7))

        if title: fig.suptitle(title, fontsize=16)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        if save: fig.savefig(filename, dpi=dpi, bbox_inches="tight")
        if disp: plt.show()
        plt.close(fig)
