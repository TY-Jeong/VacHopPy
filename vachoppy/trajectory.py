import os
import sys
import h5py
import json
import math
import copy   
import numpy as np
import matplotlib.pyplot as plt

from PIL import Image
from tqdm import tqdm
from colorama import Fore, Style
from tabulate import tabulate

from collections import defaultdict
from itertools import permutations

# For Arrow3D
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

import plotly.graph_objects as go
from vachoppy.io import Site, monitor_performance, parse_lammps

# color map for tqdm
BOLD = '\033[1m'
CYAN = '\033[36m'
MAGENTA = '\033[35m'
GREEN = '\033[92m' # Green color
RED = '\033[91m'   # Red color
RESET = '\033[0m'  # Reset to default color


class Arrow3D(FancyArrowPatch):
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
    """
    Analyzes a single molecular dynamics trajectory to trace atomic and vacancy hops.

    This class reads a trajectory from an HDF5 file, identifies which lattice
    site each atom occupies at each time step, validates hops using a
    transition state criterion, and reconstructs the movement paths of vacancies.
    It provides methods for interactive and static visualization of the results.

    Args:
        traj (str): 
            Path to the HDF5 trajectory file.
        site (Site): 
            An instance of the `Site` class containing pre-analyzed
            information about the crystal lattice sites.
        t_interval (float): 
            The time interval in picoseconds (ps) between
            each analysis step. This determines the coarse-graining level.
        verbose (bool, optional): 
            Verbosity tag. Defaults to True.

    Attributes:
        occupation (np.ndarray): 
            A 2D array of shape (num_atoms, num_steps)
            storing the index of the lattice site occupied by each atom at each step.
        trace_arrows (dict): 
            A dictionary mapping each time step to a list of
            hop events (arrows) that occurred at that step.
        vacancy_trajectory_index (dict): 
            A dictionary mapping each time step to an
            array of occupied vacancy site indices.
        vacancy_trajectory_coord_cart (dict): 
            A dictionary mapping each time step
            to the Cartesian coordinates of the occupied vacancy sites.
        unwrapped_vacancy_trajectory_coord_cart (dict):
            A dictionary mapping each time step to the Cartesian coordinates of vacancies. 
            These coordinates are unwrapped across periodic boundaries to represent the 
            true, continuous diffusion path.
        hopping_sequence (dict): 
            A dictionary mapping each time step to the
            reconstructed vacancy hop paths.
    """

    def __init__(self,
                 traj: str,
                 site,
                 t_interval: float,
                 verbose: bool = True):
        
        self._validate_traj(traj)
        self.cmap = self._custrom_cmap()
        
        self.traj = traj
        self.site = site
        self.t_interval = t_interval    # interval
        self.verbose = verbose
        
        # Read cond
        self.total_frames = None    # nsw
        self.temperature = None     # temp
        self.dt = None              # potim
        self.symbol = None
        self.num_frames = None      # nsw_cut
        self.num_steps = None       # num_step
        self.num_atoms = None       # num_atom
        self.frame_interval = None  # interval_nsw
        self.lattice_parameter = None
        self._read_cond()
        
        # Read site
        self.lattice_sites = None       # lattice_point
        self.lattice_sites_cart = None  # lattice_point_C
        self.num_lattice_sites = None   # num_lattice_point
        self.num_vacancies = None       # num_vac
        self._read_site()
        
        self.occupation = None
        self._get_occupation()          # atomic_trajectory
        
        self.trace_arrows = None
        self._get_trace_arrows()        # get_trace_arrows
        
        # Vacancy trajectory
        self.hopping_sequence = {}
        self.vacancy_trajectory_index = {}
        self.vacancy_trajectory_coord_cart = {} # vacancy_trajectory_coord_C
        self.transient_vacancy = {}
        self._get_vacancy_trajectory()  # get_vacancy_trajectory
        
        # Unwrapped trajectory
        self.unwrapped_vacancy_trajectory_coord_cart = None
        self._get_unwrapped_vacancy_trajectory()
    
    @monitor_performance
    def animate_vacancy_trajectory(self,
                                   vacancy_indices: int | list,
                                   step_init: int = 0,
                                   step_final: int = None,
                                   unwrap: bool = False,
                                   max_trace_length: int = 10,
                                   update_alpha: float = 0.8,
                                   fps: int = 5,
                                   save: bool = True,
                                   filename: str = "video.html",
                                   verbose: bool = True) -> None:
        """
        Generates an interactive 3D animation of vacancy trajectories using Plotly.

        The animation includes a play/pause button and a slider to navigate through
        time steps. The real-time position of each vacancy is shown as a circle,
        and its recent path fades out over time to indicate the direction of movement.

        Args:
            vacancy_indices (int | list): 
                Index or list of indices for the vacancies to be plotted.
            step_init (int, optional): 
                The starting step for the animation. Defaults to 0.
            step_final (int, optional): 
                The ending step for the animation. If None, animates up to the 
                last available step. Defaults to None.
            unwrap (bool, optional): 
                If True, plots the continuous, unwrapped path in a supercell view. 
                Defaults to False.
            max_trace_length (int, optional): 
                The maximum number of past path segments to display as a fading
                tail. Defaults to 15.
                **Warning**: Increasing this value significantly can lead to
                longer computation times and higher memory usage, potentially
                slowing down the animation generation and increasing the
                output file size.
            update_alpha (float, optional): 
                Fading rate for past trajectory lines. A value closer to 0 
                makes paths fade faster. Defaults to 0.75.
            fps (int, optional): 
                Frames per second for the animation playback when the 'Play'
                button is pressed. Defaults to 20.
            save (bool, optional): 
                If True, saves the animation as a standalone HTML file. If False,
                displays it directly (e.g., in a Jupyter notebook). Defaults to True.
            filename (str, optional): 
                The name of the output HTML file if save is True.
                Defaults to "video.html".
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

        for s_idx, step in enumerate(available_steps):
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
            if self.verbose: print(f"'{filename}' created.")
        else:
            fig.show()

    def plot_vacancy_trajectory(self, 
                                vacancy_indices: int | list, 
                                step_init: int = 0, 
                                step_final: int = None,
                                unwrap: bool = False,
                                alpha: float = 0.6,
                                save: bool = False,
                                filename: str = "traj.html") -> None:
        """
        Generates an interactive 3D plot of specified vacancy trajectories using Plotly.

        Args:
            vacancy_indices (int | list): The index or list of indices for the
                vacancies to be plotted. E.g., 0 for the first vacancy, [0, 1] for both.
            step_init (int, optional): The starting step for the trajectory plot. 
                Defaults to 0.
            step_final (int, optional): The ending step for the trajectory plot.
                If None, plots up to the last available step. Defaults to None.
            unwrap (bool, optional): If True, plots the continuous, unwrapped path
                in a supercell view. Defaults to False.
            alpha (float, optional): Opacity for the trajectory lines. Lower values
                make lines more transparent, causing overlapping paths to appear
                darker. Defaults to 0.6.
            save (bool, optional): If True, saves the plot as a standalone HTML file.
                If False, displays the plot directly. Defaults to False.
            filename (str, optional): The name of the output HTML file if save is True.
                Defaults to "traj.html".
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
            if self.verbose:
                print(f"'{filename}' created.")
        else:
            fig.show()
        
    def animation(self,
                  index: list | str = 'all',
                  steps: list | str = 'all',
                  vac: bool = True,
                  gif: bool = True,
                  filename: str = 'traj.gif',
                  foldername: str = 'snapshots',
                  update_alpha: float = 0.75,
                  fps: int = 5,
                  loop: int = 0,
                  dpi: int = 150,
                  legend: bool = False,
                  label: bool = False) -> None:
        """
        Generates an animation of the atomic trajectory as a GIF file.

        Args:
            index (list | str, optional): 
                Indices of atoms to display. Defaults to 'all'.
            steps (list | str, optional): 
                Time steps to include in the animation. Defaults to 'all'.
            vac (bool, optional): 
                If True, displays vacancies. Defaults to True.
            gif (bool, optional): 
                If True, creates a GIF file from snapshots. Defaults to True.
            filename (str, optional): 
                Output GIF file name. Defaults to 'traj.gif'.
            foldername (str, optional): 
                Directory to save snapshots. Defaults to 'snapshots'.
            update_alpha (float, optional): 
                Fading rate for trace arrows. Defaults to 0.75.
            fps (int, optional): 
                Frames per second for the GIF. Defaults to 5.
            loop (int, optional): 
                Number of loops for the GIF (0 for infinite). Defaults to 0.
            dpi (int, optional): 
                DPI for saved snapshots. Defaults to 150.
            legend (bool, optional): 
                If True, displays a legend for atoms. Defaults to False.
            label (bool, optional): 
                If True, displays labels for lattice sites. Defaults to False.
        """
        if not os.path.isdir(foldername):
            os.mkdir(foldername)

        atom_indices = np.arange(self.num_atoms) if index == 'all' else np.array(index)
        step_indices = np.arange(self.num_steps) if steps == 'all' else np.array(steps)
        
        # tqdm progress bar color setup
        RED = Fore.RED; BOLD = Style.BRIGHT; RESET = Style.RESET_ALL
        
        files = []
        for step in tqdm(step_indices,
                         bar_format='{l_bar}%s{bar:35}%s{r_bar}{bar:-10b}'%(Fore.GREEN, RESET),
                         ascii=False,
                         desc=f'{RED}{BOLD}progress{RESET}'):
            
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(111, projection='3d')

            self._plot_lattice(ax, label=label)

            # Plot atoms
            for i, atom_idx in enumerate(atom_indices):
                site_idx = self.occupation[atom_idx, step]
                coords = self.lattice_sites_cart[site_idx]
                ax.scatter(*coords,
                           s=100,
                           facecolor=self.cmap[atom_idx % len(self.cmap)],
                           edgecolor='k',
                           alpha=0.8,
                           label=f"Atom {atom_idx}")
            
            # Plot trace arrows with fading alpha
            alpha = 1.0
            # Iterate backwards from the current step to create a fading effect
            for i in reversed(range(step + 1)):
                if alpha < 0.1: break # Stop drawing very faint arrows
                for arrow in self.trace_arrows.get(i, []):
                    arrow_prop = dict(mutation_scale=15, arrowstyle='->',
                                      color=arrow['c'], alpha=alpha,
                                      shrinkA=0, shrinkB=0, lw=1.5)
                    disp_arrow = Arrow3D(*arrow['p'].T, **arrow_prop)
                    ax.add_artist(disp_arrow)
                alpha *= update_alpha

            # Plot vacancies
            if vac:
                # Plot true vacancies
                vac_coords = self.vacancy_trajectory_coord_cart.get(step, [])
                if len(vac_coords) > 0:
                    ax.scatter(*vac_coords.T,
                              s=120, facecolor='yellow', edgecolor='k',
                              marker='o', alpha=0.8, zorder=10, label='Vacancy')
                
                # Plot transient vacancies
                trans_vac_indices = self.transient_vacancy.get(step, [])
                if len(trans_vac_indices) > 0:
                    trans_vac_coords = self.lattice_sites_cart[trans_vac_indices]
                    ax.scatter(*trans_vac_coords.T,
                              s=120, facecolor='orange', edgecolor='k',
                              marker='s', alpha=0.8, zorder=10, label='Transient Vacancy')

            # Set title with time and step information
            time = step * self.t_interval
            time_tot = (self.num_steps - 1) * self.t_interval
            ax.set_title(f"Time: {time:.2f} ps / {time_tot:.2f} ps (Step: {step}/{self.num_steps-1})")

            if legend:
                ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1.0))

            snapshot = os.path.join(foldername, f"snapshot_{step:04d}.png")
            files.append(snapshot)
            plt.savefig(snapshot, dpi=dpi, bbox_inches='tight')
            plt.close(fig)
        
        # Create GIF from snapshots
        if gif:
            print(f"Merging {len(files)} snapshots into a GIF...")
            imgs = [Image.open(file) for file in files]
            if imgs:
                imgs[0].save(fp=filename, format='GIF', append_images=imgs[1:],
                             save_all=True, duration=int(1000 / fps), loop=loop)
                print(f"Successfully created '{filename}'.")
    
    def distance_PBC(self, 
                     coord1: np.ndarray, 
                     coord2: np.ndarray) -> float | np.ndarray:
        """Calculates the PBC-aware distance between fractional coordinates.

        This method computes the shortest distance between one or more initial
        points (coord1) and a final point (coord2), respecting the periodic
        boundary conditions of the lattice.

        Args:
            coord1 (np.ndarray): 
                Initial coordinate(s) in fractional form.
                Can be a single point (1D array) or multiple points (2D array).
            coord2 (np.ndarray): 
                Final coordinate in fractional form (1D array).

        Returns:
            float | np.ndarray: The calculated distance(s) in Cartesian units.
                Returns a single float if coord1 is 1D, or a 1D array of
                distances if coord1 is 2D.
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
    
    def _validate_traj(self, traj: str) -> None:
        """
        Validates the structure and content of the HDF5 trajectory file.

        This method checks for the correct file extension, file existence, and
        the presence of required datasets ('positions', 'forces') and metadata
        attributes ('symbol', 'nsw', 'dt', 'temperature', 'atom_counts', 'lattice').

        Args:
            traj (str): The file path to validate.

        Raises:
            ValueError: If the file extension is not '.h5' or if required
                metadata or datasets are missing.
            FileNotFoundError: If the trajectory file does not exist.
            IOError: If the file cannot be read as an HDF5 file.
        """
        if not traj.endswith('.h5'):
            raise ValueError(f"Error: Trajectory file must have a .h5 extension, but got '{traj}'.")

        if not os.path.isfile(traj):
            raise FileNotFoundError(f"Error: Input file '{traj}' not found.")
        
        try:
            with h5py.File(traj, "r") as f:
                required_datasets = ["positions", "forces"]
                for dataset in required_datasets:
                    if dataset not in f:
                        raise ValueError(f"Error: Required dataset '{dataset}' not found in '{traj}'.")

                metadata_str = f.attrs.get("metadata")
                if not metadata_str:
                    raise ValueError(f"Error: Required attribute 'metadata' not found in '{traj}'.")
                
                cond = json.loads(metadata_str)
                required_keys = ["symbol", "nsw", "dt", "temperature", "atom_counts", "lattice"]
                for key in required_keys:
                    if key not in cond:
                        raise ValueError(f"Error: Required key '{key}' not found in metadata of '{traj}'.")

        except (IOError, OSError) as e:
            raise IOError(f"Error: Failed to read '{traj}' as an HDF5 file. Reason: {e}")
          
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
        with h5py.File(self.traj, "r") as f:
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
                    f"but found '{self.symbol}' in '{self.traj}'."
                )

            if not np.all(np.abs(self.lattice_parameter - self.site.lattice_parameter) <= eps):
                raise ValueError(
                    f"Lattice parameter mismatch between site object and trajectory file '{self.traj}'."
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
    
    def _get_occupation(self) -> None:
        """
        Analyzes the atomic trajectory to determine site occupations over time,
        applying a transition state (TS) criterion to validate hops.
        This version strictly mimics the calculation flow of the original code.
        """
        with h5py.File(self.traj, 'r') as f:
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
                        
                    if np.isnan(cos_init) or np.isnan(cos_final):
                        if self.verbose:
                            print(f"WARNING: NaN in cos_init/final at step {i} in file {self.traj}")
                    
                    if cos_init > cos_final:
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

            unwrapped_coords[step] = current_unwrapped
            prev_unwrapped = current_unwrapped
            prev_indices = current_indices
        
        self.unwrapped_vacancy_trajectory_coord_cart = unwrapped_coords


class TrajectoryAnalyzer:
    """
    Analyzes and quantifies hopping statistics from Trajectory and Site objects.

    This class takes pre-processed lattice information (from a Site object) and
    a trajectory analysis (from a Trajectory object) to calculate detailed
    statistics about vacancy diffusion, such as hop counts, residence times,
    and mean squared displacement.

    Args:
        trajectory (Trajectory): 
            An initialized Trajectory object containing
            the calculated occupation and vacancy trajectory data.
        site (Site, optional): 
            An initialized Site object containing lattice structure and 
            pre-defined hopping path information. If None, the `site` object 
            associated with the `trajectory` object will be used.
            Defaults to None.
        eps (float, optional): 
            A tolerance value for floating-point comparisons when categorizing hops. 
            Defaults to 1e-3.
        verbose (bool, optional): 
            Verbosity tag. Defaults to True.

    Attributes:
        hopping_history (list[list[dict]]): 
            A list of lists, where each sublist contains 
            a time-ordered sequence of hop dictionaries for a single vacancy.
        counts (np.ndarray): 
            A 2D array of shape (num_vacancies, num_paths)
            counting the occurrences of each pre-defined hop path.
        counts_unknown (np.ndarray): 
            A 2D array counting occurrences of dynamically discovered "unknown" hop paths.
        residence_time (np.ndarray): 
            A 2D array of shape (num_vacancies, num_sites)
            storing the total time (in ps) each vacancy spent at each inequivalent site.
        msd_rand (float): 
            The Mean Squared Displacement calculated from all
            observed hops based on random walk theory.
    """
    def __init__(self,
                 trajectory,
                 site=None,
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
        
        # site(lattice) 정보
        self.path = site.path
        self.path_name = site.path_name
        self.site_name = site.site_name
        self.path_distance = np.array([p['distance'] for p in site.path])
        self.lattice_sites_info = site.lattice_sites

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
        """Prints a comprehensive summary of the hopping analysis."""
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
    """
    Analyzes vacancy encounters and calculates diffusion correlation factors.

    An "encounter" is defined as the sequence of hops made by a single atom,
    mediated by one or more vacancies, before it loses correlation with the
    initial vacancy. This class processes the hopping history from a
    TrajectoryAnalyzer to identify these encounters and compute properties like
    the mean squared displacement (MSD) and the correlation factor (f).

    Args:
        analyzer (TrajectoryAnalyzer): 
            An initialized TrajectoryAnalyzer object containing the full hopping statistics.
        use_incomplete_encounter (bool, optional): 
            If True, encounters that were still in progress at the end of the simulation are 
            included in the statistical analysis. Defaults to True.
        verbose (bool, optional): 
            Verbosity flag. Defaults to True.
    """
    def __init__(self,
                 analyzer,
                 use_incomplete_encounter: bool = True,
                 verbose: bool = True):
        self.analyzer = analyzer
        self.trajectory = analyzer.trajectory
        self.site = analyzer.site # [NEW] site 객체 직접 참조
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
        
        self.encounter_all = self.encounter_complete
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
        # [MODIFIED] Trajectory 객체의 unwrapped 데이터 직접 사용
        source = self.trajectory.unwrapped_vacancy_trajectory_coord_cart
        num_steps = self.trajectory.num_steps
        num_vac = self.trajectory.num_vacancies
        
        # 딕셔너리를 (num_steps, num_vac, 3) 모양의 배열로 변환
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
                            self.encounter_complete.append(encounter_match.copy())
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
                        self.encounter_complete.append(encounter_match.copy())
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

    def summary(self):
        """Prints a comprehensive, formatted summary of the encounter analysis."""
        print("# Encounter Analysis ###")
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
        
        print(f"# Pathwise Counts in Encounters ###")
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
