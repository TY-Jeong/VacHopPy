"""
vachoppy.einstein
=================

Provides tools for calculating diffusivity from molecular dynamics trajectories
using the Einstein relation (MSD = 6Dt).

This module contains the primary classes for Mean Squared Displacement (MSD)
analysis. It is designed to handle both individual simulations and ensembles of
simulations at various temperatures.

Main Components
---------------
- **MSDCalculator**: A class for analyzing a single trajectory file. It calculates
  the MSD over time and fits it to determine the diffusivity.
- **MSDEnsemble**: A class that manages a collection of trajectories (typically
  at different temperatures). It runs `MSDCalculator` for each, aggregates the
  results, and performs Arrhenius analysis to find the activation energy for
  diffusion.
- **Einstein**: The main factory function and user entry point. It automatically
  creates either an `MSDCalculator` or an `MSDEnsemble` object based on the
  input path (file or directory), simplifying the analysis workflow.

Typical Usage
-------------
The primary way to use this module is through the `Einstein` factory function:

.. code-block:: python

    from vachoppy import Einstein

    # Analyze a single simulation at 1000 K
    analyzer_single = Einstein(
        path_traj='TRAJ_O_1000K.h5',
        symbol='O',
        skip=1.0  # Skip first 1 ps for equilibration
    )
    analyzer_single.calculate()
    analyzer_single.summary()

    # Analyze a multi-temperature ensemble
    analyzer_ensemble = Einstein(
        path_traj='path/to/all_trajectories/',
        symbol='O',
        skip=1.0
    )
    analyzer_ensemble.calculate()
    analyzer_ensemble.plot_D() # Show Arrhenius plot
"""

from __future__ import annotations

__all__ =['MSDCalculator', 'MSDEnsemble', 'Einstein']

import os
import h5py
import json
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from tqdm.auto import tqdm
from tabulate import tabulate
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from joblib import Parallel, delayed
from typing import List, Union, Optional

from vachoppy.trajectory import TrajectoryBundle
from vachoppy.utils import monitor_performance


class MSDCalculator:
    """"Analyzes a trajectory to calculate Mean Squared Displacement and diffusivity.

    This class reads a single HDF5 trajectory file for a specific atomic species
    to calculate its Mean Squared Displacement (MSD) and subsequently its
    diffusivity via the Einstein relation.

    The primary workflow is to initialize the class and then call the `.calculate()`
    method. Results are stored as attributes and can be plotted using `.plot_msd()`
    or summarized with `.summary()`. The class is designed to be memory-efficient
    by discarding the large position array after the MSD calculation is complete.

    Args:
        path_traj (str):
            Path to the HDF5 trajectory file.
        symbol (str):
            The chemical symbol of the target diffusing species to analyze.
        skip (float, optional):
            Initial time in picoseconds (ps) to skip for thermal equilibration.
            Defaults to 0.0.
        segment_length (float | None, optional):
            Time length in picoseconds (ps) for each segment used in statistical
            averaging. If None, the entire trajectory after 'skip' is treated
            as a single segment. Defaults to None.
        start (float, optional):
            Start time in picoseconds (ps) for the linear fitting range of
            the MSD curve. Defaults to 1.0.
        end (float | None, optional):
            End time in picoseconds (ps) for the linear fitting range. If None,
            the end of the trajectory is used. Defaults to None.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.

    Attributes:
        msd (numpy.ndarray):
            The calculated Mean Squared Displacement values in Å².
        timestep (numpy.ndarray):
            The time values in picoseconds (ps) corresponding to each MSD point.
        diffusivity (float):
            The calculated diffusivity in m²/s, derived from the linear fit.
        intercept (float):
            The y-intercept of the linear fit of the MSD curve.
        temperature (float):
            Temperature in Kelvin, read from the trajectory file's metadata.
        dt (float):
            Timestep in femtoseconds (fs), read from metadata.
        nsw (int):
            Total number of steps in the trajectory, read from metadata.
        lattice (numpy.ndarray):
            The 3x3 lattice matrix, read from metadata.

    Raises:
        FileNotFoundError:
            If the input trajectory file is not found.
        ValueError:
            If the specified `symbol` is not present in the trajectory file,
            or if `segment_length` is too small.

    Examples:
        >>> msd_calc = MSDCalculator(
        ...     path_traj='TRAJ_Li_1000K.h5',
        ...     symbol='Li',
        ...     skip=10.0,
        ...     start=5.0,
        ...     end=20.0
        ... )
        >>> msd_calc.calculate()
        >>> msd_calc.summary()
        # Access results directly:
        # print(f"Calculated Diffusivity: {msd_calc.diffusivity:.3e} m^2/s")
    """
    
    def __init__(self,
                 path_traj: str,
                 symbol: str,
                 skip: float = 0.0,
                 segment_length: float | None = None,
                 start: float = 1.0,
                 end: float | None = None,
                 verbose: bool = True):
        
        self.path_traj = path_traj
        self.symbol = symbol
        self.skip = skip
        self.segment_length = segment_length
        self.start = start
        self.end = end
        self.verbose = verbose

        self.dt: float = None
        self.nsw: int = None
        self.lattice: np.ndarray = None
        self.temperature: float = None
        self._positions_unwrapped: np.ndarray = None
        
        self.msd: np.ndarray = None
        self.timestep: np.ndarray = None
        self.diffusivity: float = None
        self.intercept: float = None
        
        self._read_and_prepare()

    def _read_and_prepare(self):
        """Reads the HDF5 file, validates metadata, and loads atomic positions."""
        if not os.path.isfile(self.path_traj):
            raise FileNotFoundError(f"Input file not found: {self.path_traj}")
            
        with h5py.File(self.path_traj, 'r') as f:
            meta = json.loads(f.attrs['metadata'])
            self.dt = meta['dt']
            self.nsw = meta['nsw']
            self.lattice = np.array(meta['lattice'], dtype=np.float64)
            self.temperature = meta.get('temperature', np.nan)
            
            full_atom_counts = meta['atom_counts']
            if self.symbol not in full_atom_counts:
                raise ValueError(f"Symbol '{self.symbol}' not found in {self.path_traj}. "
                                 f"Available symbols: {list(full_atom_counts.keys())}")

            if meta.get('symbol') == self.symbol:
                self._positions_unwrapped = f['positions'][:].astype(np.float64)
            else:
                symbols_order = sorted(full_atom_counts.keys())
                start_idx = 0
                for sym in symbols_order:
                    num_atoms = full_atom_counts[sym]
                    if sym == self.symbol:
                        self._positions_unwrapped = f['positions'][:, start_idx : start_idx + num_atoms, :].astype(np.float64)
                        break
                    start_idx += num_atoms
        
        if self.end is None:
            self.end = self.nsw * self.dt / 1000.0

    def calculate(self):
        """Runs the MSD and diffusivity calculation pipeline."""
        self._calculate_msd()
        self._calculate_diffusivity()
        return self

    def _calculate_msd(self):
        """Calculates the Mean Squared Displacement (MSD) from the trajectory."""
        skip_steps = round(self.skip * 1000 / self.dt)
        total_steps = self.nsw - skip_steps
        
        if self.segment_length is None or self.segment_length * 1000 / self.dt >= total_steps:
            segments = 1
            seg_len_steps = total_steps
        else:
            seg_len_steps = round(self.segment_length * 1000 / self.dt)
            if seg_len_steps == 0: raise ValueError("`segment_length` is too small, resulting in zero steps per segment.")
            segments = total_steps // seg_len_steps
        
        if segments == 0: raise ValueError("`skip` or `segment_length` is too large, resulting in zero segments.")

        msd_segments = []
        for i in range(segments):
            start_frame = skip_steps + i * seg_len_steps
            end_frame = start_frame + seg_len_steps
            segment_pos = self._positions_unwrapped[start_frame:end_frame]
            
            displacement_frac = segment_pos - segment_pos[0]
            displacement_cart = np.dot(displacement_frac, self.lattice)
            
            squared_disp = np.sum(displacement_cart**2, axis=2)
            msd_per_atom = np.mean(squared_disp, axis=1)
            msd_segments.append(msd_per_atom)
            
        self.msd = np.mean(np.array(msd_segments), axis=0)
        self.timestep = np.arange(len(self.msd)) * self.dt / 1000.0 # in ps
        del self._positions_unwrapped # RAM Optimization

    def _calculate_diffusivity(self):
        """Calculates diffusivity by linear fitting of the MSD curve."""
        if self.msd is None: self._calculate_msd()

        start_idx = int(round(self.start * 1000 / self.dt))
        end_idx = int(round(self.end * 1000 / self.dt)) if self.end is not None else len(self.msd)
        end_idx = min(end_idx, len(self.msd))

        time_fit_fs = self.timestep[start_idx:end_idx] * 1000.0
        msd_fit = self.msd[start_idx:end_idx]

        if len(time_fit_fs) < 2:
            if self.verbose:
                print(f"Warning: Fitting range for {os.path.basename(self.path_traj)} is too small. Cannot calculate diffusivity.")
            self.diffusivity, self.intercept = np.nan, np.nan
            return

        slope, self.intercept = np.polyfit(time_fit_fs, msd_fit, 1)
        self.diffusivity = slope * (1e-5 / 6.0) # in m^2/s

    def plot_msd(self,
                 disp: bool = True,
                 save: bool = True,
                 filename: str = 'msd.png',
                 dpi: int = 300,
                 ax: Axes | None = None,
                 **kwargs) -> Line2D:
        """Plots the calculated MSD versus time.

        This method visualizes the MSD curve. If diffusivity has been calculated,
        it also overlays the linear fit and the fitting range boundaries. The plot
        can be displayed interactively, saved to a file, or drawn on an existing
        matplotlib Axes object for creating complex figures.

        Args:
            disp (bool, optional):
                If True, displays the plot interactively (`plt.show()`).
                Defaults to True.
            save (bool, optional):
                If True, saves the plot to a file. Defaults to True.
            filename (str, optional):
                The path and name of the file to save the plot.
                Defaults to 'msd.png'.
            dpi (int, optional):
                The resolution in dots per inch for the saved figure.
                Defaults to 300.
            ax (Axes | None, optional):
                A matplotlib Axes object to plot on. If None, a new figure and
                axes are created automatically. This is useful for creating
                subplots. Defaults to None.
            **kwargs:
                Additional keyword arguments to be passed directly to the
                `ax.plot()` function for the main MSD curve (e.g., `label`,
                `color`, `linewidth`).

        Returns:
            Line2D:
                The Line2D object for the plotted MSD curve, which can be used
                for further customization (e.g., adding to a legend).

        Examples:
            >>> # Simple plot, shown and saved as 'msd.png'
            >>> msd_calc.calculate()
            >>> msd_calc.plot_msd(label='My System')

            >>> # Plotting on an existing subplot figure
            >>> import matplotlib.pyplot as plt
            >>> fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
            >>> msd_calc_1.plot_msd(ax=ax1, disp=False, save=False, label='Run 1')
            >>> msd_calc_2.plot_msd(ax=ax2, disp=False, save=False, label='Run 2')
            >>> fig.suptitle("Comparison of MSD")
            >>> fig.tight_layout()
            >>> fig.savefig("msd_comparison.png")
        """
        if self.msd is None: self._calculate_msd()
        
        fig = None
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 5))
            ax.set_xlabel("Time (ps)", fontsize=13)
            ax.set_ylabel(r"MSD (Å$^2$)", fontsize=13)
        else:
            fig = ax.get_figure()
        
        line, = ax.plot(self.timestep, self.msd, **kwargs)
        
        if self.diffusivity is not None and not np.isnan(self.diffusivity):
            start_idx = int(round(self.start * 1000 / self.dt))
            end_idx = int(round(self.end * 1000 / self.dt)) if self.end is not None else len(self.msd)
            end_idx = min(end_idx, len(self.msd))
            time_fit_ps = self.timestep[start_idx:end_idx]
            fit_line = (self.diffusivity * 6 / 1e-5) * (time_fit_ps * 1000) + self.intercept
            ax.plot(time_fit_ps, fit_line, 'k--', lw=1.5)
            
            end_time = self.end if self.end is not None else self.timestep[-1]
            ax.axvline(self.start, color='grey', linestyle=':', lw=1)
            ax.axvline(end_time, color='grey', linestyle=':', lw=1)
            
        if save: 
            if fig: fig.tight_layout()
            plt.savefig(filename, dpi=dpi)
        if disp: 
            if fig: fig.tight_layout()
            plt.show()
        if fig: plt.close(fig)
        return line

    def summary(self):
        """Prints a formatted summary of the MSD analysis results."""
        if not hasattr(self, 'diffusivity') or self.diffusivity is None:
            self.calculate()
            
        header_width = 60
        title = "MSD Analysis Summary (Single)"
        padding = (header_width - len(title)) // 2
        
        print("\n" + "=" * header_width)
        print(" " * padding + title)
        print("=" * header_width)
        
        print("\n-- Input Parameters --")
        print(f"  - Trajectory File : {os.path.basename(self.path_traj)}")
        print(f"  - Target Symbol   : {self.symbol}")
        print(f"  - Temperature     : {self.temperature:.1f} K")
        print(f"  - Skip Time       : {self.skip:.2f} ps")
        print(f"  - Segment Length  : {'Full' if self.segment_length is None else f'{self.segment_length:.2f} ps'}")

        print("\n-- Fitting Range --")
        end_time = self.end if self.end is not None else self.timestep[-1]
        print(f"  - Start Time      : {self.start:.2f} ps")
        print(f"  - End Time        : {end_time:.2f} ps")
        
        print("\n-- Results --")
        if not np.isnan(self.diffusivity):
            print(f"  - Diffusivity (D) : {self.diffusivity:.3e} m^2/s")
        else:
            print("  - Diffusivity (D) : Not calculated (check warnings).")   
        print("=" * header_width + "\n")


def _run_single_msd_task(args):
    """[Parallel Worker] Initializes and runs an MSDCalculator for a single file."""
    path_traj, symbol, skip, segment_length, start, end = args
    try:
        calc = MSDCalculator(path_traj, symbol, skip, segment_length, start, end, verbose=False)
        calc.calculate()
        return calc
    except Exception as e:
        if "verbose" not in str(e):
             print(f"Warning: Failed to process {path_traj}. Error: {e}")
        return None


class MSDEnsemble:
    """Analyzes an ensemble of trajectories to find temperature-dependent diffusivity.

    This class manages multiple `MSDCalculator` instances, typically for simulations
    run at different temperatures. It calculates the diffusivity for each trajectory
    (or temperature group) in parallel and then performs an Arrhenius fit to
    determine the activation energy (Ea_D) and pre-exponential factor (D0).

    The main workflow is to initialize the class, call `.calculate()` to run the
    analysis, and then use plotting methods like `.plot_msd()` and `.plot_D()`
    or `.summary()` to view the results.

    Args:
        path_traj (str):
            Path to a directory containing HDF5 trajectory files.
        symbol (str):
            The chemical symbol of the target diffusing species to analyze.
        skip (float, optional):
            Initial time in picoseconds (ps) to skip for thermal equilibration.
            Defaults to 0.0.
        segment_length (float | list[float] | None, optional):
            Time length in ps for statistical averaging segments. Can be a single
            float for all temperatures or a list of floats corresponding to each
            temperature. If None, the entire trajectory is used. Defaults to None.
        start (float, optional):
            Start time in picoseconds (ps) for the linear fitting range of
            the MSD curve. Defaults to 1.0.
        end (float | None, optional):
            End time in picoseconds (ps) for the linear fitting range. If None,
            the end of each trajectory is used. Defaults to None.
        verbose (bool, optional):
            Verbosity flag. Defaults to True.
        **kwargs (Any):
            Additional keyword arguments passed to the `TrajectoryBundle` class
            for file discovery (e.g., `prefix`, `depth`).

    Attributes:
        D (numpy.ndarray):
            An array of the average diffusivity (m²/s) for each temperature.
        Ea_D (float):
            Activation energy (eV) from the Arrhenius fit.
        D0 (float):
            Pre-exponential factor (m²/s) from the Arrhenius fit.
        R2 (float):
            The R-squared value indicating the goodness of the Arrhenius fit.
        calculators (list[MSDCalculator]):
            A list of the successfully completed `MSDCalculator` instances.
        temperatures (numpy.ndarray):
            An array of the unique temperatures found in the ensemble.

    Raises:
        ValueError:
            If the length of `segment_length` (if provided as a list) does not
            match the number of unique temperatures found.
    """
    def __init__(self,
                 path_traj: str,
                 symbol: str,
                 skip: float = 0.0,
                 segment_length: float | list[float] | None = None,
                 start: float = 1.0,
                 end: float | None = None,
                 verbose: bool = True,
                 **kwargs):
        
        self.path_traj = path_traj
        self.symbol = symbol
        self.skip = skip
        self.segment_length = segment_length
        self.start = start
        self.end = end
        self.verbose = verbose
        self.kwargs = kwargs
        
        bundle_keys = ['prefix', 'depth', 'eps']
        bundle_kwargs = {key: kwargs[key] for key in bundle_keys if key in kwargs}
        self.bundle = TrajectoryBundle(path_traj=self.path_traj, 
                                       symbol=self.symbol, 
                                       verbose=False, 
                                       **bundle_kwargs)

        self.temperatures = np.array(self.bundle.temperatures, dtype=np.float64)
        self.all_traj_paths = [
            path for temp_paths in self.bundle.traj for path in temp_paths
        ]
        
        if isinstance(self.segment_length, (list, np.ndarray)):
            if len(self.segment_length) != len(self.temperatures):
                raise ValueError(
                    f"Length of `segment_length` ({len(self.segment_length)}) "
                    f"must match the number of temperatures ({len(self.temperatures)})."
                )
                
        self.calculators: List[MSDCalculator] = []
        self.D: np.ndarray = None
        self.Ea_D: float = None
        self.D0: float = None
        self.R2: float = None
    
    @ monitor_performance
    def calculate(self, n_jobs: int = -1, verbose=True) -> MSDEnsemble:
        """Runs the full analysis pipeline in parallel for all trajectories.

            This method orchestrates the MSD calculation for each trajectory file,
            aggregates the results by temperature, and performs an Arrhenius fit if
            at least two temperatures are present.

            Args:
                n_jobs (int, optional):
                    Number of CPU cores for parallel processing via joblib.
                    -1 uses all available cores. Defaults to -1.

            Returns:
                MSDEnsemble:
                    Returns the instance itself to allow for method chaining
                    (e.g., `ensemble.calculate().summary()`).
        """
        tasks = []
        for i, temp_group in enumerate(self.bundle.traj):
            if isinstance(self.segment_length, (list, np.ndarray)):
                seg_len_for_temp = self.segment_length[i]
            else:
                seg_len_for_temp = self.segment_length

            for path_traj in temp_group:
                tasks.append(
                    (path_traj, self.symbol, self.skip, 
                     seg_len_for_temp, self.start, self.end)
                )
        
        results = Parallel(n_jobs=n_jobs)(
            delayed(_run_single_msd_task)(task)
            for task in tqdm(tasks, 
                             desc=f'Compute MSD',
                             bar_format='{l_bar}{bar:30}{r_bar}',
                             ascii=True)
        )
        
        successful_results = [res for res in results if res is not None]
        path_order = {path: i for i, path in enumerate(self.all_traj_paths)}
        successful_results.sort(key=lambda calc: path_order.get(calc.path_traj, -1))
        self.calculators = successful_results
        
        if self.verbose: 
            print(f"\nAnalysis complete: {len(successful_results)} successful, " +
                  f"{len(results) - len(successful_results)} failed.")
        
        self._aggregate_results()
        
        if len(self.temperatures) >= 2: self._fit_arrhenius()

    def _aggregate_results(self):
        """Aggregates diffusivity results from all calculators, grouped by temperature."""
        temp_avg_D = []
        for temp in self.temperatures:
            calcs_for_temp = [calc for calc in self.calculators 
                              if np.isclose(calc.temperature, temp)]
            if not calcs_for_temp:
                temp_avg_D.append(np.nan)
                continue

            all_msds = [calc.msd for calc in calcs_for_temp]
            max_len = max(len(msd) for msd in all_msds)
            
            msd_sum = np.zeros(max_len, dtype=np.float64)
            frame_counts = np.zeros(max_len, dtype=int)

            for msd in all_msds:
                msd_len = len(msd)
                msd_sum[:msd_len] += msd
                frame_counts[:msd_len] += 1
            
            mean_msd = np.divide(msd_sum, frame_counts, where=frame_counts!=0)
            
            dt = calcs_for_temp[0].dt
            timestep = np.arange(max_len) * dt / 1000.0
            
            start_idx = int(round(self.start / (dt / 1000.0)))
            end_idx = int(round(self.end / (dt / 1000.0))) if self.end is not None else len(mean_msd)
            end_idx = min(end_idx, len(mean_msd))

            time_fit_fs = timestep[start_idx:end_idx] * 1000.0
            msd_fit = mean_msd[start_idx:end_idx]

            if len(time_fit_fs) < 2:
                temp_avg_D.append(np.nan)
                continue

            slope, _ = np.polyfit(time_fit_fs, msd_fit, 1)
            diffusivity = slope * (1e-5 / 6.0)
            temp_avg_D.append(diffusivity)
            
        self.D = np.array(temp_avg_D, dtype=np.float64)

    def _fit_arrhenius(self):
        """Performs an Arrhenius fit on the temperature-dependent diffusivity data."""
        kb = 8.61733326e-5
        temps = np.asarray(self.temperatures, dtype=np.float64)
        diffs = np.asarray(self.D, dtype=np.float64)
        valid_mask = ~np.isnan(diffs) & (diffs > 0)
        
        if np.sum(valid_mask) < 2:
            if self.verbose: print("Warning: Less than 2 valid data points for Arrhenius fit. Skipping.")
            self.Ea_D, self.D0, self.R2 = np.nan, np.nan, np.nan
            return
            
        temps_valid = temps[valid_mask]
        D_valid = diffs[valid_mask]
        
        x = 1 / temps_valid
        y = np.log(D_valid)
        slope, intercept = np.polyfit(x, y, 1)
        
        self.Ea_D = -slope * kb
        self.D0 = np.exp(intercept)
        
        y_pred = slope * x + intercept
        ss_total = np.sum((y - np.mean(y))**2)
        if ss_total < 1e-12: self.R2 = 1.0
        else: self.R2 = 1 - np.sum((y - y_pred)**2) / ss_total
        
    def plot_msd(self, 
                 disp: bool = True,
                 save: bool = True,
                 filename: str = 'msd.png',
                 dpi: int = 300) -> None:
        """Plots temperature-averaged MSD curves for each temperature.

            This method visualizes the averaged MSD curve for each unique temperature
            in the ensemble. The linear fit for each curve is also shown.

            Args:
                disp (bool, optional):
                    If True, displays the plot interactively (`plt.show()`).
                    Defaults to True.
                save (bool, optional):
                    If True, saves the plot to a file. Defaults to True.
                filename (str, optional):
                    The path and name of the file to save the plot.
                    Defaults to 'msd_ensemble.png'.
                dpi (int, optional):
                    The resolution in dots per inch for the saved figure.
                    Defaults to 300.
        """
        if not self.calculators:
            raise RuntimeError("Please call the .calculate() method first.")

        fig, ax = plt.subplots(figsize=(7, 6))
        cmap = plt.get_cmap("viridis", len(self.temperatures))

        for i, temp in enumerate(self.temperatures):
            calcs_for_temp = [
                calc for calc in self.calculators 
                if np.isclose(calc.temperature, temp)
            ]
            if not calcs_for_temp: continue

            all_msds = [calc.msd for calc in calcs_for_temp]
            max_len = max(len(msd) for msd in all_msds)
            msd_sum = np.zeros(max_len)
            frame_counts = np.zeros(max_len, dtype=int)
            
            for msd in all_msds:
                msd_len = len(msd)
                msd_sum[:msd_len] += msd
                frame_counts[:msd_len] += 1
                
            mean_msd = np.divide(msd_sum, frame_counts, where=frame_counts!=0)
            timestep = np.arange(max_len) * calcs_for_temp[0].dt / 1000.0

            ax.plot(timestep, mean_msd, color=cmap(i), label=f"{temp:.0f} K")
            
            if self.D is not None and not np.isnan(self.D[i]):
                slope = self.D[i] * 6 / 1e-5
                intercept = np.mean(
                    [c.intercept for c in calcs_for_temp 
                     if c.intercept is not None and not np.isnan(c.intercept)]
                )
                
                dt_ps = calcs_for_temp[0].dt / 1000.0
                start_idx = int(round(self.start / dt_ps))
                end_ps_val = self.end if self.end is not None else timestep[-1]
                end_idx = int(round(end_ps_val / dt_ps))
                end_idx = min(end_idx, len(timestep))

                time_fit_ps = timestep[start_idx:end_idx]
                fit_line = slope * (time_fit_ps * 1000) + intercept
                ax.plot(time_fit_ps, fit_line, 'k--', lw=1.5)
        
        end_time_ps = self.end if self.end is not None else self.calculators[0].timestep[-1]
        ax.axvline(self.start, color='grey', linestyle=':', lw=1)
        ax.axvline(end_time_ps, color='grey', linestyle=':', lw=1)

        ax.set_xlabel("Time (ps)", fontsize=13)
        ax.set_ylabel(r"MSD (Å$^2$)", fontsize=13)
        ax.legend(title="Temperature", fontsize=9)
        ax.grid(True, linestyle='--')
        plt.tight_layout()
        
        if save: plt.savefig(filename, dpi=dpi)
        if disp: plt.show()
        plt.close(fig)

    def plot_D(self, 
               disp: bool = True,
               save: bool = True,
               filename: str = 'D_atom.png',
               dpi: int = 300) -> None:
        """Plots the Arrhenius plot of diffusivity (D) vs. inverse temperature (1/T).
        
        Args:
            disp (bool, optional):
                If True, displays the plot interactively (`plt.show()`).
                Defaults to True.
            save (bool, optional):
                If True, saves the plot to a file. Defaults to True.
            filename (str, optional):
                The path and name of the file to save the plot.
                Defaults to 'D_arrhenius.png'.
            dpi (int, optional):
                The resolution for the saved figure. Defaults to 300.
        """
        if len(self.temperatures) < 2:
            print("Warning: Cannot create Arrhenius plot with less than 2 temperatures.")
            return

        fig, ax = plt.subplots(figsize=(6, 5))
        kb = 8.61733326e-5
        
        valid_mask = ~np.isnan(self.D) & (self.D > 0)
        temps_valid = self.temperatures[valid_mask]
        D_valid = self.D[valid_mask]
        
        x_points = 1000 / temps_valid
        y_points = D_valid
        
        ax.scatter(x_points, y_points, c=temps_valid, cmap='viridis', s=50, zorder=3)
        ax.set_yscale('log')
        
        if self.Ea_D is not None and not np.isnan(self.Ea_D):
            x_fit_temps = np.linspace(min(temps_valid), max(temps_valid), 100)
            y_fit = self.D0 * np.exp(-self.Ea_D / (kb * x_fit_temps))
            ax.plot(1000 / x_fit_temps, y_fit, 'k--', lw=1.5, label=f"Fit (Ea={self.Ea_D:.2f} eV)")
            
            text = f"$D_0 = {self.D0:.2e}$ m$^2$/s\n$E_a = {self.Ea_D:.2f}$ eV\n$R^2 = {self.R2:.3f}$"
            ax.text(0.05, 0.05, text, transform=ax.transAxes, 
                    va='bottom', ha='left', 
                    bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

        ax.set_xlabel("1000 / T (K$^{-1}$)", fontsize=13)
        ax.set_ylabel(r"Diffusivity, D (m$^2$/s)", fontsize=13)
        ax.grid(True, which='both', linestyle='--')
        ax.legend()
        plt.tight_layout()
        
        if save: plt.savefig(filename, dpi=dpi)
        if disp: plt.show()
        
    def summary(self) -> None:
        """Prints a formatted summary of the ensemble analysis results.

        The summary includes a table of average diffusivities per temperature
        and the results of the Arrhenius fit if applicable.
        """
        print("="*50)
        print(f"{' ' * 10}MSD Ensemble Analysis Summary")
        print("="*50)
        
        headers = ["Temp (K)", "Avg. Diffusivity (m^2/s)"]
        table_data = []
        for T, D in zip(self.temperatures, self.D):
            table_data.append([f"{T:.0f}", f"{D:.3e}" if not np.isnan(D) else "N/A"])
        
        print(tabulate(table_data, headers=headers, tablefmt="simple"))
        
        if self.Ea_D is not None and not np.isnan(self.Ea_D):
            print("\n-- Arrhenius Fit Results --")
            print(f"  - Activation Energy (Ea_D) : {self.Ea_D:.3f} eV")
            print(f"  - Pre-factor (D0)          : {self.D0:.3e} m^2/s")
            print(f"  - R-squared                : {self.R2:.4f}")
        print("="*50)
        
    def save_parameters(self, filename: str = "einstein.json") -> None:
        """Saves calculated diffusivities and Arrhenius parameters to a JSON file.

        The output JSON is self-documenting, including a 'description' key that
        explains each saved parameter. NumPy arrays are converted to lists for
        JSON compatibility.

        Args:
            filename (str, optional):
                The name of the output JSON file.
                Defaults to "einstein_ensemble.json".

        Raises:
            RuntimeError: If `.calculate()` has not been run successfully.

        Examples:
            >>> ensemble = MSDEnsemble(...)
            >>> ensemble.calculate()
            >>> ensemble.save_parameters("results.json")
        """
        if self.D is None:
            raise RuntimeError("Cannot save parameters. Please run the .calculate() method first.")

        description = {
            'temperatures': 'List of temperatures (K) at which simulations were run, (n_temperatures,)',
            'D'   : 'Temperature-dependent tracer diffusivity (m^2/s), (n_temperatures,)',
            'Ea_D': 'Activation energy (eV) from the Arrhenius fit of D.',
            'D0'  : 'Pre-exponential factor (m^2/s) from the Arrhenius fit of D.',
        }

        contents = {
            'symbol': self.symbol,
            'temperatures': self.temperatures.tolist(),
            'description': description
        }

        params_to_save = ['D', 'Ea_D', 'D0']
        json_keys = ['D', 'Ea_D', 'D0']

        for param, key in zip(params_to_save, json_keys):
            if hasattr(self, param):
                value = getattr(self, param)
                if isinstance(value, (np.ndarray, np.generic)):
                    value_list = value.tolist()
                    contents[key] = [None if np.isnan(v) else v for v in value_list] if isinstance(value_list, list) else (None if np.isnan(value_list) else value_list)
                else:
                    contents[key] = value if not (isinstance(value, float) and np.isnan(value)) else None
            else:
                contents[key] = None

        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(contents, f, indent=4)
        
        if self.verbose:
            print(f"Einstein relation parameters saved to '{filename}'")
        

def Einstein(path_traj: str,
             symbol: str,
             skip: float = 0.0,
             segment_length: float | list[float] | None = None,
             start: float = 1.0,
             end: float | None = None,
             **kwargs) -> Union[MSDEnsemble, MSDCalculator]:
    """Creates an appropriate MSD analysis object based on the input path.

    This factory function serves as the main user entry point for Einstein
    relation analysis. It intelligently selects the correct calculator class:
    - `MSDEnsemble`: If `path_traj` is a directory (for multi-temperature analysis).
    - `MSDCalculator`: If `path_traj` is a single file.

    Args:
        path_traj (str):
            Path to a directory containing HDF5 trajectory files or to a single
            HDF5 trajectory file.
        symbol (str):
            The chemical symbol of the target diffusing species to analyze.
        skip (float, optional):
            Initial time in picoseconds (ps) to skip for equilibration.
            Defaults to 0.0.
        segment_length (float | list[float] | None, optional):
            Time length in ps for statistical averaging segments. If a list, it
            must match the number of temperatures (for ensemble analysis). If
            None, the entire trajectory is used. Defaults to None.
        start (float, optional):
            Start time in picoseconds (ps) for the linear fitting range.
            Defaults to 1.0.
        end (float | None, optional):
            End time in picoseconds (ps) for the linear fitting range. If None,
            the end of the trajectory is used. Defaults to None.
        **kwargs (Any):
            Additional keyword arguments passed to `MSDEnsemble` or
            `TrajectoryBundle` (e.g., `prefix`, `depth`).

    Returns:
        Union[MSDEnsemble, MSDCalculator]:
            An initialized analysis object ready for calculation.

    Raises:
        FileNotFoundError:
            If the provided `path_traj` does not exist or is not a file/directory.
        TypeError:
            If `segment_length` is a list when `path_traj` points to a single file.

    Examples:
        >>> # Analyze an ensemble of trajectories in a directory
        >>> ensemble_analyzer = einstein_analyzer(
        ...     path_traj='path/to/trajectories/',
        ...     symbol='Li'
        ... )
        >>> ensemble_analyzer.calculate()

        >>> # Analyze a single trajectory file
        >>> single_analyzer = einstein_analyzer(
        ...     path_traj='path/to/TRAJ_Li_1000K.h5',
        ...     symbol='Li'
        ... )
        >>> single_analyzer.calculate()
    """
    p = Path(path_traj)
    if p.is_dir():
        return MSDEnsemble(path_traj, symbol, skip, segment_length, start, end, **kwargs)
    
    elif p.is_file():
        if isinstance(segment_length, (list, np.ndarray)):
            raise TypeError("For a single file analysis, `segment_length` must be a float or None.")
        return MSDCalculator(path_traj, symbol, skip, segment_length, start, end, **kwargs)
    
    else:
        raise FileNotFoundError(f"Path not found: {path_traj}")
