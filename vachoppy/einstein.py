import os
import h5py
import json
import numpy as np
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
from typing import List, Union, Optional
from joblib import Parallel, delayed, cpu_count
from pathlib import Path
from tabulate import tabulate

from vachoppy.trajectory import TrajBundle
from vachoppy.utils import monitor_performance

# ANSI color codes
BOLD = '\033[1m'
GREEN = '\033[92m'
RED = '\033[91m'
RESET = '\033[0m'

class MSDCalculator:
    """
    Analyzes a single HDF5 trajectory to calculate Mean Squared Displacement (MSD)
    and diffusivity for a specific atomic species.

    The primary workflow is to initialize the class and then call the `.calculate()`
    method. The results (MSD, diffusivity) are stored as attributes and can be
    plotted using `.plot_msd()`. This class is designed to be memory-efficient
    by discarding the large position array after the MSD calculation is complete.

    Args:
        traj_file (str): 
            Path to the HDF5 trajectory file.
        symbol (str): 
            The chemical symbol of the target diffusing species to analyze.
        skip (float, optional): 
            Initial time in picoseconds (ps) to skip for thermal equilibration.
            Defaults to 0.0.
        segment_length (Optional[float], optional): 
            The time length of each segment in picoseconds (ps) for statistical averaging.
            If None, the entire trajectory after 'skip' is treated as a single segment.
            Defaults to None.
        start (float, optional): 
            The start time in picoseconds (ps) for the linear fitting range of the MSD.
            Defaults to 1.0.
        end (Optional[float], optional): 
            The end time in picoseconds (ps) for the linear fitting range.
            If None, the end of the trajectory is used. Defaults to None.
        verbose (bool, optional): 
            Verbosity flag. Defaults to True.
    """
    def __init__(self,
                 traj_file: str,
                 symbol: str,
                 skip: float = 0.0,
                 segment_length: Optional[float] = None,
                 start: float = 1.0,
                 end: Optional[float] = None,
                 verbose: bool = True):
        
        self.traj_file = traj_file
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
        if not os.path.isfile(self.traj_file):
            raise FileNotFoundError(f"Input file not found: {self.traj_file}")
            
        with h5py.File(self.traj_file, 'r') as f:
            meta = json.loads(f.attrs['metadata'])
            self.dt = meta['dt']
            self.nsw = meta['nsw']
            self.lattice = np.array(meta['lattice'], dtype=np.float64)
            self.temperature = meta.get('temperature', np.nan)
            
            full_atom_counts = meta['atom_counts']
            if self.symbol not in full_atom_counts:
                raise ValueError(f"Symbol '{self.symbol}' not found in {self.traj_file}. "
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
                print(f"Warning: Fitting range for {os.path.basename(self.traj_file)} is too small. Cannot calculate diffusivity.")
            self.diffusivity, self.intercept = np.nan, np.nan
            return

        slope, self.intercept = np.polyfit(time_fit_fs, msd_fit, 1)
        self.diffusivity = slope * (1e-5 / 6.0) # in m^2/s

    def plot_msd(self, ax=None, **kwargs):
        """
        Plots the calculated MSD vs. time on a given matplotlib axis.

        If a linear fit was performed, the fit line and range are also plotted.

        Args:
            ax (matplotlib.axes.Axes, optional): A matplotlib axes object to plot on.
                If None, a new figure and axes are created. Defaults to None.
            **kwargs: Additional keyword arguments passed to `ax.plot()`.

        Returns:
            matplotlib.lines.Line2D: The Line2D object for the plotted MSD curve.
        """
        if self.msd is None: self._calculate_msd()
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 5))
            ax.set_xlabel("Time (ps)", fontsize=13)
            ax.set_ylabel(r"MSD (Å$^2$)", fontsize=13)
        
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
        
        return line


def _run_single_msd_task(args):
    """[Parallel Worker] Initializes and runs an MSDCalculator for a single file."""
    traj_file, symbol, skip, segment_length, start, end = args
    try:
        calc = MSDCalculator(traj_file, symbol, skip, segment_length, start, end, verbose=False)
        calc.calculate()
        return calc
    except Exception as e:
        if "verbose" not in str(e): # Avoid printing verbose flag errors
             print(f"Warning: Failed to process {traj_file}. Error: {e}")
        return None


class MSDEnsemble:
    """
    Analyzes an ensemble of HDF5 trajectories, typically at different temperatures,
    to calculate temperature-dependent diffusivity and Arrhenius parameters.
    """
    def __init__(self,
                 path: str,
                 symbol: str,
                 skip: float = 0.0,
                 segment_length: Optional[Union[float, List[float]]] = None,
                 start: float = 1.0,
                 end: Optional[float] = None,
                 verbose: bool = True,
                 **kwargs):
        
        self.path = path
        self.symbol = symbol
        self.skip = skip
        self.segment_length = segment_length
        self.start = start
        self.end = end
        self.verbose = verbose
        self.kwargs = kwargs
        
        bundle_keys = ['prefix', 'depth', 'eps']
        bundle_kwargs = {key: kwargs[key] for key in bundle_keys if key in kwargs}
        self.bundle = TrajBundle(path=self.path, symbol=self.symbol, verbose=False, **bundle_kwargs)

        self.temperatures = np.array(self.bundle.temperatures, dtype=np.float64)
        self.all_traj_paths = [path for temp_paths in self.bundle.traj for path in temp_paths]
        
        if isinstance(self.segment_length, (list, np.ndarray)):
            if len(self.segment_length) != len(self.temperatures):
                raise ValueError(f"Length of `segment_length` ({len(self.segment_length)}) "
                                 f"must match the number of temperatures ({len(self.temperatures)}).")
                
        self.calculators: List[MSDCalculator] = []
        self.diffusivities: np.ndarray = None
        self.Ea: float = None
        self.D0: float = None
        self.R2: float = None
    
    # @monitor_performance
    def calculate(self, n_jobs: int = -1):
        """
        Runs the full analysis pipeline in parallel for all trajectories.

        This method orchestrates the MSD and diffusivity calculation for each
        trajectory file found by `TrajBundle`. It then aggregates the results
        and performs an Arrhenius fit if multiple temperatures are present.

        Args:
            n_jobs (int, optional): Number of CPU cores for parallel processing.
                -1 uses all available cores. Defaults to -1.
        """
        tasks = []
        for i, temp_group in enumerate(self.bundle.traj):
            if isinstance(self.segment_length, (list, np.ndarray)):
                seg_len_for_temp = self.segment_length[i]
            else:
                seg_len_for_temp = self.segment_length

            for traj_file in temp_group:
                tasks.append(
                    (traj_file, self.symbol, self.skip, seg_len_for_temp, self.start, self.end)
                )
        
        results = Parallel(n_jobs=n_jobs)(
            delayed(_run_single_msd_task)(task)
            for task in tqdm(tasks, 
                             desc=f'{RED}{BOLD}Progress{RESET}',
                             bar_format='{l_bar}%s{bar:20}%s{r_bar}' % (GREEN, RESET),
                             ascii=False,
                             disable=not self.verbose)
        )
        
        successful_results = [res for res in results if res is not None]
        path_order = {path: i for i, path in enumerate(self.all_traj_paths)}
        successful_results.sort(key=lambda calc: path_order.get(calc.traj_file, -1))
        self.calculators = successful_results
        
        self._aggregate_results()
        
        if len(self.temperatures) >= 2:
            self._fit_arrhenius()

    def _aggregate_results(self):
        """Aggregates diffusivity results from all calculators, grouped by temperature."""
        calc_temps = np.array([calc.temperature for calc in self.calculators])
        diffs = np.array([calc.diffusivity for calc in self.calculators])
        
        temp_avg_D = []
        for T in self.temperatures:
            mask = np.isclose(calc_temps, T)
            if np.any(mask):
                temp_avg_D.append(np.nanmean(diffs[mask]))
            else:
                temp_avg_D.append(np.nan)
        self.diffusivities = np.array(temp_avg_D, dtype=np.float64)
        
    def _fit_arrhenius(self):
        """Performs an Arrhenius fit on the temperature-dependent diffusivity data."""
        kb = 8.61733326e-5
        temps = np.asarray(self.temperatures, dtype=np.float64)
        diffs = np.asarray(self.diffusivities, dtype=np.float64)
        valid_mask = ~np.isnan(diffs) & (diffs > 0)
        
        if np.sum(valid_mask) < 2:
            if self.verbose: print("Warning: Less than 2 valid data points for Arrhenius fit. Skipping.")
            self.Ea, self.D0, self.R2 = np.nan, np.nan, np.nan
            return
            
        temps_valid = temps[valid_mask]
        D_valid = diffs[valid_mask]
        
        x = 1 / temps_valid
        y = np.log(D_valid)
        slope, intercept = np.polyfit(x, y, 1)
        
        self.Ea = -slope * kb
        self.D0 = np.exp(intercept)
        
        y_pred = slope * x + intercept
        ss_total = np.sum((y - np.mean(y))**2)
        if ss_total < 1e-12: self.R2 = 1.0
        else: self.R2 = 1 - np.sum((y - y_pred)**2) / ss_total
        
    def plot_msd(self, **kwargs):
        """
        Plots the temperature-averaged Mean Squared Displacement (MSD) for each
        unique temperature in the ensemble.
        """
        if not self.calculators:
            raise RuntimeError("Please call the .calculate() method first.")

        fig, ax = plt.subplots(figsize=(7, 6))
        cmap = plt.get_cmap("viridis", len(self.temperatures))

        for i, temp in enumerate(self.temperatures):
            calcs_for_temp = [calc for calc in self.calculators if np.isclose(calc.temperature, temp)]
            if not calcs_for_temp: continue

            timestep = calcs_for_temp[0].timestep
            all_msds_for_temp = [calc.msd for calc in calcs_for_temp]
            mean_msd = np.mean(np.array(all_msds_for_temp), axis=0)

            ax.plot(timestep, mean_msd, color=cmap(i), label=f"{temp:.0f} K")
            
            if self.diffusivities is not None and not np.isnan(self.diffusivities[i]):
                slope = self.diffusivities[i] * 6 / 1e-5
                intercept = np.mean([c.intercept for c in calcs_for_temp if c.intercept is not None])
                
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
        plt.show()

    def plot_D(self, **kwargs):
        """Plots the Arrhenius plot of D vs. 1/T."""
        if len(self.temperatures) < 2:
            print("Warning: Cannot create Arrhenius plot with less than 2 temperatures.")
            return

        fig, ax = plt.subplots(figsize=(6, 5))
        kb = 8.61733326e-5
        
        valid_mask = ~np.isnan(self.diffusivities) & (self.diffusivities > 0)
        temps_valid = self.temperatures[valid_mask]
        D_valid = self.diffusivities[valid_mask]
        
        x_points = 1000 / temps_valid
        y_points = D_valid
        
        ax.scatter(x_points, y_points, c=temps_valid, cmap='viridis', s=50, zorder=3)
        ax.set_yscale('log')
        
        if self.Ea is not None and not np.isnan(self.Ea):
            x_fit_temps = np.linspace(min(temps_valid), max(temps_valid), 100)
            y_fit = self.D0 * np.exp(-self.Ea / (kb * x_fit_temps))
            ax.plot(1000 / x_fit_temps, y_fit, 'k--', lw=1.5, label=f"Fit (Ea={self.Ea:.2f} eV)")
            
            text = f"$D_0 = {self.D0:.2e}$ m$^2$/s\n$E_a = {self.Ea:.2f}$ eV\n$R^2 = {self.R2:.3f}$"
            ax.text(0.05, 0.05, text, transform=ax.transAxes, 
                    va='bottom', ha='left', 
                    bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

        ax.set_xlabel("1000 / T (K$^{-1}$)", fontsize=13)
        ax.set_ylabel(r"Diffusivity, D (m$^2$/s)", fontsize=13)
        ax.grid(True, which='both', linestyle='--')
        ax.legend()
        plt.tight_layout()
        plt.show()
        
    def summary(self):
        """Prints a summary of the ensemble analysis."""
        print("\n" + "="*60)
        print(f"{' ' * 15}MSD Ensemble Analysis Summary")
        print("="*60)
        
        headers = ["Temp (K)", "Avg. Diffusivity (m^2/s)"]
        table_data = []
        for T, D in zip(self.temperatures, self.diffusivities):
            table_data.append([f"{T:.0f}", f"{D:.3e}" if not np.isnan(D) else "N/A"])
        
        print(tabulate(table_data, headers=headers, tablefmt="simple"))
        
        if self.Ea is not None and not np.isnan(self.Ea):
            print("\n-- Arrhenius Fit Results --")
            print(f"  - Activation Energy (Ea) : {self.Ea:.3f} eV")
            print(f"  - Pre-factor (D0)        : {self.D0:.3e} m^2/s")
            print(f"  - R-squared              : {self.R2:.4f}")
        print("="*60)
        
    def save_parameters(self, 
                        filename: str = "einstein.json") -> None:
        """Saves the calculated diffusivity data and Arrhenius fit results to a JSON file.

        The output JSON includes a 'description' key that explains each parameter,
        making the file self-documenting. It automatically converts NumPy arrays
        to lists for JSON compatibility.

        Args:
            filename (str, optional): The name of the output JSON file.
                Defaults to "einstein.json".
        """
        if self.diffusivities is None:
            raise RuntimeError("Cannot save parameters. Please run the .calculate() method first.")

        description = {
            'temperatures': 'List of temperatures (K) at which simulations were run.',
            'D'   : 'Temperature-dependent tracer diffusivity (m^2/s).',
            'Ea_D': 'Activation energy (eV) from the Arrhenius fit of D.',
            'D0'  : 'Pre-exponential factor (m^2/s) from the Arrhenius fit of D.',
        }

        contents = {
            'symbol': self.symbol,
            'temperatures': self.temperatures.tolist(),
            'description': description
        }

        params_to_save = ['diffusivities', 'Ea', 'D0']
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
        

def Einstein(path: str,
             symbol: str,
             skip: float = 0.0,
             segment_length: Optional[Union[float, List[float]]] = None,
             start: float = 1.0,
             end: Optional[float] = None,
             **kwargs) -> Union[MSDEnsemble, MSDCalculator]:
    """
    Factory function for Einstein relation analysis.

    Creates an `MSDEnsemble` for a directory path or an `MSDCalculator`
    for a single file path. This function serves as the main user entry point.

    Args:
        path (str):
            Path to a directory containing trajectory files or to a single file.
        symbol (str):
            The chemical symbol of the target diffusing species to analyze.
        skip (float, optional):
            Initial time in picoseconds (ps) to skip for equilibration. Defaults to 0.0.
        segment_length (Optional[Union[float, List[float]]], optional):
            Time length of each segment in picoseconds (ps) for statistical averaging.
            If a list, must match the number of temperatures. If None, the entire
            trajectory is used as one segment. Defaults to None.
        start (float, optional):
            Start time in picoseconds (ps) for the linear fitting range. Defaults to 1.0.
        end (Optional[float], optional):
            End time in picoseconds (ps) for the linear fitting range. If None,
            the end of the trajectory/segment is used. Defaults to None.
        **kwargs:
            Additional keyword arguments passed to `TrajBundle` (e.g., `prefix`, `depth`).

    Returns:
        Union[MSDEnsemble, MSDCalculator]: An initialized analysis object.
    """
    p = Path(path)
    if p.is_dir():
        return MSDEnsemble(path, symbol, skip, segment_length, start, end, **kwargs)
    elif p.is_file():
        if isinstance(segment_length, (list, np.ndarray)):
            raise TypeError("For a single file analysis, `segment_length` must be a float or None.")
        return MSDCalculator(path, symbol, skip, segment_length, start, end, **kwargs)
    else:
        raise FileNotFoundError(f"Path not found: {path}")