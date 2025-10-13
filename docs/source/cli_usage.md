# Command-Line Interface (CLI)

In addition to its Python library, `VacHopPy` provides a powerful Command-Line Interface (CLI) to perform key tasks directly from your terminal. This allows you to run analyses quickly without writing custom scripts.

The general command structure is:
```bash
vachoppy <command> [options]
````

-----

## Available Commands

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

-----

## Getting Help for a Specific Command

Each command has its own set of specific options and arguments. You can view the detailed help for any command by adding `-h` after the command name.

For example, to see all options for the `trajectory` command, run:

```bash
vachoppy trajectory -h
```
