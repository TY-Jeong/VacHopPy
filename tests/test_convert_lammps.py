import os
import json
import h5py
import pytest
import numpy as np
from pathlib import Path
from numpy.testing import assert_allclose
from vachoppy.core import parse_lammps

@pytest.fixture
def converted_files(tmp_path):
    """
    Fixture to run parse_lammps and return paths for answer and generated files.
    Generated files are created in a temporary directory.
    """
    current_dir = Path(__file__).parent
    lammps_data = (current_dir / 'test_data' /
                   '1.convert' / '2.lammps' / 'lammps.data')
    lammps_dump = (current_dir / 'test_data' /
                   '1.convert' / '2.lammps' / 'lammps.dump')
    traj_O_answer = (current_dir / 'test_data' /
                     '1.convert' / '2.lammps' / 'TRAJ_O.h5')
    traj_Ti_answer = (current_dir / 'test_data' /
                      '1.convert' / '2.lammps' / 'TRAJ_Ti.h5')

    if not lammps_data.is_file():
        pytest.fail(f"Test input file not found: {lammps_data}")
    if not lammps_dump.is_file():
        pytest.fail(f"Test input file not found: {lammps_dump}")

    try:
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        LABEL = 'TEST'
        parse_lammps(lammps_data=str(lammps_data),
                     lammps_dump=str(lammps_dump),
                     atom_style_data='id type x y z',
                     atom_style_dump='id type x y z fx fy fz',
                     atom_symbols={1: 'Ti', 2: 'O'},
                     temperature=2100.0,
                     dt=2.0,
                     label=LABEL,
                     verbose=False)
        os.chdir(original_cwd)
    except Exception as e:
        os.chdir(original_cwd)
        pytest.fail(f"Failed to run parse_lammps(): {e}")

    traj_O_test = tmp_path / f'TRAJ_O_{LABEL}.h5'
    traj_Ti_test = tmp_path / f'TRAJ_Ti_{LABEL}.h5'

    if not traj_O_test.is_file():
        pytest.fail(f"Test output file not found: {traj_O_test}")
    if not traj_Ti_test.is_file():
        pytest.fail(f"Test output file not found: {traj_Ti_test}")

    return {
        "O": {"answer": traj_O_answer, "test": traj_O_test},
        "Ti": {"answer": traj_Ti_answer, "test": traj_Ti_test}
    }


@pytest.mark.parametrize("element", ["O", "Ti"])
def test_convert_metadata(converted_files, element):
    """Tests if the metadata in the generated file matches the answer file."""
    paths = converted_files[element]

    with h5py.File(paths["answer"], 'r') as f_answer, \
         h5py.File(paths["test"], 'r') as f_test:

        expected_metadata = json.loads(f_answer.attrs['metadata'])
        actual_metadata = json.loads(f_test.attrs['metadata'])

        # Compare nested dict separately
        expected_counts = expected_metadata.pop('atom_counts')
        actual_counts = actual_metadata.pop('atom_counts')
        assert actual_counts == expected_counts

        # Compare remaining metadata with approx for floats
        assert actual_metadata == pytest.approx(expected_metadata)


@pytest.mark.parametrize("element", ["O", "Ti"])
def test_convert_positions(converted_files, element):
    """Tests if the 'positions' dataset in the generated file matches the answer."""
    paths = converted_files[element]

    with h5py.File(paths["answer"], 'r') as f_answer, \
         h5py.File(paths["test"], 'r') as f_test:

        expected_positions = f_answer['positions'][:]
        actual_positions = f_test['positions'][:]
        assert_allclose(actual_positions, expected_positions)


@pytest.mark.parametrize("element", ["O", "Ti"])
def test_convert_forces(converted_files, element):
    """Tests if the 'forces' dataset in the generated file matches the answer."""
    paths = converted_files[element]

    with h5py.File(paths["answer"], 'r') as f_answer, \
         h5py.File(paths["test"], 'r') as f_test:

        expected_forces = f_answer['forces'][:]
        actual_forces = f_test['forces'][:]
        assert_allclose(actual_forces, expected_forces)
