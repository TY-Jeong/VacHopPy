import os
import json
import pytest
import numpy as np
from pathlib import Path
from numpy.testing import assert_allclose
from vachoppy.fingerprint import get_fingerprint


def test_fingerprint():
    """Tests the get_fingerprint function against a known answer."""
    current_dir = Path(__file__).parent
    path_traj = current_dir / 'test_data' / 'POSCAR_HfO2'
    path_answer = current_dir / 'test_data' / '6.fingerprint' / 'answer_fingerprint.npy'

    if not path_traj.is_file():
        pytest.fail(f"Test input file not found: {path_traj}")
    if not path_answer.is_file():
        pytest.fail(f"Test answer file not found: {path_answer}")
    else:
        expected_fingerprint = np.load(path_answer)

    try:
        actual_fingerprint = get_fingerprint(
            str(path_traj),
            filename=None,
            atom_pairs=[('Hf', 'Hf'), ('Hf', 'O'), ('O', 'O')],
            disp=False
        )
    except Exception as e:
        pytest.fail(f"Failed to call get_fingerprint() object: {e}")

    assert_allclose(actual_fingerprint, expected_fingerprint)