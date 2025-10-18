import os
import json
import pytest
import numpy as np
from pathlib import Path
from numpy.testing import assert_allclose
from vachoppy.vibration import Vibration


@pytest.fixture(scope='module')
def vibration_data(site_data):
    """
    Fixture to create a Vibration object, run calculations,
    and provide the answer file path.
    """
    current_dir = Path(__file__).parent
    path_traj = current_dir / 'test_data' / 'traj' / 'TRAJ_O_2000K_CUT.h5'
    path_answer = current_dir / 'test_data' / '2.vibration' / 'answer_vibration.json'

    if not path_traj.is_file():
        pytest.fail(f"Test input file not found: {path_traj}")
    if not path_answer.is_file():
        pytest.fail(f"Test answer file not found: {path_answer}")

    site_object, _ = site_data

    try:
        vib_object = Vibration(str(path_traj), site_object)
        vib_object.calculate()
    except Exception as e:
        pytest.fail(f"Failed during Vibration object initialization or calculation: {e}")
        
    return vib_object, path_answer


def test_vibration_displacement(vibration_data):
    """Tests the displacement properties (mu, sigma) of the Vibration object."""
    vib_object, path_answer = vibration_data

    with open(path_answer, 'r') as f: answer = json.load(f)

    expected_mu = answer['mu_displacements']
    expected_sigma = answer['sigma_displacements']

    actual_mu = vib_object.mu_displacements
    actual_sigma = vib_object.sigma_displacements

    assert actual_mu == pytest.approx(expected_mu)
    assert actual_sigma == pytest.approx(expected_sigma)


def test_vibration_frequency(vibration_data):
    """Tests the mean_frequency attribute of the Vibration object."""
    vib_object, path_answer = vibration_data

    with open(path_answer, 'r') as f: answer = json.load(f)

    expected_mean = answer['mean_frequency']
    actual_mean = vib_object.mean_frequency

    assert actual_mean == pytest.approx(expected_mean)
