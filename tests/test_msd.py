import os
import json
import pytest
import numpy as np
from pathlib import Path
from numpy.testing import assert_allclose
# Assuming Einstein class is imported correctly
from vachoppy.einstein import Einstein

@pytest.fixture(scope='module')
def msd_data():
    """
    Fixture to initialize the Einstein object, run calculations,
    and load answer data.
    """

    current_dir = Path(__file__).parent
    path_traj = current_dir / 'test_data' / 'traj'
    path_answer_einstein = current_dir / 'test_data' / '5.msd' / 'answer_einstein.json'
    path_answer_msd = current_dir / 'test_data' / '5.msd' / 'answer_msd.npy'
    
    if not path_traj.is_dir():
        pytest.fail(f"Test input directory not found: {path_traj}")
    if not path_answer_einstein.is_file():
        pytest.fail(f"Test answer file not found: {path_answer_einstein}")
    else:
        with open(path_answer_einstein, 'r') as f: answer_einstein = json.load(f)
    if not path_answer_msd.is_file():
        pytest.fail(f"Test answer file not found: {path_answer_msd}")
    else:
        answer_msd = np.load(path_answer_msd)

    symbol = 'O'

    try:
        einstein_object = Einstein(str(path_traj), symbol, verbose=False)
    except Exception as e:
        pytest.fail(f"Failed to initialize Einstein object: {e}")

    try:
        einstein_object.calculate()
    except Exception as e:
        pytest.fail(f"Failed to call Einstein.calculate() method: {e}")

    return einstein_object, answer_einstein, answer_msd


def test_msd_msd(msd_data):
    """Tests the msd attribute of the first internal calculator."""
    einstein_object, _, answer_msd = msd_data
    expected_value = answer_msd

    if einstein_object.calculators:
         actual_value = einstein_object.calculators[0].msd
         assert_allclose(actual_value, expected_value)
    else:
         pytest.skip("Skipping MSD test as no calculators were found in Einstein object.")


def test_msd_temperatures(msd_data):
    """Tests the temperatures attribute."""
    einstein_object, answer_einstein, _ = msd_data
    expected_value = np.array(answer_einstein['temperatures'])
    actual_value = einstein_object.temperatures
    assert_allclose(actual_value, expected_value)


def test_msd_D(msd_data):
    """Tests the D attribute (diffusivity array)."""
    einstein_object, answer_einstein, _ = msd_data
    expected_value = np.array(answer_einstein['D'])
    actual_value = einstein_object.D
    assert_allclose(actual_value, expected_value)


def test_msd_Ea_D(msd_data):
    """Tests the Ea_D attribute (activation energy)."""
    einstein_object, answer_einstein, _ = msd_data
    expected_value = answer_einstein['Ea_D']
    actual_value = float(einstein_object.Ea_D)
    assert actual_value == pytest.approx(expected_value)


def test_msd_D0(msd_data):
    """Tests the D0 attribute (pre-exponential factor)."""
    einstein_object, answer_einstein, _ = msd_data
    expected_value = answer_einstein['D0']
    actual_value = float(einstein_object.D0)
    assert actual_value == pytest.approx(expected_value)
