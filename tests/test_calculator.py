import os
import json
import pytest
import numpy as np
from pathlib import Path
from numpy.testing import assert_allclose
# Assuming Calculator class is imported correctly
from vachoppy.core import Calculator


@pytest.fixture(scope='module')
def calculator_data(site_data, tmp_path_factory):
    """Fixture to create Calculator object, run analyses, and return results."""

    current_dir = Path(__file__).parent
    path_traj = current_dir/'test_data'/'traj'
    path_neb = current_dir/'test_data'/'neb_HfO2.csv'
    path_answer_parameters = current_dir/'test_data'/'4.calculator'/'answer_parameters.json'
    path_answer_decompose = current_dir/'test_data'/'4.calculator'/'answer_decompose.json'

    if not path_traj.is_dir():
        pytest.fail(f"Test input directory not found: {path_traj}")
    if not path_neb.is_file():
        pytest.fail(f"Test input file not found: {path_neb}")
    if not path_answer_parameters.is_file():
        pytest.fail(f"Test answer file not found: {path_answer_parameters}")
    else:
        with open(path_answer_parameters, 'r') as f: answer_parameters = json.load(f)
    if not path_answer_decompose.is_file():
        pytest.fail(f"Test answer file not found: {path_answer_decompose}")
    else:
        with open(path_answer_decompose, 'r') as f: answer_decompose = json.load(f)

    site_object, _ = site_data
    temp_dir = tmp_path_factory.mktemp("calculator_run")
    temp_param_file = temp_dir / "parameters.json"

    try:
        calc_object = Calculator(str(path_traj), site_object, verbose=False)
    except Exception as e:
        pytest.fail(f"Failed to initialize Calculator object: {e}")

    try:
        calc_object.calculate(verbose=False)
    except Exception as e:
        pytest.fail(f"Failed to call Calculator.calculate() method: {e}")

    try:
        calc_object.calculate_attempt_frequency(
            neb_csv=str(path_neb),
            filename=str(temp_param_file)
        )
    except Exception as e:
        pytest.fail(f"Failed to call Calculator.calculate_attempt_frequency() method: {e}")

    try:
        calc_object.decompose_diffusivity(verbose=False)
    except Exception as e:
        pytest.fail(f"Failed to call Calculator.decompose_diffusivity() method: {e}")

    return calc_object, answer_parameters, answer_decompose


def test_calculator_temperature(calculator_data):
    """Tests the temperatures attribute."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['temperatures'])
    actual_value = calc_object.temperatures
    assert_allclose(actual_value, expected_value)


def test_calculator_num_vacancies(calculator_data):
    """Tests the num_vacancies attribute."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = answer_parameters['num_vacancies']
    actual_value = calc_object.num_vacancies
    assert actual_value == expected_value


def test_calculator_D(calculator_data):
    """Tests the D attribute (diffusivity)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['D'])
    actual_value = calc_object.D
    assert_allclose(actual_value, expected_value)


def test_calculator_D0(calculator_data):
    """Tests the D0 attribute (pre-exponential factor for D)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = answer_parameters['D0']
    actual_value = calc_object.D0
    assert actual_value == pytest.approx(expected_value)


def test_calculator_Ea_D(calculator_data):
    """Tests the Ea_D attribute (activation energy for D)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = answer_parameters['Ea_D']
    actual_value = calc_object.Ea_D
    assert actual_value == pytest.approx(expected_value)


def test_calculator_D_rand(calculator_data):
    """Tests the D_rand attribute (random walk diffusivity)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['D_rand'])
    actual_value = calc_object.D_rand
    assert_allclose(actual_value, expected_value)


def test_calculator_D_rand0(calculator_data):
    """Tests the D_rand0 attribute (pre-exponential factor for D_rand)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = answer_parameters['D_rand0']
    actual_value = calc_object.D_rand0
    assert actual_value == pytest.approx(expected_value)


def test_calculator_Ea_D_rand(calculator_data):
    """Tests the Ea_D_rand attribute (activation energy for D_rand)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = answer_parameters['Ea_D_rand']
    actual_value = calc_object.Ea_D_rand
    assert actual_value == pytest.approx(expected_value)


def test_calculator_f(calculator_data):
    """Tests the f attribute (correlation factor)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['f'])
    actual_value = calc_object.f
    assert_allclose(actual_value, expected_value)


def test_calculator_f0(calculator_data):
    """Tests the f0 attribute (pre-exponential factor for f)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = answer_parameters['f0']
    actual_value = calc_object.f0
    assert actual_value == pytest.approx(expected_value)


def test_calculator_Ea_f(calculator_data):
    """Tests the Ea_f attribute (activation energy for f)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = answer_parameters['Ea_f']
    actual_value = calc_object.Ea_f
    assert actual_value == pytest.approx(expected_value)


def test_calculator_tau(calculator_data):
    """Tests the tau attribute (residence time)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['tau'])
    actual_value = calc_object.tau
    assert_allclose(actual_value, expected_value)


def test_calculator_tau0(calculator_data):
    """Tests the tau0 attribute (pre-exponential factor for tau)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = answer_parameters['tau0']
    actual_value = calc_object.tau0
    assert actual_value == pytest.approx(expected_value)


def test_calculator_Ea_tau(calculator_data):
    """Tests the Ea_tau attribute (activation energy for tau)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = answer_parameters['Ea_tau']
    actual_value = calc_object.Ea_tau
    assert actual_value == pytest.approx(expected_value)


def test_calculator_a(calculator_data):
    """Tests the a attribute (effective hopping distance)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['a'])
    actual_value = calc_object.a
    assert_allclose(actual_value, expected_value)


def test_calculator_a_path(calculator_data):
    """Tests the a_path attribute (path-wise hopping distance)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['a_path'])
    actual_value = calc_object.a_path
    assert_allclose(actual_value, expected_value)


def test_calculator_nu(calculator_data):
    """Tests the nu attribute (effective attempt frequency)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['nu'])
    actual_value = calc_object.nu
    assert_allclose(actual_value, expected_value)


def test_calculator_nu_path(calculator_data):
    """Tests the nu_path attribute (path-wise attempt frequency)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['nu_path'])
    actual_value = calc_object.nu_path
    assert_allclose(actual_value, expected_value)


def test_calculator_z(calculator_data):
    """Tests the z attribute (effective coordination number)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['z'])
    actual_value = calc_object.z
    assert_allclose(actual_value, expected_value)


def test_calculator_z_path(calculator_data):
    """Tests the z_path attribute (path-wise coordination number)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['z_path'])
    actual_value = calc_object.z_path
    assert_allclose(actual_value, expected_value)


def test_calculator_z_mean(calculator_data):
    """Tests the z_mean attribute (mean coordination number)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['z_mean'])
    actual_value = calc_object.z_mean
    assert_allclose(actual_value, expected_value)


def test_calculator_P_site(calculator_data):
    """Tests the P_site attribute (site occupation probability)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['P_site'])
    actual_value = calc_object.P_site
    assert_allclose(actual_value, expected_value)


def test_calculator_times_site(calculator_data):
    """Tests the times_site attribute (total time spent at each site type)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['times_site'])
    actual_value = calc_object.times_site
    assert_allclose(actual_value, expected_value)


def test_calculator_counts_hop(calculator_data):
    """Tests the counts_hop attribute (hop counts per path per temperature)."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['counts_hop'])
    actual_value = calc_object.counts_hop
    assert_allclose(actual_value, expected_value)


def test_attempt_frequency_Ea_path(calculator_data):
    """Tests the Ea_path attribute from the attempt_frequency object."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['Ea_path'])
    actual_value = calc_object.attempt_frequency.Ea_path
    assert_allclose(actual_value, expected_value)


def test_attempt_frequency_m_mean(calculator_data):
    """Tests the m_mean attribute from the attempt_frequency object."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['m_mean'])
    actual_value = calc_object.attempt_frequency.m_mean
    assert_allclose(actual_value, expected_value)


def test_attempt_frequency_P_esc(calculator_data):
    """Tests the P_esc attribute from the attempt_frequency object."""
    calc_object, answer_parameters, _ = calculator_data
    expected_value = np.array(answer_parameters['P_esc'])
    actual_value = calc_object.attempt_frequency.P_esc
    assert_allclose(actual_value, expected_value)


def test_calculator_msd(calculator_data):
    """Tests the msd attribute (list of MSD arrays per temperature)."""
    calc_object, _, answer_decompose = calculator_data
    expected_value = answer_decompose['msd'] # List of lists of lists
    actual_value = calc_object.msd          # List of numpy arrays

    assert len(actual_value) == len(expected_value), "Number of temperature segments mismatch in MSD."

    for i, (actual_arr, expected_list) in enumerate(zip(actual_value, expected_value)):
        expected_arr = np.array(expected_list)
        assert_allclose(
            actual_arr,
            expected_arr,
            err_msg=f"MSD data mismatch in temperature segment {i}."
        )


def test_calculator_Dx(calculator_data):
    """Tests the Dx attribute (decomposed diffusivity in x)."""
    calc_object, _, answer_decompose = calculator_data
    expected_value = np.array(answer_decompose['Dx'])
    actual_value = calc_object.Dx
    assert_allclose(actual_value, expected_value)


def test_calculator_Dy(calculator_data):
    """Tests the Dy attribute (decomposed diffusivity in y)."""
    calc_object, _, answer_decompose = calculator_data
    expected_value = np.array(answer_decompose['Dy'])
    actual_value = calc_object.Dy
    assert_allclose(actual_value, expected_value)


def test_calculator_Dz(calculator_data):
    """Tests the Dz attribute (decomposed diffusivity in z)."""
    calc_object, _, answer_decompose = calculator_data
    expected_value = np.array(answer_decompose['Dz'])
    actual_value = calc_object.Dz
    assert_allclose(actual_value, expected_value)


def test_calculator_Ea_x(calculator_data):
    """Tests the Ea_x attribute (activation energy for Dx)."""
    calc_object, _, answer_decompose = calculator_data
    expected_value = answer_decompose['Ea_x']
    actual_value = float(calc_object.Ea_x)
    assert actual_value == pytest.approx(expected_value)


def test_calculator_Ea_y(calculator_data):
    """Tests the Ea_y attribute (activation energy for Dy)."""
    calc_object, _, answer_decompose = calculator_data
    expected_value = answer_decompose['Ea_y']
    actual_value = float(calc_object.Ea_y)
    assert actual_value == pytest.approx(expected_value)


def test_calculator_Ea_z(calculator_data):
    """Tests the Ea_z attribute (activation energy for Dz)."""
    calc_object, _, answer_decompose = calculator_data
    expected_value = answer_decompose['Ea_z']
    actual_value = float(calc_object.Ea_z)
    assert actual_value == pytest.approx(expected_value)


def test_calculator_D0_x(calculator_data):
    """Tests the D0_x attribute (pre-exponential factor for Dx)."""
    calc_object, _, answer_decompose = calculator_data
    expected_value = answer_decompose['D0_x']
    actual_value = float(calc_object.D0_x)
    assert actual_value == pytest.approx(expected_value)


def test_calculator_D0_y(calculator_data):
    """Tests the D0_y attribute (pre-exponential factor for Dy)."""
    calc_object, _, answer_decompose = calculator_data
    expected_value = answer_decompose['D0_y']
    actual_value = float(calc_object.D0_y)
    assert actual_value == pytest.approx(expected_value)


def test_calculator_D0_z(calculator_data):
    """Tests the D0_z attribute (pre-exponential factor for Dz)."""
    calc_object, _, answer_decompose = calculator_data
    expected_value = answer_decompose['D0_z']
    actual_value = float(calc_object.D0_z)
    assert actual_value == pytest.approx(expected_value)
