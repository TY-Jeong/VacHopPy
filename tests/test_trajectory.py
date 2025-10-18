import os
import json
import pytest
import numpy as np
import collections.abc
from pathlib import Path
from numpy.testing import assert_allclose
from vachoppy.trajectory import Trajectory, TrajectoryAnalyzer, Encounter


def convert_to_comparable(item):
    """Recursively converts NumPy types to standard Python types."""
    if isinstance(item, np.ndarray):
        return item.tolist()
    elif isinstance(item, (np.integer, np.floating)):
        return item.item()
    elif isinstance(item, np.bool_):
        return bool(item)
    elif isinstance(item, collections.abc.Mapping):
        return {k: convert_to_comparable(v) for k, v in item.items()}
    elif isinstance(item, collections.abc.Sequence) and not isinstance(item, (str, bytes)):
        return [convert_to_comparable(elem) for elem in item]
    else:
        return item


@pytest.fixture(scope='module')
def trajectory_data(site_data):
    """Fixture to create Trajectory, Analyzer, and Encounter objects and provide answer file paths."""
    current_dir = Path(__file__).parent
    path_traj = current_dir/'test_data'/'traj'/'TRAJ_O_1900K_CUT.h5'
    answer_trajectory = current_dir/'test_data'/'3.trajectory'/'answer_trajectory.npy'
    answer_analyzer = current_dir/'test_data'/'3.trajectory'/'answer_analyzer.json'
    answer_encounter = current_dir/'test_data'/'3.trajectory'/'answer_encounter.json'

    # Check for file existence
    if not path_traj.is_file():
        pytest.fail(f"Test input file not found: {path_traj}")
    if not answer_trajectory.is_file():
        pytest.fail(f"Test answer file not found: {answer_trajectory}")
    if not answer_analyzer.is_file():
        pytest.fail(f"Test answer file not found: {answer_analyzer}")
    if not answer_encounter.is_file():
        pytest.fail(f"Test answer file not found: {answer_encounter}")

    site_object, _ = site_data

    # Initialize objects
    try:
        trajectory_object = Trajectory(str(path_traj), site_object, t_interval=0.052)
    except Exception as e:
        pytest.fail(f"Failed to initialize Trajectory object: {e}")

    try:
        analyzer_object = TrajectoryAnalyzer(trajectory_object, verbose=False)
    except Exception as e:
        pytest.fail(f"Failed to initialize TrajectoryAnalyzer object: {e}")

    try:
        encounter_object = Encounter(analyzer_object, verbose=False)
    except Exception as e:
        pytest.fail(f"Failed to initialize Encounter object: {e}")

    # Return data in a dictionary
    return {
        'trajectory': (trajectory_object, answer_trajectory),
        'analyzer': (analyzer_object, answer_analyzer),
        'encounter': (encounter_object, answer_encounter)
    }

def test_trajectory(trajectory_data):
    """Tests the unwrapped_vacancy_trajectory_coord_cart attribute of the Trajectory object."""
    trajectory_object, answer_trajectory = trajectory_data['trajectory']

    expected_coord = np.load(answer_trajectory)
    sorted_keys = sorted(trajectory_object.unwrapped_vacancy_trajectory_coord_cart.keys())
    actual_coord = np.array([trajectory_object.unwrapped_vacancy_trajectory_coord_cart[k] for k in sorted_keys])

    assert_allclose(actual_coord, expected_coord, err_msg="Unwrapped coordinates do not match expected values.")


def test_analyzer_hopping_history(trajectory_data):
    """Tests the hopping_history attribute of the TrajectoryAnalyzer object."""
    analyzer_object, answer_analyzer = trajectory_data['analyzer']

    with open(answer_analyzer, 'r') as f: answer = json.load(f)
    expected_history = answer['hopping_history']

    actual_history = analyzer_object.hopping_history
    actual_comparable = convert_to_comparable(actual_history)

    assert len(actual_comparable) == len(expected_history), "History list lengths differ"

    for i, vacancy_history_actual in enumerate(actual_comparable):
        vacancy_history_expected = expected_history[i]
        assert len(vacancy_history_actual) == len(vacancy_history_expected), f"Vacancy {i} history lengths differ"

        for j, hop_actual in enumerate(vacancy_history_actual):
            hop_expected = vacancy_history_expected[j]

            assert hop_actual.keys() == hop_expected.keys(), f"Keys differ in hop {j} for vacancy {i}"

            for key in hop_expected.keys():
                actual_value = hop_actual[key]
                expected_value = hop_expected[key]

                if key in ['distance', 'coord_init', 'coord_final']:
                    assert actual_value == pytest.approx(expected_value), f"Value mismatch for key '{key}' in hop {j}, vacancy {i}"
                else:
                    assert actual_value == expected_value, f"Value mismatch for key '{key}' in hop {j}, vacancy {i}"


def test_analyzer_msd_rand(trajectory_data):
    """Tests the msd_rand attribute of the TrajectoryAnalyzer object."""
    analyzer_object, answer_analyzer = trajectory_data['analyzer']

    with open(answer_analyzer, 'r') as f: answer = json.load(f)
    expected_msd_rand = answer['msd_rand']

    actual_msd_rand = analyzer_object.msd_rand

    assert actual_msd_rand == pytest.approx(expected_msd_rand)


def test_analyzer_residence_time(trajectory_data):
    """Tests the residence_time attribute of the TrajectoryAnalyzer object."""
    analyzer_object, answer_analyzer = trajectory_data['analyzer']

    with open(answer_analyzer, 'r') as f: answer = json.load(f)
    expected_residence_time = np.array(answer['residence_time'])
    actual_residence_time = analyzer_object.residence_time

    assert_allclose(actual_residence_time, expected_residence_time)


def test_analyzer_counts(trajectory_data):
    """Tests the counts attribute of the TrajectoryAnalyzer object."""
    analyzer_object, answer_analyzer = trajectory_data['analyzer']

    with open(answer_analyzer, 'r') as f: answer = json.load(f)
    expected_counts = np.array(answer['counts'])
    actual_counts = analyzer_object.counts

    assert_allclose(actual_counts, expected_counts)


def test_encounter_encounter_all(trajectory_data):
    """Tests the encounter_all attribute of the Encounter object."""
    encounter_object, answer_encounter_path = trajectory_data['encounter']

    with open(answer_encounter_path, 'r') as f: answer = json.load(f)
    expected_encounter_all = answer['encounter_all']

    actual_encounter_all = encounter_object.encounter_all
    actual_comparable = convert_to_comparable(actual_encounter_all)

    assert len(actual_comparable) == len(expected_encounter_all), "Encounter list lengths differ"
    for i, encounter_actual in enumerate(actual_comparable):
        encounter_expected = expected_encounter_all[i]

        assert encounter_actual.keys() == encounter_expected.keys(), f"Keys differ in encounter {i}"

        for key in encounter_expected.keys():
            actual_value = encounter_actual[key]
            expected_value = encounter_expected[key]

            if key in ['coord_init', 'coord_final']:
                assert actual_value == pytest.approx(expected_value), f"Value mismatch for key '{key}' in encounter {i}"
            elif key == 'hopping_history':
                 assert actual_value == expected_value, f"Value mismatch for key '{key}' in encounter {i}"
            elif isinstance(expected_value, (int, float)):
                 assert actual_value == pytest.approx(expected_value), f"Value mismatch for key '{key}' in encounter {i}"
            else:
                assert actual_value == expected_value, f"Value mismatch for key '{key}' in encounter {i}"


def test_encounter_msd(trajectory_data):
    """Tests the msd attribute of the Encounter object."""
    encounter_object, answer_encounter_path = trajectory_data['encounter']

    with open(answer_encounter_path, 'r') as f: answer = json.load(f)
    expected_msd = answer['msd']

    actual_msd = float(encounter_object.msd)

    assert actual_msd == pytest.approx(expected_msd)


def test_encounter_path_count(trajectory_data):
    """Tests the path_count attribute of the Encounter object."""
    encounter_object, answer_encounter_path = trajectory_data['encounter']

    with open(answer_encounter_path, 'r') as f: answer = json.load(f)
    expected_path_count = np.array(answer['path_count'])

    actual_path_count = encounter_object.path_count

    assert_allclose(actual_path_count, expected_path_count)


def test_encounter_f_cor(trajectory_data):
    """Tests the f_cor attribute of the Encounter object."""
    encounter_object, answer_encounter_path = trajectory_data['encounter']

    with open(answer_encounter_path, 'r') as f: answer = json.load(f)
    expected_f_cor = answer['f_cor']

    actual_f_cor = float(encounter_object.f_cor)

    assert actual_f_cor == pytest.approx(expected_f_cor)
