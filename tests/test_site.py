import numpy as np
from numpy.testing import assert_allclose


def test_site_lattice_parameter(site_data):
    """Tests if the Site object's lattice parameters are correct."""
    site_object, answer_data = site_data
    expected_lattice = np.array(answer_data['lattice'])
    assert_allclose(site_object.lattice_parameter, expected_lattice)

def test_site_path_names(site_data):
    """Tests if the Site object's path_name list is correct."""
    site_object, answer_data = site_data
    assert site_object.path_name == answer_data['path_name']

def test_site_site_names(site_data):
    """Tests if the Site object's site_name list is correct."""
    site_object, answer_data = site_data
    assert site_object.site_name == answer_data['site_name']

def test_site_lattice_sites(site_data):
    """Tests if the Site object's lattice site coordinates are correct."""
    site_object, answer_data = site_data
    expected_coord = np.array(answer_data['coord'])
    actual_coord = np.array([d['coord'] for d in site_object.lattice_sites])
    assert_allclose(actual_coord, expected_coord)

