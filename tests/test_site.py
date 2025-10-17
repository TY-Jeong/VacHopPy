import numpy as np
from numpy.testing import assert_allclose
from vachoppy.core import Site


def test_site_lattice_parameter(site_data):
    """Site 객체의 격자 상수가 정확한지 테스트합니다."""
    site_object, answer_data = site_data
    expected_lattice = np.array(answer_data['lattice'])
    
    assert_allclose(site_object.lattice_parameter, expected_lattice)

def test_site_path_names(site_data):
    """Site 객체의 path_name 리스트가 정확한지 테스트합니다."""
    site_object, answer_data = site_data
    
    assert site_object.path_name == answer_data['path_name']

def test_site_site_names(site_data):
    """Site 객체의 site_name 리스트가 정확한지 테스트합니다."""
    site_object, answer_data = site_data
    
    assert site_object.site_name == answer_data['site_name']