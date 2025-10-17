import sys
import json
import pytest
import zipfile
import subprocess
from pathlib import Path
from vachoppy.core import Site


FILE_ID = "1pG8QNTUanKXMKyfQq51Wqr4EoKemEvJ2"
OUTPUT_ZIP_PATH = Path(__file__).parent / "test_data.zip"
TARGET_DATA_DIR = Path(__file__).parent / "test_data"


@pytest.fixture(scope="session", autouse=True)
def setup_test_data():
    if TARGET_DATA_DIR.exists():
        print(f"\n'{TARGET_DATA_DIR}' directory found. Skipping data setup.")
        return

    print(f"\n'{TARGET_DATA_DIR}' directory not found. Starting test data setup...")

    try:
        import gdown
    except (ImportError, ModuleNotFoundError):
        print("'gdown' package not found. Installing now...")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "gdown"])
            import gdown
            print("'gdown' installed successfully.")
        except Exception as e:
            pytest.fail(f"ERROR: Failed to install 'gdown'. Please install it manually. Reason: {e}")
            
    try:
        print(f"Downloading test data to '{OUTPUT_ZIP_PATH}'...")
        gdown.download(id=FILE_ID, output=str(OUTPUT_ZIP_PATH), quiet=False)
        print("Download complete.")
    except Exception as e:
        pytest.fail(f"ERROR: Failed to download the test data file. Reason: {e}")

    try:
        print(f"Unzipping '{OUTPUT_ZIP_PATH}'...")
        with zipfile.ZipFile(OUTPUT_ZIP_PATH, 'r') as zip_ref:
            zip_ref.extractall(OUTPUT_ZIP_PATH.parent)
        print("Unzip complete.")
    except zipfile.BadZipFile:
        pytest.fail("ERROR: Failed to unzip. The downloaded file might be corrupted.")
    finally:
        if OUTPUT_ZIP_PATH.exists():
            print("Cleaning up downloaded zip file...")
            OUTPUT_ZIP_PATH.unlink()

    if not TARGET_DATA_DIR.exists():
        pytest.fail(f"ERROR: Expected '{TARGET_DATA_DIR}' directory not found after unzipping.")
    
    print(f"Test data is ready in '{TARGET_DATA_DIR}'!")
    
    
@pytest.fixture(scope="session")
def site_data():
    """
    테스트에 필요한 Site 객체와 정답 데이터를 로드하는 공용 픽스처.
    """
    current_dir = Path(__file__).parent
    path_structure = current_dir / 'test_data' / 'POSCAR_HfO2'
    path_answer = current_dir / 'test_data' / '0.site' / 'answer_site.json'

    if not path_structure.is_file():
        pytest.fail(f"Test input file not found: {path_structure}")
    if not path_answer.is_file():
        pytest.fail(f"Test answer file not found: {path_answer}")

    site_object = Site(path_structure, 'O')
    with open(path_answer, 'r') as f:
        answer_data = json.load(f)
    
    return site_object, answer_data