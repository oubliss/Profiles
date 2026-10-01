import pathlib
from os.path import join

from profiles.conf import coef_info
import profiles.utils as utils

TEST_DIR = pathlib.Path(__file__).resolve().parent
BASE_TEST_PATH = TEST_DIR / 'data'
COEF_PATH = BASE_TEST_PATH / 'coefs'
BASELINE_PATH = BASE_TEST_PATH / 'baseline'


def get_data_file_path(file):
    return join(str(BASE_TEST_PATH), file)


# Point the package at the frozen test coefficient fixture rather than
# ~/.wxuas or the repo's top-level coefs/ (which cannot resolve this flight -
# see test/data/coefs/README.md). reset_coef_manager() matters because the
# manager caches, and conf may already have been read by an earlier import.
coef_info.USE_AZURE = "NO"
coef_info.FILE_PATH = str(COEF_PATH)
utils.reset_coef_manager()
