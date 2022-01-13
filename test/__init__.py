from os.path import join
from profiles.conf import coef_info
import pathlib

BASE_TEST_PATH = 'test/data'

def get_data_file_path(file): return join(BASE_TEST_PATH, file)

coef_info.FILE_PATH = join(pathlib.Path().resolve(), 'coefs')
