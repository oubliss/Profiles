import setuptools
import profiles

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="oucass-profiles",
    version=profiles.__version__,
    author="Jessica Blunt, Tyler Bell, Brian Greene, Gus Azevedo, Ariel Jacobs",
    author_email="cass@ou.edu",
    description="Tools to process atmospheric data collected by UAS along either vertical or horizontal lines",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/oucass/Profiles",
    packages=setuptools.find_packages(),
    install_requires = [
        'metpy>=1.1.0',
        'pymavlink',
        'netCDF4',
        'matplotlib',
        'pandas',
        'cmocean',
        'numpy>=1.16.5'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    include_package_data=True,
)
