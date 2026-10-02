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
        'metpy>=1.6',
        'pymavlink>=2.4.40',
        'netCDF4>=1.6',
        'matplotlib>=3.8',
        'pandas>=2.0',
        'scipy>=1.11',
        'xarray>=2023.1',
        'requests',
        'cmocean',
        'numpy>=2.0'
    ],
    extras_require={
        'test': ['pytest>=7.0'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    include_package_data=True,
)
