# Profiles

This package was built to handle processing data primarily from the CopterSonde, a thermodynamic and kinematic profiling UAS.

See https://oucass.github.io/Profiles for detailed documentation of the API.

## Installation

The package requires two git submodules to work properly. 

The `dronelogbook` submodule acts to interface with DroneLogBook API to download CopterSonde data and collect metadata from each flight. You will also use this to sync the registry of CopterSondes scoops that are in use to a local directory `~/.wxuas`

The `SensorCoefficients` submodule contains a CSV that is continually updated with information about calibrations and corrections for the sensors contained in each scoop.

To install everything, use the following commands

##### 1. Clone the Profiles package to the desired location
```
git clone git@github.com:oucass/Profiles.git
```

##### 2. Initiate the submodules
This will pull down the submodules that are currently compatible with the Profiles package
```
cd Profiles
git submodules init
git submodules update
```

##### 3. Install the dronelogbook submodule and sync your system with DroneLogBook to test the API
The install will also prompt you to input your DroneLogBook API key and configure the `~/.wxuas` directory
```
cd dronelogbook
python setup.py install
cd scripts
python sync_dlb.py 
cd ../../
```

##### 4. Install the profiles package
```
python setup.py install
```

##### 5. Softlink the MasterCoefList to `~/.wxuas/`
```
ln -s SensorCoefficients/MasterCoefList.csv ~/.wxuas/MasterCoefList.csv
```


### References:

https://www.atmos-meas-tech.net/13/2833/2020/ \
https://www.atmos-meas-tech.net/13/3855/2020/ \
https://www.mdpi.com/1424-8220/19/12/2720


