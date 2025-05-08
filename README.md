# README
Exposure Time Calculator (ETC) developed for Wavefront and Context Camera (WCC).

Originally written by Aaron Goldtooth

Additional Contributions from Jess Johnson, Justin Hom and Sanchit Sabhlok 

## Requirements
The package requires a number of other configuration repos to be setup to fully configure the STP or UM telescopes or instruments. These are -
1. [config_um](https://github.com/uasal/config_um) - The configuration repo for Ultramarine 3 m telescope.
2. [config_stp](https://github.com/uasal/config_stp) - The configuration repo for the Space Telescope Project.
3. [config_stp_wcc](https://github.com/uasal/config_stp_wcc) - Configuration repo for the WCC instrument. 

All the information regarding the telescope and instrument is drawn from these repos. The repos also contain supporting data files for the filters, coatings, sensors, etc...

For full functionality including stellar and galactic spectra, you will need to set up access to the [UASAL archive](https://github.com/uasal/uasal_archive) and add a path variable `$UASAL_ARCHIVE` to your environment, pointing to the location of the cloned UASAL Archive directory. 

The installation instructions for these packages can be found on their corresponding github repos. 

Additional package dependencies required for the ETC - astropy, numpy, matplotlib, scipy and synphot.

Optional dependencies - pandas

## Installation

The Exposure Time Calculator is a python package that can be installed via a download from github and then installing on your system locally. You can directly clone the [ETC github repo](https://github.com/uasal/etc_wcc) or fork the repo and clone the fork. 
```
$ git clone git@github.com:uasal/etc_wcc.git
$ cd etc_wcc 
$ pip install .
```
The ETC should now be installed on your local machine. To confirm installation, the following import on python should work - 
```
import etc_wcc
etc_wcc.__version__
```

## Example notebooks
Three demo notebooks are provided with the ETC. The `ETC_Demo.ipynb` notebook walks through the full functionality of the ETC, whereas the `Example_ETC_SNR_Calculation_SN1a.ipynb` walks through the specific science case of Type IA Supernovae embedded in host galaxies. `Throughput.ipynb` calculates the total throughput of a given setup.

## Report a Bug or Request a Feature
The ETC is currently actively maintained on github. Bug Fixes or Features can be requested by opening an issue on the github. When opening an issue, please attach a functional code segment demonstrating the behavior and a description of the desired behavior. This will help us address issues and patch the ETC in a timely manner. 