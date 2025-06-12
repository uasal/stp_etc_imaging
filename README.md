# README
Exposure Time Calculator (ETC) for estimating sensitivity of complex optical systems with commercial imaging sensors for astronomical space telescope program (STP) science.

Originally written by Aaron Goldtooth

Additional Contributions from Jess Johnson, Justin Hom and Sanchit Sabhlok 

## Requirements

The package requires  configuration repos that define a telescope or instrument system. While the ETC can function without them, this is neither recommended nor a supported mode. All the information regarding the telescope and instrument is drawn from these repos. The repos also contain supporting data files for the filters, coatings, sensors, etc...

### Telescope

1. [config_um](https://github.com/uasal/config_um) - The configuration repo for UA Ultramarine 3 m telescope concept.
2. [config_stp](https://github.com/uasal/config_stp) - The configuration repo for the Space Telescope Project 6.5m (https://arxiv.org/abs/2309.04934).

### Instrument

1. [config_stp_wcc](https://github.com/uasal/config_stp_wcc) - Configuration repo for the WCC instrument on the STP 6.5m telescope. 
2. [config_um_wcc](https://github.com/uasal/config_um_wcc) - Configuration repo for the WCC instrument on the UA Ultramarine telescope. 

### UASAL Archive
For full functionality including stellar and galactic spectra, you will need to set up access to the [UASAL archive](https://github.com/uasal/uasal_archive) and add a path variable `$UASAL_ARCHIVE` to your environment, pointing to the location of the cloned UASAL Archive directory. 

The installation instructions for these packages can be found on their corresponding GitHub repos. 

Additional package dependencies required for the ETC - astropy, numpy, matplotlib, scipy and synphot.

Optional dependencies - pandas

## Installation

The Exposure Time Calculator is a Python package that can be installed via a download from GitHub and then installing on your system locally. You can directly clone the [ETC GitHub repo](https://github.com/uasal/stp_etc_imaging) or fork the repo and clone the fork. 
```
$ git clone git@github.com:uasal/stp_etc_imaging.git
$ cd stp_etc_imaging 
$ pip install .
```
The ETC should now be installed on your local machine. To confirm installation, the following import on Python should work - 
```
import stp_etc_imaging
stp_etc_imaging.__version__
```

## Example notebooks
Three demo notebooks are provided with the ETC. The `ETC_Demo.ipynb` notebook walks through the full functionality of the ETC, whereas the `Example_ETC_SNR_Calculation_SN1a.ipynb` walks through the specific science case of Type IA Supernovae embedded in host galaxies. `Throughput.ipynb` calculates the total throughput of a given setup.

## Report a Bug or Request a Feature
The ETC is currently actively maintained on GitHub. Bug Fixes or Features can be requested by opening an issue on the GitHub. When opening an issue, please attach a functional code segment demonstrating the behavior and a description of the desired behavior. This will help us address issues and patch the ETC in a timely manner. 
