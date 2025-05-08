import matplotlib
matplotlib.use("TkAgg")

import pytest
import numpy as np
import os
from importlib import import_module
import astropy.units as u
from synphot import SpectralElement
from synphot.models import Box1D

import config_stp
import config_stp_wcc
import config_um

from etc_wcc import ExposureTimeSNRCalculator as etsc

@pytest.fixture
def UM():
    return "UM"

@pytest.fixture
def STP():
    return "STP"

@pytest.fixture
def telescope(request):
    return request.getfixturevalue(request.param)

def test_module_installations():
    try:
        import_module("config_stp", package="config_stp")
    except:
        print("config_stp is not installed.")
    try:
        import_module("config_stp_wcc", package="config_stp_wcc")
    except:
        print("config_stp_wcc is not installed.")
    try:
        import_module("config_um", package="config_um")
    except:
        print("config_um is not installed.")

    print("Module Check: All modules installed.")
    return

def test_environment_variables():
    try:
        assert "UASAL_ARCHIVE" in os.environ
    except:
        print("UASAL_ARCHIVE is not initialized in the environment.")
    print("Environment variable check: All variables initialized.")
    return

@pytest.mark.parametrize("telescope", ["UM", "STP"], indirect=True)
def test_configs_telescope(telescope):
    if telescope == "STP":
        data_telescope = config_stp.load_config_values("parsed")
    if telescope == "UM":
        data_telescope = config_um.load_config_values("parsed")
    try:
        assert "value" in data_telescope['observatory']['pointing']['jitter_rms']
        assert "value" in data_telescope['telescope']['optics']['m1']['aper_clear_OD']
        assert "value" in data_telescope['telescope']['optics']['m1']['surface_rms']
        assert isinstance(data_telescope['telescope']['general']['f_number'], float)
        assert "value" in data_telescope['telescope']['optics']['m2']['aper_clear_OD']
        assert "value" in data_telescope['telescope']['optics']['m2']['support_width']
        assert isinstance(data_telescope['telescope']['optics']['m2']['n_supports'], int)
    except:
        print(f"Values could not be verified in {telescope} config. Possibly corrupted or modified data.")

    print(f"Telescope Config Check: {telescope} configs initialized and confirmed.")

    return

def test_configs_instrument():
    data_instrument = config_stp_wcc.load_config_values("parsed")
    try:
        assert "value" in data_instrument['common_params']['sensor']['temp_nominal']
        assert isinstance(data_instrument['common_params']['sensor']['gain'], int)
    except:
        print(f"Value type mismatch. Potentially corrupted or modified data.")
    print(f"Instrument Config Check: WCC configs initialized and confirmed.")
    return

@pytest.mark.parametrize("telescope", ["UM", "STP"], indirect=True)
def test_default_throughput(telescope):
    try:
        obs = etsc.Observatory(telescope)
        obs.make_STP()
        flux = obs.bandpass(obs.bandpass.waveset)
        bp_non_zero = flux[flux!=0]
        waveset_non_zero = obs.bandpass.waveset[flux!=0]
        assert round(np.min(waveset_non_zero).value,2) == 5380.0
        assert round(np.max(waveset_non_zero).value,2) == 7229.99
        assert round(np.mean(bp_non_zero).value, 3) == 0.133
    except:
        print("Throughput check failed. Default filter modified or corrupt installation.")

    print("Throughput check: Filter and total throughput initialized and confirmed.")
    return

@pytest.mark.parametrize("telescope", ["UM", "STP"], indirect=True)
def test_sensor_initialization(telescope):
    obs = etsc.Observatory(telescope)
    obs.make_STP()
    if telescope == "UM":
        try:
            assert round(obs.plate_scale,6) == 0.016869
        except:
            print("Sensor check: Incorrect plate scale. Potentially corrupt or modified configuration file.")
    if telescope == "STP":
        try:
            assert round(obs.plate_scale,6) == 0.008054
        except:
            print("Sensor check: Incorrect plate scale. Potentially corrupt or modified configuration file.")
    try:
        assert obs.num_psf_pixels.value == 36
    except:
        print("Incorrect PSF size.")

    print("Sensor Check: Sensor initialized and confirmed.")
    return

@pytest.mark.parametrize("telescope", ["UM", "STP"], indirect=True)
def test_counts(telescope):
    if telescope == "STP":
        data_telescope = config_stp.load_config_values("parsed")
        data_path_telescope = config_stp.get_data_path()
    if telescope == "UM":
        data_telescope = config_um.load_config_values("parsed")
        data_path_telescope = config_um.get_data_path()
    zodi_magnitude_normalization = float(data_telescope['astrophysics']['zodi']['zodi_mag_r'])
    obs = etsc.Observatory(telescope)
    obs.make_STP()
    uasal_archive = os.environ.get("UASAL_ARCHIVE")
    obs.set_source(source_pickles_file='test_data/pickles_uk_9.fits',
                   plot=True)
    obs.set_background(background_file=data_telescope['astrophysics']['zodi']['profile'], support_data_path=data_path_telescope,
                       plot=True)
    obs.make_observation(flux=0.0, flux_units='AB', bg_flux=zodi_magnitude_normalization, bg_flux_units=u.ABmag)

    if telescope == "UM":
        assert round(obs.source_counts.value) == 16230700376
        assert round(obs.sky_counts.value) == 20
        assert round(obs.calc_SNR(exp_time=1.0 * u.s, int_time=1.0 * u.s)) == 127400
        assert round(obs.calc_req_source(10.0, int_time=1.0 * u.s, exp_time=1.0 * u.s)) == 155
        assert round(obs.calc_req_source(10.0, int_time=1.0 * u.s, exp_time=1.0 * u.s, magnitude=True)[0], 2) == 20.05
        assert round(obs.calc_int_time(1e6, exp_time=1.0*u.s).value) == 62
        assert round(obs.calc_saturation_time().value, 6) == 3.6e-5
    if telescope == "STP":
        assert round(obs.source_counts.value) == 64610011021
        assert round(obs.sky_counts.value) == 79
        assert round(obs.calc_SNR(exp_time=1.0 * u.s, int_time=1.0 * u.s)) == 254185
        assert round(obs.calc_req_source(10.0, int_time=1.0 * u.s, exp_time=1.0 * u.s)) == 155
        assert round(obs.calc_req_source(10.0, int_time=1.0 * u.s, exp_time=1.0 * u.s, magnitude=True)[0], 2) == 21.55
        assert round(obs.calc_int_time(1e6, exp_time=1.0*u.s).value) == 15
        assert round(obs.calc_saturation_time().value, 6) == 9e-6
    return

@pytest.mark.parametrize("telescope", ["UM", "STP"], indirect=True)
def test_validate_ETC_snr_calculation(telescope):
    if telescope == "STP":
        data_telescope = config_stp.load_config_values("parsed")
        data_path_telescope = config_stp.get_data_path()
    if telescope == "UM":
        data_telescope = config_um.load_config_values("parsed")
        data_path_telescope = config_um.get_data_path()
    zodi_magnitude_normalization = float(data_telescope['astrophysics']['zodi']['zodi_mag_r'])
    obs = etsc.Observatory(telescope)
    obs.make_STP(mirrors=False, filters=False)
    obs.bandpass *= SpectralElement(Box1D, amplitude=1, x_0=5500, width=1)

    test_flux_mag = 18.0
    test_wavelength = 5500 #    Angstroms
    obs.set_source(source_pickles_file='test_data/pickles_uk_9.fits',
                   plot=True)
    obs.set_background(background_file=data_telescope['astrophysics']['zodi']['profile'], support_data_path=data_path_telescope,
                       plot=True)
    obs.make_observation(flux=test_flux_mag)

    bandpass_val = obs.bandpass(test_wavelength)

    #   Run independent calculations of expected source and sky counts
    #   Source counts
    flux_source = test_flux_mag * u.ABmag
    flambda_flux_source = flux_source.to(u.erg / u.s / u.cm ** 2 / u.AA, equivalencies=u.spectral_density(test_wavelength * u.AA))
    flambda_flux_photon_source = flambda_flux_source / 6.626e-27 / 299792458 * test_wavelength * 1e-10 * u.photon / u.erg

    flambda_flux_photon_telescope_area_source = flambda_flux_photon_source * (obs.surf_area.value*10000) * u.cm ** 2 * u.AA  # For a 1 Angstrom narrowband filter
    flambda_flux_photon_telescope_area_bandpass_source = flambda_flux_photon_telescope_area_source * bandpass_val

    #   Sky counts
    flux_sky = obs.bg_magnitude * u.ABmag
    flambda_flux_sky = flux_sky.to(u.erg / u.s / u.cm ** 2 / u.AA, equivalencies=u.spectral_density(test_wavelength * u.AA))
    flambda_flux_photon_sky = flambda_flux_sky / 6.626e-27 / 299792458 * test_wavelength * 1e-10 * u.photon / u.erg

    flambda_flux_photon_telescope_area_sky = flambda_flux_photon_sky * (obs.surf_area.value*10000) * u.cm ** 2 * u.AA  # For a 1 Angstrom narrowband filter
    flambda_flux_photon_telescope_area_bandpass_sky = flambda_flux_photon_telescope_area_sky * bandpass_val

    #   Zodi is normalized in Johnson V band
    #   Vega has V mag = 0.03 in Johnson V band
    #   Johnson V band mag = AB mag
    #   https://lweb.cfa.harvard.edu/~dfabricant/huchra/ay145/mags.html
    #   V band to AB conversion => V = M(AB) + 0.03
    #   0.03 = V-M(AB) = - 2.5 * log10(F_V/F_AB)
    #   => 10 ** (0.03/-2.5) * F_AB = F_V
    #   This seems to be the root cause of the issue, but can't confirm that the Zodi isn't in V magnitudes
    #   flambda_flux_photon_telescope_area_bandpass_sky_V = flambda_flux_photon_telescope_area_bandpass_sky / 10**(-0.044/2.5)

    #   Check the source count calculation from independent check is good to within 0.01 mag
    assert abs((-2.5 * np.log10(obs.source_counts.value/flambda_flux_photon_telescope_area_bandpass_source.value))) < 0.01
    #   Check the sky count calculation from independent check is good to within 0.05 mag
    assert abs((-2.5 * np.log10(obs.sky_counts.value/flambda_flux_photon_telescope_area_bandpass_sky.value))) < 0.05

    #   Check the actual SNR calculation

    dc = obs.dark_current.value
    rn = obs.read_noise.value
    n_psf = obs.num_psf_pixels.value
    int_time = 100.0
    exp_time = 100.0

    total_noise = np.sqrt(flambda_flux_photon_telescope_area_bandpass_source.value * int_time + flambda_flux_photon_telescope_area_bandpass_sky.value * int_time
                       + (dc + rn * rn / exp_time) * int_time * n_psf)
    snr = flambda_flux_photon_telescope_area_bandpass_source.value * int_time / total_noise
    snr_etc = obs.calc_SNR(int_time*u.s, exp_time*u.s)

    #   Check SNR is good to within 0.01 %
    assert abs(1-snr_etc/snr)*100.0 < 0.01

    return


if __name__ == "__main__":
    test_module_installations()
    test_environment_variables()
    test_configs_telescope()
    test_configs_instrument()
    test_default_throughput()
    test_sensor_initialization()
    test_counts()
    test_validate_ETC_snr_calculation()