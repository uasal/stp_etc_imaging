import matplotlib
matplotlib.use("TkAgg")

import pytest
import numpy as np
import astropy.units as u
from synphot import SpectralElement
from synphot.models import Box1D

try:
    import config_stp
    import config_stp_wcc
    import config_um
except:
    print("Could not import configuration repos.")

try:
    import utils_config
except:
    print("Could not import utils_config. Cannot parse tomls. Exiting.")
    exit()

from stp_etc_imaging import ExposureTimeSNRCalculator as etsc

@pytest.fixture
def UM():
    return "UM"

@pytest.fixture
def STP():
    return "STP"

@pytest.fixture
def test_STP():
    return "test_STP"

@pytest.fixture
def test_UM():
    return "test_UM"

@pytest.fixture
def telescope(request):
    return request.getfixturevalue(request.param)


@pytest.mark.parametrize("telescope", ["UM", "STP", "test_STP"], indirect=True)
def test_configs_telescope(telescope):
    if telescope == "STP":
        data_telescope = config_stp.load_config_values("parsed")
    if telescope == "UM":
        data_telescope = config_um.load_config_values("parsed")
    if telescope == "test_STP":
        test_loader = utils_config.ConfigLoader('test_data/test_STP', "parsed")
        data_telescope = test_loader.load_configs()
    if telescope == "test_UM":
        test_loader = utils_config.ConfigLoader('test_data/test_UM', "parsed")
        data_telescope = test_loader.load_configs()

    assert "value" in data_telescope['observatory']['pointing']['jitter_rms'], f"data_telescope['observatory']['pointing']['jitter_rms'] could not be verified in {telescope} config"
    assert "value" in data_telescope['telescope']['optics']['m1']['aper_clear_OD'], f"data_telescope['telescope']['optics']['m1']['aper_clear_OD'] could not be verified in {telescope} config"
    assert "value" in data_telescope['telescope']['optics']['m1']['surface_rms'], f"data_telescope['telescope']['optics']['m1']['surface_rms'] could not be verified in {telescope} config"
    assert isinstance(data_telescope['telescope']['general']['f_number'], float), f"data_telescope['telescope']['general']['f_number'] could not be verified in {telescope} config"
    assert "value" in data_telescope['telescope']['optics']['m2']['aper_clear_OD'], f"data_telescope['telescope']['optics']['m2']['aper_clear_OD'] could not be verified in {telescope} config"
    assert "value" in data_telescope['telescope']['optics']['m2']['support_width'], f"data_telescope['telescope']['optics']['m2']['support_width'] could not be verified in {telescope} config"
    assert isinstance(data_telescope['telescope']['optics']['m2']['n_supports'], int), f"data_telescope['telescope']['optics']['m2']['n_supports'] could not be verified in {telescope} config"

    return

def test_configs_instrument():
    data_instrument = config_stp_wcc.load_config_values("parsed")

    assert "value" in data_instrument['common_params']['sensor']['temp_nominal'], "Value type mismatch. Potentially corrupted or modified data."
    assert isinstance(data_instrument['common_params']['sensor']['gain'], int), "Value type mismatch. Potentially corrupted or modified data."

    return

@pytest.mark.parametrize("telescope", ["UM", "STP", "test_STP", "test_UM"], indirect=True)
def test_default_throughput(telescope):
    obs = etsc.Observatory(telescope)
    obs.make_STP()
    flux = obs.bandpass(obs.bandpass.waveset)
    bp_non_zero = flux[flux!=0]
    waveset_non_zero = obs.bandpass.waveset[flux!=0]
    assert round(np.min(waveset_non_zero).value,2) == 5380.0, "Throughput check failed. Default filter modified or corrupt installation."
    assert round(np.max(waveset_non_zero).value,2) == 7229.99, "Throughput check failed. Default filter modified or corrupt installation."
    assert round(np.mean(bp_non_zero).value, 3) == 0.133, "Throughput check failed. Default filter modified or corrupt installation."

    return

@pytest.mark.parametrize("telescope", ["UM", "STP"], indirect=True)
def test_sensor_initialization(telescope):
    obs = etsc.Observatory(telescope)
    obs.make_STP()
    if telescope == "UM":
        assert round(obs.plate_scale,6) == 0.016869, "Sensor check: Incorrect plate scale. Potentially corrupt or modified configuration file."
    if telescope == "STP":
        assert round(obs.plate_scale,6) == 0.008054, "Sensor check: Incorrect plate scale. Potentially corrupt or modified configuration file."

    assert obs.num_psf_pixels.value == 36, "Incorrect PSF size."

    return

@pytest.mark.parametrize("telescope", ["UM", "test_UM", "STP", "test_STP"], indirect=True)
def test_counts(telescope):
    if telescope == "STP":
        data_telescope = config_stp.load_config_values("parsed")
        data_path_telescope = config_stp.get_data_path()
    if telescope == "UM":
        data_telescope = config_um.load_config_values("parsed")
        data_path_telescope = config_um.get_data_path()
    if telescope == "test_STP":
        test_loader = utils_config.ConfigLoader('test_data/test_STP', "parsed")
        data_telescope = test_loader.load_configs()
        data_path_telescope = config_stp.get_data_path()
    if telescope == "test_UM":
        test_loader = utils_config.ConfigLoader('test_data/test_UM', "parsed")
        data_telescope = test_loader.load_configs()
        data_path_telescope = config_um.get_data_path()

    zodi_magnitude_normalization = float(data_telescope['astrophysics']['zodi']['zodi_mag_r'])
    obs = etsc.Observatory(telescope)
    obs.make_STP()
    obs.set_source(source_pickles_file='test_data/pickles_uk_9.fits',
                   plot=True)
    obs.set_background(background_file=data_telescope['astrophysics']['zodi']['profile'], support_data_path=data_path_telescope,
                       plot=True)
    obs.make_observation(flux=0.0, flux_units='AB', bg_flux=zodi_magnitude_normalization, bg_flux_units=u.ABmag)

    #   For Version 1.0.0
    if telescope == "UM":
        assert round(obs.source_counts.value) == 16230700376
        assert round(obs.sky_counts.value) == 20
        assert round(obs.calc_SNR(exp_time=1.0 * u.s, int_time=1.0 * u.s)) == 127400
        assert round(obs.calc_req_source(10.0, int_time=1.0 * u.s, exp_time=1.0 * u.s)) == 155
        assert round(obs.calc_req_source(10.0, int_time=1.0 * u.s, exp_time=1.0 * u.s, magnitude=True)[0], 2) == 20.05
        assert round(obs.calc_int_time(1e6, exp_time=1.0*u.s).value) == 62
        assert round(obs.calc_saturation_time().value, 6) == 3.6e-5
    if telescope == "test_UM":
        assert round(obs.source_counts.value) == 16147078162
        assert round(obs.sky_counts.value) == 20
        assert round(obs.calc_SNR(exp_time=1.0 * u.s, int_time=1.0 * u.s)) == 127071
        assert round(obs.calc_req_source(10.0, int_time=1.0 * u.s, exp_time=1.0 * u.s)) == 155
        assert round(obs.calc_req_source(10.0, int_time=1.0 * u.s, exp_time=1.0 * u.s, magnitude=True)[0], 2) == 20.05
        assert round(obs.calc_int_time(1e6, exp_time=1.0*u.s).value) == 62
        assert round(obs.calc_saturation_time().value, 6) == 3.6e-5
    #   Version for config_STP 1.0.0
    if telescope == "STP" or telescope == "test_STP":
        assert round(obs.source_counts.value) == 64610011021
        assert round(obs.sky_counts.value) == 79
        assert round(obs.calc_SNR(exp_time=1.0 * u.s, int_time=1.0 * u.s)) == 254185
        assert round(obs.calc_req_source(10.0, int_time=1.0 * u.s, exp_time=1.0 * u.s)) == 155
        assert round(obs.calc_req_source(10.0, int_time=1.0 * u.s, exp_time=1.0 * u.s, magnitude=True)[0], 2) == 21.55
        assert round(obs.calc_int_time(1e6, exp_time=1.0*u.s).value) == 15
        assert round(obs.calc_saturation_time().value, 6) == 9e-6
    return

@pytest.mark.parametrize("telescope", ["UM", "STP", "test_UM", "test_STP"], indirect=True)
def test_validate_ETC_snr_calculation(telescope):
    if telescope == "STP":
        data_telescope = config_stp.load_config_values("parsed")
        data_path_telescope = config_stp.get_data_path()
    if telescope == "UM":
        data_telescope = config_um.load_config_values("parsed")
        data_path_telescope = config_um.get_data_path()
    if telescope == "test_STP":
        test_loader = utils_config.ConfigLoader('test_data/test_STP', "parsed")
        data_telescope = test_loader.load_configs()
        data_path_telescope = config_stp.get_data_path()
    if telescope == "test_UM":
        test_loader = utils_config.ConfigLoader('test_data/test_UM', "parsed")
        data_telescope = test_loader.load_configs()
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

