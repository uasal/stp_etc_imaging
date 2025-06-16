"""
ExposureTimeSNRCalculator.py

This code facilitates the calculation of both exposure times/snr for a given snr/exposure time respectively.

The code contained here is derived from the HST Exposure Time calculator from the STSci website. 
Link here: https://etc.stsci.edu/etc/results/ACS.im.1813528/

The code uses the synphot package developed by STSci. 
Link here: https://synphot.readthedocs.io/en/latest/
"""

from astropy import units as u
from synphot import units, SourceSpectrum, SpectralElement, Observation, Empirical1D
from synphot.models import BlackBodyNorm1D, GaussianFlux1D, Box1D
import numpy as np
from numpy import sqrt
from scipy.interpolate import interp1d
import pandas as pd
from math import ceil, floor, log10

### Defining Classes

try:
    import config_stp
except:
    print("config_stp could not be imported. Please check installation.")
try:
    import config_um
except:
    print("config_um could not be imported. Please check installation.")
try:
    import config_stp_wcc
except:
    print("config_stp could not be imported. Please check installation.")

import utils_config
data_stp_wcc = config_stp_wcc.load_config_values()

#   One line function prepends the support path to variable s2 if s1 is provided, else returns s2 as is.
prepend_if_not_none = lambda s1, s2: f"{s1}{s2}" if s1 is not None else s2

####################################
class Observatory:
    
    attribute = "space observatory"
    
    def __init__(self, telescope_name):

        self.name = telescope_name
        self.diameter_primary = None
        self.surf_area = None
        self.num_mirrors = 0
        self.focal_len = None
        self.f_num = None
        self.psf_diameter = None
        self.gain = None
        self.dark_current = None
        self.read_noise = None
        self.num_pixels = None
        self.num_psf_pixels = None
        self.psf_area = None
        self.pixel_size = None
        self.sensor_area = None
        self.sensor_temp = None
        self.plate_scale = None
        self.sky_counts = 0
        self.well_depth = None
        self.bandpass = SpectralElement(Box1D, amplitude=1, x_0=6000, width=7000)
        self.qe_curves = []
        self.filters = []
        self.qe_wpeak = None
        self.source_spectrum = None
        self.source_name = None
        self.source_z = None
        self.source_counts = None
        self.background_spectrum = None
        self.background_name = None
        self.bg_surface_brightness = None   #   The surface brightness of the background source
        self.bg_magnitude = None        #   The AB magnitude of the background source (Requires PSF area to calculate)
        self.zero_mag_counts = None     #   Source counts for an AB mag = 0 PSF


    # Add Quantum Efficiency properties
    def add_qe_curve(self, qe_fits_file, wave_unit='angstrom', num_curves=1, plot=False):
                
        self.qe_curves.append(f"{qe_fits_file} x {num_curves}")
        bp = SpectralElement.from_file(qe_fits_file, wave_unit=wave_unit)
        for num in range(num_curves-1):
            bp *= SpectralElement.from_file(qe_fits_file, wave_unit=wave_unit)
        self.bandpass *= bp

        if plot==True:
            bp.plot()
    
    # Add Sensor
    def add_sensor(self, num_curves=1
                   , gain_setting= 100  # (0.1 dB)
                   , sensor_temp = 0 * u.Celsius
                   , sensor_area = None
                   , sensor_pixel_size = None
                   , gain = None
                   , dark_current = None
                   , read_noise = None
                   , well_depth = None
                   , sensor_toml = None
                   , support_data_path=None
                   , plot=False
                   , plot_title="Sensor"
                   ):
        
        if sensor_temp['unit'] != u.Celsius:
            sensor_temp.to(u.Celsius, equivalencies=u.temperature())

        wave_unit = 'nm'

        # physical settings
        if sensor_area is None:
            sensor_area = 36*u.mm * 24*u.mm
        if sensor_pixel_size is None:
            sensor_pixel_size = 3.76*(u.um/u.pix)

        #   plate scale [arcsecond/pixel] = pixel size / (telescope diameter * f#)  * ( 206265 * arcseconds/radian)
        self.plate_scale = (sensor_pixel_size.value * 1e-6 /self.diameter_primary['value'] / self.f_num * 206265)

        # gain(gain-setting)

        if sensor_toml is not None:

            gain_cols = ['gain_setting', 'gain']
            gain = self.get_interpolated_value(prepend_if_not_none(support_data_path, sensor_toml['gain_curve']), gain_setting, gain_cols) * (u.electron / u.ct)

            dark_current_cols = ['sensor_temperature', 'dark_current']
            dark_current = self.get_interpolated_value(prepend_if_not_none(support_data_path, sensor_toml['dark_current']), sensor_temp['value'], dark_current_cols) * (u.electron / (u.s * u.pix))

            read_noise_cols = ['gain_setting', 'read_noise']
            read_noise = self.get_interpolated_value(prepend_if_not_none(support_data_path, sensor_toml['read_noise']), gain_setting, read_noise_cols) * sqrt(1.0 * u.electron / u.pix)

            well_depth_cols = ['gain_setting', 'well_depth']
            well_depth = self.get_interpolated_value(prepend_if_not_none(support_data_path, sensor_toml['well_depth']), gain_setting, well_depth_cols) * (u.electron / u.pix)

        sensor_qe = prepend_if_not_none(support_data_path, sensor_toml['qe'])
        self.qe_curves.append(f"{sensor_qe} x {num_curves}")

        self.sensor_area = sensor_area
        self.pixel_size = sensor_pixel_size

        sensor_num_pixels = sensor_area / (sensor_pixel_size**2)
        self.num_pixels = sensor_num_pixels.to('pix2').value * (u.pix)

        self.gain = gain
        self.dark_current = dark_current
        self.read_noise = read_noise
        self.well_depth = well_depth
        self.sensor_temp = sensor_temp
        if sensor_qe[-5:] == '.fits':
            sensor_qe = sensor_qe[:-5] + '.csv'
        bp = SpectralElement.from_file(sensor_qe, wave_unit=wave_unit)

        for num in range(num_curves-1):
            bp *= SpectralElement.from_file(sensor_qe, wave_unit=wave_unit)
        self.bandpass *= bp

        if plot==True:
            bp.plot(title=plot_title)
        
    # Add mirrors
    def add_mirror(self, mirror_qe_fits_file, num_curves=1 ,wave_unit='nm', support_data_path=None, plot=False,plot_title="Mirror"):

        self.num_mirrors += num_curves
        mirror_qe_fits_file = prepend_if_not_none(support_data_path,mirror_qe_fits_file)
        self.qe_curves.append(f"{mirror_qe_fits_file} x {num_curves}")
        
        bp = SpectralElement.from_file(mirror_qe_fits_file, wave_unit=wave_unit)
        for num in range(num_curves-1):
            bp *= SpectralElement.from_file(mirror_qe_fits_file, wave_unit=wave_unit)
        self.bandpass *= bp
        
        if plot==True:
            bp.plot(title=plot_title)

        return


    # Add filter properties
    def add_filter(self, filter_fits_file, wave_unit='angstrom', num_curves=1, support_data_path=None, plot=False, plot_title="Filter"):
        
        filter_fits_file = prepend_if_not_none(support_data_path, filter_fits_file)
        self.filters.append(f"{filter_fits_file} x {num_curves}")
        
        bp = SpectralElement.from_file(filter_fits_file)

        for num in range(num_curves-1):
            bp *= SpectralElement.from_file(filter_fits_file)

        self.bandpass *= bp
        
        if plot==True:
            bp.plot(title=plot_title)
    
    # Calculate the mean PSF based on the mean wavelength of combined Spectral Elements
    # OR from a user-selected wavelength
    def calc_PSF(self, wavelength=None, approx_type='sq', verbose=False):
        
        if not wavelength == None:
            psf_diameter = 2*1.22*wavelength*self.f_num
        else:
            wavelength = self.bandpass.wpeak()
            self.qe_wpeak = wavelength
            psf_diameter = (2*1.22*wavelength*self.f_num)
            
        self.psf_diameter = psf_diameter.to('um')
        
        # determine total number of pixels with the psf
        n = ceil((self.psf_diameter / self.pixel_size).value)
        
        if approx_type.lower() in ['s', 'sq', 'square']:  # Gives result of nxn square
            self.num_psf_pixels = n**2 * (u.pix)
        elif approx_type.lower() in ['c', 'circ', 'circle', 'circular']:  # Gives result of all pixels with centers within a circle of radius n
            if n%2 == 0:  # n even
                L = int(n/2)
                sequence = [floor(0.5 + 0.5*sqrt(n**2 - (2*y-1)**2)) for y in list(range(1,L+1))]
                self.num_psf_pixels = 4 * sum(sequence) * (u.pix)
            elif n%2 == 1:  # n odd
                L = int((n-1)/2)
                sequence = [floor(1 + 0.5*sqrt(n**2 - 4*(y**2))) for y in list(range(1,L+1))]
                self.num_psf_pixels = 1 + 4 * sum(sequence) * (u.pix)

        #   Total PSF Area = Area of 1 pixel (arcsecond * arcsecond)/pixel * num of pixels (pixels)
        #   Gives area in square arcseconds IF self.plate_scale is in arcseconds/pixel
        self.psf_area = (self.plate_scale*self.plate_scale) * self.num_psf_pixels

        #   Calculate the total magnitude of the background spectrum given the PSF area (in square arcseconds)
        self.calculate_bg_normalization_magnitude()

        if verbose:
            print(f"Number of PSF pixels: {self.num_psf_pixels}")
        return psf_diameter, self.num_psf_pixels
        
    # NOTE: This function is "complete", meaning you would only need to add observational parameters
    # to make an observation ('set_source', 'set_background', and 'make_observation').

    def make_STP(self, plot=False
                 , filters=True
                 , sensor=True
                 , mirrors=True
                 , bg_sb = None
                 , custom_toml_dir = None
                 , support_data_dir = None
                 # , left=4000 ,right=8000  # for plotting (units in angstroms)
                 ):

        if self.name == 'STP':
            data_telescope = config_stp.load_config_values("parsed")
            support_data_telescope = config_stp.get_data_path()

        if self.name == 'UM':
            data_telescope = config_um.load_config_values("parsed")
            support_data_telescope = config_um.get_data_path()

        if self.name == 'test_STP':
            test_loader = utils_config.ConfigLoader('tests/test_data/test_STP', "parsed")
            data_telescope = test_loader.load_configs()
            support_data_telescope = config_stp.get_data_path()

        if self.name == 'test_UM':
            test_loader = utils_config.ConfigLoader('tests/test_data/test_UM', "parsed")
            data_telescope = test_loader.load_configs()
            support_data_telescope = config_um.get_data_path()
        if self.name == 'custom_toml':
            if custom_toml_dir is None:
                print("Directory for custom toml files must be provided. Exiting...")
                exit()
            if support_data_dir is None:
                print("Support data directory for custom configuration must be provided. Exiting...")
                exit()
            custom_loader = utils_config.ConfigLoader(custom_toml_dir, "parsed")
            data_telescope = custom_loader.load_configs()
            support_data_telescope = support_data_dir



        self.data_observatory = data_telescope
        self.data_wcc = config_stp_wcc.load_config_values("parsed")

        #   Setup Paths

        self.support_data_observatory = support_data_telescope if support_data_telescope.endswith('/') else support_data_telescope + '/'
        support_data_wcc = config_stp_wcc.get_data_path()
        self.jitter_rms         = self.data_observatory['observatory']['pointing']['jitter_rms']
        self.diameter_primary = self.data_observatory['telescope']['optics']['m1']['aper_clear_OD']
        self.surf_area = np.pi * (0.5 * self.diameter_primary['value']) ** 2
        self.num_mirrors = 0
        self.rms_surf           = self.data_observatory['telescope']['optics']['m1']['surface_rms']
        self.f_num = self.data_observatory['telescope']['general']['f_number']  # telescope_focal_len / telescope_diameter
        self.focal_len = self.f_num * self.diameter_primary['value']
        self.diameter_secondary = self.data_observatory['telescope']['optics']['m2']['aper_clear_OD']
        self.support_width      = self.data_observatory['telescope']['optics']['m2']['support_width']
        self.n_supports         = self.data_observatory['telescope']['optics']['m2']['n_supports']

        # account for light obscuration by support structures and secondary mirror
        if self.data_observatory['observatory']['description'] == "STP Spacecraft Baseline Configuration":
            self.surf_area          = self.surf_area - (np.pi * (0.5*self.diameter_secondary['value'])**2) - (self.n_supports*self.support_width['value']*(self.diameter_primary['value']-self.diameter_secondary['value']))

        #   Fixed units here, when passing this area to countrate(area=<>) the default assumption is cm^2
        self.surf_area *= (u.meter * u.meter)
        self.sensor_temp        = self.data_wcc['common_params']['sensor']['temp_nominal']
        self.gain_setting = self.data_wcc['common_params']['sensor']['gain']
        if bg_sb is None:
            self.bg_surface_brightness = float(self.data_observatory['astrophysics']['zodi']['zodi_mag_r'])
        else:
            self.bg_surface_brightness = bg_sb

        if sensor:
            # Add STP sensor
            if plot:
                print("Adding Sensor")
            self.add_sensor( num_curves=1
                           , gain_setting= self.gain_setting # (0.1 dB)
                           , sensor_temp = self.sensor_temp
                           , sensor_area = None
                           , sensor_pixel_size = None
                           , gain = None
                           , dark_current = None
                           , read_noise = None
                           , well_depth = None
                           , sensor_toml = self.data_wcc['common_params']['sensor']
                           , support_data_path=support_data_wcc
                           , plot=plot
                           , plot_title = "ZWO_ASI6200MM"
                           )

        # Add mirrors

        if mirrors:
            if plot:
                print("Adding Mirrors")
            self.add_mirror(self.data_observatory['telescope']['optics']['m1']['coating_refl'], num_curves=1, support_data_path=self.support_data_observatory, plot=plot, plot_title="M1")   # M1
            self.add_mirror(self.data_observatory['telescope']['optics']['m2']['coating_refl'], num_curves=1, support_data_path=self.support_data_observatory, plot=plot, plot_title="M2")   # M2
            self.add_mirror(self.data_observatory['telescope']['optics']['m3']['coating_refl'], num_curves=1, support_data_path=self.support_data_observatory, plot=plot, plot_title="M3")   # M3
            self.add_mirror(self.data_observatory['telescope']['optics']['m4']['coating_refl'], num_curves=1, support_data_path=self.support_data_observatory, plot=plot, plot_title="M4")   # M4


        if filters:
            if plot:
                print("Adding Filters")
            # Add filter (Option for not using a filter)
            self.add_filter(filter_fits_file=self.data_wcc['WCC_ETC']['filters']['filter'], wave_unit='angstrom'
                            , num_curves=1, support_data_path=support_data_wcc, plot=plot, plot_title=self.data_wcc['WCC_ETC']['filters']['filter'])

        # Calculate PSF
        if sensor:
            self.calc_PSF(wavelength=None, approx_type='sq')

    # Create source for observation
    def set_source(self, source_pickles_file, source_z=0, support_data_path=None, plot=False):

        """
        Source: https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/pickles-atlas.html
        Source Spectrum: Pickles Stellar Atlas
        filename        sptype    T_eff
        --------------  --------  -------
        pickles_uk_1	O5V	      39810.7
        pickles_uk_2	O9V	      35481.4
        pickles_uk_3	B0V	      28183.8
        pickles_uk_4	B1V	      22387.2
        pickles_uk_5	B3V	      19054.6
        pickles_uk_6	B5-7V	  14125.4
        pickles_uk_7	B8V	      11749.0
        pickles_uk_9	A0V	      9549.93
        pickles_uk_10	A2V	      8912.51
        pickles_uk_11	A3V	      8790.23
        pickles_uk_12	A5V	      8491.80
        pickles_uk_14	F0V	      7211.08
        pickles_uk_15	F2V	      6776.42
        pickles_uk_16	F5V	      6531.31
        pickles_uk_20	F8V	      6039.48
        pickles_uk_23	G0V	      5807.64
        pickles_uk_26	G2V	      5636.38 ***
        pickles_uk_27	G5V	      5584.70
        pickles_uk_30	G8V	      5333.35
        pickles_uk_31	K0V	      5188.00
        pickles_uk_33	K2V	      4886.52
        pickles_uk_36	K5V	      4187.94
        pickles_uk_37	K7V	      3999.45
        pickles_uk_38	M0V	      3801.89
        pickles_uk_40	M2V	      3548.13
        pickles_uk_43	M4V	      3111.72 ***
        pickles_uk_44	M5V	      2951.21
        pickles_uk_46	B2IV	  19952.6
        pickles_uk_47	B6IV	  12589.3
        pickles_uk_48	A0IV	  9727.47
        pickles_uk_49	A4-7IV	  7943.28
        pickles_uk_50	F0-2IV	  7030.72
        pickles_uk_51	F5IV	  6561.45
        pickles_uk_52	F8IV	  6151.77
        pickles_uk_53	G0IV	  5929.25
        pickles_uk_54	G2IV	  5688.53
        pickles_uk_55	G5IV	  5597.57
        pickles_uk_56	G8IV	  5308.84
        pickles_uk_57	K0IV	  5011.87
        pickles_uk_58	K1IV	  4786.30
        pickles_uk_59	K3IV	  4570.88
        pickles_uk_60	O8III	  31622.8
        pickles_uk_61	B1-2III	  19952.6
        pickles_uk_63	B5III	  14791.1
        pickles_uk_64	B9III	  11091.8
        pickles_uk_65	A0III	  9571.94
        pickles_uk_67	A5III	  8452.79
        pickles_uk_69	F0III	  7585.78
        pickles_uk_71	F5III	  6531.31
        pickles_uk_72	G0III	  5610.48
        pickles_uk_73	G5III	  5164.16
        pickles_uk_76	G8III	  5011.87
        pickles_uk_78	K0III	  4852.89
        pickles_uk_87	K3III	  4365.16
        pickles_uk_93	K5III	  4008.67
        pickles_uk_95	M0III	  3819.44
        pickles_uk_100	M5III	  3419.79
        pickles_uk_105	M10III	  2500.35
        pickles_uk_106	B2II	  15995.6
        pickles_uk_107	B5II	  12589.3
        pickles_uk_108	F0II	  7943.28
        pickles_uk_109	F2II	  7328.25
        pickles_uk_110	G5II	  5248.07
        pickles_uk_111	K0-1II	  5011.87
        pickles_uk_112	K3-4II	  4255.98
        pickles_uk_113	M3II	  3411.93
        pickles_uk_114	B0I	      26001.6
        pickles_uk_117	B5I	      13396.8
        pickles_uk_118	B8I	      11194.4
        pickles_uk_119	A0I	      9727.47
        pickles_uk_121	F0I	      7691.30
        pickles_uk_122	F5I	      6637.43
        pickles_uk_123	F8I	      6095.37
        pickles_uk_124	G0I	      5508.08
        pickles_uk_126	G5I	      5046.61
        pickles_uk_127	G8I	      4591.98
        pickles_uk_128	K2I	      4255.98
        pickles_uk_130	K4I	      3990.25
        pickles_uk_131	M2I	      3451.44
        """

        pickles_file_path = prepend_if_not_none(support_data_path,source_pickles_file)
        self.source_spectrum = SourceSpectrum.from_file(pickles_file_path)

        self.source_spectrum.z = source_z
        self.source_name = source_pickles_file
        self.source_z = source_z

        if plot==True:
            self.source_spectrum.plot()

    def create_supernova_spectrum(self, epoch, source_z = 0.0, support_data_path=None, synphot_spectrum=True, plot=False):
        """
        Read spectral data from a file and return wavelength and flux arrays for a given epoch.

        Parameters:
        filename : str
            Path to the file containing spectral data.
        epoch : float
            The epoch value to extract the corresponding spectrum.

        Returns:
        tuple
            A tuple containing wavelength and flux arrays.
        """
        template_file = self.data_observatory['astrophysics']['sources']['supernova_templates']
        if support_data_path is not None:
            template_file = prepend_if_not_none(support_data_path, self.data_observatory['astrophysics']['sources']['supernova_templates'])

        data = np.loadtxt(template_file)

        # Extract relevant rows for the given epoch
        mask = data[:, 0] == epoch
        filtered_data = data[mask]

        if filtered_data.size == 0:
            raise ValueError(f"Epoch {epoch} not found in the file.")

        wavelength = filtered_data[:, 1] * u.AA  # Wavelength in Angstroms
        flux = filtered_data[:, 2] * units.FLAM  # Flux in flambda units

        if synphot_spectrum:

            self.source_spectrum = SourceSpectrum(Empirical1D, points=wavelength, lookup_table=flux, z=source_z)
            self.source_spectrum.z = source_z
            self.source_z = source_z
            self.source_name = self.data_observatory['astrophysics']['sources']['supernova_templates']

            if plot == True:
                self.source_spectrum.plot()
            return
        else:
            return wavelength, flux

    def set_background(self, background_file, support_data_path=None,  plot=False):

        # Create background for observation
        # info from: https://etc.stsci.edu/etcstatic/users_guide/1_ref_9_background.html
        if support_data_path is not None:
            background_file = prepend_if_not_none(support_data_path, background_file)
        self.background_spectrum = SourceSpectrum.from_file(background_file)
        self.background_name = background_file

        if plot==True:
            self.background_spectrum.plot(#left=left, right=right
                                         )

    def make_observation(self, flux=0, flux_units=u.ABmag, bg_flux=None, bg_flux_units=u.ABmag, plot=False):
        # Create observation using source and return countrate

        #   Calculate Background Normalization magnitude
        if bg_flux is None:
            if self.bg_magnitude is None:
                try:
                    self.calculate_bg_normalization_magnitude()
                except:
                    print("Error: Could not calculate Background magnitude.")
                    exit()
            bg_flux = self.bg_magnitude

        # Add Source
        self.source_spectrum.z = self.source_z  # make sure source spectra has proper redshift

        if flux_units in ['vega', units.VEGAMAG]:
            vega = SourceSpectrum.from_vega()  # For unit conversion
            normalization_units = flux * units.VEGAMAG
            sp_rn = self.source_spectrum.normalize(normalization_units
                                                      , self.bandpass
                                                      , vegaspec=vega
                                                      # , force='taper'
                                                      , force='extrap'
                                                     )

        elif flux_units in ['AB', 'ABmag', 'AB mag', 'AB magnitude', u.ABmag]:
            normalization_units = flux * u.ABmag
            sp_rn = self.source_spectrum.normalize(normalization_units
                                                      , self.bandpass
                                                      # , vegaspec=vega
                                                      # , force='taper'
                                                      , force='extrap'
                                                     )
        else:
            raise NotImplementedError("User-defined source flux units not currently implemented")

        sp_obs = Observation(sp_rn, self.bandpass, force='extrap')

        if plot==True:
            sp_obs.plot(title='Source')

        # Get countrate for observation
        source_counts = sp_obs.countrate(area=self.surf_area) * u.electron/u.ct
        self.source_counts = source_counts    #    Makes the counts in units of e/s

        # Add Background
        # Background needs to be normalized in the Johnson V band
        johnson_v_passband =  SpectralElement.from_filter('johnson_v')
        bg_rn = self.background_spectrum.normalize( bg_flux * bg_flux_units
                                                      , johnson_v_passband
                                                      #, vegaspec=vega
                                                      # , force='taper'
                                                      , force='extrap'
                                                     )

        bg_obs = Observation(bg_rn, self.bandpass, force='extrap')

        if plot==True:
            bg_obs.plot(title='Background')

        # Get countrate for observation
        background_counts = bg_obs.countrate(area=self.surf_area)* u.electron/u.ct
        self.sky_counts = background_counts

        return source_counts, background_counts

    def calc_saturation_time(self):
        # Calculates the time for any one pixel on a sensor to completely fill it's well-depth.
        # NOTE: Requires setting of source and/or background and performing 'make_observation'
        return self.well_depth / ( (self.source_counts / self.num_psf_pixels) + (self.sky_counts / self.num_psf_pixels) )

    def calc_SNR(self, int_time, exp_time):
        # Calculate SNR for a given total integration time with set frame exposure times
        # NOTE: Requires setting of source and/or background and performing 'make_observation'

        total_noise = sqrt( self.source_counts*int_time + self.sky_counts*int_time
                              + ( self.dark_current + self.read_noise*self.read_noise/exp_time)*int_time*self.num_psf_pixels)
        snr = self.source_counts*int_time / total_noise

        return snr.value
    
    def calc_int_time(self, snr, exp_time):
        # Calculate total integration time to achieve a given snr and frame exposure time
        # NOTE: Requires setting of source and/or background and performing 'make_observation'
        
        snr = snr * sqrt(1.0 * u.ct)  # to ensure units match
        A = ((self.source_counts/snr)**2) * exp_time
        B = self.source_counts*exp_time + self.sky_counts*exp_time + ( self.dark_current*exp_time + self.read_noise*self.read_noise)*self.num_psf_pixels
        int_time = B/A

        return int_time * u.electron / u.ct

    def calc_zero_mag(self, flux=0, plot=False):
        zero_mag_counts, sky_counts = self.make_observation(flux=flux, flux_units=u.ABmag
                                                            , bg_flux=self.bg_magnitude, bg_flux_units=u.ABmag
                                                            , plot=plot)
        self.zero_mag_counts = zero_mag_counts
        return

    def calc_req_source(self, snr, int_time, exp_time, magnitude=False):
        # Calculate required Source Counts for given exposure time and snr
        # NOTE: can be converted to magnitude if 'magnitude=True' is given.
        # NOTE: Requires setting of source and/or background first
        
        ### FOLLOW UP NOTE (SS 3/13/25) - Rejigged the zero_mag call to run it only once per ETC construction unless called explicitly.

        if self.zero_mag_counts is None:
            self.calc_zero_mag(flux=0,plot=False)
        snr = snr * sqrt(1.0 * u.ct)  # to ensure units match

        A = int_time/snr/snr
        B = -1
        C = -1 * (self.sky_counts + self.num_psf_pixels * ( self.dark_current + self.read_noise*self.read_noise/exp_time))
        req_source_counts = (-B/2/A.value)*(1+sqrt(1-(4*A*C/B/B).value))

        if magnitude==False:
            return req_source_counts
        
        elif magnitude==True:            
            mag = -2.5 * log10(req_source_counts/self.zero_mag_counts.value)
            return mag, self.filters


    def calculate_bg_normalization_magnitude(self):
        """
        Convert the Background Surface Brightness into the total magnitude given the PSF area (in arcseconds squared)
        The area needs to be in square arcseconds since this the typical definition of Surface Brightness is in units
        of magnitudes per arcseconds^2
        :return: None
        """
        self.bg_magnitude = self.bg_surface_brightness - 2.5 * np.log10(self.psf_area.value)
        return

    def set_gain(self, gain):
        # Set gain
        self.gain = gain

    def set_dark_current(self, dark_current):
        # Set dark current
        self.dark_current = dark_current

    def set_read_noise(self, read_noise):
        # Set read noise
        self.read_noise = read_noise

    def set_pixels(self, num_pixels):
        # Set number of pixels
        self.num_pixels = num_pixels

    def set_psf_pixels(self, num_psf_pixels):
        # Set number of psf pixels
        self.num_psf_pixels = num_psf_pixels

    def set_sky_counts(self, sky_counts):
        # Set sky counts
        self.sky_counts = sky_counts

    def get_info(self):
        # Get info as strings
        for key, value in self.__dict__.items():
            print(f'***__{key}={value}', '\n')

    def as_array(self):
        # Get info as array
        return np.array(list(self.__dict__.items()), dtype=object)


    def as_df(self):
        # Get info as pandas data frame
        df = pd.DataFrame(np.transpose(np.array(list(self.__dict__.items()), dtype=object)))
        df.columns = df.iloc[0]
        df = df[1:]
        return df

    def get_interpolated_value(self,input_file, interpolation_xval, col_headers):
        """
        Interpolate the given file columns to get the value at interpolation_xval
        :param input_file: File to use x columns and y columns on
        :param interpolation_xval: Value that the interpolation function takes as argument
        :param col_headers: Names of column headers as a list
        :return: The value of the interpolated function at interpolation_xval
        """

        df = pd.read_csv(input_file, skiprows=1, names=[col_headers[0], col_headers[1]])
        xlist = df[col_headers[0]]
        ylist = df[col_headers[1]]
        interp = interp1d(xlist, ylist)

        return interp(interpolation_xval)