SLSTR FM02 Spectral Response README
-----------------------------------

        Tim Nightingale (tim.nightingale@stfc.ac.uk) 28/05/2015



This note accompanies the nine spectral response files:

        SLSTR_FM02_S1_20150122.nc
	SLSTR_FM02_S2_20150122.nc
	SLSTR_FM02_S3_20150122.nc
        SLSTR_FM02_S4_20150122.nc
	SLSTR_FM02_S5_20150122.nc
	SLSTR_FM02_S6_20150122.nc
        SLSTR_FM02_S7_20150122.nc
	SLSTR_FM02_S8_20150122.nc
	SLSTR_FM02_S9_20150122.nc

These files contain the draft measured spectral responses of the nine standard channels S1 - S9 in the focal plane of the SLSTR FM02 instrument, to fly on the Sentinel 3A platform. The responses of the fire channels F1 and F2 are identical to the responses of standard channels S7 and S8. The responses are derived from measurements taken during a spectral calibration campaign conducted at the STFC Rutherford Appleton Laboratory (RAL) in July and August 2014 (www.stfc.ac.uk).

The files are encoded in netCDF-4 and follow the Climate and Forecast (CF) metadata convention (cfconventions.org). The basic structures of the files are identical. As an example, the structure of file SLSTR_FM02_S1_20150122 is shown below:

netcdf SLSTR_FM02_S1_20150122 {
dimensions:
	n_data = 4001 ;
variables:
	double wavelength(n_data) ;
		wavelength:long_name = "wavelength" ;
		wavelength:units = "µm" ;
	double response(n_data) ;
		response:long_name = "spectral response" ;
		response:units = "1" ;
		response:coordinates = "wavenumber" ;
		response:ancillary_variables = "uncertainty_a uncertainty_b" ;
		response:_FillValue = -1.e+37 ;
		response:scale_factor = 0.000774030680430831 ;
		response:add_offset = 0. ;
	double uncertainty_a(n_data) ;
		uncertainty_a:long_name = "spectral response uncertainty (type A)" ;
		uncertainty_a:units = "1" ;
		uncertainty_a:coordinates = "wavenumber" ;
		uncertainty_a:_FillValue = -1.e+37 ;
		uncertainty_a:scale_factor = 0.000774030680430831 ;
		uncertainty_a:add_offset = 0. ;
	double uncertainty_b(n_data) ;
		uncertainty_b:long_name = "spectral response uncertainty (type B)" ;
		uncertainty_b:units = "1" ;
		uncertainty_b:coordinates = "wavenumber" ;
		uncertainty_b:_FillValue = -1.e+37 ;
		uncertainty_b:scale_factor = 0.000774030680430831 ;
		uncertainty_b:add_offset = 0. ;
	double relative_variance_a(n_data) ;
		relative_variance_a:long_name = "spectral response relative variance (type A)" ;
		relative_variance_a:units = "1" ;
		relative_variance_a:coordinates = "wavenumber" ;
		relative_variance_a:_FillValue = -1.e+37 ;
	double relative_variance_b(n_data) ;
		relative_variance_b:long_name = "spectral response relative variance (type B)" ;
		relative_variance_b:units = "1" ;
		relative_variance_b:coordinates = "wavenumber" ;
		relative_variance_b:_FillValue = -1.e+37 ;

// global attributes:
		:title = "Resampled Sentinel-3A SLSTR S1 unpolarised FPA response at 87 K" ;
		:created = "2015-01-22T17:20:15Z" ;
		:Conventions = "CF-1.6" ;
		:institution = "STFC Rutherford Appleton Laboratory" ;
		:source = "Bruker IFS 120A spectrometer measurements" ;
		:history = "2015-01-22T17:20:15Z: IDL> fpa_resample" ;
		:comment = "" ;
		:data_sources = 
"S3A_0802_01_S1_T087_P999_FTS_obs_20140718T155657Z.nc
    S3A_0802_01_S1_T087_P999_FTS_obs_20140718T155657Z.0.dat
    S3A_0801_01_S1_T087_P999_FTS_obr_20140718T155257Z.0.dat
    S3A_0803_01_S1_T087_P999_FTS_obr_20140718T160040Z.0.dat
    S3A_1400_01_SV_T999_P999_FTS_tta_20141217T080941Z.nc
        S3A_1404_02_SV_T999_P999_FTS_sga_20141217T100938Z.0.dat
        S3A_1405_02_SV_T999_P999_FTS_sga_20141217T100948Z.0.dat
        S3A_1402_02_SV_T999_P999_FTS_sgb_20141217T083821Z.0.dat
        S3A_1407_02_SV_T999_P999_FTS_sgb_20141217T130020Z.0.dat
        S3A_1401_02_SV_T999_P999_FTS_rfa_20141217T080941Z.0.dat
        S3A_1408_02_SV_T999_P999_FTS_rfa_20141217T132904Z.0.dat
        S3A_1403_02_SV_T999_P999_FTS_rfb_20141217T100913Z.0.dat
        S3A_1406_02_SV_T999_P999_FTS_rfb_20141217T101008Z.0.dat
    diode_182_20130625.nc" ;
}

The file contains a number of global attributes containing basic information about the file, including an indented list “data_sources” tracing the contributing measurement files, and a set of variables describing the spectral response:

1) Variable “wavelength” contains the wavelength axis in units of microns.

2) Variable “response” contains an estimate, in raw measurement units, of the spectral responsivity of the channel to unpolarised light. Multiplication of the raw response values by the scaling factor “response:scale_factor” normalises the response to a peak value of 1.0. The draft calibration report “S3-RP-RAL-SL-102 SLSTR B FPA Spectral Calibration - Draft 0.2.pdf” contains an discussion of the derivation of the spectral responsivities from the direct measurements.

3) Variables “uncertainty_a” and “uncertainty_b” contain estimates, in raw measurement units, of the one-sigma type A and type B uncertainties associated with the measured spectral responsivity. Multiplication of the raw uncertainty values by the scaling factors “uncertainty_a:scale_factor” and “uncertainty_b:scale_factor” normalises them to the uncertainties associated with responses with a peak value of 1.0. Type A uncertainties are derived from the statistics of repeated measurements (e.g. detector noise). Type B uncertainties are those derived by other than statistical methods (e.g. uncertainties in the positions of the measurement baselines). See e.g. physics.nist.gov/cuu/Uncertainty/basic.html. The total estimated measurement uncertainty is the root sum square (RSS) of the two quantities. The draft calibration report contains a discussion of the derivations of the uncertainties.

4) Variables “relative_variance_a” and “relative_variance_b” are derived values: relative_variance_a = (uncertainty_a / response)^2 and relative_variance_b = (uncertainty_b / response)^2.
