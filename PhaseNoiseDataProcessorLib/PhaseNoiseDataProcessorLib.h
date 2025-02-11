#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "band_data.h"
#include "fft.h"
#include "filter.h"

using namespace std;

#define MAX_DATA_SIZE 268435456
#define MAX_PACKET_SIZE 268435456
#define MAX_MEAS 65536
#define SMPL_FILENAME "samples.bin"
#define OUTPUT_FILENAME "output.csv"

namespace PhaseNoiseDataProcessor
{
	extern std::vector<double> meas_buffer[2]; //!< Sample buffers for captured data from a single measurement
	extern std::vector<double> smpl_buffer[2]; //!< Sample buffers for a single sample window
	extern std::vector<double> fft_re_buffer; //!< FFT real output buffer
	extern std::vector<double> fft_im_buffer; //!< FFT imaginary output buffer
	extern std::vector<double> output_re_buffer; //!< Buffer to store the sum of real values of cross correlations / auto correlations
	extern std::vector<double> output_im_buffer; //!< Buffer to store the sum of imaginary values of cross correlations / auto correlations
	extern std::vector<double> output_x_buffer; //!< Buffer to store x-axis frequency output values
	extern std::vector<double> output_y_buffer; //!< Buffer to store y-axis power spectral density output values

	extern uint32_t meas_count; //!< Number of measurements that have been processed since processing was started

	extern BandData band_data; //!< Generated data descibing the properties of each measurement band

	struct Config
	{
		bool cross_corr = true; //!< True for cross correlation, false for single channel
		double smpl_rate = 524288.0; //!< Sample rate of the captured data
		uint32_t num_meas = 1; //!< Number of measurements to process
		uint32_t data_size = 4194304; //!< Number of samples captured per channel for a single measurement
		uint32_t packet_size = 524288; //!< Minimum size of the FFT window
		uint32_t num_bands = 1; //!< Number of frequency bands to use
		uint32_t band_mult = 8; //!< Frequency multiplier for subsequent bands (must be a power of 2)
		uint32_t band_shift = 2; //!< Shifts each band down in multiples of band_mult by the specified number (prioritises more correlations over smaller RBW for higher frequencies)
		double conversion_fac = 0.0025; //!< Maximum input voltage reference (bipolar)
		bool window_en = false; //!< True to apply flat top window to FFT input data (NOT YET IMPLEMENTED), false for no windowing (rectangular)
	};

	extern Config config; //!< Measurement configuration settings

	//! Convert raw input sample byte data to voltage values and store to meas buffer for the selected channel
	 
	//! Byte data should be in the format of consecutive 16-bit half-words each representing one sample.
	//! Samples are converted to a voltage in the range of -(*conversion_fac*) to +(*conversion_fac*) using an offset binary representation (not 2's complement).
	void convert_samples(int channel, const void* data_bytes);

	//! Calculate and store x-axis frequency output values
	
	//! Values in Hz are calculated based on configuration settings and generated band data.
	//! Should be called after selected processing has been started.
	void calc_x_values();

	//! Calculate and store y-axis power spectral density output values
	
	//! Values in dBV/Hz are calculated based on the combined data of all processed measurements since processing was started.
	//! Should be called after measurement data are pushed for real-time monitoring of phase noise.
	void calc_y_values(uint32_t num_meas);
	void calc_y_values();

	//! Convert y-axis power spectral density output values from dBV/Hz to dBc/Hz based on the specified calibration settings
	
	//! *rf_input_level* is the level of the carrier signal at the output of the device under test in dBV.
	//! *cal_input_level* is the level of the calibration tone applied at the output of the device under test in dBV.
	//! *cal_output_level* is the level of the calibration tone measured at the input of the sample capture device in dBV.
	void dbv_to_dbc(double rf_input_level, double cal_input_level, double cal_output_level);

	//! Read sample data from file for a single measurement and store to data buffers 
	
	//! Byte data should be in the format of consecutive 16-bit half-words each representing one sample.
	//! Samples are converted to a voltage in the range of -(*conversion_fac*) to +(*conversion_fac*) using an offset binary representation (not 2's complement).
	bool read_meas_from_file(char* filename, uint32_t meas_num);

	//! Retrieve a sample packet containing a single data set from meaurement data buffers and store in sample buffers (only used for fixed window size processing)
	
	//! Setting *smpl_rate_div* to a number greater than 1 allows the sample rate to be effectively reduced by skipping samples.
	//! (e.g. setting *smpl_rate_div* to 8 will retrieve the first in every 8 samples)
	void retrieve_samples(uint32_t set_num, uint32_t smpl_rate_div);

	//! Apply a digital filter to sample buffers (only used for fixed window size processing)
	
	//! This is a 6th-order digital Butterworth filter with a cutoff at (*sampl_rate* / *band_mult*).
	//! It allows for under-sampling of the data (such as with the retrieve_samples() function) with reduced aliasing artefacts.
	void filter_samples(AntiAliasingFilter& filter_0, AntiAliasingFilter& filter_1);

	//! Overwrite selected packet in measurement data buffers with the contents of sample buffers (only used for fixed window size processing)
	void overwrite_samples(uint32_t packet_num);

	//! Set all values in measurement data buffers to 0
	void clear_meas_buffer();

	//! Set all values in output buffers to 0
	void clear_output_buffer();

	//! Setup processing functions based on the specified parameters
	void setup(bool cross_corr, double smpl_rate, uint32_t num_meas, uint32_t data_size, uint32_t min_window_size, uint32_t num_bands, uint32_t band_mult, uint32_t band_shift, double conversion_fac);
	void setup(Config new_config);

	//! Return vector containing x-axis frequency output values
	const std::vector<double>& get_x_values();

	//! Return vector containing y-axis power spectral density output values
	const std::vector<double>& get_y_values();

	//! Write output data points to file
	
	//! The file is in in CSV format and contains the properties of each frequency band followed by the list of power spectral density output values against frequency.
	bool write_output_to_file(char* filename);

	//! Process sample data from file to produce a set of output data points using a fixed window size
	
	//!Frequency bands are generated with a fixed window size and variable sampling rate (faster).
	void process_fixed_window_size();

	//! Process sample data from file to produce a set of output data points using a fixed sampling rate
	
	//! Frequency bands are generated with a fixed sampling rate and variable window size (more accurate).
	void process_fixed_smpl_rate();

	//! Start fixed window size processing where sample data is pushed in real time
	void start_fixed_window_size();

	//! Push sample data for a single measurement to be processed with a fixed window size
	bool push_data_fixed_window_size(double* data_0, double* data_1);

	//! Start fixed sample rate processing where sample data is pushed in real time
	void start_fixed_smpl_rate();

	//! Push sample data for a single measurement to be processed with a fixed sample rate
	bool push_data_fixed_smpl_rate(double* data_0, double* data_1);

	//! Stop real time processing and release data buffers
	void stop();
}


