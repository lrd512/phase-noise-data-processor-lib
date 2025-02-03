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

namespace ResidualNoiseDataProcessor
{
	extern std::vector<double> data_set[2]; //sample buffers for captured data from a single measurement
	extern std::vector<double> smpl_buffer[2]; //sample buffers for a single sample window
	extern std::vector<double> fft_re_buffer; //FFT real output buffer
	extern std::vector<double> fft_im_buffer; //FFT imaginary output buffer
	extern std::vector<double> output_re_buffer; //buffer to store the sum of real values of cross correlations / auto correlations
	extern std::vector<double> output_im_buffer; //buffer to store the sum of imaginary values of cross correlations / auto correlations
	extern std::vector<double> output_x_buffer; //buffer to store x-axis frequency output values
	extern std::vector<double> output_y_buffer; //buffer to store y-axis power spectral density output values

	extern uint32_t meas_count;

	extern BandData band_data;

	struct Config
	{
		bool cross_corr = true; //true for cross correlation, false for single channel
		double smpl_rate = 524288.0; //sample rate of the captured data
		uint32_t num_meas = 1; //number of measurements to process
		uint32_t data_size = 4194304; //number of samples captured per channel for a single measurement
		uint32_t packet_size = 524288; //minimum size of the FFT window
		uint32_t num_bands = 1; //number of frequency bands to use
		uint32_t band_mult = 8; //frequency multiplier for subsequent bands (must be a power of 2)
		uint32_t band_shift = 2; //shifts each band down in multiples of band_mult by the specified number (prioritises more correlations over smaller RBW for higher frequencies)
		double conversion_fac = 0.0025; //maximum input voltage reference (bipolar)
		bool window_en = false; //true to apply flat top window to FFT input data, false for no windowing (rectangular)
	};

	extern Config config;

	//Convert raw input sample byte data to voltage values and store to data buffers
	void convert_samples(int set_num, const void* data_bytes);

	//Calculate and store x-axis frequency output values
	void calc_x_values();

	//Calculate and store y-axis power spectral density output values
	void calc_y_values(uint32_t num_meas);
	void calc_y_values();

	//convert y-axis power spectral density output values from dBV/Hz to dBc/Hz based on the specified calibration settings
	void dbv_to_dbc(double rf_input_level, double cal_input_level, double cal_output_level);

	//Read sample data from file for a single measurement and store to data buffers 
	bool read_meas_from_file(char* filename, uint32_t meas_num);

	//Retrieve a sample packet from data buffers and store in sample buffers
	void retrieve_samples(uint32_t set_num, uint32_t smpl_rate_div);

	//Apply a digital filter to sample buffers
	void filter_samples(AntiAliasingFilter& filter_0, AntiAliasingFilter& filter_1);

	//Overwrite selected packet in data buffers with the contents of sample buffers
	void overwrite_samples(uint32_t packet_num);

	//Set all values in data buffers to 0
	void clear_data_set();

	//Set all values in output buffers to 0
	void clear_output_buffer();

	//setup processing functions based on the specified parameters
	void setup(bool cross_corr, double smpl_rate, uint32_t num_meas, uint32_t data_size, uint32_t min_window_size, uint32_t num_bands, uint32_t band_mult, uint32_t band_shift, double conversion_fac);
	void setup(Config new_config);

	//return vector containing x-axis frequency output values
	const std::vector<double>& get_x_values();

	//return vector containing y-axis power spectral density output values
	const std::vector<double>& get_y_values();

	//Write output data points to file
	bool write_output_to_file(char* filename);

	//Process input data from file to produce a set of output data points
	//Each frequency band uses a fixed window size with variable sampling rates (faster)
	void process_fixed_window_size();

	//Process input data from file to produce a set of output data points
	//Each frequency band uses a fixed sampling with variable window sizes (more accurate)
	void process_fixed_smpl_rate();

	//Start fixed window size processing where data is pushed in real time
	void start_fixed_window_size();

	//Push a data set to be processed with a fixed window size
	bool push_data_fixed_window_size(double* data_0, double* data_1);

	//Start fixed sample rate processing where data is pushed in real time
	void start_fixed_smpl_rate();

	//Push a data set to be processed with a fixed sample rate
	bool push_data_fixed_smpl_rate(double* data_0, double* data_1);

	//Stop real time processing and release data buffers
	void stop();
}


