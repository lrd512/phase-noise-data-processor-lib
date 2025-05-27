#define _USE_MATH_DEFINES
#include <cmath>

#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <string>

#include <PhaseNoiseDataProcessorLib.h>

using namespace std;

#define MAX_DATA_SIZE 268435456
#define MAX_WINDOW_SIZE 268435456
#define MAX_MEAS 65536
#define OUTPUT_FILENAME "output.csv"

static bool cross_corr_en = true; //True for cross correlation, false for single channel
static double smpl_rate = 524288.0; //Sample rate of the captured data
static uint32_t num_meas = 1; //Number of measurements to process
static uint32_t data_size = 4194304; //Number of samples captured per channel for a single measurement (must be a power of 2)
static uint32_t min_window_size = 524288; //Minimum size of the FFT window (must be a power of 2)
static uint32_t num_bands = 1; //Number of frequency bands to use
static uint32_t band_mult = 8; //Frequency multiplier for subsequent bands (must be a power of 2)
static uint32_t band_shift = 2; //Shifts each band down in multiples of band_mult by the specified number (prioritises more correlations over smaller RBW for higher frequencies)
static double conversion_fac = 0.0025; //Maximum input voltage reference (bipolar)

//Parse input arguments to obtain measurement settings
bool parseArguments(char* argv[])
{
	//Measurement mode
	string mode_s = argv[1];
	try
	{
		if ((mode_s == "c") || (mode_s == "C")) cross_corr_en = true;
		else if ((mode_s == "s") || (mode_s == "S")) cross_corr_en = false;
		else throw(mode_s);
	}
	catch (string arg)
	{
		cout << "Invalid Mode argument: " << arg;
		return false;
	}

	//Sampling rate
	const char* smpl_rate_s = argv[2];
	try
	{
		if (smpl_rate_s == nullptr) throw((string) "");
		
		double val = strtod(smpl_rate_s, NULL);
		if (val == 0.0) throw((string)smpl_rate_s);
		else smpl_rate = val;
	}
	catch (string arg)
	{
		cout << "SmplRate argument out of range: " << arg;
		return false;
	}

	//Number of measurements
	const char* num_meas_s = argv[3];
	try
	{
		if (num_meas_s == nullptr) throw((string) "");
		
		unsigned long val = strtoul(num_meas_s, NULL, 0);
		if ((val <= 0) || (val > MAX_MEAS)) throw((string)num_meas_s);
		else num_meas = val;
	}
	catch (string arg)
	{
		cout << "NumMeas argument out of range: " << arg;
		return false;
	}

	//Data size
	const char* data_size_s = argv[4];
	try
	{
		if (data_size_s == nullptr) throw((string) "");
		
		unsigned long val = strtoul(data_size_s, NULL, 0);
		if ((val <= 0) || (val > MAX_DATA_SIZE)) throw((string)data_size_s);
		else data_size = val;
	}
	catch (string arg)
	{
		cout << "DataSize argument out of range: " << arg;
		return false;
	}

	//Minimum window size
	const char* min_window_size_s = argv[5];
	try
	{
		if (min_window_size_s == nullptr) throw((string) "");
		
		unsigned long val = strtoul(min_window_size_s, NULL, 0);
		if ((val <= 0) || (val > MAX_WINDOW_SIZE) || (val > data_size)) throw((string)min_window_size_s);
		else min_window_size = val;
	}
	catch (string arg)
	{
		cout << "WindowSize argument out of range: " << arg;
		return false;
	}

	//Number of bands
	const char* num_bands_s = argv[6];
	try
	{
		if (num_bands_s == nullptr) throw((string) "");
		
		unsigned long val = strtoul(num_bands_s, NULL, 0);
		if ((val <= 0) || (val > UINT16_MAX)) throw((string)num_bands_s);
		else num_bands = val;
	}
	catch (string arg)
	{
		cout << "NumBands argument out of range: " << arg;
		return false;
	}

	//Band multiple
	const char* band_mult_s = argv[7];
	try
	{
		if (band_mult_s == nullptr) throw((string) "");
		
		unsigned long val = strtoul(band_mult_s, NULL, 0);
		if ((val <= 0) || (val > UINT16_MAX)) throw((string)band_mult_s);
		else band_mult = val;
	}
	catch (string arg)
	{
		cout << "BandMult argument out of range: " << arg;
		return false;
	}

	//Band shift
	const char* band_shift_s = argv[8];
	try
	{
		if (band_shift_s == nullptr) throw((string) "");
		
		unsigned long val = strtoul(band_shift_s, NULL, 0);
		if ((val < 0) || (val > UINT16_MAX)) throw((string)band_shift_s);
		else band_shift = val;
	}
	catch (string arg)
	{
		cout << "BandShift argument out of range: " << arg;
		return false;
	}

	//Conversion factor
	const char* conversion_fac_s = argv[9];
	try
	{
		if (conversion_fac_s == nullptr) throw((string) "");
		
		double val = strtod(conversion_fac_s, NULL);
		if (val == 0.0) throw((string)conversion_fac_s);
		else conversion_fac = val;
	}
	catch (string arg)
	{
		cout << "Invalid ConvFac argument: " << arg;
		return false;
	}

	return true;
}

int main(int argc, char* argv[])
{
	//Parse the command line arguments and start the progrram if successful
	if (!parseArguments(argv)) return 0;
	cout << "Started.\n\r\n\r";

	//Set up processor according to arguments
	PhaseNoiseDataProcessor::setup(cross_corr_en, smpl_rate, num_meas, data_size, min_window_size, num_bands, band_mult, band_shift, conversion_fac);

	//Read data samples from 'samples.bin' and process
	PhaseNoiseDataProcessor::process_fixed_smpl_rate();

	//Write the processed measurement to 'output.csv'
	cout << "Writing to file...\n\r\n\r";
	if (!PhaseNoiseDataProcessor::write_output_to_file(OUTPUT_FILENAME)) cout << "Error: could not open output file";
	else cout << "Finished!";

	return 0;
}
