#define _USE_MATH_DEFINES
#include <cmath>

#include <string>
#include <vector>
#include <time.h>
#include "PhaseNoiseDataProcessorLib.h"

//paths to configuration and calibration files
string CONFIG_FILENAME = "config.txt";
string CALIB_FILENAME = "calib.txt";

//buffers for storing raw capture data
std::vector<int16_t> data_buffer_A;
std::vector<int16_t> data_buffer_B;

//buffers for storing converted capture data
std::vector<double> smpl_buffer_A;
std::vector<double> smpl_buffer_B;

//calibration values
bool calib_en = false; //perform calibration to measure dBc/Hz if true, measure absoulte value in dBV/Hz if false
double rf_input_level = 0.0; //RF signal level at DUT input in dBV
double cal_input_level = 0.0; //calibration tone level at DUT input in dBV
double cal_output_level = 0.0; //calibration tone level at anti-aliasing filter output in dBV

//read configuration data from config.txt file
bool readConfigFile()
{
	ifstream config_file(CONFIG_FILENAME);
	string line;
	string var_name;
	string var_value;
	size_t pos;

	PhaseNoiseDataProcessor::Config c;

	if (config_file.is_open())
	{
		while (getline(config_file, line))
		{
			pos = line.find("=");
			var_name = line.substr(0, pos);
			var_value = line.substr(pos + 1);

			if (var_name == "cross_corr")
			{
				if (var_value == "true") c.cross_corr = true;
				else c.cross_corr = false;
			}
			else if (var_name == "smpl_rate") c.smpl_rate = strtod(var_value.data(), NULL);
			else if (var_name == "num_meas") c.num_meas = strtoul(var_value.data(), NULL, 0);
			else if (var_name == "data_size") c.data_size = strtoul(var_value.data(), NULL, 0);
			else if (var_name == "min_window_size") c.min_window_size = strtoul(var_value.data(), NULL, 0);
			else if (var_name == "num_bands") c.num_bands = strtoul(var_value.data(), NULL, 0);
			else if (var_name == "band_mult") c.band_mult = strtoul(var_value.data(), NULL, 0);
			else if (var_name == "band_shift") c.band_shift = strtoul(var_value.data(), NULL, 0);
			else if (var_name == "window_en")
			{
				if (var_value == "true") c.window_en = true;
				else c.window_en = false;
			}
			else if (var_name == "calib_en")
			{
				if (var_value == "true") calib_en = true;
				else calib_en = false;
			}
		}

		c.conversion_fac = 1.0;

		PhaseNoiseDataProcessor::setup(c);

		return true;
	}
	else printf("Error: could not open config file\n\n");

	return false;
}

//read calibration data from calib.txt file
bool readCalibFile()
{
	ifstream calib_file(CALIB_FILENAME);
	string line;
	string var_name;
	string var_value;
	size_t pos;

	if (calib_file.is_open())
	{
		while (getline(calib_file, line))
		{
			pos = line.find("=");
			var_name = line.substr(0, pos);
			var_value = line.substr(pos + 1);

			if (var_name == "rf_input_level") rf_input_level = strtod(var_value.data(), NULL);
			else if (var_name == "cal_input_level") cal_input_level = strtod(var_value.data(), NULL);
			else if (var_name == "cal_output_level") cal_output_level = strtod(var_value.data(), NULL);
		}

		return true;
	}
	else printf("Error: could not open calibration file\n\n");

	return false;
}

void displayMeasurement(uint32_t meas)
{
	PhaseNoiseDataProcessor::calc_y_values();
	if (calib_en) PhaseNoiseDataProcessor::dbv_to_dbc(rf_input_level, cal_input_level, cal_output_level);

	system("cls");

	printf("Measurement %d / %d\n\n", meas, PhaseNoiseDataProcessor::config.num_meas);

	for (int b = 0; b < PhaseNoiseDataProcessor::band_data.num_bands; b++)
	{
		uint32_t index = (b == 0) ? 1 : PhaseNoiseDataProcessor::band_data.bands[b].output_start_pos;
		printf("\r                                                     \r");
		if (calib_en) printf("                - %f dBc/Hz\r%f Hz\n", PhaseNoiseDataProcessor::output_y_buffer[index], PhaseNoiseDataProcessor::output_x_buffer[index]);
		else printf("                - %f dBV/Hz\r%f Hz\n", PhaseNoiseDataProcessor::output_y_buffer[index], PhaseNoiseDataProcessor::output_x_buffer[index]);
	}
}

bool openDevice()
{
	//
	//	Code for opening communication with capture device should go here.
	//  Return true if successful.
	//
}

void setupDevice(uint32_t smpl_rate, uint32_t num_smpls)
{
	//
	// Code for setting up capture device should go here.
	//
}

void closeDevice()
{
	//
	//	Code for closing communication with capture device should go here.
	//
}

void startCapture(uint32_t num_smpls)
{
	//
	// This function should tell the capture device to begin capturing num_smpls samples to be stored in data_buffer_A (and data_buffer_B in cross-correlation mode)
	// Should be non-blocking.
	//
}

void waitForCapture()
{
	//
	// This function should block until the requested capture has finished.
	// In case the device does not transfer captured samples directly into data buffers, this function should also be responsible for retrieving data from the device
	//	and storing it into data_buffer_A and data_buffer_B
	//
}

//setup channels and data buffers to perform residual noise measurement
void setupMeasurement()
{
	uint32_t num_smpls = PhaseNoiseDataProcessor::config.data_size;

	data_buffer_A.resize(num_smpls);
	data_buffer_B.resize(num_smpls);
	smpl_buffer_A.resize(num_smpls);
	smpl_buffer_B.resize(num_smpls);

	setupDevice(PhaseNoiseDataProcessor::config.smpl_rate, num_smpls);
}

//perform measurement and data processing based on configuration data
void startMeasurement()
{
	uint32_t num_meas = PhaseNoiseDataProcessor::config.num_meas;
	uint32_t num_smpls = PhaseNoiseDataProcessor::config.data_size;

	PhaseNoiseDataProcessor::start_fixed_smpl_rate();

	PhaseNoiseDataProcessor::calc_x_values();

	startCapture(num_smpls);

	for (int meas = 0; meas < num_meas; meas++)
	{
		displayMeasurement(meas);

		waitForCapture();

		convertSamples(CAPTURE_RANGE_MV);

		startCapture(num_smpls);
		
		PhaseNoiseDataProcessor::push_data_fixed_smpl_rate(smpl_buffer_A.data(), smpl_buffer_B.data());

	}

	displayMeasurement(num_meas);

	printf("\n\n");

	PhaseNoiseDataProcessor::stop();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int32_t main(void)
{
	printf("Real-Time Phase Noise Measurement\n");

	printf("\n\nOpening the device...\n");

	if (!openDevice())
	{
		printf("Unable to open device\n");
		while (!_kbhit());
		exit(99);
	}

	printf("Device opened successfully.\n\n");

	time_t start_t, stop_t;
	double meas_time;

	if (!readConfigFile())
	{
		closeDevice();
		printf("Unable to read config file\n");
		while (!_kbhit());
		exit(99);
	}

	if (calib_en)
	{
		if (!readCalibFile())
		{
			closeDevice();
			printf("Unable to read calib file\n");
			while (!_kbhit());
			exit(99);
		}
	}

	time(&start_t);

	setupMeasurement();
	startMeasurement();

	time(&stop_t);

	int ch = ' ';

	while (true)
	{
		printf("Writing to file...\n\n");
		if (!PhaseNoiseDataProcessor::write_output_to_file(OUTPUT_FILENAME))
		{
			printf("Error: could not open output file\n\n");
			printf("Retry? (Y/N)\n\n");

			ch = toupper(_getch());
			if (ch != 'Y') break;
		}
		else
		{
			printf("Finished!\n\n");
			break;
		}
	}

	meas_time = difftime(stop_t, start_t);
	int min = (int)(meas_time / 60.0);
	int sec = (int)meas_time % 60;

	printf("Measurement time: %d min %d sec\n\n", min, sec);

	printf("Press X to exit...");

	while (true)
	{
		ch = toupper(_getch());
		if (ch == 'X') break;
	}

	closeDevice();

	return 0;
}
