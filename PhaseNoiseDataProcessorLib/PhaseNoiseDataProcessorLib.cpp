#include "PhaseNoiseDataProcessorLib.h"

namespace PhaseNoiseDataProcessor
{
	std::vector<double> meas_buffer[2];
	std::vector<double> smpl_buffer[2];
	std::vector<double> fft_re_buffer;
	std::vector<double> fft_im_buffer;
	std::vector<double> output_re_buffer;
	std::vector<double> output_im_buffer;
	std::vector<double> output_x_buffer;
	std::vector<double> output_y_buffer;

	uint32_t meas_count = 0;

	BandData band_data;

	Config config;

	void convert_samples(int channel, const void* data_bytes)
	{
		double sum = 0.0;
		uint16_t* data_2bytes = (uint16_t*)data_bytes;

		for (int i = 0; i < config.data_size; i++)
		{
			meas_buffer[channel].at(i) = (((double)data_2bytes[i] + INT16_MIN) / -INT16_MIN) * config.conversion_fac;
		}
	}

	void calc_x_values()
	{
		int pos_cout = 0;
		double rbw;

		for (int b = 0; b < band_data.num_bands; b++)
		{
			rbw = config.smpl_rate / (band_data.bands[b].window_size * band_data.bands[b].smpl_rate_div);
			for (int i = band_data.bands[b].dft_start_pos; i <= band_data.bands[b].dft_end_pos; i++)
			{
				output_x_buffer.at(pos_cout) = i * rbw;
				pos_cout += 1;
			}
		}
	}

	void calc_y_values(uint32_t num_meas)
	{
		int pos_count = 0;
		double fac, mag, rbw;

		for (int b = 0; b < band_data.num_bands; b++)
		{
			//calculate resolution bandwidth of current band
			rbw = config.smpl_rate / (band_data.bands[b].window_size * band_data.bands[b].smpl_rate_div);

			for (int i = band_data.bands[b].dft_start_pos; i <= band_data.bands[b].dft_end_pos; i++)
			{
				if (config.cross_corr) mag = sqrt(pow(output_re_buffer[pos_count], 2) + pow(output_im_buffer[pos_count], 2));
				else mag = abs(output_re_buffer[pos_count]);
				mag /= ((double)num_meas * (double)band_data.bands[b].num_sets * (double)band_data.bands[b].window_size * (double)band_data.bands[b].window_size);
				output_y_buffer.at(pos_count) = 10 * log10(mag / (2 * rbw));
				pos_count += 1;
			}
		}
	}
	void calc_y_values()
	{
		calc_y_values(meas_count);
	}

	void dbv_to_dbc(double rf_input_level, double cal_input_level, double cal_output_level)
	{
		for (int i = 0; i < output_y_buffer.size(); i++) output_y_buffer[i] += ((cal_input_level - rf_input_level) - cal_output_level - 6.0);
	}

	bool read_meas_from_file(char* filename, uint32_t meas_num)
	{
		ifstream file(filename, ios::in | ios::binary);
		if (file.is_open())
		{
			streampos pos = (config.cross_corr) ? (meas_num * config.data_size * 4) : (meas_num * config.data_size * 2);
			file.seekg(pos);

			std::vector<char> data_bytes(config.data_size * 4);

			file.read(data_bytes.data(), config.data_size * 2);
			convert_samples(0, data_bytes.data());

			if (config.cross_corr)
			{
				file.read(data_bytes.data(), config.data_size * 2);
				convert_samples(1, data_bytes.data());
			}

			file.close();

			return true;
		}

		return false;
	}

	bool write_output_to_file(char* filename)
	{
		ofstream file(filename, ios::out | ios::trunc);
		if (file.is_open())
		{
			string line;

			for (int b = 0; b < band_data.num_bands; b++)
			{
				file << "Band," << to_string(b + 1) << "\n";
				file << "num_sets," << to_string(band_data.bands[b].num_sets) << "\n";
				file << "smpl_rate_div," << to_string(band_data.bands[b].smpl_rate_div) << "\n";
				file << "output_start_pos," << to_string(band_data.bands[b].output_start_pos) << "\n";
				file << "window_size," << to_string(band_data.bands[b].window_size) << "\n";
				file << "dft_size," << to_string(band_data.bands[b].dft_size) << "\n";
				file << "dft_start_pos," << to_string(band_data.bands[b].dft_start_pos) << "\n";
				file << "dft_end_pos," << to_string(band_data.bands[b].dft_end_pos) << "\n\n";
			}

			for (int i = 0; i < output_x_buffer.size(); i++)
			{
				line = to_string(output_x_buffer[i]) + ',' + to_string(output_y_buffer[i]) + ',' + "\n";
				file << line;
			}

			file.close();

			return true;
		}

		return false;
	}

	void retrieve_samples(uint32_t set_num, uint32_t smpl_rate_div)
	{
		uint32_t offset = smpl_buffer[0].size() * set_num * smpl_rate_div;
		double* data_buffer = meas_buffer[0].data() + offset;
		uint32_t pos = 0;

		for (int i = 0; i < smpl_buffer[0].size(); i++)
		{
			smpl_buffer[0].at(i) = data_buffer[pos];
			pos += smpl_rate_div;
		}

		if (config.cross_corr)
		{
			offset = smpl_buffer[1].size() * set_num * smpl_rate_div;
			data_buffer = meas_buffer[1].data() + offset;
			pos = 0;

			for (int i = 0; i < smpl_buffer[1].size(); i++)
			{
				smpl_buffer[1].at(i) = data_buffer[pos];
				pos += smpl_rate_div;
			}
		}
	}

	void filter_samples(AntiAliasingFilter& filter_0, AntiAliasingFilter& filter_1)
	{
		for (int i = 0; i < smpl_buffer[0].size(); i++) smpl_buffer[0].at(i) = filter_0.pass(smpl_buffer[0][i]);
		for (int i = 0; i < smpl_buffer[1].size(); i++) smpl_buffer[1].at(i) = filter_1.pass(smpl_buffer[1][i]);
	}

	void overwrite_samples(uint32_t packet_num)
	{
		uint32_t offset = smpl_buffer[0].size() * packet_num;
		double* data_buffer = meas_buffer[0].data() + offset;
		for (int i = 0; i < smpl_buffer[0].size(); i++) data_buffer[i] = smpl_buffer[0][i];

		if (config.cross_corr)
		{
			offset = smpl_buffer[1].size() * packet_num;
			data_buffer = meas_buffer[1].data() + offset;
			for (int i = 0; i < smpl_buffer[1].size(); i++) data_buffer[i] = smpl_buffer[1][i];
		}
	}

	void clear_meas_buffer()
	{
		for (int i = 0; i < meas_buffer[0].size(); i++) meas_buffer[0].at(i) = 0.0;
		for (int i = 0; i < meas_buffer[1].size(); i++) meas_buffer[1].at(i) = 0.0;
	}

	void clear_output_buffer()
	{
		for (int i = 0; i < output_re_buffer.size(); i++) output_re_buffer.at(i) = 0.0;
		for (int i = 0; i < output_im_buffer.size(); i++) output_im_buffer.at(i) = 0.0;
	}

	void process_fixed_window_size()
	{
		start_fixed_window_size();

		for (int m = 0; m < config.num_meas; m++)
		{
			cout << "Measurement " << to_string(m + 1) << ": ";

			if (!read_meas_from_file(SMPL_FILENAME, m)) cout << "Error: could not open input file";
			if (config.cross_corr) push_data_fixed_window_size(meas_buffer[0].data(), meas_buffer[1].data());
			else push_data_fixed_window_size(meas_buffer[0].data(), nullptr);

			cout << "done\n\r";
		}

		stop();

		meas_count = config.num_meas;

		calc_x_values();
		calc_y_values();
	}

	void setup(bool cross_corr, double smpl_rate, uint32_t num_meas, uint32_t data_size, uint32_t min_window_size, uint32_t num_bands, uint32_t band_mult, uint32_t band_shift, double conversion_fac)
	{
		config.cross_corr = cross_corr;
		config.smpl_rate = smpl_rate;
		config.num_meas = num_meas;
		config.data_size = data_size;
		config.packet_size = min_window_size;
		config.num_bands = num_bands;
		config.band_mult = band_mult;
		config.band_shift = band_shift;
		config.conversion_fac = conversion_fac;
	}

	void setup(Config new_config)
	{
		config = new_config;
	}

	const std::vector<double>& get_x_values()
	{
		return output_x_buffer;
	}

	const std::vector<double>& get_y_values()
	{
		return output_y_buffer;
	}

	void process_fixed_smpl_rate()
	{
		meas_buffer[0].resize(config.data_size);
		if (config.cross_corr) meas_buffer[1].resize(config.data_size);

		start_fixed_smpl_rate();

		for (int m = 0; m < config.num_meas; m++)
		{
			cout << "Measurement " << to_string(m + 1) << ": ";

			if (!read_meas_from_file(SMPL_FILENAME, m)) cout << "Error: could not open input file";
			if (config.cross_corr) push_data_fixed_smpl_rate(meas_buffer[0].data(), meas_buffer[1].data());
			else push_data_fixed_smpl_rate(meas_buffer[0].data(), nullptr);

			cout << "done\n\r";
		}

		stop();

		meas_count = config.num_meas;

		calc_x_values();
		calc_y_values();
	}

	void start_fixed_window_size()
	{
		band_data = BandData(config.num_bands, config.band_mult, config.data_size, config.packet_size, config.band_shift, false);

		meas_buffer[0].resize(config.data_size);
		if (config.cross_corr) meas_buffer[1].resize(config.data_size);
		smpl_buffer[0].resize(config.packet_size);
		if (config.cross_corr) smpl_buffer[1].resize(config.packet_size);
		fft_re_buffer.resize(config.packet_size);
		fft_im_buffer.resize(config.packet_size);
		output_re_buffer.resize(band_data.output_size);
		if (config.cross_corr) output_im_buffer.resize(band_data.output_size);

		output_x_buffer.resize(band_data.output_size);
		output_y_buffer.resize(band_data.output_size);

		set_data_size(config.packet_size);

		clear_output_buffer();

		meas_count = 0;
	}

	bool push_data_fixed_window_size(double* data_0, double* data_1)
	{
		for (int i = 0; i < meas_buffer[0].size(); i++) meas_buffer[0].at(i) = data_0[i];
		for (int i = 0; i < meas_buffer[1].size(); i++) meas_buffer[1].at(i) = data_1[i];

		for (int b = (band_data.num_bands - 1); b >= 0; b--)
		{
			AntiAliasingFilter filter_0(config.band_mult);
			AntiAliasingFilter filter_1(config.band_mult);

			for (int s = 0; s < band_data.bands[b].num_sets; s++)
			{
				if (b == (band_data.num_bands - 1)) retrieve_samples(s, band_data.bands[b].smpl_rate_div);
				else retrieve_samples(s, band_data.band_mult);

				if (config.cross_corr)
				{
					pfft_dual_real(smpl_buffer[0].data(), smpl_buffer[1].data(), fft_re_buffer.data(), fft_im_buffer.data(), band_data.bands[b].dft_end_pos + 1);
					cross_corr_dual_real(fft_re_buffer.data(), fft_im_buffer.data(), output_re_buffer.data() + band_data.bands[b].output_start_pos,
						output_im_buffer.data() + band_data.bands[b].output_start_pos, band_data.bands[b].dft_start_pos, band_data.bands[b].dft_end_pos);
				}
				else
				{
					rpfft(smpl_buffer[0].data(), fft_re_buffer.data(), fft_im_buffer.data(), band_data.bands[b].dft_end_pos + 1);
					auto_corr_real(fft_re_buffer.data(), fft_im_buffer.data(), output_re_buffer.data() + band_data.bands[b].output_start_pos,
						band_data.bands[b].dft_start_pos, band_data.bands[b].dft_end_pos);
				}

				filter_samples(filter_0, filter_1);
				overwrite_samples(s);
			}
		}

		meas_count++;
		if (meas_count >= config.num_meas) return true;
		return false;
	}

	void start_fixed_smpl_rate()
	{
		band_data = BandData(config.num_bands, config.band_mult, config.data_size, config.packet_size, config.band_shift, true);

		fft_re_buffer.resize(config.data_size);
		fft_im_buffer.resize(config.data_size);
		output_re_buffer.resize(band_data.output_size);
		if (config.cross_corr) output_im_buffer.resize(band_data.output_size);

		output_x_buffer.resize(band_data.output_size);
		output_y_buffer.resize(band_data.output_size);

		clear_output_buffer();

		meas_count = 0;
	}

	bool push_data_fixed_smpl_rate(double* data_0, double* data_1)
	{
		for (int b = (band_data.num_bands - 1); b >= 0; b--)
		{
			set_data_size(band_data.bands[b].window_size);

			for (int s = 0; s < band_data.bands[b].num_sets; s++)
			{
				uint32_t data_offset = s * band_data.bands[b].window_size;

				if (config.cross_corr)
				{
					pfft_dual_real(data_0 + data_offset, data_1 + data_offset, fft_re_buffer.data(), fft_im_buffer.data(), band_data.bands[b].dft_end_pos + 1);
					cross_corr_dual_real(fft_re_buffer.data(), fft_im_buffer.data(), output_re_buffer.data() + band_data.bands[b].output_start_pos,
						output_im_buffer.data() + band_data.bands[b].output_start_pos, band_data.bands[b].dft_start_pos, band_data.bands[b].dft_end_pos);
				}
				else
				{
					rpfft(data_0 + data_offset, fft_re_buffer.data(), fft_im_buffer.data(), band_data.bands[b].dft_end_pos + 1);
					auto_corr_real(fft_re_buffer.data(), fft_im_buffer.data(), output_re_buffer.data() + band_data.bands[b].output_start_pos,
						band_data.bands[b].dft_start_pos, band_data.bands[b].dft_end_pos);
				}
			}
		}

		meas_count++;
		if(meas_count >= config.num_meas) return true;
		return false;
	}

	void stop()
	{
		meas_buffer[0].resize(0);
		if (config.cross_corr) meas_buffer[1].resize(0);
		smpl_buffer[0].resize(0);
		if (config.cross_corr) smpl_buffer[1].resize(0);
		fft_re_buffer.resize(0);
		fft_im_buffer.resize(0);
	}
}

