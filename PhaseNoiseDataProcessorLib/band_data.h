#pragma once

#include <vector>

//! Data structure for a single phase noise measurement band
class Band
{
public:

	uint32_t num_sets = 0; //!< Number of separate data sets that a single measurement is split into (equal to the number of correlation/autocorrelations performed)
	uint16_t smpl_rate_div = 0; //!< Ratio of under-sampling relative to the sampling rate of the captured data
	uint32_t output_start_pos = 0; //!< Start position in the output data buffers of where output data for this band should be written
	uint32_t window_size = 0; //!< The sampling window size for a single data set
	uint32_t dft_size = 0; //!< Size of the output of the Fourier transform (only positive frequencies)
	uint32_t dft_start_pos = 0; //!< Start position of useful output data from Fourier transforms calculated for this band
	uint32_t dft_end_pos = 0; //!< End position of useful output data from Fourier transforms calculated for this band
};


//! Generates data for a set of phase noise measurement bands based on the supplied parameters
class BandData
{
public:

	uint32_t smpl_data_size = 0; //!< Number of samples per channel in a single measurement data capture (must be a power of 2)
	uint16_t band_mult = 1; //!< Frequency multiplier for subsequent bands (must be a power of 2)
	uint16_t num_bands = 0; //!< Number of frequency bands to use
	uint32_t min_window_size = 0; //!< FFT window size to use for the highest frequency band (must be a power of 2)
	uint16_t band_shift = 0; //!< Shifts each band down in multiples of band_mult by the specified number (prioritises more correlations over smaller RBW for higher frequencies)
	bool fixed_smpl_rate = true; //!< True if each band uses a fixed sample rate (more accurate), false if each band uses a fixed window size (faster)

	uint32_t output_size = 0; //!< Number of output data points that will be generated based on current settings
	std::vector<Band> bands; //!< Vector array of band data structures

	BandData() {};

	//! Instatiates the BandData object and generates the Bands
	BandData(uint16_t num_bands, uint16_t band_mult, uint32_t smpl_data_size, uint32_t min_window_size, uint16_t band_shift, bool fixed_smpl_rate)
	{
		this->smpl_data_size = smpl_data_size;
		this->band_mult = band_mult;
		this->num_bands = num_bands;
		this->min_window_size = min_window_size;
		this->band_shift = band_shift;
		this->fixed_smpl_rate = fixed_smpl_rate;
		

		uint16_t smpl_rate_div = 1;
		uint32_t num_sets = smpl_data_size / min_window_size;

		this->bands.resize(num_bands);

		for (int b = (this->num_bands - 1); b >= 0; b--)
		{
			this->bands[b].smpl_rate_div = smpl_rate_div;
			this->bands[b].num_sets = num_sets;

			if(!fixed_smpl_rate) smpl_rate_div *= this->band_mult;
			num_sets /= this->band_mult;
		}

		uint32_t pos_count = 0;
		uint16_t shift = 1;
		for (int i = 0; i < band_shift; i++) shift *= this->band_mult;

		for (int b = 0; b < this->num_bands; b++)
		{
			this->bands[b].window_size = fixed_smpl_rate ? (smpl_data_size / this->bands[b].num_sets) : min_window_size;
			this->bands[b].dft_size = (this->bands[b].window_size / 2) + 1;

			this->bands[b].output_start_pos = pos_count;
			this->bands[b].dft_start_pos = (b == 0) ? 0 : (min_window_size / (2 * band_mult * shift));
			this->bands[b].dft_end_pos = (b == (this->num_bands - 1)) ? (min_window_size / 2) : ((min_window_size / (2 * shift)) - 1);
			pos_count += (this->bands[b].dft_end_pos - this->bands[b].dft_start_pos + 1);
		}

		this->output_size = pos_count;
	}
};
