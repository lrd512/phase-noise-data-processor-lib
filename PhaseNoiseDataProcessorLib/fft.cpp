#include "fft.h"

uint32_t g_data_size = 1;
TrigTable g_trig = TrigTable();

void set_data_size(uint32_t size)
{
	g_data_size = size;
	g_trig = TrigTable(size);
}

void dft(double* x_re, double* x_im, double* out_re, double* out_im, uint32_t k_start, uint32_t k_stop)
{
	uint32_t pos;
	double cos_fac, sin_fac;

	for(uint32_t k = k_start; k <= k_stop; k++)
	{
		pos = 0;

		for(uint32_t n = 0; n < g_data_size; n++)
		{
			cos_fac = g_trig.cos_lookup(pos);
			sin_fac = g_trig.sin_lookup(pos);

			out_re[k] += (x_re[n] * cos_fac) - (x_im[n] * sin_fac);
			out_im[k] += (x_re[n] * sin_fac) + (x_im[n] * cos_fac);

			pos += k;
		}
	}
}

void rdft(double* x_re, double* out_re, double* out_im, uint32_t k_start, uint32_t k_stop)
{
	uint32_t pos;
	double cos_fac, sin_fac;

	for (uint32_t k = k_start; k <= k_stop; k++)
	{
		pos = 0;

		for (uint32_t n = 0; n < g_data_size; n++)
		{
			cos_fac = g_trig.cos_lookup(pos);
			sin_fac = g_trig.sin_lookup(pos);

			out_re[k] += (x_re[n] * cos_fac);
			out_im[k] += (x_re[n] * sin_fac);

			pos += k;
		}
	}
}

void fft(double* x_re, double* x_im, double* out_re, double* out_im, uint32_t data_size, uint32_t div)
{
	if (data_size == 1)
	{
		out_re[0] = x_re[0];
		out_im[0] = x_im[0];
	}
	else
	{
		uint32_t half_data_size = data_size >> 1;
		uint32_t div_2 = div << 1;
		fft(x_re, x_im, out_re, out_im, half_data_size, div_2);
		fft(x_re + div, x_im + div, out_re + half_data_size, out_im + half_data_size, half_data_size, div_2);

		double p_re, p_im, q_re, q_im;
		double cos_fac, sin_fac;
		uint32_t pos = 0;
		uint32_t k2 = half_data_size;

		for (uint32_t k = 0; k < half_data_size; k++)
		{
			cos_fac = g_trig.cos_lookup(pos);
			sin_fac = g_trig.sin_lookup(pos);

			p_re = out_re[k];
			p_im = out_im[k];

			q_re = (cos_fac * out_re[k2]) + (sin_fac * out_im[k2]);
			q_im = (cos_fac * out_im[k2]) - (sin_fac * out_re[k2]);

			out_re[k] = p_re + q_re;
			out_re[k2] = p_re - q_re;
			out_im[k] = p_im + q_im;
			out_im[k2] = p_im - q_im;

			pos += div;
			k2++;
		}
	}
}
void fft(double* x_re, double* x_im, double* out_re, double* out_im) { fft(x_re, x_im, out_re, out_im, g_data_size, 1); };

void pfft(double* x_re, double* x_im, double* out_re, double* out_im, uint32_t k_num, uint32_t data_size, uint32_t div)
{
	uint32_t half_data_size = data_size >> 1;
	uint32_t div_2 = div << 1;

	if (k_num > (half_data_size >> 1))
	{
		fft(x_re, x_im, out_re, out_im, half_data_size, div_2);
		fft(x_re + div, x_im + div, out_re + half_data_size, out_im + half_data_size, half_data_size, div_2);
	}
	else
	{
		pfft(x_re, x_im, out_re, out_im, k_num, half_data_size, div_2);
		pfft(x_re + div, x_im + div, out_re + half_data_size, out_im + half_data_size, k_num, half_data_size, div_2);
	}

	double p_re, p_im, q_re, q_im;
	double cos_fac, sin_fac;
	uint32_t pos = 0;
	uint32_t k2 = half_data_size;

	for (uint32_t k = 0; k < k_num; k++)
	{
		cos_fac = g_trig.cos_lookup(pos);
		sin_fac = g_trig.sin_lookup(pos);

		p_re = out_re[k];
		p_im = out_im[k];

		q_re = (cos_fac * out_re[k2]) + (sin_fac * out_im[k2]);
		q_im = (cos_fac * out_im[k2]) - (sin_fac * out_re[k2]);

		out_re[k] = p_re + q_re;
		out_im[k] = p_im + q_im;

		pos += div;
		k2++;
	}
}
void pfft(double* x_re, double* x_im, double* out_re, double* out_im, uint32_t k_num)
{
	if(k_num >= (g_data_size >> 1)) fft(x_re, x_im, out_re, out_im, g_data_size, 1);
	else pfft(x_re, x_im, out_re, out_im, k_num, g_data_size, 1);
}

void pfft_dual_real(double* x_re, double* x_im, double* out_re, double* out_im, uint32_t k_num, uint32_t data_size, uint32_t div)
{
	uint32_t half_data_size = data_size >> 1;
	uint32_t div_2 = div << 1;

	if (k_num > (half_data_size >> 2))
	{
		fft(x_re, x_im, out_re, out_im, half_data_size, div_2);
		fft(x_re + div, x_im + div, out_re + half_data_size, out_im + half_data_size, half_data_size, div_2);
	}
	else
	{
		pfft_dual_real(x_re, x_im, out_re, out_im, k_num, half_data_size, div_2);
		pfft_dual_real(x_re + div, x_im + div, out_re + half_data_size, out_im + half_data_size, k_num, half_data_size, div_2);
	}

	double p1_re, p1_im, q1_re, q1_im, p2_re, p2_im, q2_re, q2_im;
	uint32_t pos = 0;
	uint32_t pos2 = half_data_size - div;
	uint32_t k2 = half_data_size;
	uint32_t k3 = data_size - 1;
	uint32_t k4 = half_data_size - 1;
	double cos_fac = 1;
	double sin_fac = 0;

	for (uint32_t k = 0; k < k_num; k++)
	{
		p1_re = out_re[k];
		p1_im = out_im[k];

		q1_re = (cos_fac * out_re[k2]) + (sin_fac * out_im[k2]);
		q1_im = (cos_fac * out_im[k2]) - (sin_fac * out_re[k2]);

		out_re[k] = p1_re + q1_re;
		out_im[k] = p1_im + q1_im;

		pos += div;
		cos_fac = g_trig.cos_lookup(pos);
		sin_fac = g_trig.sin_lookup(pos);

		p2_re = out_re[k4];
		p2_im = out_im[k4];

		q2_re = (cos_fac * out_re[k3]) - (sin_fac * out_im[k3]);
		q2_im = (cos_fac * out_im[k3]) + (sin_fac * out_re[k3]);

		out_re[k3] = p2_re + q2_re;
		out_im[k3] = p2_im + q2_im;

		k2++;
		k3--;
		k4--;
	}
}
void pfft_dual_real(double* x_re, double* x_im, double* out_re, double* out_im, uint32_t k_num)
{
	if (k_num >= (g_data_size >> 1)) fft(x_re, x_im, out_re, out_im, g_data_size, 1);
	else pfft_dual_real(x_re, x_im, out_re, out_im, k_num, g_data_size, 1);
}

void rfft(double* x_re, double* out_re, double* out_im, uint32_t data_size, uint32_t div)
{
	if (data_size == 1)
	{
		out_re[0] = x_re[0];
		out_im[0] = 0;
	}
	else
	{
		uint32_t half_data_size = data_size >> 1;
		uint32_t div_2 = div << 1;
		rfft(x_re, out_re, out_im, half_data_size, div_2);
		rfft(x_re + div, out_re + half_data_size, out_im + half_data_size, half_data_size, div_2);

		double p_re, p_im, q_re, q_im;
		double cos_fac, sin_fac;
		uint32_t pos = 0;
		uint32_t k2 = half_data_size;

		for (uint32_t k = 0; k < half_data_size; k++)
		{
			cos_fac = g_trig.cos_lookup(pos);
			sin_fac = g_trig.sin_lookup(pos);

			p_re = out_re[k];
			p_im = out_im[k];

			q_re = (cos_fac * out_re[k2]) + (sin_fac * out_im[k2]);
			q_im = (cos_fac * out_im[k2]) - (sin_fac * out_re[k2]);

			out_re[k] = p_re + q_re;
			out_re[k2] = p_re - q_re;
			out_im[k] = p_im + q_im;
			out_im[k2] = p_im - q_im;

			pos += div;
			k2++;
		}
	}
}
void rfft(double* x_re, double* out_re, double* out_im) { rfft(x_re, out_re, out_im, g_data_size, 1); };

void rpfft(double* x_re, double* out_re, double* out_im, uint32_t k_num, uint32_t data_size, uint32_t div)
{
	uint32_t half_data_size = data_size >> 1;
	uint32_t div_2 = div << 1;

	if (k_num > (half_data_size >> 1))
	{
		rfft(x_re, out_re, out_im, half_data_size, div_2);
		rfft(x_re + div, out_re + half_data_size, out_im + half_data_size, half_data_size, div_2);
	}
	else
	{
		rpfft(x_re,  out_re, out_im, k_num, half_data_size, div_2);
		rpfft(x_re + div, out_re + half_data_size, out_im + half_data_size, k_num, half_data_size, div_2);
	}

	double p_re, p_im, q_re, q_im;
	double cos_fac, sin_fac;
	uint32_t pos = 0;
	uint32_t k2 = half_data_size;

	for (uint32_t k = 0; k < k_num; k++)
	{
		cos_fac = g_trig.cos_lookup(pos);
		sin_fac = g_trig.sin_lookup(pos);

		p_re = out_re[k];
		p_im = out_im[k];

		q_re = (cos_fac * out_re[k2]) + (sin_fac * out_im[k2]);
		q_im = (cos_fac * out_im[k2]) - (sin_fac * out_re[k2]);

		out_re[k] = p_re + q_re;
		out_im[k] = p_im + q_im;

		pos += div;
		k2++;
	}
}
void rpfft(double* x_re, double* out_re, double* out_im, uint32_t k_num)
{
	if (k_num >= (g_data_size >> 1)) rfft(x_re, out_re, out_im, g_data_size, 1);
	else rpfft(x_re, out_re, out_im, k_num, g_data_size, 1);
}

void cross_corr_dual_real(double* x_re, double* x_im, double* out_re, double* out_im, uint32_t k_start, uint32_t k_stop)
{
	double f_re, f_im, gc_re, gc_im;
	uint32_t output_pos = 0;
	uint32_t k_start_t = k_start;

	if (k_start_t == 0)
	{
		f_re = x_re[0];
		gc_re = x_im[0];

		out_re[0] += f_re * gc_re;
		k_start_t = 1;
		output_pos++;
	}

	for (uint32_t k = k_start_t; k <= k_stop; k++)
	{
		f_re = x_re[k] + x_re[g_data_size - k];
		f_im = x_im[k] - x_im[g_data_size - k];
		gc_re = x_im[k] + x_im[g_data_size - k];
		gc_im = x_re[k] - x_re[g_data_size - k];

		out_re[output_pos] += (f_re * gc_re) - (f_im * gc_im);
		out_im[output_pos] += (f_re * gc_im) + (f_im * gc_re);

		output_pos++;
	}
}
void cross_corr_dual_real(double* x_re, double* x_im, double* out_re, double* out_im)
{
	cross_corr_dual_real(x_re, x_im, out_re, out_im, 0, g_data_size >> 1);
}

void cross_corr(double* f_re, double* f_im, double* g_re, double* g_im, double* out_re, double* out_im, uint32_t k_start, uint32_t k_stop)
{
	uint32_t output_pos = 0;

	for (uint32_t k = k_start; k <= k_stop; k++)
	{
		out_re[output_pos] += (f_re[k] * g_re[k]) + (f_im[k] * g_im[k]);
		out_im[output_pos] += (f_im[k] * g_re[k]) - (f_re[k] * g_im[k]);

		output_pos++;
	}
}
void cross_corr(double* f_re, double* f_im, double* g_re, double* g_im, double* out_re, double* out_im)
{
	cross_corr(f_re, f_im, g_re, g_im, out_re, out_im, 0, g_data_size - 1);
}

void auto_corr_real(double* x_re, double* x_im, double* out_re, uint32_t k_start, uint32_t k_stop)
{
	uint32_t output_pos = 0;
	uint32_t k_start_t = k_start;

	if (k_start_t == 0)
	{
		out_re[output_pos] += (x_re[0] * x_re[0]) + (x_im[0] * x_im[0]);
		k_start_t = 1;
		output_pos++;
	}

	for (uint32_t k = k_start_t; k <= k_stop; k++)
	{
		out_re[output_pos] += ((x_re[k] * x_re[k]) + (x_im[k] * x_im[k])) * 4;
		output_pos++;
	}
}

void auto_corr_real(double* x_re, double* x_im, double* out_re)
{
	auto_corr_real(x_re, x_im, out_re, 0, g_data_size >> 1);
}
