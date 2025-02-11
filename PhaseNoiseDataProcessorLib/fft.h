#pragma once

#include <stdint.h>
#include "trig_table.h"

extern uint32_t g_data_size;
extern TrigTable g_trig;

//! Sets the size of the FFT window to be used and precalculates sine/cosine values

//! Must be set before any transform functions are performed.
void set_data_size(uint32_t size);

//! Calculate the DFT for an array of real and imaginary data points

//! Only output values in the range *k_start* to *k_stop* are calculated.
void dft(double* x_re, double* x_im, double* out_re, double* out_im, uint32_t k_start, uint32_t k_stop);

//! Calculate the DFT for an array of real data points

//! Only output values in the range *k_start* to *k_stop* are calculated.
void rdft(double* x_re, double* out_re, double* out_im, uint32_t k_start, uint32_t k_stop);

//! Calculate the FFT for an array of real and imaginary data points
void fft(double* x_re, double* x_im, double* out_re, double* out_im);
void fft(double* x_re, double* x_im, double* out_re, double* out_im, uint32_t data_size, uint32_t div);

//! Calculate the partial FFT for an array of real and imaginary data points

//! Only the first *k_num* output values are calculated.
void pfft(double* x_re, double* x_im, double* out_re, double* out_im, uint32_t k_num);
void pfft(double* x_re, double* x_im, double* out_re, double* out_im, uint32_t k_num, uint32_t data_size, uint32_t div);

//! Calculate the partial FFT for two arrays of real data points to be used as the input to a dual-real FFT cross-correlation function

//! //! Only the first *k_num* output values are calculated.
void pfft_dual_real(double* x_re1, double* x_re2, double* out_re, double* out_im, uint32_t k_num);
void pfft_dual_real(double* x_re1, double* x_re2, double* out_re, double* out_im, uint32_t k_num, uint32_t data_size, uint32_t div);

//! Calculate the FFT for an array of real data points
void rfft(double* x_re, double* out_re, double* out_im);
void rfft(double* x_re, double* out_re, double* out_im, uint32_t data_size, uint32_t div);

//! Calculate the partial FFT for an array of real data points

//! Only the first *k_num* output values are calculated.
void rpfft(double* x_re, double* out_re, double* out_im, uint32_t k_num);
void rpfft(double* x_re, double* out_re, double* out_im, uint32_t k_num, uint32_t data_size, uint32_t div);

//! Performs the cross correlation of complex data sets *f* and *g*

//! Only output values in the range *k_start* to *k_stop* are calculated.
void cross_corr(double* f_re, double* f_im, double* g_re, double* g_im, double* out_re, double* out_im, uint32_t k_start, uint32_t k_stop);
void cross_corr(double* f_re, double* f_im, double* g_re, double* g_im, double* out_re, double* out_im);

//! Performs the cross correlation of an FFT of two sets of real data

//! The input data is the FFT generated where one of the real data sets is used as the real part of the FFT input, and the other real data set is used as the imaginary part of the FFT input.
//! Only output values in the range *k_start* to *k_stop* are calculated.
//! Unlike the standard cross_corr() function, the output values are the sum of the positive and negative frequency components (i.e. multiplied by 2).
void cross_corr_dual_real(double* x_re, double* x_im, double* out_re, double* out_im, uint32_t k_start, uint32_t k_stop);
void cross_corr_dual_real(double* x_re, double* x_im, double* out_re, double* out_im);

//! Performs the auto correlation of the result of an FFT of real data

//! Only output values in the range *k_start* to *k_stop* are calculated.
//! The output values are the sum of the positive and negative frequency components (i.e. multiplied by 2).
void auto_corr_real(double* x_re, double* x_im, double* out_re, uint32_t k_start, uint32_t k_stop);
void auto_corr_real(double* x_re, double* x_im, double* out_re);
