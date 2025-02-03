#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <stdint.h>

//Digital biquadratic filter
class BiquadFilter
{
    protected:
        
        double input_coeff[3] = { 0, 0, 0 };
        double feedback_coeff[2] = { 0, 0 };

        double input_buffer[3] = { 0, 0, 0 };
		double feedback_buffer[2] = { 0, 0 };
		
		double output_value = 0;

    public:

        BiquadFilter() {};
        BiquadFilter(int cutoff_div, double q)
        {
            generateCoefficients(cutoff_div, q);
        }

        void generateCoefficients(int cutoff_div, double q)
        {
			double k_wa = 1.0 / tanf(M_PI / cutoff_div);
			double k_wa_sqr = powf(k_wa, 2);

			double fac = 1.0 / (k_wa_sqr + (k_wa / q) + 1.0);

            input_coeff[0] = fac;
            input_coeff[1] = 2 * fac;
            input_coeff[2] = fac;

            feedback_coeff[0] = (2.0 * k_wa_sqr - 2.0) * fac;
            feedback_coeff[1] = -(k_wa_sqr - (k_wa / q) + 1.0) * fac;
        }
		
		double pass(double input)
		{
			input_buffer[2] = input_buffer[1];
			input_buffer[1] = input_buffer[0];
			input_buffer[0] = input;
			
			feedback_buffer[1] = feedback_buffer[0];
			feedback_buffer[0] = output_value;
			output_value = (input_buffer[0] * input_coeff[0]) + (input_buffer[1] * input_coeff[1]) + (input_buffer[2] * input_coeff[2]) +
									(feedback_buffer[0] * feedback_coeff[0]) + (feedback_buffer[1] * feedback_coeff[1]);
									
			return output();
		}
		
		void reset()
		{
			input_buffer[2] = 0;
			input_buffer[1] = 0;
			input_buffer[0] = 0;
			
			feedback_buffer[1] = 0;
			feedback_buffer[0] = 0;
			
			output_value = 0;
		}
		
		double output() { return output_value; };
};

//6th order Butterworth digital filter
class AntiAliasingFilter
{
    protected:
        const double Q1 = 1.931851989228;
		const double Q2 = 0.707106562373;
		const double Q3 = 0.517637997114;

        BiquadFilter filt1;
        BiquadFilter filt2;
        BiquadFilter filt3;

		double output_value = 0;

    public:
        
        AntiAliasingFilter(int cutoff_div)
        {
            filt1.generateCoefficients(cutoff_div, Q1);
            filt2.generateCoefficients(cutoff_div, Q2);
            filt3.generateCoefficients(cutoff_div, Q3);
        }

        float pass(double input)
        {
            output_value = filt3.pass(filt2.pass(filt1.pass(input)));
            return output();
        }
		
		void reset()
		{
			filt1.reset();
			filt2.reset();
			filt3.reset();
			
			output_value = 0;
		}

		double output() { return output_value; };
};
