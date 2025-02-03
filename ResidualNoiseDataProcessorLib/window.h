#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>

namespace Window
{
	class Window
	{
	protected:

		size_t size = 0;
		double rbw_factor;
		std::vector<double> data;

		virtual double rbwFactor() = 0;
		virtual void generate() = 0;

	public:

		Window(size_t size)
		{
			this->size = size;
			data.resize(size);
			this->rbw_factor = rbwFactor();
			generate();
		}
	};

	class FlatTop : Window
	{
	protected:

		double rbwFactor() override { return 3.8193596; };

		void generate()
		{
			double a0 = 0.21557895;
			double a1 = 0.41663158;
			double a2 = 0.277263158;
			double a3 = 0.083578947;
			double a4 = 0.006947368;
			double fac = M_PI / size;

			for (int i = 0; i < data.size(); i++)
			{
				data[i] = a0 - (a1 * cos(2 * i * fac)) + (a2 * cos(4 * i * fac)) - (a3 * cos(6 * i * fac)) + (a4 * cos(8 * i * fac));
			}
		}

	public:

		FlatTop(size_t size) : Window(size) {};
	};
}
