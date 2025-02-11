#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <stdint.h>
#include <vector>

//! Generates a lookup table of precalculated sine/cosine values
class TrigTable
{
private:
	
	uint32_t mod_mask;
	uint32_t cos_shift;
	std::vector<double> sin_table;

public:
	TrigTable(uint32_t size)
	{
		mod_mask = size - 1;
		cos_shift = size >> 2;
		sin_table.resize(size);
		for (int i = 0; i < size; i++) sin_table[i] = sin((((double) i / size)) * 2 * M_PI);
	}
	TrigTable()
	{
		mod_mask = 0x0;
		cos_shift = 0;
		sin_table.resize(1);
		sin_table[0] = 0.0;
	}
	//! Look up sine value at the given position
	double sin_lookup(uint64_t pos)
	{
		uint32_t pos_m = (uint32_t)(pos & mod_mask);
		return sin_table[pos_m];
	}

	//! Look up cosine value at the given position
	double cos_lookup(uint64_t pos)
	{
		uint32_t pos_m = (uint32_t)((pos + cos_shift) & mod_mask);
		return sin_table[pos_m];
	}
};

 