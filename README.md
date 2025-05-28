# phase-noise-data-processor
This is a C++ library of efficient digital signal processing algorithms for phase noise measurements.

## Overview

This library contains a number of functions designed to enable the efficient calculation of absolute/residual phase noise from a set (or sets) of time-domain voltage samples. Phase noise can be calculated from either single-channel input data or dual-channel cross-correlated data for improved measurement noise floor.

Specifically, this library aims to reduce the time taken to obtain useful phase noise measurements of ultra-low phase noise devices where a large number of cross-correlated measurements may be required to achieve an acceptable measurement noise floor. This is done by splitting the output data into separate frequency bands each with a different ratio of resolution bandwidth to number of correlations in a given time period. This enables measurements with high internal noise floor suppression at far-out frequencies whilst still maintaining a narrow resolution bandwidth close-to-carrier.

## Using the Library

This is provided as a ready-built static C++ library (.lib file) and can be found in the Build directory.

Alternatively, you can build the library yourself either directly from the source code or using the included Visual Studio project file.

## Documentation

Documentation for the library is generated using Doxygen (https://www.doxygen.nl/) and a reference manual containing descriptions of available functions and classes can be found in the Documentation directory.

This library was originally produced as part of my PhD Thesis "Ultra-Low Phase Noise 100 MHz Crystal Oscillator & Residual Phase Noise Measurement Systems" (https://etheses.whiterose.ac.uk/id/eprint/36102/) and further information on the motivation for this work, the principles behind the techniques used, including detailed descriptions of algorithms, can be found in Chapter 3. 

## Example Code

Two example programs showing how the library can be used to calculate phase noise are provided in the Examples directory.

"residual_noise_data_processor_app.cpp" - This is a simple command line interface for generating phase noise data from a set of time-domain input samples. The input should stored as consecutive 16-bit offset binary in a "samples.bin" file to produce a CSV of power spectral density versus offset frequency.

"real_time_measurement.cpp" - This demonstrates how the library can be used to perform phase noise measurements in real-time using an appropriate capture device.
