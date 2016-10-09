#ifndef FOURIER_HPP
#define FOURIER_HPP


#include <complex>
#include <vector>

typedef std::vector<std::complex<double>> signal_struct;
void fourier_nr(const signal_struct &data_input,signal_struct &data_output ,int nn, int isign);
void fourier_ct(const signal_struct &data_input,signal_struct &data_output);
void fourier_cto(const signal_struct &data_input,signal_struct &data_output);

#endif
