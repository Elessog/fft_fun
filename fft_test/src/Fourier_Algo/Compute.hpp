#ifndef COMPUTE_HPP
#define COMPUTE_HPP

#include "Utilities/Parser.hpp"

class ComputeFourier
{
public:
  ComputeFourier(const signal_struct &signal,double frequency);
  ~ComputeFourier(){};

  void applyFourier();
private:
  signal_struct signal_input;
  signal_struct signal_output;
  int nn;
  double sampling_frequency;

  void select_fftnr();
  void select_fftct();
  void select_fftcto();
  void fftshift(std::vector<double> &data);
};

#endif
