#ifndef SIGNALS_HPP
#define SIGNALS_HPP

#include <vector>
#include <string>

class Signal
{
public:
  Signal(double duration,double frequency,char *filename);
  ~Signal(){};


  void create_signal();
  void save_file();

private:
  std::vector<double> signal_data;
  double sampling_frequency;
  double duration;
  std::string output_name;
  int n;

  void add_cos(double frequency, double amplitude, double phase = 0);
  void add_sin(double frequency, double amplitude, double phase = 0);
  void add_rectangle(double duration, double amplitude, double time_middle);
  void select_sin();
  void select_rectangle();
  void select_cos();
};

#endif
