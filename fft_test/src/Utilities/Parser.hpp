#ifndef PARSER_HPP
#define PARSER_HPP

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <complex>
typedef std::vector<std::complex<double>> signal_struct;

class Parser{
public:
  Parser(char *filename);

  ~Parser(){};

  bool load_file();
  bool parse();
  signal_struct get_data();
  double get_frequency();

private:
  std::string filename;
  std::ifstream file_to_read;
  int number_lines;
  signal_struct data;
  double sampling_frequency;
};


#endif
