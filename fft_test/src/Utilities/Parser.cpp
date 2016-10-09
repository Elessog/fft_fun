#include "Utilities/Parser.hpp"
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <complex>

Parser::Parser(char *filename_):filename(filename_)
{


}

bool Parser::load_file(){
  struct stat buffer;
  if (! (stat(filename.c_str(),&buffer)==0) ){
    std::cerr << "File does not exist !" << std::endl;
    return false;
  }
  file_to_read.open(filename);
  if (file_to_read.fail()){
    std::cerr << "File could not be opened" << std::endl;
    return false;
  }
  return true;
}

bool Parser::parse(){
  std::string line;
  std::regex first_line_regex ("FE[ ]*=[ ]*([0-9]*[.]?[0-9]+)");
  std::regex data_line_regex ("[ ]*([-+]?[0-9]*[.]?[0-9]+)[ ]*,?[ ]*([-+]?[0-9]*[.]?[0-9]+)?[ ]*");
  //get number of line for proper allocation of vector
  std::size_t data_line_count = 0;
  std::smatch sm;

  while(std::getline(file_to_read,line))
      ++data_line_count;

  file_to_read.clear();
  file_to_read.seekg(0,std::ios::beg);

  std::getline(file_to_read,line);

  if (data_line_count<2){
    file_to_read.close();
    std::cerr<<"Error File is empty"<<std::endl;
    return false;
  }

  //putting the real number of line available (if file well written)
  --data_line_count;

	if(std::regex_match( line, sm,first_line_regex ))
	{
  	sampling_frequency = stof(sm[1].str());
    std::cout<< "Frequency of sample is : "<<sampling_frequency << " Hz\n" << "Number of data lines : "<< data_line_count <<std::endl;
	}
  else
  {
    file_to_read.close();
    std::cerr<<"Error File is not to correct forme: \n\t \"FE[ ]*=[ ]*([0-9]*[.])?[0-9]+\" ex : FE = 1000.50"<<std::endl;
    return false;
  }

  data.reserve(data_line_count+10);

  while(std::getline(file_to_read,line))
  {
    double real,im;
    im = 0;
    if (std::regex_match(line, sm,data_line_regex ))
    {
      real = stof(sm[1].str());

      if (sm[2].str().size()>0)
      {
        im = stof(sm[2].str());
      }
    }
    else
      break;
    data.push_back(std::complex<double>(real,im));
  }

  std::cout<<"Data count expected : "<<data_line_count<<" Data count gotten : "<< data.size()<<std::endl;
  return true;
}

signal_struct Parser::get_data()
{
 return data;
}

double Parser::get_frequency()
{
  return sampling_frequency;
}
