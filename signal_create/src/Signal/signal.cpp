#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "Signal/signal.hpp"

#define PI	M_PI	/* pi to machine precision, defined in math.h */
#define TWOPI	(2.0*PI)

Signal::Signal(double duration_,double frequency,char *filename):
  sampling_frequency(frequency),duration(duration_),output_name(filename)
{
  n = sampling_frequency*duration;
  duration = (double) n /sampling_frequency;
  signal_data.assign(n,0);
}


void Signal::add_cos(double frequency, double amplitude, double phase)
{
  for(int i = 0; i<n; ++i)
  {
    signal_data[i] += amplitude*cos((((double)i)/sampling_frequency)*TWOPI*frequency+phase);
  }
}

void Signal::add_sin(double frequency, double amplitude, double phase)
{
  for(int i = 0; i<n; ++i)
  {
    signal_data[i] += amplitude*sin((((double)i)/sampling_frequency)*TWOPI*frequency+phase);
  }
}

void Signal::add_rectangle(double duration_, double amplitude, double time_middle)
{
  for(int i = (int) (time_middle-duration_/2)*sampling_frequency; i< (int) (time_middle+duration_/2)*sampling_frequency; ++i)
  {
    signal_data[i] += amplitude;
  }
}

void Signal::create_signal()
{

  int i=0;

  while(i!=4)
  {
    std::cout << "Add to signal :\n\t1)\tcosinus <frequency> <amplitude> <phase>\n"
                 "\t2)\tsinus <frequency> <amplitude> <phase>\n"
                 "\t3)\trectangle <duration> <amplitude> <time>\n"
                 "\t4)\tExit and save signal to textfile\n"<< std::endl;

    std::cout << "Enter choice :";
    std::cin  >> i;
    std::cout<<std::endl;

    switch(i)
    {
      case 1: select_cos(); break;
      case 2: select_sin(); break;
      case 3: select_rectangle(); break;
      case 4: break;
      default: std::cout << "No correct value chosen please retry ... \n"<< std::endl;
    }


  }
  std::cout << "Creation of signal finished saving now ..."<<std::endl;
}

void Signal::select_sin()
{
  double frequency,amplitude,phase;
  std::cout << "Select Frequency : ";
  std::cin >> frequency;
  std::cout << "Select Amplitude : ";
  std::cin >> amplitude;
  std::cout << "Select Phase : ";
  std::cin >> phase;
  std::cout<<"Computing ..."<< std::endl;
  add_sin(frequency,amplitude,phase);
  std::cout<<"Added to signal !"<< std::endl;
}


void Signal::select_cos()
{
  double frequency,amplitude,phase;
  std::cout << "Select Frequency : ";
  std::cin >> frequency;
  std::cout << "Select Amplitude : ";
  std::cin >> amplitude;
  std::cout << "Select Phase : ";
  std::cin >> phase;
  std::cout<<"Computing ..."<< std::endl;
  add_cos(frequency,amplitude,phase);
  std::cout<<"Added to signal !"<< std::endl;
}


void Signal::select_rectangle()
{
  double duration_,amplitude,time_;
  std::cout << "Select Duration : ";
  std::cin >> duration_;
  std::cout << "Select Amplitude : ";
  std::cin >> amplitude;
  std::cout << "Select Time place : ";
  std::cin >> time_;
  std::cout<<"Computing ..."<< std::endl;
  add_rectangle(duration_,amplitude,time_);
  std::cout<<"Added to signal !"<< std::endl;
}

void Signal::save_file(){
  std::ofstream file;
  file.open(output_name);
  file << "FE = " << sampling_frequency << "\n";

  std::stringstream temp;
  temp << std::setprecision(10) << std::fixed;
  for(auto it=signal_data.cbegin(); it!=signal_data.cend(); ++it)
     temp << (*it) << "\n";
  file << temp.str();
  file.close();
}
