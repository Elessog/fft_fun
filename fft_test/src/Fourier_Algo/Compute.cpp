#include "Fourier_Algo/Compute.hpp"
#include "Fourier_Algo/Fourier.hpp"
#include "Utilities/gnuplot_i.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

ComputeFourier::ComputeFourier(const signal_struct &signal,double frequency):
signal_input(signal),sampling_frequency(frequency)
{
  //zero padding
  int add_n_zero = ceil(log(signal_input.size())/log(2));
  nn = add_n_zero;
  signal_input.insert(signal_input.end(),pow(2,add_n_zero)-signal_input.size(),0);

}


void ComputeFourier::applyFourier(){
  int i=0;
  std::cout << "Apply :\n\t1)\tApply FFT/IFFT NR\n"
               "\t2)\tApply FFT Cooley-Tukey\n"
               "\t3)\tApply FFT Cooley-Tukey Optimized"<< std::endl;

  std::cout << "Enter choice :";
  std::cin  >> i;
  std::cout<<std::endl;

  switch(i)
  {
    case 1: select_fftnr(); break;
    case 2: select_fftct(); break;
    case 3: select_fftcto(); break;
    default: std::cout << "No correct value chosen please retry ... \n"<< std::endl;
  }

}

void ComputeFourier::select_fftnr()
{

  int isign;
  std::cout << "Select FFT 1 IFFT -1 : ";
  std::cin >> isign;

  fourier_nr(signal_input,signal_output ,nn, isign);

  std::vector<double> abs_fft;
  std::vector<double> frequencies;

  abs_fft.reserve(signal_output.size());
  frequencies.reserve(signal_output.size());

  int n = 0;
  for (auto it=signal_output.cbegin(); it!=signal_output.cend(); ++it)
  {
       abs_fft.push_back(abs(*it));
       frequencies.push_back((double(n)/double(signal_output.size()))*sampling_frequency);
       ++n;
  }

  Gnuplot g1("FFT");

  g1.reset_all();
  g1.set_grid();
  g1.set_style("lines").plot_xy(frequencies,abs_fft,"user-defined points 2d");
  wait_for_key();
}

void ComputeFourier::select_fftct()
{

  fourier_ct(signal_input,signal_output);

  std::vector<double> abs_fft;
  std::vector<double> frequencies;

  abs_fft.reserve(signal_output.size());
  frequencies.reserve(signal_output.size());

  int n = 0;
  for (auto it=signal_output.cbegin(); it!=signal_output.cend(); ++it)
  {
       abs_fft.push_back(abs(*it));
       frequencies.push_back((double(n)/double(signal_output.size()))*sampling_frequency-sampling_frequency/2.0);
       ++n;
  }


  fftshift(abs_fft);

  Gnuplot g1("FFT");

  g1.reset_all();
  g1.set_grid();
  g1.set_style("lines").plot_xy(frequencies,abs_fft,"user-defined points 2d");
  wait_for_key();
}

void ComputeFourier::select_fftcto()
{

  fourier_cto(signal_input,signal_output);

  std::vector<double> abs_fft;
  std::vector<double> frequencies;

  abs_fft.reserve(signal_output.size());
  frequencies.reserve(signal_output.size());

  int n = 0;
  for (auto it=signal_output.cbegin(); it!=signal_output.cend(); ++it)
  {
       abs_fft.push_back(abs(*it));
       frequencies.push_back((double(n)/double(signal_output.size()))*sampling_frequency-sampling_frequency/2.0);
       ++n;
  }
  fftshift(abs_fft);

  Gnuplot g1("FFT");

  g1.reset_all();
  g1.set_grid();
  g1.set_style("lines").plot_xy(frequencies,abs_fft,"user-defined points 2d");
  wait_for_key();
}

void ComputeFourier::fftshift(std::vector<double> &data)
{
    int k = 0;
    int c = (int) floor((float)data.size()/2);
    // For odd and for even numbers of element use different algorithm
    if (data.size() % 2 == 0)
    {
        for (k = 0; k < c; k++){
          double tmp = data[k];
          data[k] = data[k+c];
          data[k+c] = tmp;
        }
    }
    else
    {
        double tmp = data[0];
        for (k = 0; k < c; k++)
        {
            data[k] = data[c + k + 1];
            data[c + k + 1] = data[k + 1];
        }
        data[c] = tmp;
    }
}

/*void ifftshift(complex *data, int count)
{
    int k = 0;
    int c = (int) floor((float)count/2);
    if (count % 2 == 0)
    {
        for (k = 0; k < c; k++)
            swap(&data[k], &data[k+c]);
    }
    else
    {
        complex tmp = data[count - 1];
        for (k = c-1; k >= 0; k--)
        {
            data[c + k + 1] = data[k];
            data[k] = data[c + k];
        }
        data[c] = tmp;
    }
}*/
