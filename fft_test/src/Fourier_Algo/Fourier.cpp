/************************************************
* FFT code from the book Numerical Recipes in C *
* Visit www.nr.com for the licence.             *
************************************************/

// The following line must be defined before including math.h to correctly define M_PI
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <vector>
#include <iostream>
#include <valarray>

#include "Fourier_Algo/Fourier.hpp"
#include "Utilities/Timer.h"

#define PI	M_PI	/* pi to machine precision, defined in math.h */
#define TWOPI	(2.0*PI)


typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

/*
 FFT/IFFT routine. (see pages 507-508 of Numerical Recipes in C)

 Inputs:
	data[] : array of complex* data points of size 2*NFFT+1.
		data[0] is unused,
		* the n'th complex number x(n), for 0 <= n <= length(x)-1, is stored as:
			data[2*n+1] = real(x(n))
			data[2*n+2] = imag(x(n))
		if length(Nx) < NFFT, the remainder of the array must be padded with zeros

	nn : FFT order NFFT. This MUST be a power of 2 and >= length(x).
	isign:  if set to 1,
				computes the forward FFT
			if set to -1,
				computes Inverse FFT - in this case the output values have
				to be manually normalized by multiplying with 1/NFFT.
 Outputs:
	data[] : The FFT or IFFT results are stored in data, overwriting the input.
*/

void fourier_NR(std::vector<double> &data,int nn, int isign)
{

    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2)
    {
	    if (j > i)
      {
	      tempr = data[j];     data[j] = data[i];     data[i] = tempr;
	      tempr = data[j+1];   data[j+1] = data[i+1]; data[i+1] = tempr;
	    }
	    m = n >> 1;
	    while (m >= 2 && j > m)
      {
	      j -= m;
	      m >>= 1;
  	  }
	    j += m;
    }
    mmax = 2;
    while (n > mmax)
    {
	    istep = 2*mmax;
	    theta = TWOPI/(isign*mmax);
	    wtemp = sin(0.5*theta);
	    wpr = -2.0*wtemp*wtemp;
    	wpi = sin(theta);
	    wr = 1.0;
    	wi = 0.0;
	    for (m = 1; m < mmax; m += 2)
      {
	      for (i = m; i <= n; i += istep)
        {
		      j =i + mmax;
		      tempr = wr*data[j]   - wi*data[j+1];
		      tempi = wr*data[j+1] + wi*data[j];
		      data[j]   = data[i]   - tempr;
		      data[j+1] = data[i+1] - tempi;
		      data[i] += tempr;
		      data[i+1] += tempi;
	      }
	    wr = (wtemp = wr)*wpr - wi*wpi + wr;
	    wi = wi*wpr + wtemp*wpi + wi;
    	}
	    mmax = istep;
    }
}


void fourier_nr(const signal_struct &data_input,signal_struct &data_output ,int nn, int isign){
  std::vector<double> data;
  Timer timer;
  data.reserve(data_input.size()*2+10);
  data.push_back(0);//due to algorithm
  for(auto iter_data=data_input.cbegin(); iter_data != data_input.cend(); ++iter_data)
  {
    data.push_back((*iter_data).real());
    data.push_back((*iter_data).imag());
  }

  timer.start();
  fourier_NR(data, nn, isign);
  timer.stop();
  std::cout << "Time passed : " << timer.timePassed() << std::endl;

  data_output.clear();
  data_output = data_input;
  int n = 0;
  for(auto iter_data=data_output.begin();iter_data!=data_output.end();++iter_data)
  {
    *(iter_data) = std::complex<double>(data[2*n+1],data[2*n+2]);
    ++n;
  }

}




// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft_ct(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;

    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];

    // conquer
    fft_ct(even);
    fft_ct(odd);

    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

// Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
// Better optimized but less intuitive
void fft_cto(CArray &x)
{
	// DFT
	unsigned int N = x.size(), k = N, n;
	double thetaT = 3.14159265358979323846264338328L / N;
	Complex phiT = Complex(cos(thetaT), sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				Complex t = x[a] - x[b];
				x[a] += x[b];
				x[b] = t * T;
			}
			T *= phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			Complex t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}
	//// Normalize (This section make it not working correctly)
	//Complex f = 1.0 / sqrt(N);
	//for (unsigned int i = 0; i < N; i++)
	//	x[i] *= f;
}

// inverse fft (in-place)
/*void ifft(CArray& x)
{
    // conjugate the complex numbers
    x = x.apply(std::conj);

    // forward fft
    fft( x );

    // conjugate the complex numbers again
    x = x.apply(std::conj);

    // scale the numbers
    x /= x.size();
}*/


void fourier_ct(const signal_struct &data_input,signal_struct &data_output){
  CArray data(data_input.data(),data_input.size());
  Timer timer;

  timer.start();
  fft_ct(data);
  timer.stop();
  std::cout << "Time passed : " << timer.timePassed() << std::endl;

  data_output.assign(std::begin(data),std::end(data));
}

void fourier_cto(const signal_struct &data_input,signal_struct &data_output){
  CArray data(data_input.data(),data_input.size());
  Timer timer;

  timer.start();
  fft_cto(data);
  timer.stop();
  std::cout << "Time passed : " << timer.timePassed() << std::endl;

  data_output.assign(std::begin(data),std::end(data));
}
