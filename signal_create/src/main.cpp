#include <iostream>
#include <cstdlib>

#include "Signal/signal.hpp"

int main(int argc, char **argv)
{

    if (argc==2)
    {
       double sampling_frequency,duration;
       std::cout << "Select Sampling Frequency : ";
       std::cin >> sampling_frequency;
       std::cout << "Select Duration of signal : ";
       std::cin >> duration;
       std::cout << "Now creating signal of duration "<< duration << " s and with of sampling frequency of "<< sampling_frequency << " Hz"<< std::endl;
       Signal signal(duration,sampling_frequency,argv[1]);
       signal.create_signal();
       signal.save_file();


    }
    else{
      std::cerr << "usage is \n" << argv[0] <<" output_name" << std::endl;
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
