#include <iostream>
#include <cstdlib>

#include "Utilities/Parser.hpp"
#include "Fourier_Algo/Compute.hpp"

int main(int argc, char **argv)
{

    if (argc==2)
    {
       Parser parser(argv[1]);

       if (!parser.load_file())
          return EXIT_FAILURE;

       if (!parser.parse())
          return EXIT_FAILURE;

       ComputeFourier compute_fourier(parser.get_data(),parser.get_frequency());
       compute_fourier.applyFourier();



    }
    else{
      std::cerr << "usage is \n" << argv[0] <<" <input_file>" << std::endl;
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
