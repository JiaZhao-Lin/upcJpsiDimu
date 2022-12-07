#include "uncertainties/uncertainties/impl.hpp"
#include "uncertainties/uncertainties/io.hpp"
#include "uncertainties/uncertainties/ureal.hpp"
#include "uncertainties/uncertainties/math.hpp"

void Test()
{
    uncertainties::udouble x(2, 1), y(2, 1);
    uncertainties::udouble a = x/y;
    uncertainties::udouble b = sqrt(x);
    std::cout << a << ", " << b << "\n";
    // x = {3,2};
    // std::cout << a.s() <<"\n";
    // std::cout << a << ", " << b << "\n";
}