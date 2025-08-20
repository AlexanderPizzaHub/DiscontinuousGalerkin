#include "settings.hpp"

namespace settings
{
    namespace hyperbolic
    {
        // 算例不好，解是多项式，高阶格式的数值误差会dominate，看起来精度反而会掉
        double G(double x)
        {
            // Define the function G(x) here
            return -2.0 * x;
            //return 1.0;
            //return 0.0;
            //return 2.0;
        }
        
        double S(double x)
        {
            // Define the function S(x) here
            return 2.0 * x + 2.0 * x * x * x;
            //return 2.0*x - x*x;
            //return 1.0;
            //return std::pow(x,5);
        }

        double Truth(double x)
        {
            return std::exp(-x*x) + x*x;
            //return std::exp(x) + x*x;
            //return 1.0 + x;
            //return 1.0 + std::pow(x,6)/6.0;
        }
    }
}