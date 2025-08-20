#include "numerics.hpp"
#include "femmesh.hpp"
#include "femspace.hpp"
#include "femapplication.hpp"

#include "test.hpp"

#include "../eigen-3.4.0/Eigen/Dense"

#include <iostream>

int main(int argc, char *argv[])
{
     test::TestQuadrature();
    // test::TestJacobiCoef();
    // test::TestMesh();
    // test::TestRefFE();
    // test::TestElemBC();
    
    //for (int NxRefine = 2; NxRefine <= 8; NxRefine++)
    //{
        // for (int order = 1; order <= 10; order++)
        // {
        //     //double l2error = test::TestHyperbolic(order, std::pow(2, NxRefine));
        //     double l2error = test::TestHyperbolic(order, 1);
        //     std::cout << "Order: " << order << ", Nx: " << 1<< ", L2 Error: " << l2error << std::endl;
        // }
    //}
    //double l2error = test::TestHyperbolic(2, 1);
    //std::cout << "L2 Error: " << l2error << std::endl;

    //test::TestQuadrature(5, -1.0, 1.0);

    return 0;
}