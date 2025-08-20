#include "numerics.hpp"
#include "femmesh.hpp"
#include "femspace.hpp"
#include "femapplication.hpp"

#include "test.hpp"

#include "../eigen-3.4.0/Eigen/Dense"

#include <iostream>

int main(int argc, char *argv[])
{
    // test::TestQuadrature();
    // test::TestQuadrature(5, -1.0, 1.0);
    // test::TestJacobiCoef();
    // test::TestMesh();
    // test::TestRefFE();
    // test::TestElemBC();

    // for(int order= 1; order <= 10; order++)
    // {
    //     test::TestApproximation(order, 0.0, 1.0);
    // }
    // test::TestProjection(1, 0.0, 1.0);

    for (int NxRefine = 0; NxRefine <= 3; NxRefine++)
    {
        for (int order = 1; order <= 10; order++)
        {
            double l2error = test::TestHyperbolic(order, std::pow(2, NxRefine));
            std::cout << "Order: " << order << ", Nx: " << std::pow(2, NxRefine) << ", L2 Error: " << l2error << std::endl;
        }
    }

    return 0;
}