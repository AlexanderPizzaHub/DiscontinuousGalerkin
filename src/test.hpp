/*
test examples for modules
*/

#pragma once

#include "numerics.hpp"
#include "femspace.hpp"
#include "femmesh.hpp"
#include "settings.hpp"
#include "femapplication.hpp"

#include "../eigen-3.4.0/Eigen/Dense"

#include <iostream>
#include <chrono>

namespace test
{
    using namespace Eigen;

    void TestQuadrature();

    void TestJacobiCoef();

    void TestMesh();

    void TestRefFE();

    void TestApproximation(int order, double xl, double xr);

    void TestProjection(int order, double xl, double xr);

    void TestElemBC();
    double TestQuadrature(int order, double xl, double xr);

    double TestHyperbolic(int order, int Nx);

    void TestShakhov1D1V(int order);

}