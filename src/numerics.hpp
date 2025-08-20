/*
参考单元上的：
1. 求积式：求积点，权重
2. Jacobi多项式系数计算，插值点（多项式根+端点计算）
*/

#pragma once 

#include "../eigen-3.4.0/Eigen/Dense"

#include <iostream>
#include <cmath>

namespace numerics
{
    using namespace Eigen;
    void ComputeGaussQuadrature(int order,double xl, double xr, VectorXd& node_coords, VectorXd& node_weights);

    void ComputeJacobiQuadrature(int order, double alpha, double beta, VectorXd& node_coords, VectorXd& node_weights);

    void ComputeGaussLobattoQuadrature(int order, VectorXd& node_coords, VectorXd& node_weights);

    void ComputeJacobiCoefficients(int order, double alpha, double beta, MatrixXd& coeffs);

    double ComputeConditionNumber(const Eigen::MatrixXd& A);
}