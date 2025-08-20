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

    MatrixXd ComputeJacobiMatrix(int n, double alpha, double beta);

    double gamma_function(double x);

    // void ComputeJacobiQuadrature(int order, double alpha, double beta, VectorXd& node_coords, VectorXd& node_weights);

    void ComputeGaussLobattoQuadrature(int order, VectorXd& node_coords, VectorXd& node_weights);

    void ComputeJacobiCoefficients(int order, double alpha, double beta, MatrixXd& coeffs);

    void ComputeJacobiCoefficients_deepseek(int order, double alpha, double beta, MatrixXd& coeffs);

    double ComputeConditionNumber(const Eigen::MatrixXd& A);

    void Trapezoid(const std::vector<MatrixXd>& f,const std::vector<double> x, MatrixXd& result);

    void Trapezoid(const std::vector<MatrixXd>& f, const std::vector<double> scale, const std::vector<double> x, MatrixXd& result);

    void Trapezoid(const std::vector<double>& f, const std::vector<double> x, double& result);
}