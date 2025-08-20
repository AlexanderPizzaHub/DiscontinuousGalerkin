/*
编写实际算例。输入：settings.hpp

1. 稳态双曲方程application
femmesh -> femspace -> femapplication
*/

#pragma once

#include "femmesh.hpp"
#include "femspace.hpp"
#include "numerics.hpp"
#include "settings.hpp"

#include "../eigen-3.4.0/Eigen/Dense"

#include <string>

namespace femapplication
{
    using namespace Eigen;

    class Hyperbolic
    {
        public:
            Hyperbolic(femspace::FEMSpace& femspace, double left_bc_value, double right_bc_value, double velocity,double delta);
            ~Hyperbolic() = default;

            void Prepare();
            // 求解
            void Solve();

            void OutputResults(const std::string& filename);

            std::vector<VectorXd> GetSolution() const { return solution_; }

        private:
            femspace::FEMSpace* femspace_; // FEM空间
            double left_bc_value_; // 左侧边界值
            double right_bc_value_; // 右侧边界值
            double velocity_; // 速度
            double delta_;

            std::vector<MatrixXd> system_lhs_matrix; // [Nx, Np, Np] 系统矩阵
            std::vector<VectorXd> system_rhs_vector; // [Nx, Np] 系

            std::vector<VectorXd> solution_; // [Nx, Np] 解向量
    };
}