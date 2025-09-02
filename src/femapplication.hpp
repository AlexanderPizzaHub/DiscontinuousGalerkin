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
        Hyperbolic(femspace::FEMSpace &femspace, double left_bc_value, double right_bc_value, double velocity, double delta);
        ~Hyperbolic() = default;

        void Prepare();

        void UpdateBCValues(double left_bc_value, double right_bc_value);
        // 求解
        void Solve();

        void Solve(const MatrixXd G, const MatrixXd S, MatrixXd& solution);

        void OutputResults(const std::string &filename);

        std::vector<VectorXd> GetSolution() const { return solution_; }

    private:
        femspace::FEMSpace *femspace_; // FEM空间
        double left_bc_value_;         // 左侧边界值
        double right_bc_value_;        // 右侧边界值
        double velocity_;              // 速度
        double delta_;

        std::vector<MatrixXd> system_lhs_matrix; // [Nx, Np, Np] 系统左端矩阵
        std::vector<VectorXd> system_rhs_vector; // [Nx, Np] 系统右端项

        std::vector<VectorXd> solution_; // [Nx, Np] 解向量
    };

    class Shakhov1D1V
    {
    public:
        Shakhov1D1V(std::vector<double> v, femspace::FEMSpace& femspace);
        ~Shakhov1D1V() = default;

        // solve both f1 and f2
        void Solve(double tolerance = 1e-6, int max_iter = 1000);

        void Initialize();

        std::vector<double> dv_vec_;

        // getters
        std::vector<MatrixXd> GetF1() const { return f1_; }
        std::vector<MatrixXd> GetF2() const { return f2_; }
        MatrixXd GetMacroDensity() const { return macro_density_; }
        MatrixXd GetMacroVelocity() const { return macro_velocity_; }
        MatrixXd GetMacroTemperature() const { return macro_temperature_; }
        MatrixXd GetMacroEnergyFlux() const { return macro_energy_flux_; }

    private:
        int Nv_; // number of velocity nodal points
        double v_left_, v_right_;

        int Nx_, Np_; // Nx: number of elements, Np: number of polynomial order + 1

        femspace::FEMSpace *femspace_; // FEM空间

        std::vector<Hyperbolic*>  hyperbolic_; // [Nv] 双曲方程

        std::vector<double> v_vec_;
        

        // VDF function. For each v, each row: Np quadrature points, each column: Nx elements
        std::vector<MatrixXd> f1_; // [Nv, Nx, Np] 
        std::vector<MatrixXd> f2_; // [Nv, Nx, Np]

        // VDF source, S and G. S与速度有关，G与速度无关
        std::vector<MatrixXd> S1_; // [Nv, Nx, Np]
        std::vector<MatrixXd> S2_; // [Nv, Nx, Np]

        MatrixXd G_; // [Nx, Np] G与速度无关

        // VDF boundary, left and right
        VectorXd f1_bc_left_;  // [Nv,1]
        VectorXd f1_bc_right_; // [Nv,1]
        VectorXd f2_bc_left_;  // [Nv,1]
        VectorXd f2_bc_right_; // [Nv,1]

        // VDF boundary value
        VectorXd f1_left_values_;  // [Nv,1]
        VectorXd f1_right_values_; // [Nv,1]

        // macro variables on cell center
        MatrixXd macro_density_;     // [Nx, Np]
        MatrixXd macro_velocity_;    // [Nx, Np]
        MatrixXd macro_temperature_; // [Nx, Np]
        MatrixXd macro_energy_flux_; // [Nx, Np]

        // macro density on boundary
        double bc_density_left_;
        double bc_density_right_;

        // macro temperature on boundary
        double bc_temperature_left_;
        double bc_temperature_right_;

        // compute macro
        void UpdateMacroVariables();

        // compute VDF boundary value
        void UpdateF1BoundaryValues();

        // boundary condition update
        void UpdateBCDensity();
        void UpdateBCVDF();

        // velocity quadrature by dv_vec
        void VelocityTrapezoid(const std::vector<MatrixXd>& f, MatrixXd& result);
        void VelocityTrapezoid(const std::vector<MatrixXd>& f, const std::vector<double> scale,  MatrixXd& result);

        // compute source
        void ComputeSourceFs();

        void ReNormalize();

        // TEMP VARIABLES 
        //MatrixXd solution_old;
    };
}