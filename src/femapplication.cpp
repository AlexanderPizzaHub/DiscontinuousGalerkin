#include "femapplication.hpp"

namespace femapplication
{
    using namespace Eigen;

    Hyperbolic::Hyperbolic(femspace::FEMSpace &femspace, double left_bc_value, double right_bc_value, double velocity, double delta)
        : femspace_(&femspace), left_bc_value_(left_bc_value), right_bc_value_(right_bc_value), velocity_(velocity), delta_(delta)
    {
        // Initialize solution vector
        int Nx = femspace.GetMesh()->GetNumCells();
        solution_.resize(Nx);
        system_lhs_matrix.resize(Nx);
        system_rhs_vector.resize(Nx);

        for (int i = 0; i < Nx; ++i)
        {
            int Np = femspace.GetNp();
            solution_[i] = VectorXd::Zero(Np);
            system_lhs_matrix[i] = MatrixXd::Zero(Np, Np);
            system_rhs_vector[i] = VectorXd::Zero(Np);
        }
    }

    void Hyperbolic::Prepare()
    {
        femspace_->VolumeMassIntegrator();
        femspace_->VolumeStiffIntegrator();
        femspace_->FaceMassIntegrator();

        femspace_->BCChange(left_bc_value_, right_bc_value_);
        femspace_->BCUpdate();

        int Nx = femspace_->GetNx();
        std::vector<MatrixXd> &mass = femspace_->GetMass();
        std::vector<MatrixXd> &stiff = femspace_->GetStiff();
        MatrixXd face_mass_right = femspace_->GetFaceMassRight();
        for (int i = 0; i < Nx; i++)
        {
            system_lhs_matrix[i] = velocity_ * face_mass_right - velocity_ * stiff[i] - delta_ * mass[i];
            //std::cout << "Elemenet "<< i << " face mass right:\n" << face_mass_right << std::endl;
        }
    }

    void Hyperbolic::Solve()
    {
        int Nx = femspace_->GetNx();
        int Np = femspace_->GetNp();

        std::vector<VectorXd> solution_old;
        solution_old.resize(Nx);
        for (int i = 0; i < Nx; i++)
        {
            solution_old[i] = solution_[i];
        }

        std::vector<VectorXd> G_values = femspace_->Interpolate(settings::hyperbolic::G);
        std::vector<VectorXd> S_values = femspace_->Interpolate(settings::hyperbolic::S);

        // 获取参考元的相关数据，用于后续计算
        MatrixXd face_mass_left = femspace_->GetFaceMassLeft();
        VectorXd lagrange_left_values = femspace_->GetLagrangeLeftValues();
        std::vector<MatrixXd> &mass = femspace_->GetMass();

        double residual = 99999.0;
        int iter = 0;
        VectorXd tmpvec = VectorXd::Zero(Np);
        VectorXd tmpvec2 = VectorXd::Zero(Np);
        while (iter < 100000 && residual > 1e-10)
        {
            //std::cout << "Iteration: " << iter << std::endl;
            for (int i = 0; i < Nx; i++)
            {
                tmpvec = (solution_old[i].array() * (G_values[i].array() - delta_) + S_values[i].array());
                //system_rhs_vector[i] = mass[i].array() * (solution_old[i].array() * (G_values[i].array() - delta_) + S_values[i].array())
                system_rhs_vector[i] = mass[i] * tmpvec;

                if (i == 0)
                {
                    system_rhs_vector[i].array() += velocity_ * left_bc_value_ * lagrange_left_values.array();
                }
                else
                {
                    system_rhs_vector[i].array() += velocity_ * (face_mass_left * solution_[i - 1]).array();
                }

                //solution_[i] = system_lhs_matrix[i].colPivHouseholderQr().solve(system_rhs_vector[i]);
                tmpvec2 = system_lhs_matrix[i].colPivHouseholderQr().solve(system_rhs_vector[i]);
                //std::cout << "tmpvec " << i << ": " << tmpvec2.transpose() << std::endl;
                for(int l=0;l<Np;l++)
                {
                    solution_[i](l) = tmpvec2(Np-l-1);
                    //std::cout << "!!" << i << ": " << solution_[i](l) << std::endl; 
                }
                //std::cout << "Element " << i << " system matrix:\n" << system_lhs_matrix[i] << std::endl;
                //std::cout << "Elemenet " << i << " stiff :\n" << femspace_->GetStiff()[i] << std::endl;
                //std::cout << "Element " << i << " mass:\n" << femspace_->GetMass()[i] << std::endl;
                //std::cout << "Element " << i << " face mass left:\n" << face_mass_left << std::endl;
                //std::cout << "Element " << i << " system rhs:\n" << system_rhs_vector[i].transpose() << std::endl;
                //std::cout << "Element " << i << " solution:\n" << solution_[i].transpose() << std::endl;

                /*
                总结bug：一堆奇奇怪怪的矩阵向量转置出问题
                */
            }
            // 计算残差
            residual = 0.0;
            for (int i = 0; i < Nx; i++)
            {
                residual += (solution_[i] - solution_old[i]).norm();
                solution_old[i] = solution_[i];
            }
            iter++;
            //std::cout << "Iteration: " << iter << ", Residual: " << residual << std::endl;
        }
    }

    void Hyperbolic::OutputResults(const std::string &filename)
    {
        return;
    }
}