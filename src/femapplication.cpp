#include "femapplication.hpp"

namespace femapplication
{
    using namespace Eigen;

#pragma region "Hyperbolic Application"
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
        MatrixXd face_mass_left = femspace_->GetFaceMassLeft();
        for (int i = 0; i < Nx; i++)
        {
            if(velocity_>=0.0)
            {
                system_lhs_matrix[i] = velocity_ * face_mass_right - velocity_ * stiff[i] - delta_ * mass[i];
            }
            else 
            {
                system_lhs_matrix[i] = - velocity_ * face_mass_left - velocity_ * stiff[i] - delta_ * mass[i];
            }
            
        }
    }

    void Hyperbolic::UpdateBCValues(double left_bc_value, double right_bc_value)
    {
        left_bc_value_ = left_bc_value;
        right_bc_value_ = right_bc_value;
    }

    void Hyperbolic::Solve()
    {
        /*
        总结bug：
        1. 一堆奇奇怪怪的矩阵向量转置出问题
        2. 误差统一放到体心上计算
        3. 单元左右面通量的计算有问题，已debug。
        4. 注意Np和order的区别。elem和refelem的order全部传成np了。
        */

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
        VectorXd lagrange_left_values = femspace_->GetLagrangeLeftValues();
        VectorXd lagrange_right_values = femspace_->GetLagrangeRightValues();
        std::vector<MatrixXd> &mass = femspace_->GetMass();

        double residual = 99999.0;
        int iter = 0;
        VectorXd tmpvec = VectorXd::Zero(Np);
        while (iter < 100000 && residual > 1e-10)
        {
            for (int i = 0; i < Nx; i++)
            {
                tmpvec = (solution_old[i].array() * (G_values[i].array() - delta_) + S_values[i].array());
                system_rhs_vector[i] = mass[i] * tmpvec;

                if (i == 0)
                {
                    system_rhs_vector[i].array() += velocity_ * left_bc_value_ * lagrange_left_values.array();
                }
                else
                {
                    system_rhs_vector[i].array() += velocity_ * (solution_[i - 1].dot(lagrange_right_values)) * lagrange_left_values.array();
                }

                solution_[i] = system_lhs_matrix[i].colPivHouseholderQr().solve(system_rhs_vector[i]);
            }

            // 计算残差
            residual = 0.0;
            for (int i = 0; i < Nx; i++)
            {
                residual += (solution_[i] - solution_old[i]).norm();
                solution_old[i] = solution_[i];
            }
            iter++;
        }
    }

    void Hyperbolic::Solve(const MatrixXd G, const MatrixXd S, MatrixXd &solution)
    {
        int Nx = femspace_->GetNx();
        int Np = femspace_->GetNp();

        MatrixXd solution_old = MatrixXd::Zero(Nx, Np);
        // solution = MatrixXd::Zero(Nx, Np);
        solution.setZero();

        // 获取参考元的相关数据，用于后续计算
        VectorXd lagrange_left_values = femspace_->GetLagrangeLeftValues();
        VectorXd lagrange_right_values = femspace_->GetLagrangeRightValues();
        std::vector<MatrixXd> &mass = femspace_->GetMass();

        double residual = 99999.0;
        int iter = 0;
        VectorXd tmpvec = VectorXd::Zero(Np);
        VectorXd G_row = VectorXd::Zero(Np);
        VectorXd S_row = VectorXd::Zero(Np);
        //std::cout << "Solving hyperbolic equation with velocity: " << velocity_ << std::endl;
        //std::cout << "G: " << G << std::endl;
        //std::cout << "S: " << S << std::endl;
        while (iter < 10000 && residual > 1e-8)
        {
            if (velocity_ >= 0.0)
            {
                for (int i = 0; i < Nx; i++)
                {
                    G_row = G.row(i).transpose();
                    S_row = S.row(i).transpose();
                    tmpvec = (solution_old.row(i).transpose().array() * (G_row.array() - delta_) + S_row.array());
                    // std::cout << "tmpvec: " << tmpvec.transpose() << std::endl;
                    // std::cout << "G row " << i << ": " << G_row << std::endl;
                    // std::cout << "S row " << i << ": " << S_row << std::endl;
                    // std::cout << "solution_old row " << i << ": " << solution_old.row(i) << std::endl;
                    // std::cout << "tmpvec: " << tmpvec.transpose() << std::endl;
                    system_rhs_vector[i] = mass[i] * tmpvec;
                    // std::cout << "system_rhs_vector[" << i << "]: " << system_rhs_vector[i].transpose() << std::endl;

                    if (i == 0)
                    {
                        system_rhs_vector[i].array() += velocity_ * left_bc_value_ * lagrange_left_values.array();
                    }
                    else
                    {
                        system_rhs_vector[i].array() += velocity_ * (solution.row(i - 1).dot(lagrange_right_values)) * lagrange_left_values.array();
                        // std::cout << "solution row " << i - 1 << ": " << solution.row(i - 1) << std::endl;
                    }
                    // std::cout << "system_rhs_vector post[" << i << "]: " << system_rhs_vector[i].transpose() << std::endl;
                    solution.row(i) = system_lhs_matrix[i].colPivHouseholderQr().solve(system_rhs_vector[i]);
                }
            }
            else
            {
                for (int i = Nx - 1; i >= 0; i--)
                {
                    // std::cout << "Processing row: " << i << std::endl;
                    G_row = G.row(i);
                    S_row = S.row(i);
                    // std::cout << "!!!" << std::endl;
                    // std::cout << "G row " << i << ": " << G_row << std::endl;
                    // std::cout << "solution_old row " << i << ": " << solution_old.row(i) << std::endl;
                    tmpvec = (solution_old.row(i).transpose().array() * (G_row.array() - delta_) + S_row.array());
                    // std::cout << "tmpvec: " << tmpvec.transpose() << std::endl;
                    system_rhs_vector[i] = mass[i] * tmpvec;
                    // std::cout << "tmpvec: " << tmpvec.transpose() << std::endl;

                    if (i == Nx - 1)
                    {
                        system_rhs_vector[i].array() -= velocity_ * right_bc_value_ * lagrange_right_values.array();
                    }
                    else
                    {
                        system_rhs_vector[i].array() -= velocity_ * (solution.row(i + 1).dot(lagrange_left_values)) * lagrange_right_values.array();
                    }
                    // std::cout << "system_rhs_vector[" << i << "]: " << system_rhs_vector[i].transpose() << std::endl;

                    solution.row(i) = system_lhs_matrix[i].colPivHouseholderQr().solve(system_rhs_vector[i]);
                }
            }

            // 计算残差
            residual = 0.0;
            residual = (solution - solution_old).norm();
            solution_old = solution;
            iter++;
        }
        if (iter >= 10000)
        {
            std::cerr << "Warning: Hyperbolic solver did not converge within 100000 iterations." << std::endl;
            std::cerr << "Last residual: " << residual << std::endl;
            std::cerr << "Last solution: " << solution.transpose() << std::endl;
        }
        //std::cout << "Iteration: " << iter << ", Residual: " << residual << std::endl;
        //std::cout << "Solution: " << solution.transpose() << std::endl;
    }

    void Hyperbolic::OutputResults(const std::string &filename)
    {
        return;
    }
#pragma endregion

#pragma region "Shakhov 1D 1V Application"
    Shakhov1D1V::Shakhov1D1V(std::vector<double> v, femspace::FEMSpace &femspace)
        : v_vec_(v), femspace_(&femspace)
    {
        // std::cout << "@@@" <<std::endl;
        //  assert v strictly increasing
        Nv_ = v.size();
        Nx_ = femspace_->GetNx();
        Np_ = femspace_->GetNp();

        f1_.resize(Nv_);
        f2_.resize(Nv_);

        S1_.resize(Nv_);
        S2_.resize(Nv_);

        hyperbolic_.resize(Nv_);
        // std::cout << "!!!" <<std::endl;

        for (int i = 0; i < Nv_; ++i)
        {
            f1_[i] = MatrixXd::Zero(Nx_, Np_);
            f2_[i] = MatrixXd::Zero(Nx_, Np_);

            S1_[i] = MatrixXd::Zero(Nx_, Np_);
            S2_[i] = MatrixXd::Zero(Nx_, Np_);

            hyperbolic_[i] = new Hyperbolic(femspace, 0.0, 0.0, v_vec_[i], settings::shakhov::delta_eff);
            // std::cout << "$$$" <<std::endl;
        }

        G_ = MatrixXd::Zero(femspace_->GetNx(), femspace_->GetNp());

        f1_bc_left_ = VectorXd::Zero(Nv_);
        f1_bc_right_ = VectorXd::Zero(Nv_);
        f2_bc_left_ = VectorXd::Zero(Nv_);
        f2_bc_right_ = VectorXd::Zero(Nv_);

        f1_left_values_ = VectorXd::Zero(Nv_);
        f1_right_values_ = VectorXd::Zero(Nv_);

        macro_density_ = MatrixXd::Zero(Nx_, Np_);
        macro_velocity_ = MatrixXd::Zero(Nx_, Np_);
        macro_temperature_ = MatrixXd::Zero(Nx_, Np_);
        macro_energy_flux_ = MatrixXd::Zero(Nx_, Np_);
        std::cout << "Shakhov1D1V initialized with Nv: " << Nv_ << ", Nx: " << Nx_ << ", Np: " << Np_ << std::endl;
    }

    void Shakhov1D1V::UpdateMacroVariables()
    {
        MatrixXd tmp_mat = MatrixXd::Zero(Nx_, Np_); // inter mat

        VectorXd tmp_vec = VectorXd::Zero(Nv_); // inter vec
        double tmp_var = 0.0;                   // inter var

        double T, q;
        // mass
       // for(int vi=0; vi < Nv_; vi++)
       // {
            //std::cout << "f1_[" << vi << "]: " << f1_[vi] << std::endl;
       // }
        numerics::Trapezoid(f1_, v_vec_, macro_density_);

        // velocity
        numerics::Trapezoid(f1_, v_vec_, v_vec_, tmp_mat);
        macro_velocity_ = tmp_mat.array() / macro_density_.array();

        // temperature and energy flux
        std::vector<MatrixXd> therm_rhs;
        std::vector<MatrixXd> flux_rhs;
        therm_rhs.resize(Nv_);
        flux_rhs.resize(Nv_);
        for (int i = 0; i < Nv_; ++i)
        {
            therm_rhs[i] = MatrixXd::Zero(Nx_, Np_);
            flux_rhs[i] = MatrixXd::Zero(Nx_, Np_);
        }

        double v;
        MatrixXd therm_speed = MatrixXd::Zero(Nx_, Np_);
        for (int vi = 0; vi < Nv_; vi++)
        {
            v = v_vec_[vi];
            therm_speed = v - macro_velocity_.array();
            therm_rhs[vi] = therm_speed.array().square() * f1_[vi].array() + f2_[vi].array();
            flux_rhs[vi] = therm_speed.array() * therm_rhs[vi].array();
        }
        numerics::Trapezoid(therm_rhs, v_vec_, tmp_mat);
        macro_temperature_ = 2.0 * tmp_mat.array() / (3.0 * macro_density_.array());
        numerics::Trapezoid(flux_rhs, v_vec_, macro_energy_flux_);

        std::cout << "in macro updating"    << std::endl;
        std::cout << "Macro density: " << macro_density_.transpose() << std::endl;
        std::cout << "Macro velocity: " << macro_velocity_.transpose() << std::endl;
        std::cout << "Macro temperature: " << macro_temperature_.transpose() << std::endl;
        std::cout << "Macro energy flux: " << macro_energy_flux_.transpose() << std::endl;
        std::cout << "macro update finished" << std::endl;
    }

    void Shakhov1D1V::UpdateF1BoundaryValues()
    {
        femspace::REFFE *ref_elem = femspace_->GetRefElem();
        VectorXd lagrange_left_values = ref_elem->lagrange_left_values_;
        VectorXd lagrange_right_values = ref_elem->lagrange_right_values_;

        VectorXd f1_left_cell = VectorXd::Zero(Np_);
        VectorXd f1_right_cell = VectorXd::Zero(Np_);

        for (int vi = 0; vi < Nv_; vi++)
        {
            f1_left_cell = f1_[vi].row(0).transpose();
            f1_right_cell = f1_[vi].row(Nx_ - 1).transpose();

            f1_left_values_(vi) = f1_left_cell.dot(lagrange_left_values);
            f1_right_values_(vi) = f1_right_cell.dot(lagrange_right_values);
        }
    }

    void Shakhov1D1V::UpdateBCDensity()
    {
        // fully diffuse wall
        double T_wall_l = settings::shakhov::T_wall_l;
        double T_wall_r = settings::shakhov::T_wall_r;
        double u_wall_l = settings::shakhov::u_wall_l;
        double u_wall_r = settings::shakhov::u_wall_r;

        double left_in = 0.0, left_out = 0.0;
        double right_in = 0.0, right_out = 0.0;

        std::vector<double> left_in_rhs, left_out_rhs, right_in_rhs, right_out_rhs;
        left_in_rhs.resize(Nv_);
        left_out_rhs.resize(Nv_);
        right_in_rhs.resize(Nv_);
        right_out_rhs.resize(Nv_);
        for (int vi = 0; vi < Nv_; vi++)
        {
            double v = v_vec_[vi];
            left_in_rhs[vi] = 0.0;
            left_out_rhs[vi] = 0.0;
            right_in_rhs[vi] = 0.0;
            right_out_rhs[vi] = 0.0;

            if (v >= settings::shakhov::u_wall_l)
            {
                left_in_rhs[vi] = (v - u_wall_l) * std::sqrt(1.0 / (M_PI * T_wall_l)) * std::exp(-(v - u_wall_l) * (v - u_wall_l) / T_wall_l);
            }
            if (v <= settings::shakhov::u_wall_l)
            {
                // left_out_rhs[vi] = (v - u_wall_l) * f1_(0, vi);
                left_out_rhs[vi] = (v - u_wall_l) * f1_left_values_(vi);
            }

            if (v <= settings::shakhov::u_wall_r)
            {
                right_in_rhs[vi] = (v - u_wall_r) * std::sqrt(1.0 / (M_PI * T_wall_r)) * std::exp(-(v - u_wall_r) * (v - u_wall_r) / T_wall_r);
            }
            if (v >= settings::shakhov::u_wall_r)
            {
                // right_out += (v - u_wall_r) * f1_(2 * Nx_ - 1, vi) * dv_;
                right_out_rhs[vi] = (v - u_wall_r) * f1_right_values_(vi);
            }
        }
        numerics::Trapezoid(left_in_rhs, v_vec_, left_in);
        numerics::Trapezoid(left_out_rhs, v_vec_, left_out);
        numerics::Trapezoid(right_in_rhs, v_vec_, right_in);
        numerics::Trapezoid(right_out_rhs, v_vec_, right_out);

        bc_density_left_ = -left_out / left_in;
        bc_density_right_ = -right_out / right_in;
        bc_temperature_left_ = T_wall_l;
        bc_temperature_right_ = T_wall_r;
        // std::cout << "bc_density_left: " << bc_density_left_ << ", bc_temperature_left: " << T_wall_l << std::endl;
    }

    void Shakhov1D1V::UpdateBCVDF()
    {
        double v;
        double u_l = settings::shakhov::u_wall_l;
        double u_r = settings::shakhov::u_wall_r;

        for (int vi = 0; vi < Nv_; vi++)
        {
            v = v_vec_[vi];

            f1_bc_left_(vi) = bc_density_left_ / (std::sqrt(M_PI * bc_temperature_left_)) * std::exp(-(v - u_l) * (v - u_l) / bc_temperature_left_);

            f1_bc_right_(vi) = bc_density_right_ / (std::sqrt(M_PI * bc_temperature_right_)) * std::exp(-(v - u_r) * (v - u_r) / bc_temperature_right_);

            f2_bc_left_(vi) = f1_bc_left_(vi) * bc_temperature_left_;
            f2_bc_right_(vi) = f1_bc_right_(vi) * bc_temperature_right_;

            /*
            f1_bc_left_(vi) = bc_density_left_ / (std::sqrt(M_PI * constants::T_wall_l)) * std::exp(-(v - u_l) * (v - u_l) / constants::T_wall_l);

            f1_bc_right_(vi) = bc_density_right_ / (std::sqrt(M_PI * constants::T_wall_r)) * std::exp(-(v - u_r) * (v - u_r) / constants::T_wall_r);


            f2_bc_left_(vi) = f1_bc_left_(vi) * constants::T_wall_l;
            f2_bc_right_(vi) = f1_bc_right_(vi) * constants::T_wall_r;
            */

            // std::cout << "f1_bc_left: " << f1_bc_left_(vi) << ", f1_bc_right: " << f1_bc_right_(vi) << std::endl;
            // std::cout << "f2_bc_left: " << f2_bc_left_(vi) << ", f2_bc_right: " << f2_bc_right_(vi) << std::endl;
        }
    }

    void Shakhov1D1V::ComputeSourceFs()
    {
        MatrixXd nu = MatrixXd::Zero(Nx_, Np_);

        nu = settings::shakhov::delta_0 * macro_density_.array() * macro_temperature_.array().pow(1.0 - settings::shakhov::omega);
        G_.array() = -nu;

        MatrixXd therm_v = MatrixXd::Zero(Nx_, Np_);

        MatrixXd f1s = MatrixXd::Zero(Nx_, Np_);
        MatrixXd f2s = MatrixXd::Zero(Nx_, Np_);

        MatrixXd fs_tmp = MatrixXd::Zero(Nx_, Np_);
        MatrixXd energy_flux_tmp = MatrixXd::Zero(Nx_, Np_);

        MatrixXd testtmp = MatrixXd::Zero(Nx_, Np_);
        for (int vi = 0; vi < Nv_; vi++)
        {
            therm_v = v_vec_[vi] - macro_velocity_.array();
            // std::cout << v_vec[vi] << std::endl;

            // fs_tmp = macro_density_.array() / (M_PI * macro_temperature_.array()).sqrt() *
            //  (-(therm_v.array().square() / macro_temperature_.array()).exp

            fs_tmp = macro_density_.array() / ((M_PI * macro_temperature_.array()).sqrt()) * ((-therm_v.array().square() / macro_temperature_.array()).exp());

            // std::cout << "fs_tmp: " << fs_tmp << std::endl;
            energy_flux_tmp = 0.8 * (1.0 - settings::shakhov::Pr) * macro_energy_flux_.array() * therm_v.array() / (macro_density_.array() * macro_temperature_.array() * macro_temperature_.array());

            f1s = fs_tmp.array() * (1.0 + energy_flux_tmp.array() * ((therm_v.array().square()) / macro_temperature_.array() - 1.5));

            f2s = fs_tmp.array() * macro_temperature_.array() * (1.0 + energy_flux_tmp.array() * (therm_v.array().square() / macro_temperature_.array() - 0.5));

            S1_[vi] = nu.array() * f1s.array();
            S2_[vi] = nu.array() * f2s.array();
            //std::cout << "f1_source_S: " << S1_[vi] << std::endl;
            //std::cout << "f2_source_S: " << S2_[vi] << std::endl;
        }

        // std::cout << "f2_source_S: " << f2_source_S_.row(xi)<< std::endl;
        //std::cout << "Source_G: " << G_ << std::endl;
        // std::cout << "f2_source_G: " << f2_source_G_.row(xi) << std::endl;
    }

    void Shakhov1D1V::Initialize()
    {
        bc_density_left_ = 0.0;
        bc_density_right_ = 0.0;
        bc_temperature_left_ = 0.0;
        bc_temperature_right_ = 0.0;

        double v, value;
        for (int vi = 0; vi < Nv_; vi++)
        {
            v = v_vec_[vi];
            value = std::sqrt(1.0 / M_PI) * std::exp(-v * v);
            // value = numerics::MaxwellianDistribution(v, 1.5*1.5);
            // f1_[vi] = MatrixXd::Constant(Nx_, Np_, value);
            // f2_[vi] = MatrixXd::Constant(Nx_, Np_, value);
            f1_[vi].array() = value;
            f2_[vi].array() = value;
        }
    }

    void Shakhov1D1V::ReNormalize()
    {
        using namespace settings::shakhov;
        //std::cout << "in normalization " << std::endl;
        double f1_sum = 0.0;
        for (int vi = 0; vi < Nv_; vi++)
        {
            f1_sum += f1_[vi].sum();
           // std::cout << "fi_[" << vi << "] : " << f1_[vi] << std::endl;
            //std::cout << "f1_[" << vi << "] sum: " << f1_[vi].sum() << std::endl;
        }
        //std::cout << "f1_sum: " << f1_sum << std::endl;
        f1_sum = f1_sum / Nv_ * 2.0 * vm / (Nx_ * Np_) * (xr-xl);
        

        for (int vi = 0; vi < Nv_; vi++)
        {
            f1_[vi] /= (f1_sum);
            f2_[vi] /= f1_sum;
            
            //std::cout << "f1_[" << vi << "] sum: " << f1_[vi].sum() << std::endl;
        }
    }

    void Shakhov1D1V::Solve(double tolerance, int max_iter)
    {
        // Initialize VDF
        Initialize();

        MatrixXd temperature_old = MatrixXd::Zero(Nx_, Np_);
        MatrixXd density_old = MatrixXd::Zero(Nx_, Np_);
        MatrixXd velocity_old = MatrixXd::Zero(Nx_, Np_);
        MatrixXd energy_flux_old = MatrixXd::Zero(Nx_, Np_);

        UpdateMacroVariables();

        temperature_old = macro_temperature_;
        density_old = macro_density_;
        velocity_old = macro_velocity_;
        energy_flux_old = macro_energy_flux_;

        // Update boundary conditions
        UpdateF1BoundaryValues();
        UpdateBCDensity();
        UpdateBCVDF();

        double residual = 999.0;
        double resT, resn, resu, resq = 0.0;
        double MASS;

        std::vector<double> mass_hist;

        // numerics::Trapezoid(macro_density_,MASS);
        // MASS = macro_density_.sum() * dx_ / 2.0;
        double volume = 1.0;
        MASS = macro_density_.sum() / volume;
        mass_hist.push_back(MASS);
        int iter = 0;

        // Prepare hyperbolic solver
        for (int vi = 0; vi < Nv_; ++vi)
        {
            hyperbolic_[vi]->Prepare();
        }

        while (residual > tolerance && iter < max_iter)
        {
            MASS = macro_density_.sum() / volume;
            // std::cout << "Mass: " << MASS << std::endl;
            mass_hist.push_back(MASS);

            if (iter % settings::shakhov::output_interval == 0)
            {
                if (settings::shakhov::printinfo)
                {
                    std::cout << "Iteration: " << iter << std::endl;
                    std::cout << "Macro density: " << macro_density_.transpose() << std::endl;
                    std::cout << "macro temperature: " << macro_temperature_.transpose() << std::endl;
                    std::cout << "macro velocity: " << macro_velocity_.transpose() << std::endl;
                    MASS = macro_density_.sum() / volume;
                    std::cout << "Mass: " << MASS << std::endl;

                    //std::cout << "Residual in shakhov: " << residual << std::endl;
                    //std::cout << "resT: " << resT << ", resn: "
                    //          << resn << ", resu: " << resu << ", resq: " << resq << std::endl;
                    std::cout << "mean temperature: " << macro_temperature_.mean() << std::endl;
                    //std::cout << "min temperature: " << macro_temperature_.minCoeff() << std::endl;
                    //std::cout << "max temperature: " << macro_temperature_.maxCoeff() << std::endl;
                }
            }

            // Compute macro variables
            UpdateMacroVariables();

            // Update boundary conditions
            UpdateF1BoundaryValues();
            UpdateBCDensity();
            UpdateBCVDF();

            // Compute source terms
            ComputeSourceFs();

            // Solve f1 and f2 using DG method
            for (int vi = 0; vi < Nv_; ++vi)
            {
                hyperbolic_[vi]->UpdateBCValues(f1_bc_left_(vi), f1_bc_right_(vi));
                hyperbolic_[vi]->Solve(G_, S1_[vi], f1_[vi]);

                hyperbolic_[vi]->UpdateBCValues(f2_bc_left_(vi), f2_bc_right_(vi));
                hyperbolic_[vi]->Solve(G_, S2_[vi], f2_[vi]);
            }

            ReNormalize();

            resT = (macro_temperature_ - temperature_old).norm() / macro_temperature_.norm();
            resn = (macro_density_ - density_old).norm() / macro_density_.norm();
            resu = (macro_velocity_ - velocity_old).norm(); // 速度趋于0所以不用相对误差
            resq = (macro_energy_flux_ - energy_flux_old).norm() / macro_energy_flux_.norm();
            residual = std::max({resT, resn, resq});
            // residual = (resT + resn + resq) / 3.0;

            temperature_old = macro_temperature_;
            density_old = macro_density_;
            velocity_old = macro_velocity_;
            energy_flux_old = macro_energy_flux_;
            iter++;
        }
        // numerics::WriteVec(filename + "mass_hist.csv", mass_hist);
        std::cout << "final num iter: " << iter << std::endl;
    }

}

#pragma endregion