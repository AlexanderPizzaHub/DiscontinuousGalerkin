#include "femspace.hpp"

namespace femspace
{
    using namespace Eigen;

#pragma region FE
    FE::FE(int Np, double xl, double xr)
        : Np_(Np), xl_(xl), xr_(xr), dx_(xr - xl)
    {
        // Initialize the quadrature points and weights
        VectorXd node_weights;
        // numerics::ComputeGaussLobattoQuadrature(Np_, node_coords, node_weights);
        numerics::ComputeGaussQuadrature(Np_, xl, xr, node_coords_, node_weights);

        node_values_ = VectorXd::Zero(Np_);
    }
#pragma endregion

#pragma region REFFE
    REFFE::REFFE(int order)
    {
        order_ = order;
        Np_ = order + 1; // 插值点数

        // 求积点
        numerics::ComputeGaussQuadrature(order, -1.0, 1.0, node_coords_, node_weights_);

        // 多项式系数
        ComputeLegendreCoefficients();

        // 预计算矩阵
        ComputeVandermonde();
        ComputeDx();
        ComputeBCExtrap();
    }

    void REFFE::ComputeLegendreCoefficients()
    {
        // 计算Legendre多项式系数
        legendre_coefs_ = MatrixXd::Zero(Np_, Np_);
        numerics::ComputeJacobiCoefficients(order_, 0.0, 0.0, legendre_coefs_);
    }

    void REFFE::ComputeVandermonde()
    {
        Vandermonde_ = MatrixXd::Zero(Np_, Np_);

        for (int i = 0; i < Np_; ++i)
        {
            for (int k = 0; k < Np_; ++k)
            {
                double x = node_coords_(i);
                // double xr = 1.0;
                for (int j = 0; j < Np_; ++j)
                {
                    Vandermonde_(i, j) += legendre_coefs_(k, j) * std::pow(x, k);
                }
            }
        }
        Vandinv_ = Vandermonde_.inverse();
    }

    void REFFE::ComputeDx()
    {
        MatrixXd Vandx = MatrixXd::Zero(Np_, Np_);
        MatrixXd Jacobi11coef = MatrixXd::Zero(Np_, Np_);
        numerics::ComputeJacobiCoefficients(order_, 1.0, 1.0, Jacobi11coef);

        for (int j = 1; j < Np_; j++)
        {
            for (int i = 0; i < Np_; i++)
            {
                double x = node_coords_(i);
                for (int k = 0; k < j; k++)
                {
                    Vandx(i, j) += Jacobi11coef(k, j - 1) * std::pow(x, k);
                }
                Vandx(i, j) *= std::sqrt(j * (j + 1.0));
            }
        }
        /*
        MatrixXd Vand2 = MatrixXd::Zero(Np_, Np_);

        for (int i = 0; i < Np_; i++)
        {
            for (int j=0;j< Np_;j++)
            {
                for(int k=1;k<=j;k++)
                {
                    Vand2(i,j) += k * legendre_coefs_(k,j) * std::pow(node_coords_(i), k-1);
                }
            }
        }*/
        // 直接用系数计算，用于测试 --> 已测试正确
        Dx_ = Vandx * Vandinv_;
    }

    void REFFE::ComputeMass()
    {
        MatrixXd VVT = Vandermonde_ * Vandermonde_.transpose();
        mass_ = VVT.inverse();
    }

    void REFFE::ComputeStiff()
    {
        if (std::abs(mass_(0, 0)) < 0.00000001)
        {
            // 如果未计算质量矩阵，则先计算质量矩阵
            ComputeMass();
        }
        stiff_ = (mass_ * Dx_).transpose();
    }

    void REFFE::ComputeBCExtrap()
    {
        VectorXd left_values = VectorXd::Zero(Np_);
        VectorXd right_values = VectorXd::Zero(Np_);

        for (int i = 0; i < Np_; i++)
        {
            for (int k = 0; k <= i; k++)
            {
                left_values(i) += legendre_coefs_(k, i) * std::pow(-1.0, k);
                right_values(i) += legendre_coefs_(k, i);
            }
        }
        lagrange_left_values_ = Vandinv_.transpose() * left_values;
        lagrange_right_values_ = Vandinv_.transpose() * right_values;
    }

    void REFFE::ComputeLagrangeExtrapolation(double xnew, FE fe, VectorXd& lagrange_values)
    {
        // 计算外插点的Lagrange值
        // 将xnew映射到[-1, 1]区间
        double x_mapped = (xnew - fe.xl_) / fe.dx_ * 2.0 - 1.0; // [-1, 1] 区间映射

        VectorXd values = VectorXd::Zero(Np_);

        for (int i = 0; i < Np_; i++)
        {
            for (int k = 0; k <= i; k++)
            {
                values(i) += legendre_coefs_(k, i) * std::pow(x_mapped, k);
            }
        }
        lagrange_values = VectorXd::Zero(Np_);
        lagrange_values = Vandinv_.transpose() * values;
    }

    void REFFE::ComputeFuncExtrapolation(double xnew, FE fe, VectorXd& func_values, double& fnew)
    {
        VectorXd lagrange_values;
        ComputeLagrangeExtrapolation(xnew, fe, lagrange_values);
        fnew = 0.0;
        for (int i = 0; i < Np_; i++)
        {
            fnew += lagrange_values(i) * func_values(i);
        }
    }

    void REFFE::MapToActualFE(FE &fe)
    {
        // 将参考元映射到实际单元
        double xl = fe.xl_;
        double dx = fe.dx_;
        fe.Np_ = Np_;

        fe.node_coords_ = VectorXd::Zero(Np_);
        fe.node_values_ = VectorXd::Zero(Np_);
        for (int i = 0; i < Np_; i++)
        {
            fe.node_coords_(i) = 0.5 * dx * (node_coords_(i) + 1.0) + xl;
        }
    }

#pragma endregion

#pragma region FEMSpace
    FEMSpace::FEMSpace(femmesh::FEMMesh &mesh, int order)
    {
        mesh_ = &mesh;
        order_ = order;
        Np_ = order + 1; // 插值点数

        // ref_elem_ = new REFFE(order);

        Nx_ = mesh.GetNumCells(); // 单元数目
        // elems_.resize(Nx_);
        mass_.resize(Nx_);
        stiff_.resize(Nx_);
        rhs_.resize(Nx_);
        elem_bc_values_.resize(Nx_);
        quad_weights_.resize(Nx_);

        face_mass_left_ = MatrixXd::Zero(Np_, Np_);
        face_mass_right_ = MatrixXd::Zero(Np_, Np_);

        InitRefElem();

        double xl, xr;
        MatrixXi elem2vert_indices = mesh.GetElem2VertIndices();
        VectorXd vert_coords = mesh.GetVertCoords();
        for (int i = 0; i < Nx_; i++)
        {
            mass_[i] = MatrixXd::Zero(Np_, Np_);
            stiff_[i] = MatrixXd::Zero(Np_, Np_);
            rhs_[i] = VectorXd::Zero(Np_);
            elem_bc_values_[i] = VectorXd::Zero(2);
            xl = vert_coords(elem2vert_indices(i, 0));
            xr = vert_coords(elem2vert_indices(i, 1));
            elems_.push_back(FE(Np_, xl, xr)); // 初始化单元

            quad_weights_[i] = ref_elem_->node_weights_.array() * elems_[i].dx_ * 0.5; // 求积权重
        }
    }

    void FEMSpace::InitRefElem()
    {
        // 初始化参考元
        ref_elem_ = new REFFE(order_);
    }

    void FEMSpace::VolumeMassIntegrator()
    {
        ref_elem_->ComputeMass();
        MatrixXd &refmass = ref_elem_->mass_;
        for (int i = 0; i < Nx_; i++)
        {
            mass_[i] = elems_[i].dx_ * 0.5 * refmass;
        }
    }

    void FEMSpace::VolumeStiffIntegrator()
    {
        ref_elem_->ComputeStiff();
        MatrixXd &refstiff = ref_elem_->stiff_;
        for (int i = 0; i < Nx_; i++)
        {
            stiff_[i] = refstiff;
        }
    }

    void FEMSpace::FaceMassIntegrator()
    {
        VectorXd &lagrange_left_values = ref_elem_->lagrange_left_values_;
        VectorXd &lagrange_right_values = ref_elem_->lagrange_right_values_;
        //std::cout << "Lagrange left values: " << lagrange_left_values.transpose() << std::endl;
        //std::cout << "Lagrange right values: " << lagrange_right_values.transpose() << std::endl;
        for (int i = 0; i < Np_; i++)
        {
            for (int j = 0; j < Np_; j++)
            {
                face_mass_left_(i, j) = lagrange_left_values(i) * lagrange_left_values(j);
                face_mass_right_(i, j) = lagrange_right_values(i) * lagrange_right_values(j);
            }
        }
    }

    void FEMSpace::ComputeBCFuncValues()
    {
        for (int i = 0; i < Nx_; i++)
        {
            double fl = 0.0, fr = 0.0;
            for (int j = 0; j < Np_; j++)
            {
                fl += ref_elem_->lagrange_left_values_(j) * elems_[i].node_values_(j);
                fr += ref_elem_->lagrange_right_values_(j) * elems_[i].node_values_(j);
            }
            elem_bc_values_[i](0) = fl; // 左侧边界值
            elem_bc_values_[i](1) = fr; // 右侧边界值
        } // 有chache优化空间
    }

    void FEMSpace::ComputeFlux(double flux_in, double flux_out)
    {
        // not implemented yet
    }

    void FEMSpace::BCInit()
    {
        left_bc_value_ = 0.0;
        right_bc_value_ = 0.0;
    }

    void FEMSpace::BCChange(double left_bc_value, double right_bc_value)
    {
        left_bc_value_ = left_bc_value;
        right_bc_value_ = right_bc_value;
    }

    void FEMSpace::BCUpdate()
    {
        elem_bc_values_[0](0) = left_bc_value_;
        elem_bc_values_[Nx_ - 1](1) = right_bc_value_;
    }

    /*
    void FEMSpace::Interpolate(double (*func)(double x))
    {
        for(int i=0;i<Nx_;i++)
        {
            for(int j=0;j<Np_;j++)
            {
                elems_[i].node_values_(j) = func(elems_[i].node_coords_(j));
            }
        }
    }*/

    std::vector<VectorXd> FEMSpace::Interpolate(double (*func)(double x))
    {
        std::vector<VectorXd> func_values(Nx_);

        for (int i = 0; i < Nx_; i++)
        {
            func_values[i] = VectorXd::Zero(Np_);
            for (int j = 0; j < Np_; j++)
            {
                func_values[i](j) = func(elems_[i].node_coords_(j));
            }
        }
        return func_values;
    }

    void FEMSpace::Interpolate(double (*func)(double x), double xl, double xr, VectorXd &func_values)
    {
        func_values.resize(Np_);
        FE fe(Np_, xl, xr);
        ref_elem_->MapToActualFE(fe);
        VectorXd &node_coords = fe.node_coords_;
        for (int j = 0; j < Np_; j++)
        {
            func_values(j) = func(node_coords(j));
        }
    }

    void FEMSpace::Quaduature(const VectorXd &func_values, int elem_idx, double &result)
    {
        result = 0.0;
        for (int j = 0; j < Np_; j++)
        {
            result += func_values(j) * quad_weights_[elem_idx](j);
        }
    }

#pragma endregion

}