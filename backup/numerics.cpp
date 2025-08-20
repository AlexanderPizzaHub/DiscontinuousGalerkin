#include "numerics.hpp"

namespace numerics
{
    using namespace Eigen;

    // void ComputeGaussQuadrature(int order, double xl, double xr, VectorXd &node_coords, VectorXd &node_weights)
    // {
    //     // TEST CORRECT
    //     int Np = order + 1; // Gauss-Legendre求积点数
    //     node_coords.resize(Np);
    //     node_weights.resize(Np);

    //     VectorXd y, y0, Lp;
    //     y = VectorXd::Zero(Np);
    //     y0 = VectorXd::Zero(Np);
    //     Lp = VectorXd::Zero(Np);

    //     MatrixXd L = MatrixXd::Zero(Np, Np + 1);

    //     for (int i = 0; i < Np; i++)
    //     {
    //         y[i] = std::cos((2.0 * i + 1.0) * 4.0 * std::atan(1.0) / (2.0 * (Np - 1.0) + 2.0)) +
    //                0.27 / Np * std::sin(4.0 * std::atan(1.0) * (-1.0 + i * 2.0 / (Np - 1.0)) * (Np - 1.0) / (Np + 1));
    //         y0[i] = 2.0;
    //     }

    //     // Compute the zeros of the Np+1 Legendre Polynomial using recursion and Newton's method
    //     while (true)
    //     {
    //         for (int i = 0; i < Np; ++i)
    //         {
    //             L(i, 0) = 1.0;
    //             L(i, 1) = y[i];
    //         }

    //         for (int K = 2; K < Np; ++K)
    //         {
    //             for (int i = 0; i < Np; ++i)
    //             {
    //                 L(i, K) = ((2.0 * K - 1.0) * y[i] * L(i, K - 1) - (K - 1) * L(i, K - 2)) / K;
    //             }
    //         }

    //         for (int i = 0; i < Np; ++i)
    //         {
    //             Lp[i] = (Np + 1) * (L(i, Np - 1) - y[i] * L(i, Np)) / (1.0 - y[i] * y[i]);
    //         }

    //         y0 = y;
    //         for (int i = 0; i < Np; ++i)
    //         {
    //             y[i] = y0[i] - L(i, Np) / Lp[i];
    //         }

    //         double max_diff = 0.0;
    //         for (int i = 0; i < Np; ++i)
    //         {
    //             max_diff = std::max(max_diff, std::abs(y[i] - y0[i]));
    //         }

    //         if (max_diff < 1.0e-13)
    //         {
    //             break;
    //         }
    //     }

    //     // Linear map from [-1, 1] to [a, b]
    //     for (int i = 0; i < Np; ++i)
    //     {
    //         node_coords[i] = (xl * (1.0 - y[i]) + xr * (1.0 + y[i])) / 2.0;
    //         node_weights[i] = (Np + 1) * (Np + 1) * (xr - xl) / ((1.0 - y[i] * y[i]) * Lp[i] * Lp[i]) / (Np * Np);
    //     }
    // }

    void ComputeGaussQuadrature(int order, double xl, double xr, VectorXd &node_coords, VectorXd &node_weights) {
        if (order <= 0) {
            throw std::invalid_argument("order must be positive");
        }
        int n = order + 1; // Gauss-Legendre求积点数
    
        // 创建对称三对角矩阵
        MatrixXd T = MatrixXd::Zero(n, n);
        for (int i = 1; i < n; ++i) {
            double d = i / std::sqrt(4.0 * i * i - 1.0);
            T(i, i-1) = d;
            T(i-1, i) = d;
        }
    
        // 计算特征值和特征向量
        SelfAdjointEigenSolver<MatrixXd> solver(T);
        VectorXd eigenvalues = solver.eigenvalues();
        MatrixXd eigenvectors = solver.eigenvectors();
    
        // 求积点为特征值
        node_coords = (eigenvalues.array())* (xr - xl) / 2.0 + (xr + xl) / 2.0;
    
        // 权重为第一行特征向量的平方乘以2
        node_weights = (2.0 * eigenvectors.row(0).array().square())* (xr - xl) / 2.0;
    }

    void ComputeJacobiQuadrature(int order, double alpha, double beta, VectorXd &node_coords, VectorXd &node_weights)
    {
        // TEST NOT CORRECT

        int Np = order + 1;

        if (Np == 1)
        {
            node_coords.resize(1);
            node_weights.resize(1);
            node_coords(0) = (alpha - beta) / (alpha + beta + 2.0);
            return;
        }

        MatrixXd J = MatrixXd::Zero(Np, Np);
        VectorXd h1 = 2.0 * VectorXd::LinSpaced(Np, 0, Np - 1).array() + alpha + beta;
        // std::cout << h1.transpose() << std::endl;

        for (int i = 0; i < Np; ++i)
        {
            if (i < Np - 1)
            {
                // double a = 2.0 / ((h1(i) + 2.0) * std::sqrt((i + 1.0) * (i + 1.0 + alpha + beta) * (i + 1.0 + alpha) * (i + 1.0 + beta) / ((h1(i) + 1.0) * (h1(i) + 3.0))));
                double a = 2.0 / (h1(i + 1) + 2.0) * std::sqrt((i + 1.0) * (i + 1.0 + alpha + beta) * (i + 1.0 + alpha) * (i + 1.0 + beta) / ((h1(i + 1) + 1.0) * (h1(i + 1) + 3.0)));
                J(i, i + 1) = a;
                J(i + 1, i) = a;
            }
            // J(i, i) = -0.5 * (alpha * alpha - beta * beta) / ((h1(i) + 2.0) * h1(i));
            J(i, i) = -(alpha * alpha - beta * beta) / ((h1(i) + 2.0) * h1(i));
        }

        if (std::abs(alpha + beta) < 10 * std::numeric_limits<double>::epsilon())
        {
            J(0, 0) = 0.0;
        }
        //std::cout << J << std::endl;

        // Compute quadrature by eigenvalue solve
        SelfAdjointEigenSolver<MatrixXd> solver(J);
        VectorXd x = solver.eigenvalues();
        MatrixXd V = solver.eigenvectors();

        node_coords = x;
        node_weights = V.row(0).array().square() * std::pow(2.0, alpha + beta + 1.0) / (alpha + beta + 1.0) * std::tgamma(alpha + 1.0) * std::tgamma(beta + 1.0) / std::tgamma(alpha + beta + 1.0);
        return;
    }

    void ComputeJacobiCoefficients(int order, double alpha, double beta, MatrixXd &coeffs)
    {
        // TEST CORRECT
        int Np = order + 1;

        coeffs = MatrixXd::Zero(Np, Np);
        double tmp = 0.0;

        coeffs(0, 0) = std::sqrt(
            std::pow(2.0, -alpha - beta - 1.0) * std::tgamma(alpha + beta + 2.0) / (std::tgamma(beta + 1.0) * std::tgamma(alpha + 1.0)));
        if (Np == 1)
        {
            return;
        }

        //coeffs(0, 0) = std::sqrt(
        //    std::pow(2.0, -alpha - beta - 1.0) * std::tgamma(alpha + beta + 2.0) / (std::tgamma(beta + 1.0) * std::tgamma(alpha + 1.0)));

        tmp = 0.5 * coeffs(0, 0) * std::sqrt((alpha + beta + 3.0) / ((alpha + 1.0) * (beta + 1.0)));
        //std::cout << tmp << std::endl;
        coeffs(0, 1) = tmp * (alpha - beta);
        coeffs(1, 1) = tmp * (alpha + beta + 2.0);

        if (Np == 2)
        {
            return;
        }

        double an, anprev, bnprev;

        anprev = 2.0 / (2.0 + alpha + beta) * std::sqrt((1.0 + alpha + beta) * (1.0 + alpha) * (1.0 + beta) / ((1.0 + alpha + beta) * (3.0 + alpha + beta)));

        for (int n = 2; n < Np; n++)
        {
            bnprev = -(alpha * alpha - beta * beta) / ((2.0 * n - 2.0 + alpha + beta) * (2.0 * n + alpha + beta));

            an = 2.0 / (2.0 * n + alpha + beta) * std::sqrt(n * (n + alpha + beta) * (n + alpha) * (n + beta) / ((2.0 * n - 1.0 + alpha + beta) * (2.0 * n + 1.0 + alpha + beta)));

            coeffs(0, n) = -bnprev / an * coeffs(0, n - 1) - anprev / an * coeffs(0, n - 2);

            for (int r = 1; r <= n - 2; r++)
            {
                coeffs(r, n) = coeffs(r - 1, n - 1) / an - bnprev / an * coeffs(r, n - 1) - anprev / an * coeffs(r, n - 2);
            }

            coeffs(n - 1, n) = coeffs(n - 2, n - 1) / an - bnprev / an * coeffs(n - 1, n - 1);
            coeffs(n, n) = coeffs(n - 1, n - 1) / an;

            anprev = an;
        }
        //std::cout << coeffs << std::endl;
        return;
    }

    double ComputeConditionNumber(const Eigen::MatrixXd& A) {
        // 使用奇异值分解
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
        // 获取奇异值
        double max_singular = svd.singularValues()(0); // 最大奇异值
        double min_singular = svd.singularValues()(svd.singularValues().size() - 1); // 最小奇异值
    
        // 如果最小奇异值为 0，条件数为无穷大
        if (min_singular == 0) {
            return std::numeric_limits<double>::infinity();
        }
    
        // 返回条件数
        return max_singular / min_singular;
    }

}