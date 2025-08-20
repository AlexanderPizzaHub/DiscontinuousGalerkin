#include "numerics.hpp"

namespace numerics
{
    using namespace Eigen;

    void ComputeGaussQuadrature(int order, double xl, double xr, VectorXd &node_coords, VectorXd &node_weights)
    {
        if (order <= 0)
        {
            throw std::invalid_argument("order must be positive");
        }
        int n = order + 1; // Gauss-Legendre求积点数

        // 创建对称三对角矩阵
        MatrixXd T = MatrixXd::Zero(n, n);
        for (int i = 1; i < n; ++i)
        {
            double d = i / std::sqrt(4.0 * i * i - 1.0);
            T(i, i - 1) = d;
            T(i - 1, i) = d;
        }

        // 计算特征值和特征向量
        SelfAdjointEigenSolver<MatrixXd> solver(T);
        VectorXd eigenvalues = solver.eigenvalues();
        MatrixXd eigenvectors = solver.eigenvectors();

        // 求积点为特征值
        node_coords = (eigenvalues.array()) * (xr - xl) / 2.0 + (xr + xl) / 2.0;

        // 权重为第一行特征向量的平方乘以2
        node_weights = (2.0 * eigenvectors.row(0).array().square()) * (xr - xl) / 2.0;
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

        // coeffs(0, 0) = std::sqrt(
        //     std::pow(2.0, -alpha - beta - 1.0) * std::tgamma(alpha + beta + 2.0) / (std::tgamma(beta + 1.0) * std::tgamma(alpha + 1.0)));

        tmp = 0.5 * coeffs(0, 0) * std::sqrt((alpha + beta + 3.0) / ((alpha + 1.0) * (beta + 1.0)));
        // std::cout << tmp << std::endl;
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
         // std::cout << coeffs << std::endl;
        return;
    }

    void ComputeJacobiCoefficients_deepseek(int order, double alpha, double beta, MatrixXd &coeffs)
    {
        // 计算Jacobi多项式的系数矩阵
        // 返回的矩阵大小为 (n+1) x (n+1)，第j列存储j阶Jacobi多项式的系数
            if (order <= 0)
            {
                throw std::invalid_argument("order must be positive");
            }
            if (alpha <= -1.0 || beta <= -1.0)
            {
                throw std::invalid_argument("alpha and beta must be > -1");
            }
            int n = order + 1;
   
            coeffs = MatrixXd::Zero(n, n);

            // 0阶Jacobi多项式: P_0^{(α,β)}(x) = 1
            coeffs(0, 0) = 1.0;

            if (order >= 1)
            {
                // 1阶Jacobi多项式: P_1^{(α,β)}(x) = (α - β + (α + β + 2)x)/2
                coeffs(0, 1) = (alpha - beta) / 2.0;
                coeffs(1, 1) = (alpha + beta + 2.0) / 2.0;
                if(order == 1){return;}
            }

            
            // 使用递推关系计算高阶Jacobi多项式系数
            for (int j = 2; j < n; ++j)
            {
                // 递推系数
                double a1 = 2.0 * j * (j + alpha + beta) * (2.0 * j + alpha + beta - 2.0);
                double a2 = (2.0 * j + alpha + beta - 1.0) * (alpha * alpha - beta * beta);
                double a3 = (2.0 * j + alpha + beta - 2.0) * (2.0 * j + alpha + beta - 1.0) * (2.0 * j + alpha + beta);
                double a4 = 2.0 * (j + alpha - 1.0) * (j + beta - 1.0) * (2.0 * j + alpha + beta);

                // 递推关系:
                // P_j^{(α,β)}(x) = [(a2 + a3*x) * P_{j-1}^{(α,β)}(x) - a4 * P_{j-2}^{(α,β)}(x)] / a1

                // 计算 (a2 + a3*x) * P_{j-1}^{(α,β)}(x)
                VectorXd term1 = VectorXd::Zero(j + 1);
                for (int i = 0; i < j; ++i)
                {
                    term1(i) += (a2 / a1) * coeffs(i, j - 1);     // a2 * P_{j-1}
                    term1(i + 1) += (a3 / a1) * coeffs(i, j - 1); // a3 * x * P_{j-1}
                }

                // 计算 a4 * P_{j-2}^{(α,β)}(x)
                VectorXd term2 = VectorXd::Zero(j + 1);
                for (int i = 0; i < j - 1; ++i)
                {
                    term2(i) += (a4 / a1) * coeffs(i, j - 2);
                }

                // 组合得到P_j^{(α,β)}(x)
                for (int i = 0; i <= j; ++i)
                {
                    coeffs(i, j) = term1(i) - term2(i);
                }
            }
            return;
        }

        double ComputeConditionNumber(const Eigen::MatrixXd &A)
        {
            // 使用奇异值分解
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

            // 获取奇异值
            double max_singular = svd.singularValues()(0);                               // 最大奇异值
            double min_singular = svd.singularValues()(svd.singularValues().size() - 1); // 最小奇异值

            // 如果最小奇异值为 0，条件数为无穷大
            if (min_singular == 0)
            {
                return std::numeric_limits<double>::infinity();
            }

            // 返回条件数
            return max_singular / min_singular;
        }

        MatrixXd ComputeJacobiMatrix(int n, double alpha, double beta)
        {
            MatrixXd J = MatrixXd::Zero(n, n);

            for (int i = 0; i < n; ++i)
            {
                // 对角元素
                double a = (beta * beta - alpha * alpha) /
                           ((2.0 * i + alpha + beta) * (2.0 * i + alpha + beta + 2.0));
                J(i, i) = a;

                // 非对角元素
                if (i < n - 1)
                {
                    double b = 2.0 / (2.0 * i + alpha + beta + 2.0) *
                               std::sqrt((i + 1.0) * (i + alpha + 1.0) *
                                         (i + beta + 1.0) * (i + alpha + beta + 1.0) /
                                         ((2.0 * i + alpha + beta + 1.0) *
                                          (2.0 * i + alpha + beta + 3.0)));
                    J(i, i + 1) = b;
                    J(i + 1, i) = b;
                }
            }

            return J;
        }

        // 计算Gamma函数近似值
        double gamma_function(double x)
        {
            // Lanczos近似公式
            const double g = 7.0;
            const double p[] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                                771.32342877765313, -176.61502916214059, 12.507343278686905,
                                -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

            if (x < 0.5)
            {
                return M_PI / (std::sin(M_PI * x) * gamma_function(1.0 - x));
            }

            x -= 1.0;
            double a = p[0];
            double t = x + g + 0.5;

            for (int i = 1; i < 9; ++i)
            {
                a += p[i] / (x + i);
            }

            return std::sqrt(2.0 * M_PI) * std::pow(t, x + 0.5) * std::exp(-t) * a;
        }

        // 计算Jacobi求积点和权重
        // void ComputeJacobiQuadrature(int order, double alpha, double beta, VectorXd &node_coords, VectorXd &node_weights)
        // {
        //     if (order <= 0)
        //     {
        //         throw std::invalid_argument("n must be positive");
        //     }
        //     if (alpha <= -1.0 || beta <= -1.0)
        //     {
        //         throw std::invalid_argument("alpha and beta must be > -1");
        //     }

        //     int n = order + 1; // Jacobi求积点数
        //     // 计算Jacobi矩阵
        //     MatrixXd J = ComputeJacobiMatrix(n, alpha, beta);

        //     // 计算特征值和特征向量
        //     SelfAdjointEigenSolver<MatrixXd> solver(J);
        //     VectorXd eigenvalues = solver.eigenvalues();
        //     MatrixXd eigenvectors = solver.eigenvectors();

        //     // 求积点为特征值
        //     node_coords = eigenvalues;

        //     // 计算权重
        //     node_weights.resize(n);

        //     // 计算归一化常数
        //     double mu0 = std::pow(2.0, alpha + beta + 1.0) *
        //                  gamma_function(alpha + 1.0) *
        //                  gamma_function(beta + 1.0) /
        //                  gamma_function(alpha + beta + 2.0);

        //     // 权重计算
        //     for (int i = 0; i < n; ++i)
        //     {
        //         double phi0 = eigenvectors(0, i);
        //         node_weights(i) = mu0 * phi0 * phi0;
        //     }
        // }
    }