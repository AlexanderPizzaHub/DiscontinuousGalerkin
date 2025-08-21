#include "test.hpp"

namespace test
{
    using namespace Eigen;
    void TestQuadrature()
    {
        // Example test for Gauss Quadrature
        int order = 2;    // Example order
        double xl = -1.0; // Left boundary
        double xr = 1.0;  // Right boundary
        VectorXd node_coords1, node_weights1;
        VectorXd node_coords2, node_weights2;

        numerics::ComputeGaussQuadrature(order, xl, xr, node_coords1, node_weights1);
        // numerics::ComputeJacobiQuadrature(order, 0.0, 0.0, node_coords2, node_weights2);

        std::cout << "Node Coordinates: " << node_coords1.transpose() << std::endl;
        std::cout << "Node Weights: " << node_weights1.transpose() << std::endl;

        // std::cout << "Jacobi Node Coordinates: " << node_coords2.transpose() << std::endl;
        // std::cout << "Jacobi Node Weights: " << node_weights2.transpose()
        //<< std::endl;
    }

    void TestJacobiCoef()
    {
        int order = 2;
        MatrixXd J;

        numerics::ComputeJacobiCoefficients(order, 0.0, 0.0, J);

        std::cout << "Jacobi Coefficients: " << std::endl;
        std::cout << J << std::endl;
    }

    void TestMesh()
    {
        std::vector<double> vertex_coords;
        std::vector<int> vert_indices;
        std::vector<int> cell_indices;
        std::vector<std::vector<int>> elem2vert_indices;

        int Nx = 10;
        for (int i = 0; i < Nx + 1; i++)
        {
            vertex_coords.push_back(i * 0.1); // Example coordinates
            vert_indices.push_back(i);
        }
        for (int i = 0; i < Nx; i++)
        {
            cell_indices.push_back(i);
            elem2vert_indices.push_back({i, i + 1}); // Example connectivity
        }
        femmesh::FEMMesh mesh(vertex_coords, vert_indices, cell_indices, elem2vert_indices);

        // MatrixXd coords = mesh.GetVertCoords();
        // std::cout << "Vertex Coordinates: " << std::endl;
        // std::cout << coords << std::endl;
    }

    void TestRefFE()
    {
        std::vector<double> vertex_coords;
        std::vector<int> vert_indices;
        std::vector<int> cell_indices;
        std::vector<std::vector<int>> elem2vert_indices;

        int Nx = 35;
        for (int i = 0; i < Nx + 1; i++)
        {
            vertex_coords.push_back(i * 0.1); // Example coordinates
            vert_indices.push_back(i);
        }
        for (int i = 0; i < Nx; i++)
        {
            cell_indices.push_back(i);
            elem2vert_indices.push_back({i, i + 1}); // Example connectivity
        }
        femmesh::FEMMesh mesh(vertex_coords, vert_indices, cell_indices, elem2vert_indices);

        int order = 1; // Example order
        femspace::REFFE ref_elem(order);
        ref_elem.ComputeMass();
        ref_elem.ComputeStiff();
        std::cout << "Vandermonde Matrix: " << std::endl;
        std::cout << ref_elem.Vandermonde_ << std::endl;

        std::cout << "Vandermonde determinant: " << std::endl;
        // double cond = numerics::ComputeConditionNumber(ref_elem.Vandermonde_);
        std::cout << std::abs(ref_elem.Vandermonde_.determinant()) << std::endl;

        std::cout << "Mass Matrix: " << std::endl;
        std::cout << ref_elem.mass_ << std::endl;
        std::cout << "mass sum" << std::endl;
        std::cout << ref_elem.mass_.sum() << std::endl;
        std::cout << "mass determinant" << std::endl;
        std::cout << std::abs(ref_elem.mass_.determinant()) << std::endl;
        std::cout << "mass condition number: " << std::endl;
        double cond = numerics::ComputeConditionNumber(ref_elem.mass_);
        std::cout << cond << std::endl;

        std::cout << "Dx Matrix: " << std::endl;
        std::cout << ref_elem.Dx_ << std::endl;
        // 刚度阵有问题 --> 已debug
        std::cout << "Stiffness Matrix: " << std::endl;
        std::cout << ref_elem.stiff_ << std::endl;

        // std::cout << "skew"<<std::endl;
        // MatrixXd skew = ref_elem.stiff_ + ref_elem.stiff_.transpose();
        // std::cout << skew << std::endl;
        std::cout << "legendre coefficients: " << std::endl;
        std::cout << ref_elem.legendre_coefs_ << std::endl;

        std::cout << "coords: " << std::endl;
        std::cout << ref_elem.node_coords_.transpose() << std::endl;
    }

    double f(double x)
    {
        return x * x;
    };

    void TestElemBC()
    {
        std::vector<double> vertex_coords;
        std::vector<int> vert_indices;
        std::vector<int> cell_indices;
        std::vector<std::vector<int>> elem2vert_indices;

        int Nx = 35;
        for (int i = 0; i < Nx + 1; i++)
        {
            vertex_coords.push_back(i * 0.1); // Example coordinates
            vert_indices.push_back(i);
        }
        for (int i = 0; i < Nx; i++)
        {
            cell_indices.push_back(i);
            elem2vert_indices.push_back({i, i + 1}); // Example connectivity
        }
        std::cout << "!!! " << std::endl;
        femmesh::FEMMesh mesh(vertex_coords, vert_indices, cell_indices, elem2vert_indices);
        std::cout << "mesh set" << std::endl;

        int order = 4; // Example order
        femspace::FEMSpace fem_space(mesh, order);
        std::cout << "space set" << std::endl;

        std::vector<VectorXd> inter = fem_space.Interpolate(f);

        std::cout << "interpolated" << std::endl;

        fem_space.ComputeBCFuncValues();

        std::cout << "BC values computed" << std::endl;

        auto nx = fem_space.GetNx();
        auto elem_bc_values = fem_space.GetElemBCValues();

        for (int i = 0; i < nx; i++)
        {
            std::cout << "Element " << i << " BC Values: " << elem_bc_values[i].transpose() << std::endl;
        }
    }

    double TestQuadrature(int order, double xl, double xr)
    {
        std::cout << "testing quadrature order: " << order << std::endl;
        VectorXd node_coords, node_weights;
        // double xl = -1.0; // Left boundary
        // double xr = 1.0;  // Right boundary
        numerics::ComputeGaussQuadrature(order, xl, xr, node_coords, node_weights);

        for (int poly_order = 0; poly_order <= 2 * order + 5; poly_order++)
        {
            double integral = 0.0;
            for (int i = 0; i < node_coords.size(); i++)
            {
                integral += node_weights[i] * std::pow(node_coords[i], poly_order);
            }
            std::cout << "Integral of x^" << poly_order << " is: " << integral << " ";
            std::cout << "Exact value: " << (std::pow(xr, poly_order + 1) - std::pow(xl, poly_order + 1)) / (poly_order + 1) << " ";
            // std::cout << "Error: " << std::abs(integral - (std::pow(xr, poly_order + 1) - std::pow(xl, poly_order + 1)) / (poly_order + 1)) << std::endl;
            std::cout << " Relative error: " << std::abs(integral - (std::pow(xr, poly_order + 1) - std::pow(xl, poly_order + 1)) / (poly_order + 1)) / std::abs((std::pow(xr, poly_order + 1) - std::pow(xl, poly_order + 1)) / (poly_order + 1)) << std::endl;
        }
    }

    void TestApproximation(int order, double xl, double xr)
    {
        femspace::REFFE ref_elem(order);
        femspace::FE fe(order, xl, xr);
        ref_elem.MapToActualFE(fe);

        VectorXd node_coords = fe.node_coords_;
        VectorXd node_weights = fe.dx_ * 0.5 * ref_elem.node_weights_;

        VectorXd func_values = VectorXd::Zero(order + 1);
        for (int i = 0; i < order + 1; i++)
        {
            // func_values(i) = std::exp(node_coords(i));
            func_values(i) = std::sin(M_PI * node_coords(i));
        }

        // compute error on fine mesh
        femspace::REFFE ref_elem_fine(order * 2);
        femspace::FE fe_fine(order * 2, xl, xr);
        ref_elem_fine.MapToActualFE(fe_fine);
        VectorXd fine_node_coords = fe_fine.node_coords_;

        VectorXd fine_func_values = VectorXd::Zero(order * 2 + 1);
        VectorXd extrapolated_values = VectorXd::Zero(order * 2 + 1);
        for (int i = 0; i < order * 2 + 1; i++)
        {
            fine_func_values(i) = std::sin(M_PI * fine_node_coords(i));
            ref_elem.ComputeFuncExtrapolation(fine_node_coords(i), fe, func_values, extrapolated_values(i));
        }
        double l2_error = (extrapolated_values - fine_func_values).norm() / std::sqrt(order * 2 + 1);
        std::cout << "L2 Error of extrapolation: " << l2_error << std::endl;
    }

    void TestProjection(int order, double xl, double xr)
    {
        femspace::REFFE ref_elem(order);
        femspace::FE fe(order, xl, xr);
        ref_elem.ComputeMass();
        ref_elem.ComputeStiff();

        ref_elem.MapToActualFE(fe);

        VectorXd node_coords, node_weights;
        node_coords = fe.node_coords_;
        node_weights = fe.dx_ * 0.5 * ref_elem.node_weights_;

        VectorXd func_values = VectorXd::Zero(order + 1);
        for (int i = 0; i < order + 1; i++)
        {
            // func_values(i) = std::exp(node_coords(i));
            func_values(i) = std::sin(M_PI * node_coords(i));
        }

        MatrixXd ProjectionMatrix = 0.5 * fe.dx_ * ref_elem.mass_;

        VectorXd rhs = VectorXd::Zero(order + 1);
        for (int i = 0; i < order + 1; i++)
        {
            rhs(i) = node_weights(i) * func_values(i);
        }
        VectorXd coeffs = ProjectionMatrix.colPivHouseholderQr().solve(rhs);
        std::cout << "Projection Coefficients: " << coeffs.transpose() << std::endl;
        std::cout << "Exact Coefficients: " << func_values.transpose() << std::endl;

        double l2_error = (coeffs - func_values).norm() / std::sqrt(order + 1);
        std::cout << "L2 Error: " << l2_error << std::endl;
    }

    double TestHyperbolic(int order, int Nx)
    {
        std::vector<double> vertex_coords;
        std::vector<int> vert_indices;
        std::vector<int> cell_indices;
        std::vector<std::vector<int>> elem2vert_indices;

        for (int i = 0; i < Nx + 1; i++)
        {
            vertex_coords.push_back(i * (1.0 / Nx)); // Example coordinates
            vert_indices.push_back(i);
        }

        // vertex_coords = {0.0, 0.05, 0.1, 0.2, 0.3,0.5, 0.7, 0.8, 0.9, 0.95, 1.0};
        for (int i = 0; i < Nx; i++)
        {
            cell_indices.push_back(i);
            elem2vert_indices.push_back({i, i + 1}); // Example connectivity
        }

        femmesh::FEMMesh mesh(vertex_coords, vert_indices, cell_indices, elem2vert_indices);

        femspace::FEMSpace fem_space(mesh, order);

        std::vector<VectorXd> G_values = fem_space.Interpolate(settings::hyperbolic::G);
        std::vector<VectorXd> S_values = fem_space.Interpolate(settings::hyperbolic::S);

        MatrixXd G = MatrixXd::Zero(Nx, order + 1);
        MatrixXd S = MatrixXd::Zero(Nx, order + 1);
        for (int i = 0; i < Nx; i++)
        {
            G.row(i) = G_values[i].transpose();
            S.row(i) = S_values[i].transpose();
        }
        MatrixXd solutionmat = MatrixXd::Zero(Nx, order + 1);

        

        femapplication::Hyperbolic hyperbolic_solver(fem_space, 1.0, 1.0, 1.0, 1.0);

        hyperbolic_solver.Prepare();

        hyperbolic_solver.Solve(G,S,solutionmat);
        

        //std::vector<VectorXd> solution;// = hyperbolic_solver.GetSolution();

        VectorXd approx_solution = VectorXd::Zero(Nx);
        VectorXd exact_solution = VectorXd::Zero(Nx);
        VectorXd misfit = VectorXd::Zero(order + 1);

        // 函数值统一插值到区间中心计算
        double center_x;
        for (int i = 0; i < Nx; i++)
        {
            femspace::FE fe = fem_space.GetElems()[i];
            femspace::REFFE *ref_elem = fem_space.GetRefElem();

            center_x = (fe.xl_ + fe.xr_) / 2.0;
            VectorXd solutionrow = solutionmat.row(i);
            ref_elem->ComputeFuncExtrapolation(center_x, fe, solutionrow,approx_solution(i));

            //std::cout << "Element solution at " << center_x << ": " << solution[i].transpose()<< std::endl;
            //std::cout << "Corresponding coords: "<< fe.node_coords_.transpose() << std::endl;
            
            exact_solution(i) = settings::hyperbolic::Truth(center_x);
        }

        double l2_error = (approx_solution - exact_solution).norm() / std::sqrt(Nx);
        //std::cout << "Approximation: " << std::endl;
        // std::cout << approx_solution.transpose() << std::endl;
        // std::cout << "Exact Solution: " << std::endl;
        // std::cout << exact_solution.transpose() << std::endl;

        return l2_error;
    }

    void TestShakhov1D1V(int order)
    {
        std::vector<double> vertex_coords;
        std::vector<int> vert_indices;
        std::vector<int> cell_indices;
        std::vector<std::vector<int>> elem2vert_indices;

        int Nx = settings::shakhov::Nx; // Number of elements
        double xl = settings::shakhov::xl; // Left boundary
        double xr = settings::shakhov::xr; // Right boundary

        int Nv = settings::shakhov::Nv; // Number of velocity points
        double v_left = -settings::shakhov::vm; // Left velocity boundary
        double v_right = settings::shakhov::vm;

        std::vector<double> v_vec(Nv);

        for (int i = 0; i < Nx + 1; i++)
        {
            vertex_coords.push_back(i * ((xr-xl) / Nx)); // Example coordinates
            vert_indices.push_back(i);
        }

        for (int i = 0; i <= Nv; i++)
        {
            v_vec[i] = v_left + i * (v_right - v_left) / (Nv);
        }   

        for (int i = 0; i < Nx; i++)
        {
            cell_indices.push_back(i);
            elem2vert_indices.push_back({i, i + 1}); // Example connectivity
        }
        femmesh::FEMMesh mesh(vertex_coords, vert_indices, cell_indices, elem2vert_indices);

        femspace::FEMSpace fem_space(mesh, order);

        femapplication::Shakhov1D1V shakhov_solver(v_vec,fem_space);

        shakhov_solver.Initialize();
        // std::cout << "shakhov solver initialized" << std::endl;
        shakhov_solver.Solve(-1.0,10);
        // std::cout << "shakhov solver solved" << std::endl;

        /*
        现问题：出现震荡，没有收敛。
        考虑：把旧版shakhov的第一步density输出来，和现在的结果对比。看看是哪里出了问题
        
        */
    }

}