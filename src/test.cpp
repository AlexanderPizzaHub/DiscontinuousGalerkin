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
        numerics::ComputeJacobiQuadrature(order, 0.0, 0.0, node_coords2, node_weights2);

        std::cout << "Node Coordinates: " << node_coords1.transpose() << std::endl;
        std::cout << "Node Weights: " << node_weights1.transpose() << std::endl;

        std::cout << "Jacobi Node Coordinates: " << node_coords2.transpose() << std::endl;
        std::cout << "Jacobi Node Weights: " << node_weights2.transpose()
                  << std::endl;
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
        for(int i=0;i<Nx+1;i++)
        {
            vertex_coords.push_back(i * 0.1); // Example coordinates
            vert_indices.push_back(i);
        }
        for(int i=0;i<Nx;i++)
        {
            cell_indices.push_back(i);
            elem2vert_indices.push_back({i, i+1}); // Example connectivity
        }
        femmesh::FEMMesh mesh(vertex_coords, vert_indices, cell_indices, elem2vert_indices);

       //MatrixXd coords = mesh.GetVertCoords();
        //std::cout << "Vertex Coordinates: " << std::endl;
        //std::cout << coords << std::endl;

    }

    void TestRefFE()
    {
        std::vector<double> vertex_coords;
        std::vector<int> vert_indices;
        std::vector<int> cell_indices;
        std::vector<std::vector<int>> elem2vert_indices;

        int Nx = 35;
        for(int i=0;i<Nx+1;i++)
        {
            vertex_coords.push_back(i * 0.1); // Example coordinates
            vert_indices.push_back(i);
        }
        for(int i=0;i<Nx;i++)
        {
            cell_indices.push_back(i);
            elem2vert_indices.push_back({i, i+1}); // Example connectivity
        }
        femmesh::FEMMesh mesh(vertex_coords, vert_indices, cell_indices, elem2vert_indices);

        int order =1; // Example order
        femspace::REFFE ref_elem(order);
        ref_elem.ComputeMass();
        ref_elem.ComputeStiff();
        std::cout << "Vandermonde Matrix: " << std::endl;
        std::cout << ref_elem.Vandermonde_ << std::endl;

        std::cout << "Vandermonde determinant: "   <<std::endl;
        //double cond = numerics::ComputeConditionNumber(ref_elem.Vandermonde_);
        std::cout << std::abs(ref_elem.Vandermonde_.determinant()) << std::endl;
        

        std::cout << "Mass Matrix: " << std::endl;
        std::cout << ref_elem.mass_ << std::endl;
        std::cout << "mass sum" << std::endl;
        std::cout << ref_elem.mass_.sum() << std::endl;
        std::cout << "mass determinant" << std::endl;
        std::cout << std::abs(ref_elem.mass_.determinant()) << std::endl;
        std::cout <<"mass condition number: " << std::endl;
        double cond = numerics::ComputeConditionNumber(ref_elem.mass_);
        std::cout << cond << std::endl;
        
       
        
        std::cout << "Dx Matrix: " << std::endl;
        std::cout << ref_elem.Dx_ << std::endl;
        // 刚度阵有问题 --> 已debug
        std::cout << "Stiffness Matrix: " << std::endl;
        std::cout << ref_elem.stiff_ << std::endl;
        
        //std::cout << "skew"<<std::endl;
        //MatrixXd skew = ref_elem.stiff_ + ref_elem.stiff_.transpose();
        //std::cout << skew << std::endl;
        std::cout << "legendre coefficients: " << std::endl;
        std::cout << ref_elem.legendre_coefs_ << std::endl;

        std::cout << "coords: " << std::endl;
        std::cout << ref_elem.node_coords_.transpose() << std::endl;
    }

    double f(double x)
        {
            return x*x;
        };

    void TestElemBC()
    {
        std::vector<double> vertex_coords;
        std::vector<int> vert_indices;
        std::vector<int> cell_indices;
        std::vector<std::vector<int>> elem2vert_indices;

        int Nx = 35;
        for(int i=0;i<Nx+1;i++)
        {
            vertex_coords.push_back(i * 0.1); // Example coordinates
            vert_indices.push_back(i);
        }
        for(int i=0;i<Nx;i++)
        {
            cell_indices.push_back(i);
            elem2vert_indices.push_back({i, i+1}); // Example connectivity
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

        for(int i=0;i<nx;i++)
        {
            std::cout << "Element " << i << " BC Values: " << elem_bc_values[i].transpose() << std::endl;
        }

    }

    double TestQuadrature(int order,double xl, double xr)
    {
        std::cout << "testing quadrature order: " << order << std::endl;
        VectorXd node_coords, node_weights;
        //double xl = -1.0; // Left boundary
        //double xr = 1.0;  // Right boundary
        numerics::ComputeGaussQuadrature(order, xl, xr, node_coords, node_weights);

        for(int poly_order = 0; poly_order <= 2*order+5; poly_order++)
        {
            double integral = 0.0;
            for(int i = 0; i < node_coords.size(); i++)
            {
                integral += node_weights[i] * std::pow(node_coords[i], poly_order);
            }
            std::cout << "Integral of x^" << poly_order<< " is: " << integral << " ";
            std::cout << "Exact value: " << (std::pow(xr, poly_order + 1) - std::pow(xl, poly_order + 1)) / (poly_order + 1) << " ";
            //std::cout << "Error: " << std::abs(integral - (std::pow(xr, poly_order + 1) - std::pow(xl, poly_order + 1)) / (poly_order + 1)) << std::endl;
            std::cout <<" Relative error: " << std::abs(integral - (std::pow(xr, poly_order + 1) - std::pow(xl, poly_order + 1)) / (poly_order + 1)) / std::abs((std::pow(xr, poly_order + 1) - std::pow(xl, poly_order + 1)) / (poly_order + 1)) << std::endl;
        }
    }

    double TestHyperbolic(int order, int Nx)
    {
        std::vector<double> vertex_coords;
        std::vector<int> vert_indices;
        std::vector<int> cell_indices;
        std::vector<std::vector<int>> elem2vert_indices;
        
        for(int i=0;i<Nx+1;i++)
        {
            vertex_coords.push_back(i * (1.0/Nx)); // Example coordinates
            vert_indices.push_back(i);
        }
            
        //vertex_coords = {0.0, 0.05, 0.1, 0.2, 0.3,0.5, 0.7, 0.8, 0.9, 0.95, 1.0};
        for(int i=0;i<Nx;i++)
        {
            cell_indices.push_back(i);
            elem2vert_indices.push_back({i, i+1}); // Example connectivity
        }
        femmesh::FEMMesh mesh(vertex_coords, vert_indices, cell_indices, elem2vert_indices);

        femspace::FEMSpace fem_space(mesh, order);
        //std::cout << "space set" << std::endl; 

        femapplication::Hyperbolic hyperbolic_solver(fem_space, 1.0, 1.0, 1.0, 1.0);
        //std::cout << "solver set" << std::endl; 

        hyperbolic_solver.Prepare();
        //std::cout << "prepared " << std::endl; 
        hyperbolic_solver.Solve();

        std::vector<VectorXd> solution = hyperbolic_solver.GetSolution();
        // for(int i=0;i<solution.size();i++)
        // {
        //     std::cout << "Solution at element " << i << ": " << solution[i].transpose() << std::endl;
        // }
        /*
        VectorXd globalcoords = VectorXd::Zero((Nx)*(order+1));
        for(int i=0;i<Nx;i++)
        {
            for(int j=0;j<order+1;j++)
            {
                globalcoords(i*(order+1)+j) = fem_space.GetElems()[i].node_coords_(j);
            }
        }
        VectorXd globalsolutions = VectorXd::Zero((Nx)*(order+1));
        for(int i=0;i<Nx;i++)
        {
            for(int j=0;j<order+1;j++)
            {
                globalsolutions(i*(order+1)+j) = solution[i](j);
            }
        }
        VectorXd exact_solution = VectorXd::Zero((Nx)*(order+1));
        for(int i=0;i<Nx*(order+1);i++)
        {
            exact_solution(i) = settings::hyperbolic::Truth(globalcoords(i));
        }

        double l2_error = (globalsolutions - exact_solution).norm() / std::sqrt((Nx)*(order+1));
        */
       VectorXd exact_solution = VectorXd::Zero(order+1);
       VectorXd residual = VectorXd::Zero(order+1);
       auto mass = fem_space.GetMass();
       double l2_error = 0.0;
       double elem_error = 0.0;
       for(int i=0;i<Nx;i++)
       {
            elem_error = 0.0;
            for(int j=0;j<order+1;j++)
            {
                exact_solution(j) = settings::hyperbolic::Truth(fem_space.GetElems()[i].node_coords_(j));
            }
            residual = solution[i] - exact_solution;
            elem_error = residual.dot(mass[i]*residual);
            l2_error += elem_error;
       }
        l2_error = std::sqrt(l2_error);
        //std::cout << "L2 Error: " << l2_error << std::endl;

        return l2_error;

    }
    
}