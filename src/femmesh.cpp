#include "femmesh.hpp"

namespace femmesh
{
    using namespace Eigen;

    FEMMesh::FEMMesh(std::vector<double> vertex_coords,
                     std::vector<int> vert_indices, 
                     std::vector<int> cell_indices,
                     std::vector<std::vector<int>> elem2vert_indices)
    {
        Nx_ = cell_indices.size(); // Number of cells
        // Convert input vectors to Eigen types

        vert_coords_.resize(Nx_ + 1,1); // Nx+1 vertices for 1D mesh
        vert_indices_.resize(Nx_ + 1);
        cell_indices_.resize(Nx_);
        elem2vert_indices_.resize(Nx_, 2); // Assuming 1D elements with 2 vertices each
        
        for(int i=0;i<=Nx_;i++)
        {
            vert_coords_(i,0) = vertex_coords[i];
            vert_indices_(i) = vert_indices[i];
        }

        for(int i=0 ;i<Nx_;i++)
        {
            cell_indices_[i] = cell_indices[i];
            elem2vert_indices_(i,0) = elem2vert_indices[i][0];
            elem2vert_indices_(i,1) = elem2vert_indices[i][1];
        }

        // Compute geometry information
        ComputeGeometry();
    }

    void FEMMesh::ComputeGeometry()
    {
        // Compute areas and normals for each element
        elem_areas_.resize(Nx_);
        face_normals_.resize(Nx_, 2); // Assuming 2D elements
        double xl, xr;
        for (int i = 0; i < Nx_; ++i) {
            xl = vert_coords_(elem2vert_indices_(i, 0),0);
            xr = vert_coords_(elem2vert_indices_(i, 1),0);
            elem_areas_(i) = std::abs(xr-xl); // Length of the edge in 1D case
            
            // Compute normal (perpendicular vector)
            face_normals_(i, 0) = -1.0; // left
            face_normals_(i, 1) = 1.0 ; // right
        }
    }
}