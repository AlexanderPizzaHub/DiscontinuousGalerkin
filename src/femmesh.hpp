/*
输入：顶点坐标，顶点索引，单元索引，顶点-单元邻接
1. 生成邻接信息
2. 计算几何信息
*/

#pragma once 

#include "numerics.hpp"

#include "../eigen-3.4.0/Eigen/Dense"

#include <vector>

namespace femmesh
{
    using namespace Eigen;

    class FEMMesh
    {
        public:
            FEMMesh(std::vector<double> vertex_coords,
                    std::vector<int> vert_indices, 
                    std::vector<int> cell_indices,
                    std::vector<std::vector<int>> elem2vert_indices);
            ~FEMMesh() = default;

            VectorXi GetVertIndices() const { return vert_indices_; }
            VectorXi GetCellIndices() const { return cell_indices_; }
            MatrixXd GetVertCoords() const { return vert_coords_; }
            MatrixXi GetElem2VertIndices() const { return elem2vert_indices_; }
            int GetNumCells() const { return Nx_; }
            int GetNumVertices() const { return Nx_ + 1; }


        private:
            /*
            直接数据
            */
            int Nx_; // 单元数目

            MatrixXd vert_coords_; // [Nx+1,1] 顶点坐标 

            VectorXi vert_indices_; //[Nx+1] 顶点索引 
            VectorXi cell_indices_; // [Nx] 单元索引 
            MatrixXi elem2vert_indices_; // [Nx, 2] 一维下，单元-顶点邻接 = 单元-面邻接

            /*
            间接数据
            */
            VectorXd elem_areas_; // [Nx, 1] 单元面积

            MatrixXd face_normals_; // [Nx, 2] 面法向量

            //void ComputeAdajacency();
            void ComputeGeometry();
            
    };
}