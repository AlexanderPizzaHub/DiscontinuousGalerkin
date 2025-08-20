/*
FE 有限元
输入：左右侧顶点坐标
1. Lagrange插值点坐标 // Gauss Lobatto
2. 插值函数值

REFFE 参考元
输入：阶数
1. 根据fem阶数，生成求积阶数，求积点，求积式
2. 生成多项式系数，插值点坐标（Gauss-Lobatto）
3. 根据插值点坐标，预计算矩阵V，D

FEMSpace
输入：femmesh，fem阶数
1. 初始化
1.1 网格格式转换
1.2 申请空间

2. 积分核函数 通过网格信息，求取各个积分
2.1 VolumeMassIntegrator()
2.2 VolumeStiffIntegrator()

3. 通量核函数
3.1 计算变量在边界的值
3.15 通量计算（application层面完成）
3.2 f(u+,u-) 

4. 边界条件核函数 // 本应在fempatch，但由于暂时处理简单边界，不拆出
4.1 边界条件初始化(bcinit)
4.2 边界更改(bcchange)
4.2 边界条件作用(bcupdate)

5. 函数计算
5.1 在求积点上计算函数
*/

#pragma once 

#include "femmesh.hpp"
#include "numerics.hpp"

#include "../eigen-3.4.0/Eigen/Dense"

#include <vector>

namespace femspace
{
    using namespace Eigen;

    class FE
    {
        public:
            FE(int order, double xl, double xr);
            ~FE() = default;

            int Np_; // 插值点数 >= 1
            int order_;

            // 几何
            const double xl_; // 左侧顶点坐标
            const double xr_; // 右侧顶点坐标
            const double dx_; // 网格大小 = xr - xl

            // 插值点
            VectorXd node_coords_; // [Np] 插值点坐标
            VectorXd node_values_; // [Np] 插值函数值
    };

    class REFFE
    {
        // default interval: [-1, 1]
        public:
            REFFE(int order);
            ~REFFE() = default;

            int order_; // 阶数
            int Np_; // 插值点数 = order + 1

            // 求积点
            VectorXd node_coords_; // [Np] Gauss-Lobatto 求积点/插值点坐标
            VectorXd node_weights_; // [Np] Gauss-Lobatto 求积权重

            // 多项式系数
            MatrixXd legendre_coefs_; // [Np, Np] Legendre多项式 系数矩阵

            // 预计算矩阵
            MatrixXd Vandermonde_; // [Np, Np] 插值矩阵
            MatrixXd Vandinv_; // [Np, Np] Vandermonde矩阵逆
            MatrixXd Dx_; // [Np, Np] 微分矩阵, uh' = Dx * uh
            // 外插点
            VectorXd lagrange_left_values_; // [Nx] 左侧插值点值
            VectorXd lagrange_right_values_; // [Nx] 右侧插值点值

            // 质量，刚度矩阵
            MatrixXd mass_, stiff_; // [Np, Np] 质量矩阵，刚度矩阵

            void ComputeMass();
            void ComputeStiff();
            void ComputeBCExtrap();
            void ComputeLagrangeExtrapolation(double xnew, FE fe, VectorXd& lagrange_values);
            void ComputeFuncExtrapolation(double xnew, FE fe, VectorXd& func_values, double& fnew);

            void MapToActualFE(FE& fe);

        private:
            void ComputeLegendreCoefficients();

            void ComputeVandermonde();

            void ComputeDx();
    };

    class FEMSpace
    {
        public:
            FEMSpace(femmesh::FEMMesh& mesh, int order); // Np = order - 1;
            ~FEMSpace() = default;

            // 积分核函数
            void VolumeMassIntegrator();
            void VolumeStiffIntegrator();
            void FaceMassIntegrator(); 

            // 通量核函数
            void ComputeBCFuncValues();
            void ComputeFlux(double flux_in, double flux_out);

            // 边界条件核函数
            void BCInit();
            void BCChange(double left_bc_value, double right_bc_value);
            void BCUpdate();

            // 插值核函数
            //void Interpolate(double (*func)(double x));
            void Interpolate(double (*func)(double x), double xl, double xr, VectorXd& func_values);
            std::vector<VectorXd> Interpolate(double (*func)(double x));

            // 求积式
            //void Quaduature(const VectorXd& func_values,double xl, double xr, double& result);
            void Quaduature(const VectorXd& func_values, int elem_idx, double& result);

            // getters
            femmesh::FEMMesh* GetMesh() const { return mesh_; }
            int GetOrder() const { return order_; }
            int GetNp() const { return Np_; }
            int GetNx() const { return Nx_; }
            std::vector<FE>& GetElems() { return elems_; }
            REFFE* GetRefElem() { return ref_elem_; }
            std::vector<VectorXd>& GetElemBCValues() { return elem_bc_values_; }
            std::vector<MatrixXd>& GetMass() { return mass_; }
            std::vector<MatrixXd>& GetStiff() { return stiff_; }
            MatrixXd GetFaceMassLeft() const { return face_mass_left_; }
            MatrixXd GetFaceMassRight() const { return face_mass_right_; }
            std::vector<VectorXd>& GetQuadWeights() { return quad_weights_; }
            VectorXd GetLagrangeLeftValues() const { return ref_elem_->lagrange_left_values_; }
            VectorXd GetLagrangeRightValues() const { return ref_elem_->lagrange_right_values_; }

        private:
            // 网格
            femmesh::FEMMesh* mesh_; // 网格

            int Nx_; // 单元数目

            // 矩阵与右端项
            int Np_; 
            int order_; 
            
            // 参考元
            REFFE* ref_elem_; // 参考元

            // 全局系统
            std::vector<FE> elems_; // [Nx] 单元
            std::vector<VectorXd> elem_bc_values_; // [Nx, 2] 单元边界值
            std::vector<MatrixXd> mass_;   // [Nx, Np, Np] 质量矩阵
            std::vector<MatrixXd> stiff_;  // [Nx, Np, Np] 刚度矩阵
            std::vector<VectorXd> rhs_; // [Nx, Np] 右端项
            MatrixXd face_mass_left_; // [Np, Np]  // 左侧面质量矩阵
            MatrixXd face_mass_right_; // [Np, Np] // 右侧面质量矩阵
            
            // 物理边界条件
            double left_bc_value_; // 左侧边界值
            double right_bc_value_; // 右侧边界值

            // 求积权重
            std::vector<VectorXd> quad_weights_; // [Nx, Np] 求积权重

            // 初始化参考元
            void InitRefElem();
    };
}

